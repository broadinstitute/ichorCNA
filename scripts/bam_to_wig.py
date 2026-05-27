#!/usr/bin/env python3
"""
Stream SAM/BAM records from stdin and output a fixed-step WIG read-count file.
Drop-in replacement for HMMcopy readCounter when no BAM index is available
(e.g. when streaming from S3 via a named pipe or process substitution).

Usage:
    samtools view -T ref.fa file.cram chr1 chr2 ... \
        | bam_to_wig.py -w 50000 -q 20 -c chr1,chr2,...

Output:
    WIG file to stdout (variableStep per chromosome, one value per bin).
"""
import sys
import argparse
from collections import defaultdict, OrderedDict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-w", "--window", type=int, default=1_000_000,
                   help="Bin size in bp (default: 1000000)")
    p.add_argument("-q", "--quality", type=int, default=0,
                   help="Minimum mapping quality (default: 0)")
    p.add_argument("-c", "--chromosomes", default=None,
                   help="Comma-separated list of chromosomes to output. "
                        "Others are still counted but not emitted. "
                        "If omitted, all chromosomes in the header are output.")
    return p.parse_args()


def main():
    args = parse_args()
    keep_chrs = set(args.chromosomes.split(",")) if args.chromosomes else None

    chr_lengths = OrderedDict()   # chromosome → length (from SAM header)
    counts = defaultdict(int)     # (chrom, bin_start) → read count

    # SAM flag bits to skip: unmapped, secondary, QC fail, duplicate, supplementary
    SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

    total_read = 0
    kept = 0
    current_chr = None
    import time
    t0 = time.time()

    for line in sys.stdin:
        if line.startswith("@"):
            # Parse @SQ header lines for chromosome lengths
            if line.startswith("@SQ"):
                fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
                chrom = fields.get("SN", "")
                length = int(fields.get("LN", 0))
                if chrom and length:
                    chr_lengths[chrom] = length
            continue

        parts = line.split("\t", 12)
        if len(parts) < 5:
            continue

        flag = int(parts[1])
        if flag & SKIP_FLAGS:
            continue

        mapq = int(parts[4])
        if mapq < args.quality:
            continue

        chrom = parts[2]
        if chrom == "*":
            continue
        if keep_chrs and chrom not in keep_chrs:
            continue

        # New chromosome → print a completion line then announce the new one
        if chrom != current_chr:
            if current_chr is not None:
                elapsed = time.time() - t0
                sys.stderr.write(
                    f"\r  [{elapsed:6.0f}s] finished {current_chr:<6} | "
                    f"{total_read/1e6:6.1f}M reads total | "
                    f"{kept/1e6:5.1f}M kept | "
                    f"{total_read/elapsed/1e6:.2f}M reads/s\n"
                )
                sys.stderr.flush()
            current_chr = chrom
            sys.stderr.write(f"  [      ] starting {current_chr} ...\r")
            sys.stderr.flush()

        total_read += 1

        # Rolling update every 1M reads within a chromosome
        if total_read % 1_000_000 == 0:
            elapsed = time.time() - t0
            sys.stderr.write(
                f"\r  [{elapsed:6.0f}s] {current_chr:<6} "
                f"{total_read/1e6:6.1f}M reads | "
                f"{kept/1e6:5.1f}M kept | "
                f"{total_read/elapsed/1e6:.2f}M reads/s   "
            )
            sys.stderr.flush()

        pos = int(parts[3])  # 1-based leftmost position
        bin_start = ((pos - 1) // args.window) * args.window + 1
        counts[(chrom, bin_start)] += 1
        kept += 1

    # ── Emit WIG ──────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    sys.stderr.write(
        f"\n  Done: {total_read/1e6:.1f}M reads in {elapsed:.0f}s "
        f"({total_read/elapsed/1e6:.2f}M reads/s) | {kept} bins populated\n"
    )
    # Iterate chromosomes in header order so ichorCNA gets a consistent sort.
    emit_chrs = [c for c in chr_lengths if keep_chrs is None or c in keep_chrs]
    # Also include any chromosomes seen in reads but not in the header
    for chrom, _ in counts:
        if chrom not in chr_lengths:
            emit_chrs.append(chrom)

    for chrom in emit_chrs:
        length = chr_lengths.get(chrom, 0)
        # Collect all bins for this chromosome and sort them
        bins = sorted(b for (c, b) in counts if c == chrom)
        if not bins:
            continue
        sys.stdout.write(f"variableStep chrom={chrom} span={args.window}\n")
        for bin_start in bins:
            sys.stdout.write(f"{bin_start}\t{counts[(chrom, bin_start)]}\n")


if __name__ == "__main__":
    main()
