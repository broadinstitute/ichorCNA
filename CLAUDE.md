# ichorCNA — Developer & Agent Continuation Guide

This file documents the environment setup, custom scripts, and current pipeline state for running ichorCNA on cfDNA CRAM files streamed directly from S3.  It is written for both human developers and AI coding agents (Claude, Copilot, etc.) continuing this work.

---

## Project Overview

We extended ichorCNA to:
1. **Stream 100× cfDNA CRAMs directly from S3** — no local download of ~600 GB files.
2. **Downsample on the fly** (e.g. 100×→1× with `-s 0.01`) so ichorCNA's ULP-WGS model applies.
3. **Replace `readCounter`** (requires a BAM index) with a custom stdin-reading WIG generator (`scripts/bam_to_wig.py`).
4. **Auto-select GC/map reference WIGs and Panel-of-Normals** based on the requested bin size.

The two analysis modes we target:

| Mode | Downsample | Bin size | PoN |
|------|-----------|---------|-----|
| Full depth (debug/QC) | none | 50 kb | none |
| 1× equivalent | 0.01 | 500 kb | bundled hg38 HD-ULP PoN |

---

## Environment

### Host

| Item | Value |
|------|-------|
| OS | Ubuntu (Linux x86_64) |
| Ref genome | `/proj/ref/Homo_sapiens_assembly38.fasta` (hg38, `.fai` + `.dict` present) |
| Output storage | `/data/Runs/` — 2 TB NVMe, ~674 GB free.  **Do NOT write to `/` (only 14 GB free).** |

### Conda environment — `ichorCNA`

Located at `/home/ubuntu/miniconda3/envs/ichorCNA`.  Activate with:

```bash
conda activate ichorCNA
# or prefix commands with:
conda run -n ichorCNA <command>
```

**Key packages installed:**

| Package | Notes |
|---------|-------|
| `r-base 4.5.3` | R runtime |
| `r-optparse` | CLI parsing for `runIchorCNA.R` |
| `bioconductor-hmmcopy` | Core HMM library |
| `bioconductor-genomicranges` | GRanges support |
| `bioconductor-genomeinfodb` | Seqinfo / genome builds |
| `r-plyr` | Data manipulation (ichorCNA dep) |
| `ichorCNA` | Installed from local repo via `remotes::install_local("/proj/src/ichorCNA")` |
| `hmmcopy` | Provides `readCounter` binary (not used in pipeline, but installed) |

**To reinstall ichorCNA package after code changes:**

```bash
conda run -n ichorCNA Rscript -e 'remotes::install_local("/proj/src/ichorCNA", upgrade="never", force=TRUE)'
```

### samtools

`/home/ubuntu/miniconda3/bin/samtools` — version 1.22, compiled with htslib **S3 support** (`S3=yes`).  This is what allows `samtools view s3://...` without pre-downloading.

### AWS credentials

CRAMs are in a private S3 bucket.  Credentials **must** be injected before running the pipeline.  The pipeline script does this automatically, but if running samtools manually:

```bash
# Inject into current shell (exports vars to all child processes)
eval "$(aws configure export-credentials --format env)"

# Now samtools can access S3:
samtools view -T /proj/ref/Homo_sapiens_assembly38.fasta \
    s3://bucket/path/sample.cram chr1 | head
```

> **Critical**: Use `--format env` (not `--format env-no-export`).  The `env-no-export` variant does NOT propagate credentials to child processes (subshells, conda run, python pipes, etc.).

There is also an alias `st` in the shell that handles single samtools calls, but **do not use it in scripts** — it only sets credentials for one call.

---

## New Files Added

### `run_pipeline.sh`

Main pipeline orchestrator.  Streams each CRAM from S3, generates a WIG file, then runs ichorCNA.

**Usage:**

```bash
bash run_pipeline.sh [sample_list] [downsample_frac] [outdir] [bin_size] [normal_panel]
```

| Argument | Default | Notes |
|----------|---------|-------|
| `sample_list` | `samples.txt` (next to script) | S3 CRAM URIs, one per line |
| `downsample_frac` | `1.0` | `0.01` = keep 1% of reads (100×→1×) |
| `outdir` | `/data/Runs/ichorCNA` | Must be on `/data/` not `/` |
| `bin_size` | `500000` | Supported: 10000 / 50000 / 500000 / 1000000 |
| `normal_panel` | auto | Path to `.rds` PoN, or `NULL`.  Auto-selected for 500 kb and 1 Mb bins |

**Examples:**

```bash
# 1× equivalent run with PoN (recommended for cfDNA CNA calling)
bash run_pipeline.sh /path/to/ref.fasta /data/Runs/ichorCNA_1x_500kb samples.txt 0.01 500000

# Full-depth run without PoN, finer bins (QC / debugging)
bash run_pipeline.sh /path/to/ref.fasta /data/Runs/ichorCNA_fullDepth samples.txt 1.0 50000

# Single sample test
echo "s3://bucket/.../sample.cram" > /tmp/test_sample.txt
bash run_pipeline.sh /path/to/ref.fasta /data/Runs/ichorCNA_test /tmp/test_sample.txt 0.01 500000
```

**Output layout:**

```
outdir/
  readDepth/<sample>.bin<size>.wig      ← read count WIG (Step 1 output)
  ichorCNA/<sample>/<sample>.*.pdf      ← genome-wide CNA plots
  ichorCNA/<sample>/<sample>.seg.txt    ← segment calls
  ichorCNA/<sample>/<sample>.params.txt ← estimated ploidy / tumour fraction
  logs/<sample>.log                     ← combined samtools + bam_to_wig + R output
```

**How streaming works:**

```
S3 CRAM
  └─ samtools view -T ref.fa [[-s frac]] s3://...cram chr1 chr2 ...
        │  (CRAI index fetches only requested chromosomes — no full download)
        │ SAM text on stdout
        ▼
     bam_to_wig.py  (reads stdin, emits WIG to stdout)
        │
        ▼
     WIG file
        ▼
     runIchorCNA.R
```

### `scripts/bam_to_wig.py`

Reads SAM text from stdin, counts reads per fixed-size genomic bin, writes WIG to stdout.  Replaces `readCounter` (which requires a BAM/CRAM index file and cannot read from a pipe).

```bash
samtools view ... | python3 -u scripts/bam_to_wig.py \
    -w 500000 -q 20 -c chr1,chr2,...,chrX
```

> **Always use `python3 -u`** (unbuffered) so progress lines appear in the log in real time.  The pipeline script does this automatically.

**Progress output** (written to stderr → appears in `.log`):

```
  [      ] starting chr1 ...
  [  42s] chr1    5.0M reads |   4.8M kept | 0.12M reads/s
  [  84s] finished chr1   |   9.6M reads total |   9.2M kept | 0.11M reads/s
  ...
  Done: 230.4M reads in 1820s (0.13M reads/s) | 1140 bins populated
```

**Flag filtering** — skips: unmapped (0x4), secondary (0x100), QC fail (0x200), duplicate (0x400), supplementary (0x800).

---

## Current Status (as of session end)

| Task | Status |
|------|--------|
| Conda env created + all R/bioc deps installed | ✅ Done |
| `run_pipeline.sh` created and debugged | ✅ Done |
| `scripts/bam_to_wig.py` created | ✅ Done |
| AWS credential injection working | ✅ Done |
| Fork `broadinstitute/ichorCNA` → `doron-st/ichorCNA` | ✅ Done |
| Push `run_pipeline.sh`, `bam_to_wig.py`, `CLAUDE.md` | ✅ Done |
| Test run on first sample (1×/500 kb/PoN) | 🔄 Was running; may need restart |
| Validate ichorCNA R step output | ❌ Not yet done |
| Full run on all 10 samples | ❌ Not yet done |

---

## Resuming the Pipeline

### Check if a run is still active

```bash
ps aux | grep run_pipeline
# or
screen -ls       # if launched in a screen session
tmux ls          # if launched in tmux
```

### Restart a single-sample test

```bash
echo "s3://your-bucket/path/to/sample.cram" > /tmp/test_sample.txt

bash /proj/src/ichorCNA/run_pipeline.sh \
    /proj/ref/Homo_sapiens_assembly38.fasta \
    /data/Runs/ichorCNA_1x_500kb \
    /tmp/test_sample.txt \
    0.01 \
    500000
```

### Monitor progress

```bash
# Follow log for the test sample (progress updates every 1M reads)
tail -f /data/Runs/ichorCNA_1x_500kb/logs/<sample-name>.log

# Check WIG was created
ls -lh /data/Runs/ichorCNA_1x_500kb/readDepth/
```

### Run all 10 samples

After the test sample validates successfully:

```bash
# Provide your own sample list — DO NOT commit files with internal S3 paths
bash /proj/src/ichorCNA/run_pipeline.sh \
    /path/to/ref.fasta \
    /path/to/outdir \
    /path/to/your/samples.txt \
    0.01 \
    500000
```

> Each sample takes roughly 30–60 minutes (streaming ~1% of a 600 GB CRAM at ~100 MB/s S3 throughput).

---

## Known Issues & Decisions

### Why not use `readCounter`?

`readCounter` (from HMMcopy) requires a BAM/CRAM **index** to be present locally.  When streaming from S3, there is no local index file.  `bam_to_wig.py` reads SAM text directly from stdin and requires no index.

### Why 0.01 downsampling (100×→1×)?

ichorCNA was designed for ultra-low-pass (0.1× coverage) WGS.  The bundled Panel-of-Normals and model parameters are calibrated for that depth.  At 100×, variance estimates and normalisation break down.  Downsampling to ~1× restores the ULP-WGS statistical regime.

> **Note:** `samtools view -s` downsampling is performed *after* CRAM decompression.  It does NOT reduce S3 bandwidth — the full set of reads for the requested chromosomes is still transferred.  If bandwidth cost matters, consider generating pre-downsampled BAMs upstream.

### Why 500 kb bins with PoN?

The bundled hg38 PoN files are only available at 500 kb and 1 Mb resolution.  50 kb bins would require a custom PoN built from matched normal samples.  500 kb bins with the bundled HD-ULP PoN is the most straightforward starting point.

### Sample list — do not commit

Your sample list contains internal S3 URIs.  Keep it outside the repo (e.g. `/proj/ichorCNA/samples/`) or ensure it matches the `*_samples.txt` pattern in `.gitignore`.  Provide one S3 CRAM URI per line.

---

## Installation from Scratch

If you need to rebuild the environment on a new machine:

```bash
# 1. Create conda environment
conda create -n ichorCNA -c conda-forge -c bioconda \
    r-base=4.5 r-optparse bioconductor-hmmcopy \
    bioconductor-genomicranges bioconductor-genomeinfodb \
    r-plyr hmmcopy -y

# 2. Install ichorCNA R package from this repo
conda run -n ichorCNA Rscript -e \
    'remotes::install_local("/proj/src/ichorCNA", upgrade="never")'

# 3. Verify samtools has S3 support
samtools --version | grep htslib
# should see "S3 support: yes" or similar; if not, install via conda:
# conda install -c bioconda samtools

# 4. Configure AWS credentials (once)
aws configure
# or use instance role / environment variables

# 5. Verify S3 access
eval "$(aws configure export-credentials --format env)"
samtools view -c -T /proj/ref/Homo_sapiens_assembly38.fasta \
    s3://your-bucket/path/sample.cram chr22
```

---

## Repository Layout (additions)

```
ichorCNA/
├── CLAUDE.md                  ← this file
├── run_pipeline.sh            ← main S3-streaming pipeline script (NEW)
├── scripts/
│   ├── bam_to_wig.py          ← SAM-stdin → WIG converter (NEW)
│   └── runIchorCNA.R          ← original ichorCNA R entry point (unchanged)
├── inst/extdata/
│   ├── gc_hg38_*.wig          ← GC reference WIGs (10kb/50kb/500kb/1000kb)
│   ├── map_hg38_*.wig         ← mappability reference WIGs
│   ├── HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds  ← bundled PoN
│   └── HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds
└── R/                         ← ichorCNA R library (unchanged)
```
