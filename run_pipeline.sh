#!/usr/bin/env bash
# ichorCNA pipeline — streams CRAMs from S3, no local copy needed
#
# Usage: bash run_pipeline.sh <ref> <outdir> [sample_list] [downsample_frac] [bin_size] [normal_panel]
#
#   ref              path to reference FASTA (must have .fai index)   [REQUIRED]
#   outdir           output root directory                             [REQUIRED]
#   sample_list      S3 CRAM URIs, one per line  (default: samples.txt next to script)
#   downsample_frac  fraction kept by samtools -s; e.g. 0.01 = 100x→1x
#                    (default: 1.0 — no downsampling)
#   bin_size         genomic bin size in bp      (default: 500000)
#   normal_panel     path to .rds PoN file, or NULL to skip
#                    (default: bundled hg38 PoN matching bin_size)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 2 ]]; then
    echo "Usage: bash run_pipeline.sh <ref> <outdir> [sample_list] [downsample_frac] [bin_size] [normal_panel]" >&2
    exit 1
fi

REF="${1}"
OUTDIR="${2}"
SAMPLE_LIST="${3:-${SCRIPT_DIR}/samples.txt}"
DOWNSAMPLE_FRAC="${4:-1.0}"
BIN_SIZE="${5:-500000}"
# Derive GC / map / PoN paths from bin size automatically
case "${BIN_SIZE}" in
  10000)   BIN_LABEL="10kb"  ;;
  50000)   BIN_LABEL="50kb"  ;;
  500000)  BIN_LABEL="500kb" ;;
  1000000) BIN_LABEL="1000kb";;
  *)       echo "ERROR: unsupported bin size ${BIN_SIZE}. Choose 10000/50000/500000/1000000" >&2; exit 1 ;;
esac

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO=${SCRIPT_DIR}
CONDA_ENV=ichorCNA

# ── Read-counting settings ───────────────────────────────────────────────────
# Only primary chromosomes — CRAI index lets samtools fetch these directly from
# S3 without downloading the entire file
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
CHRS_CSV="${CHRS// /,}"  # comma-separated list
QUAL=20
NCPU=$(nproc 2>/dev/null || echo 4)
# Parallelize across chromosomes — S3 single-stream bandwidth is the bottleneck,
# so fan out one samtools view per chromosome (each opens its own S3 connection)
# and combine the per-chrom WIGs at the end. -@ per samtools is small because
# the work is I/O-bound, not CPU-bound.
SAMTOOLS_THREADS=2

BAM_TO_WIG=${REPO}/scripts/bam_to_wig.py

# ── ichorCNA settings (hg38 / UCSC) ──────────────────────────────────────────
GC_WIG=${REPO}/inst/extdata/gc_hg38_${BIN_LABEL}.wig
MAP_WIG=${REPO}/inst/extdata/map_hg38_${BIN_LABEL}.wig
CENTROMERE=${REPO}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
# PoN: use bundled hg38 PoN if bin size matches (500kb or 1Mb), else NULL
if [[ -z "${6:-}" ]]; then
  case "${BIN_SIZE}" in
    500000)  NORMAL_PANEL=${REPO}/inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds ;;
    1000000) NORMAL_PANEL=${REPO}/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds  ;;
    *)       NORMAL_PANEL=NULL ;;
  esac
else
  NORMAL_PANEL="${6}"
fi

ICHORCNA_SCRIPT=${REPO}/scripts/runIchorCNA.R

# ── Resolve binaries from conda env ──────────────────────────────────────────
RSCRIPT=$(conda run -n "${CONDA_ENV}" which Rscript 2>/dev/null)
PYTHON=$(conda run -n "${CONDA_ENV}" which python3 2>/dev/null)

# ── Inject AWS credentials into environment for samtools S3 access ────────────
# --format env outputs "export KEY=VALUE" lines, making them available to all
# child processes (samtools background jobs, conda run subprocesses, etc.)
eval "$(aws configure export-credentials --format env)"

mkdir -p "${OUTDIR}/readDepth" "${OUTDIR}/ichorCNA" "${OUTDIR}/logs"

DOWNSAMPLE_LABEL=$(awk -v f="${DOWNSAMPLE_FRAC}" 'BEGIN{if(f+0>=1.0) print "none"; else print f}')
echo "=== ichorCNA pipeline started: $(date) ==="
echo "    Samples     : ${SAMPLE_LIST}"
echo "    Output      : ${OUTDIR}"
echo "    Bin size    : ${BIN_SIZE} (${BIN_LABEL})"
echo "    Downsample  : ${DOWNSAMPLE_LABEL}"
echo "    PoN         : ${NORMAL_PANEL}"
echo ""

# ── Per-sample loop ───────────────────────────────────────────────────────────
while IFS= read -r CRAM_S3 || [[ -n "${CRAM_S3}" ]]; do
    # Skip blank lines and comments
    [[ -z "${CRAM_S3}" || "${CRAM_S3}" =~ ^[[:space:]]*# ]] && continue

    SAMPLE=$(basename "${CRAM_S3}" .cram)
    WIG=${OUTDIR}/readDepth/${SAMPLE}.bin${BIN_SIZE}.wig
    SAMPLE_OUTDIR=${OUTDIR}/ichorCNA/${SAMPLE}
    LOG=${OUTDIR}/logs/${SAMPLE}.log

    mkdir -p "${SAMPLE_OUTDIR}"
    echo "[$(date '+%F %T')] ── ${SAMPLE}" | tee -a "${LOG}"

    # ── Step 1: CRAM → WIG via per-chromosome streaming pipes ───────────────
    # Bottleneck is S3 single-stream bandwidth, not CPU. Fan out one
    # samtools view per chromosome (each = its own S3 connection) so total
    # wall is ~max(per-chrom time) instead of sum. Each pipe writes a
    # per-chrom partial WIG; we cat them in genome order at the end.
    # -s downsamples reads (e.g. 0.01 = keep 1% → 100x becomes ~1x).
    SUBSAMPLE_FLAG=$(awk -v f="${DOWNSAMPLE_FRAC}" 'BEGIN{if(f+0<1.0) print "-s"" "f}')
    PARTIAL_DIR=${OUTDIR}/readDepth/${SAMPLE}.parts
    mkdir -p "${PARTIAL_DIR}"
    rm -f "${PARTIAL_DIR}"/*.wig "${PARTIAL_DIR}"/*.log

    echo "[$(date '+%F %T')]   launching ${SAMTOOLS_THREADS}-thread samtools per chrom (parallel)" | tee -a "${LOG}"
    for CHR in ${CHRS}; do
        (
            samtools view -@ "${SAMTOOLS_THREADS}" ${SUBSAMPLE_FLAG} -T "${REF}" "${CRAM_S3}" "${CHR}" 2>> "${PARTIAL_DIR}/${CHR}.log" \
                | ${PYTHON} -u "${BAM_TO_WIG}" -w "${BIN_SIZE}" -q "${QUAL}" -c "${CHR}" --fai "${REF}.fai" \
                > "${PARTIAL_DIR}/${CHR}.wig" 2>> "${PARTIAL_DIR}/${CHR}.log"
        ) &
    done
    wait

    : > "${WIG}"
    for CHR in ${CHRS}; do
        if [[ -s "${PARTIAL_DIR}/${CHR}.wig" ]]; then
            cat "${PARTIAL_DIR}/${CHR}.wig" >> "${WIG}"
        else
            echo "WARNING: empty WIG for ${CHR} — see ${PARTIAL_DIR}/${CHR}.log" | tee -a "${LOG}"
        fi
    done
    for CHR in ${CHRS}; do
        if [[ -s "${PARTIAL_DIR}/${CHR}.log" ]]; then
            echo "─── ${CHR} ───" >> "${LOG}"
            cat "${PARTIAL_DIR}/${CHR}.log" >> "${LOG}"
        fi
    done

    echo "[$(date '+%F %T')]   WIG written: ${WIG}" | tee -a "${LOG}"

    # ── Step 2: ichorCNA ─────────────────────────────────────────────────────
    # Notes for 100x cfDNA:
    #   - lambda is auto-estimated from data variance (appropriate for high depth)
    #   - normalPanel=NULL because the bundled PoN is calibrated for 0.1x ULP
    #   - 50kb bins give ~60k data points per sample (well within R memory)
    conda run --no-capture-output -n "${CONDA_ENV}" \
        Rscript "${ICHORCNA_SCRIPT}" \
            --id "${SAMPLE}" \
            --libdir "${REPO}" \
            --WIG "${WIG}" \
            --gcWig "${GC_WIG}" \
            --mapWig "${MAP_WIG}" \
            --normalPanel "${NORMAL_PANEL}" \
            --genomeBuild hg38 \
            --genomeStyle UCSC \
            --chrs 'c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")' \
            --chrTrain 'c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")' \
            --chrNormalize 'c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")' \
            --centromere "${CENTROMERE}" \
            --ploidy "c(2,3)" \
            --normal "c(0.5,0.6,0.7,0.8,0.9)" \
            --maxCN 5 \
            --includeHOMD FALSE \
            --scStates "c(1,3)" \
            --estimateNormal TRUE \
            --estimatePloidy TRUE \
            --estimateScPrevalence TRUE \
            --txnE 0.9999 \
            --txnStrength 10000 \
            --minMapScore 0.75 \
            --plotFileType pdf \
            --plotYLim "c(-2,4)" \
            --outDir "${SAMPLE_OUTDIR}/" \
        >> "${LOG}" 2>&1

    echo "[$(date '+%F %T')]   ichorCNA done: ${SAMPLE_OUTDIR}" | tee -a "${LOG}"
    echo ""

done < "${SAMPLE_LIST}"

echo "=== All samples complete: $(date) ==="
