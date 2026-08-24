#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 RAW_R1.fastq.gz RAW_R2.fastq.gz OUTDIR THREADS" >&2
    exit 2
fi

RAW_R1=$1
RAW_R2=$2
OUTDIR=$3
THREADS=$4

if [ ! -f "$RAW_R1" ] || [ ! -f "$RAW_R2" ]; then
    echo "Input FASTQ missing" >&2
    exit 2
fi

case "$THREADS" in
    ''|*[!0-9]*) echo "THREADS must be a positive integer" >&2; exit 2 ;;
esac
if [ "$THREADS" -lt 1 ]; then
    echo "THREADS must be >= 1" >&2
    exit 2
fi

if ! command -v fastp >/dev/null 2>&1; then
    echo "fastp not found" >&2
    exit 2
fi
if ! command -v spades.py >/dev/null 2>&1; then
    echo "spades.py not found" >&2
    exit 2
fi

FASTP_VERSION=$(fastp --version 2>&1 | head -n 1)
SPADES_VERSION=$(spades.py --version 2>&1 | head -n 1)

case "$FASTP_VERSION" in
    *"1.3.6"*) ;;
    *) echo "Expected fastp v1.3.6; observed: $FASTP_VERSION" >&2; exit 3 ;;
esac
case "$SPADES_VERSION" in
    *"4.3.0"*) ;;
    *) echo "Expected SPAdes v4.3.0; observed: $SPADES_VERSION" >&2; exit 3 ;;
esac

if [ -e "$OUTDIR" ]; then
    echo "Output path already exists: $OUTDIR" >&2
    exit 2
fi

mkdir -p "$OUTDIR/raw_checksums" "$OUTDIR/reads" "$OUTDIR/logs"

sha256sum "$RAW_R1" "$RAW_R2" > "$OUTDIR/raw_checksums/raw_fastq.sha256"
printf '%s\n' "$FASTP_VERSION" > "$OUTDIR/logs/fastp.version.txt"
printf '%s\n' "$SPADES_VERSION" > "$OUTDIR/logs/spades.version.txt"

CLEAN_R1="$OUTDIR/reads/clean_R1.fastq.gz"
CLEAN_R2="$OUTDIR/reads/clean_R2.fastq.gz"
UNPAIRED_R1="$OUTDIR/reads/unpaired_R1.fastq.gz"
UNPAIRED_R2="$OUTDIR/reads/unpaired_R2.fastq.gz"

fastp \
  --in1 "$RAW_R1" \
  --in2 "$RAW_R2" \
  --out1 "$CLEAN_R1" \
  --out2 "$CLEAN_R2" \
  --unpaired1 "$UNPAIRED_R1" \
  --unpaired2 "$UNPAIRED_R2" \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required 30 \
  --thread 4 \
  --json "$OUTDIR/fastp.json" \
  --html "$OUTDIR/fastp.html" \
  2> "$OUTDIR/logs/fastp.stderr.log"

sha256sum \
  "$CLEAN_R1" "$CLEAN_R2" \
  "$UNPAIRED_R1" "$UNPAIRED_R2" \
  > "$OUTDIR/reads/processed_fastq.sha256"

spades.py \
  --isolate \
  -1 "$CLEAN_R1" \
  -2 "$CLEAN_R2" \
  -t "$THREADS" \
  -o "$OUTDIR/spades" \
  > "$OUTDIR/logs/spades.stdout.log" \
  2> "$OUTDIR/logs/spades.stderr.log"

REQUIRED=(
  "$OUTDIR/spades/assembly_graph_with_scaffolds.gfa"
  "$OUTDIR/spades/contigs.paths"
  "$OUTDIR/spades/scaffolds.paths"
  "$OUTDIR/spades/contigs.fasta"
  "$OUTDIR/spades/scaffolds.fasta"
  "$OUTDIR/spades/params.txt"
  "$OUTDIR/spades/spades.log"
)

for path in "${REQUIRED[@]}"; do
    if [ ! -f "$path" ]; then
        echo "Required SPAdes output missing: $path" >&2
        exit 4
    fi
done

sha256sum "${REQUIRED[@]}" > "$OUTDIR/spades/canonical_outputs.sha256"

echo "G4 assembly workflow complete: $OUTDIR"
