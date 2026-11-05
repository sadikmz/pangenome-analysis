#!/bin/bash -l
# ---------------------------------------------------------------------------
# minimap2_syri.sh
# Run pairwise minimap2 alignments + syri + plotsr on a sequential genome list.
#
# Usage:
#   bash minimap2_syri.sh [--no-filter] <genome_fofn> <threads> <output_prefix>
#
# Arguments:
#   genome_fofn    File of genome paths, one .fna per line (# lines ignored).
#                  Comparisons are made sequentially: line[i] (ref) vs line[i+1] (qry).
#   threads        Number of CPU threads for minimap2/samtools/syri.
#   output_prefix  Prefix for the final plotsr output PDF.
#
# Options:
#   --no-filter    Skip SV/indel filtering; pass raw syri.out to plotsr.
#                  Default: filter out SNP/INS/DEL/HDR/CPG/CPL/TDM lines.
#
# Intermediate outputs kept:
#   <qry>_<ref>.bam               minimap2 BAM
#   <qry>_<ref>_syri.out          syri structural rearrangements
#   <qry>_<ref>_syri.filtered.out filtered syri (no SNP/INS/DEL/HDR/CPG/CPL/TDM)
#                                  (only written when --no-filter is not set)
# ---------------------------------------------------------------------------

set -euo pipefail

# ---------------------------------------------------------------------------
# Parse optional flags
# ---------------------------------------------------------------------------
FILTER=true
while [[ "${1:-}" == --* ]]; do
    case "$1" in
        --no-filter) FILTER=false; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

FOFN="${1:?Usage: $0 [--no-filter] <genome_fofn> <threads> <output_prefix>}"
THREADS="${2:?Usage: $0 [--no-filter] <genome_fofn> <threads> <output_prefix>}"
FINAL_OUT="${3:?Usage: $0 [--no-filter] <genome_fofn> <threads> <output_prefix>}"

#module load minimap2/2.26-GCCcore-12.3.0
#module load samtools/1.18-GCC-12.3.0

# Read genome list (skip blank lines and comments)
GENOMES=()
while IFS= read -r line; do
    [[ "$line" =~ ^\s*# || -z "${line// }" ]] && continue
    GENOMES+=("$line")
done < "$FOFN"

n=${#GENOMES[@]}
if [[ $n -lt 2 ]]; then
    echo "Error: need at least 2 genomes in $FOFN (found $n)" >&2
    exit 1
fi

echo "Processing $n genome(s) -> $((n-1)) pairwise comparison(s)"
echo "Threads: $THREADS | Output prefix: $FINAL_OUT | Filter: $FILTER"
echo ""

# ---------------------------------------------------------------------------
# Step 1: minimap2 + syri + filter for each consecutive pair
# ---------------------------------------------------------------------------
filtered_outputs=()

for (( i=0; i<n-1; i++ )); do
    ref_path="${GENOMES[$i]}"
    qry_path="${GENOMES[$((i+1))]}"

    ref=$(basename "$ref_path" .fna)
    qry=$(basename "$qry_path" .fna)
    out="${qry}_${ref}"

    echo "[${qry} vs ${ref}] Aligning with minimap2 ..."
    minimap2 -ax asm5 -t "$THREADS" -k 19 --eqx \
        "$ref_path" "$qry_path" \
        | samtools view -Sb -F 4 \
            -e 'qname =~ "^chr_" && rname =~ "^chr_"' \
            -@"$THREADS" \
        | samtools sort -@"$THREADS" -o "${out}.bam"

    echo "[${qry} vs ${ref}] Running syri ..."
    syri -c "${out}.bam" \
         -r "$ref_path" \
         -q "$qry_path" \
         -F B \
         --prefix "${out}_" \
         --nc "$THREADS" \
         --hdrseq

    if [[ "$FILTER" == true ]]; then
        echo "[${qry} vs ${ref}] Filtering syri output (SNP/INS/DEL/HDR/CPG/CPL/TDM) ..."
        grep -v -E "CPG|CPL|DEL|HDR|INS|SNP|TDM" "${out}_syri.out"  \
            > "${out}_syri.filtered.out"
        filtered_outputs+=("${out}_syri.filtered.out")
    else
        echo "[${qry} vs ${ref}] Skipping filter (--no-filter set) ..."
        # grep 'chr.*chr' "${out}_syri.out" > "${out}_syri.chr.out"
        filtered_outputs+=("${out}_syri.out")
    fi
    echo ""
done

# ---------------------------------------------------------------------------
# Step 2: plotsr
# ---------------------------------------------------------------------------
echo "Running plotsr ..."

plotsr_args=()
for f in "${filtered_outputs[@]}"; do
    plotsr_args+=(--sr "$f")
done

plotsr \
    "${plotsr_args[@]}" \
    --genomes genomes.txt \
    --cfg base.cfg \
    --itx \
    -R \
    -H 6 \
    -W 7 \
    -d 1000 \
    -o "${FINAL_OUT}.pdf"

echo "Done. Output: ${FINAL_OUT}.pdf"
