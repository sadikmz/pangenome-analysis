#!/usr/bin/env bash
# ============================================================================
# 02_run_orthofinder.sh
# ---------------------------------------------------------------------------
# Infer orthogroups across the Aspergillus proteome panel with OrthoFinder.
#
# OrthoFinder pipeline:
#   1. DIAMOND all-vs-all protein search (every proteome vs every proteome)
#   2. Normalised bit-score graph + MCL clustering  -> orthogroups
#   3. Per-orthogroup MSA (MUSCLE) + gene trees (FastTree)
#   4. Species tree from single-copy orthogroups (STAG/STRIDE), rooted
#
# Key outputs (under proteomes/OrthoFinder/Results_<date>/):
#   Orthogroups/Orthogroups.GeneCount.tsv   gene counts per species per OG
#   Orthogroups/Orthogroups.tsv             member gene IDs per OG
#   Orthogroups/Orthogroups_UnassignedGenes.tsv  singletons (species-specific)
#   Comparative_Genomics_Statistics/        summary stats
#   Species_Tree/SpeciesTree_rooted.txt     species tree (Newick)
#
# Environment: conda env with orthofinder, diamond, mcl, muscle, fasttree.
# Usage:
#   bash 02_run_orthofinder.sh <proteome_dir> [threads]
# ============================================================================
set -euo pipefail

PROTEOME_DIR="${1:-../data/proteomes}"
THREADS="${2:-8}"

# -S diamond : sequence search program
# -M msa     : infer gene trees from multiple sequence alignments (more accurate
#              than the default dendroblast); -A/-T pick the MSA and tree tools
# -t         : parallel threads for search; -a : parallel analysis threads
orthofinder \
    -f "${PROTEOME_DIR}" \
    -t "${THREADS}" -a 2 \
    -S diamond \
    -M msa -A muscle -T fasttree

echo "OrthoFinder finished. Results in ${PROTEOME_DIR}/OrthoFinder/Results_*/"
