# Aspergillus Pangenome & Comparative Genomics

Pan-gene / pangenome and comparative–evolutionary genomics analysis of six
*Aspergillus* reference proteomes, from RefSeq download through orthogroup
inference, pangenome partitioning, accumulation-curve modelling, and
phylogenomics.

## Species panel (RefSeq reference assemblies)

| Label        | Species                          | Accession        | Proteins |
|--------------|----------------------------------|------------------|----------|
| A_fumigatus  | *A. fumigatus* Af293             | GCF_000002655.1  | 9,823    |
| A_flavus     | *A. flavus*                      | GCF_009017415.1  | 13,707   |
| A_niger      | *A. niger*                       | GCF_000002855.4  | 14,058   |
| A_oryzae     | *A. oryzae* RIB40                | GCF_000184455.2  | 12,074   |
| A_nidulans   | *A. nidulans* FGSC A4            | GCF_000011425.1  | 10,453   |
| A_terreus    | *A. terreus* NIH2624             | GCF_000149615.1  | 10,401   |

## Directory layout

```
pangenome_project/
├── README.md                     this file
├── REPORT.md                     narrative results + interpretation
├── scripts/
│   ├── 01_download_proteomes.py  fetch RefSeq proteomes from NCBI Datasets
│   ├── 02_run_orthofinder.sh     OrthoFinder orthogroup inference
│   ├── 03_pangenome_analysis.py  pangenome partitioning + accumulation curves
│   └── 04_phylogenomics.py       gene-content NJ tree + species tree render
├── data/
│   ├── proteomes/*.faa           input proteomes
│   └── Orthogroups.txt           OrthoFinder MCL orthogroups
├── results/
│   ├── tables/*.csv              all matrices and summary tables
│   └── figures/*.png             all figures
└── logs/orthofinder.log          OrthoFinder run log
```

## Reproduce

```bash
# 1. Environment (conda / mamba)
mamba create -n pangenome -c bioconda -c conda-forge \
    orthofinder diamond mcl fastme muscle fasttree \
    biopython dendropy ete3 seaborn pandas scipy
conda activate pangenome

# 2. Run the pipeline in order
python scripts/01_download_proteomes.py
bash   scripts/02_run_orthofinder.sh data/proteomes 8
cp proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt data/
python scripts/03_pangenome_analysis.py
# copy species tree when OrthoFinder finishes its gene-tree step:
cp proteomes/OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt data/
python scripts/04_phylogenomics.py
```

## Method notes

- **Orthogroups** are inferred by OrthoFinder (DIAMOND all-vs-all → normalised
  bit-score graph → MCL clustering).
- The **pangenome matrix** treats every OrthoFinder-unassigned protein as a
  species-specific singleton, so every gene is in exactly one orthogroup and
  the pangenome is complete.
- **Categories**: core = all N species; soft-core = N−1; accessory
  (dispensable/shell) = 2…N−2; species-specific (cloud) = 1.
- **Accumulation curves** average over 500 random species orderings; the
  pangenome curve is fit to Heaps' law `pan(N) = κ·N^γ` (γ>0 ⇒ open pangenome).
- **Gene-content tree**: neighbour-joining on Jaccard distance between species'
  orthogroup presence/absence vectors (repertoire-level signal).
- **Species tree**: OrthoFinder STAG/STRIDE tree from single-copy orthologue
  alignments (sequence-level signal).
