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
workflows/orthogroup_pangenome/
├── README.md                     this file
├── scripts/
│   ├── 01_download_asmProt.py    download genomes, proteins, GFF3 and CDS for an accession list
│   ├── 01_download_proteomes.py  fetch the fixed Aspergillus proteome panel
│   ├── 02_run_orthofinder.sh     OrthoFinder orthogroup inference
│   ├── 03_pangenome_analysis.py  pangenome partitioning + accumulation curves
│   ├── 04_phylogenomics.py       gene-content NJ tree + species tree render
│   └── 05_species_tree.py        species tree from sampled single-copy orthogroups
└── data/                         inputs, not in Git: proteomes/*.faa, Orthogroups.txt

output/orthogroup_pangenome/<dataset>_<YYYYMMDD>/    run folder, not in Git
```

Each run gets its own folder, holding `RUN_INFO.txt`, `tables/`, `trees/`, `figures/`, `reports/` and `logs/`.

## Reproduce

Run from the repository root.

```bash
conda env create -f environment.yml
conda activate PangenomeAnalysis

W=workflows/orthogroup_pangenome
DATA=$W/data
RUN=output/orthogroup_pangenome/aspergillus_6species_$(date -u +%Y%m%d)

python $W/scripts/01_download_proteomes.py
bash   $W/scripts/02_run_orthofinder.sh $DATA/proteomes 8
cp $DATA/proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt $DATA/
python $W/scripts/03_pangenome_analysis.py --data-dir $DATA --run-dir $RUN
# copy the species tree when OrthoFinder finishes its gene-tree step:
cp $DATA/proteomes/OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt $DATA/
python $W/scripts/04_phylogenomics.py --data-dir $DATA --run-dir $RUN
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
