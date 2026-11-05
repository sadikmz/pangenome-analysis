# synteny_syri

Pairwise whole-genome alignment and synteny along a list of genomes:
minimap2, then samtools, SyRI and plotsr.

## Environment

```bash
conda env create -f envs/synteny.yml
conda activate PangenomeSynteny
```

## Run

Run from the run's `logs/` folder, because the scripts write to the current folder.
`minimap2_syri.sh` also reads `genomes.txt` and `base.cfg` from there, for plotsr.

```bash
RUN=output/synteny_syri/<dataset>_<YYYYMMDD>
mkdir -p "$RUN/logs" && cd "$RUN/logs"

# Optional: give sequences consistent names first
python <repo>/workflows/synteny_syri/rename_seq.py input.fasta output.fasta --name-map old_new.tsv

# Line i of the list is the reference for line i+1
bash <repo>/workflows/synteny_syri/minimap2_syri.sh [--no-filter] <genome_fofn> <threads> <output_prefix>
```
