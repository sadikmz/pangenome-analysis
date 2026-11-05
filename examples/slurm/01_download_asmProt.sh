#!/bin/bash -l
#SBATCH --job-name=download 
#SBATCH --partition=med
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=1          
#SBATCH --cpus-per-task=40   
#SBATCH --output=job.%J.%N.out
#SBATCH --error=job.%J.%N.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=END,FAIL

# =========================
# Environment setup
# =========================

# Load Python module (if required by your cluster)
module load python/3.11.9

# Activate Conda environment
HOME="/path/to/home"
miniconda="$HOME/miniconda3"
source $HOME/miniconda3/etc/profile.d/conda.sh

# conda activate env_name

# =========================
# Input files or params
# =========================   

REPO="/path/to/pangenome-analysis"
DATADIR="/path/to/data"
fgram_dir="$DATADIR/fgram"

# =========================
# Run
# =========================
start_time=$(date +%s)

# 01_download_asmProt.py

python "$REPO/workflows/orthogroup_pangenome/scripts/01_download_asmProt.py" "$REPO/examples/fgram_acc.txt" --prot-dir "$fgram_dir/proteomes" --genome-dir "$fgram_dir/genomes" --dl-dir "$fgram_dir/dl" 2>&1

end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$((runtime % 60))
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

