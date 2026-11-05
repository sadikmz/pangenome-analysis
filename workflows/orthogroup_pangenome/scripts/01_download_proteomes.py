#!/usr/bin/env python3
"""
01_download_proteomes.py
========================
Download RefSeq reference proteomes for a panel of *Aspergillus* species from
NCBI Datasets, and extract the protein FASTA (one longest-isoform protein per
protein-coding gene) for each assembly.

These proteomes are the input to OrthoFinder (see 02_run_orthofinder.sh).

Panel (RefSeq reference assemblies, resolved 2026-07):
    A_fumigatus  GCF_000002655.1   Aspergillus fumigatus Af293
    A_flavus     GCF_009017415.1   Aspergillus flavus
    A_niger      GCF_000002855.4   Aspergillus niger
    A_oryzae     GCF_000184455.2   Aspergillus oryzae RIB40
    A_nidulans   GCF_000011425.1   Aspergillus nidulans FGSC A4
    A_terreus    GCF_000149615.1   Aspergillus terreus NIH2624

Requires only the Python standard library (urllib, zipfile).
Usage:
    python 01_download_proteomes.py
Outputs:
    data/proteomes/<species>.faa   one protein FASTA per species
"""
import urllib.request
import os
import zipfile

# ----------------------------------------------------------------------------
# Species panel -> RefSeq assembly accession.
# To change the panel, resolve accessions via the NCBI Datasets taxon endpoint:
#   https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/<name>/dataset_report
#   ?filters.reference_only=true&filters.assembly_source=refseq
# ----------------------------------------------------------------------------
ACCESSIONS = {
    "A_fumigatus": "GCF_000002655.1",
    "A_flavus":    "GCF_009017415.1",
    "A_niger":     "GCF_000002855.4",
    "A_oryzae":    "GCF_000184455.2",
    "A_nidulans":  "GCF_000011425.1",
    "A_terreus":   "GCF_000149615.1",
}

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "proteomes")
DL_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "_downloads")


def download_assembly(accession: str, dest_zip: str) -> None:
    """Fetch the PROT_FASTA package for one assembly from NCBI Datasets."""
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        f"{accession}/download?include_annotation_type=PROT_FASTA"
    )
    urllib.request.urlretrieve(url, dest_zip)


def extract_protein_faa(dest_zip: str, out_faa: str) -> int:
    """Pull the single protein.faa out of the NCBI Datasets zip. Returns #seqs."""
    with zipfile.ZipFile(dest_zip) as zf:
        members = [n for n in zf.namelist() if n.endswith("protein.faa")]
        if len(members) != 1:
            raise RuntimeError(f"expected 1 protein.faa in {dest_zip}, found {members}")
        with zf.open(members[0]) as src, open(out_faa, "wb") as out:
            data = src.read()
            out.write(data)
    return data.count(b">")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(DL_DIR, exist_ok=True)
    for name, acc in ACCESSIONS.items():
        zip_path = os.path.join(DL_DIR, f"{name}.zip")
        faa_path = os.path.join(OUT_DIR, f"{name}.faa")
        download_assembly(acc, zip_path)
        n = extract_protein_faa(zip_path, faa_path)
        print(f"{name:12s} {acc:18s} {n:6d} proteins -> {faa_path}")


if __name__ == "__main__":
    main()
