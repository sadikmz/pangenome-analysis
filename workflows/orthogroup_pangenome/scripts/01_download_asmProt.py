#!/usr/bin/env python3
"""
01_download_asmProt.py
======================
Download genome FASTA, protein FASTA, GFF3 annotation and CDS FASTA
from NCBI Datasets using a FOFN (File Of File Names) of assembly accessions.

Taxon names are resolved automatically from the NCBI Datasets API
so no manual labelling is required in the input file.

FOFN format (one accession per line; lines starting with # are ignored):
    GCF_000002655.1
    GCF_009017415.1
    GCF_000002855.4

Usage:
    python 01_download_asmProt.py accessions.fofn
    python 01_download_asmProt.py accessions.fofn \\
        --prot-dir data/proteomes --genome-dir data/genomes \\
        --gff-dir data/annotations --cds-dir data/cds --dl-dir data/_dl

Outputs:
    <genome-dir>/<accession>.fna      genome FASTA per accession
    <prot-dir>/<accession>.faa        protein FASTA per accession
    <gff-dir>/<accession>.gff         GFF3 annotation per accession
    <cds-dir>/<accession>_cds.fna     CDS FASTA per accession
"""

from __future__ import annotations

import argparse
import json
import os
import urllib.request
import urllib.error
import zipfile
from pathlib import Path


NCBI_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


# ---------------------------------------------------------------------------
# FOFN parsing
# ---------------------------------------------------------------------------

def load_fofn(fofn_path: Path) -> list[str]:
    """Read assembly accessions from a FOFN file.

    Rules:
      - One accession per line.
      - Lines starting with '#' are treated as comments and skipped.
      - Blank lines are skipped.

    Returns a list of accession strings in file order.
    """
    accessions: list[str] = []
    with open(fofn_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                accessions.append(line)
    return accessions


# ---------------------------------------------------------------------------
# NCBI Datasets helpers
# ---------------------------------------------------------------------------

def fetch_taxon_name(accession: str) -> str:
    """Return the organism_name for *accession* via the NCBI Datasets API.

    Endpoint: /genome/accession/{accession}/dataset_report
    JSON path: reports[0].organism.organism_name
    """
    url = f"{NCBI_BASE}/genome/accession/{accession}/dataset_report"
    try:
        with urllib.request.urlopen(url) as resp:
            payload = json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"NCBI API returned HTTP {exc.code} for accession {accession!r}"
        ) from exc

    reports = payload.get("reports", [])
    if not reports:
        raise RuntimeError(
            f"NCBI returned no report for accession {accession!r}. "
            "Check that the accession is valid."
        )
    org_name = reports[0].get("organism", {}).get("organism_name")
    if not org_name:
        raise RuntimeError(
            f"organism_name missing in NCBI response for {accession!r}"
        )
    strain_name = reports[0].get("organism", {}).get("infraspecific_names", {}).get("strain", "")
    strain_primary = strain_name.split(";")[0].strip()  # e.g. "PH-1; NRRL 31084" -> "PH-1"
    if strain_primary and strain_primary not in org_name:
        org_name = f"{org_name} {strain_primary}"
    return org_name


def _safe_filename(name: str) -> str:
    """Convert an organism name to a filesystem-safe string."""
    for ch in (" ", "/", "\\", ":", "*", "?", '"', "<", ">", "|", "(", ")"):
        name = name.replace(ch, "_")
    return name.strip("_")


def download_assembly(accession: str, dest_zip: Path) -> None:
    """Fetch GENOME_FASTA, PROT_FASTA, GENOME_GFF and CDS_FASTA for one assembly."""
    url = (
        f"{NCBI_BASE}/genome/accession/{accession}/download"
        "?include_annotation_type=GENOME_FASTA"
        "&include_annotation_type=PROT_FASTA"
        "&include_annotation_type=GENOME_GFF"
        "&include_annotation_type=CDS_FASTA"
    )
    urllib.request.urlretrieve(url, dest_zip)


def extract_genome_fna(dest_zip: Path, out_fna: Path) -> int | None:
    """Pull the genomic FASTA (.fna) from the NCBI Datasets zip into *out_fna*.

    Returns the number of sequences (contigs/chromosomes), or None if absent.
    """
    with zipfile.ZipFile(dest_zip) as zf:
        members = [
            n for n in zf.namelist()
            if n.endswith(".fna") and "cds_from_genomic" not in n
        ]
        if not members:
            return None
        with open(out_fna, "wb") as out:
            for member in members:
                with zf.open(member) as src:
                    data = src.read()
                    out.write(data)
    return data.count(b">")


def extract_protein_faa(dest_zip: Path, out_faa: Path) -> int | None:
    """Pull protein.faa from the NCBI Datasets zip into *out_faa*.

    Returns the number of sequences extracted, or None if the zip
    contains no protein.faa (e.g. unannotated assembly).
    """
    with zipfile.ZipFile(dest_zip) as zf:
        members = [n for n in zf.namelist() if n.endswith("protein.faa")]
        if not members:
            return None
        if len(members) > 1:
            raise RuntimeError(
                f"Expected 1 protein.faa in {dest_zip}, found {members}"
            )
        with zf.open(members[0]) as src, open(out_faa, "wb") as out:
            data = src.read()
            out.write(data)
    return data.count(b">")


def extract_gff(dest_zip: Path, out_gff: Path) -> int | None:
    """Pull the GFF3 annotation file from the NCBI Datasets zip into *out_gff*.

    Returns the number of feature lines (non-comment), or None if absent.
    """
    with zipfile.ZipFile(dest_zip) as zf:
        members = [
            n for n in zf.namelist()
            if n.endswith(".gff") or n.endswith(".gff3")
        ]
        if not members:
            return None
        with open(out_gff, "wb") as out:
            for member in members:
                with zf.open(member) as src:
                    data = src.read()
                    out.write(data)
    return sum(1 for ln in data.split(b"\n") if ln and not ln.startswith(b"#"))


def extract_cds_fna(dest_zip: Path, out_cds: Path) -> int | None:
    """Pull cds_from_genomic.fna from the NCBI Datasets zip into *out_cds*.

    Returns the number of CDS sequences, or None if absent.
    """
    with zipfile.ZipFile(dest_zip) as zf:
        members = [
            n for n in zf.namelist()
            if "cds_from_genomic" in n and n.endswith(".fna")
        ]
        if not members:
            return None
        with open(out_cds, "wb") as out:
            for member in members:
                with zf.open(member) as src:
                    data = src.read()
                    out.write(data)
    return data.count(b">")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download genome, protein, GFF and CDS files from NCBI Datasets\n"
            "for a FOFN of assembly accessions.\n"
            "Taxon names are resolved automatically from the NCBI Datasets API."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "fofn",
        type=Path,
        help="FOFN file: one NCBI assembly accession per line (# lines are comments)",
    )
    parser.add_argument(
        "--genome-dir",
        type=Path,
        default=Path("data/genomes"),
        metavar="DIR",
        help="Output directory for .fna files (default: data/genomes)",
    )
    parser.add_argument(
        "--prot-dir",
        type=Path,
        default=Path("data/proteomes"),
        metavar="DIR",
        help="Output directory for .faa files (default: data/proteomes)",
    )
    parser.add_argument(
        "--gff-dir",
        type=Path,
        default=Path("data/annotations"),
        metavar="DIR",
        help="Output directory for .gff files (default: data/annotations)",
    )
    parser.add_argument(
        "--cds-dir",
        type=Path,
        default=Path("data/cds"),
        metavar="DIR",
        help="Output directory for _cds.fna files (default: data/cds)",
    )
    parser.add_argument(
        "--dl-dir",
        type=Path,
        default=Path("data/_downloads"),
        metavar="DIR",
        help="Directory for intermediate .zip files (default: data/_downloads)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    genome_dir: Path = args.genome_dir
    prot_dir: Path = args.prot_dir
    gff_dir: Path = args.gff_dir
    cds_dir: Path = args.cds_dir
    dl_dir: Path = args.dl_dir

    for d in (genome_dir, prot_dir, gff_dir, cds_dir, dl_dir):
        d.mkdir(parents=True, exist_ok=True)

    accessions = load_fofn(args.fofn)
    if not accessions:
        print("No accessions found in FOFN file. Nothing to do.")
        return

    print(f"Processing {len(accessions)} accession(s) from {args.fofn}\n")

    taxon_map: dict[str, list[str]] = {}  # taxon_name (raw NCBI) -> [accession, ...]

    for acc in accessions:
        print(f"[{acc}] Resolving taxon name ... ", end="", flush=True)
        taxon_name = fetch_taxon_name(acc)
        taxon_map.setdefault(taxon_name, []).append(acc)
        print(taxon_name)

        zip_path = dl_dir / f"{acc}.zip"
        fna_path = genome_dir / f"{acc}.fna"
        faa_path = prot_dir / f"{acc}.faa"
        gff_path = gff_dir / f"{acc}.gff"
        cds_path = cds_dir / f"{acc}_cds.fna"

        print(f"[{acc}] Downloading ... ", end="", flush=True)
        download_assembly(acc, zip_path)
        print("done")

        n_genome = extract_genome_fna(zip_path, fna_path)
        if n_genome is None:
            print(f"[{acc}]   genome : SKIP (no .fna in zip)")
        else:
            print(f"[{acc}]   genome : {n_genome:6d} sequences -> {fna_path}")

        n_prot = extract_protein_faa(zip_path, faa_path)
        if n_prot is None:
            print(f"[{acc}]   protein: SKIP (no protein.faa in zip)")
        else:
            print(f"[{acc}]   protein: {n_prot:6d} sequences -> {faa_path}")

        n_gff = extract_gff(zip_path, gff_path)
        if n_gff is None:
            print(f"[{acc}]   gff    : SKIP (no .gff in zip)")
        else:
            print(f"[{acc}]   gff    : {n_gff:6d} features  -> {gff_path}")

        n_cds = extract_cds_fna(zip_path, cds_path)
        if n_cds is None:
            print(f"[{acc}]   cds    : SKIP (no cds_from_genomic.fna in zip)")
        else:
            print(f"[{acc}]   cds    : {n_cds:6d} sequences -> {cds_path}")

        print()

    print("Taxon name mapping:")
    for taxon, accs in taxon_map.items():
        print(f"  {taxon} -> {accs}")


if __name__ == "__main__":
    main()
