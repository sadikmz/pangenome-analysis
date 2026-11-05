#!/usr/bin/env python3
"""
Rename sequence identifiers in a FASTA file and optionally apply the same
renaming to the corresponding GFF3 annotation file.

Naming rules (applied in priority order):
  1. unplaced_<number>              -> unplaced_<number>
  2. Scaffold_<number>[_string]     -> Scaffold_<number>[_string]
  3. scaffold_<number>[_string]     -> scaffold_<number>[_string]
  4. Contig_<number>[_string]       -> Contig_<number>[_string]
  5. contig_<number>[_string]       -> contig_<number>[_string]
  6. [Ss]upercont<number>.<number>  -> [Ss]upercont<number>.<number>
  7. Chromosome_<number>            -> chromosome_<number>  (only if no scaffold/contig in header)
  8. Chromosome <number>            -> chromosome_<number>  (only if no scaffold/contig in header)
  9. chromosome: <Roman>            -> chromosome_<arabic>  (only if no scaffold/contig in header)
 10. [Mm]itochondri(on|a|al)       -> <seqid>_mito
 11. Fallback                       -> original sequence ID (first whitespace-delimited token)

Usage:
    python rename_seq.py input.fasta output.fasta
    python rename_seq.py input.fasta output.fasta --name-map old_new.tsv
    python rename_seq.py input.fasta output.fasta --gff input.gff --gff-out output.gff
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Compiled patterns (order defines priority)
# ---------------------------------------------------------------------------

_PRIORITY_PATTERNS: list[re.Pattern] = [
    re.compile(r'(unplaced_\d+)'),
    re.compile(r'([Ss]caffold_\d+(?:_[^\s,;]+)?)'),
    re.compile(r'([Cc]ontig_\d+(?:_[^\s,;]+)?)'),
    re.compile(r'([Ss]upercont\d+\.\d+)'),
]

_CHROMOSOME_PATTERN        = re.compile(r'[Cc]hromosome[_ ](\d+)')
_ROMAN_CHROMOSOME_PATTERN  = re.compile(r'[Cc]hromosome:\s+([IVXivx]+)')
_SCAFFOLD_CONTIG_RE        = re.compile(r'[Ss]caffold|[Cc]ontig')
_MITO_PATTERN              = re.compile(r'[Mm]itochondri')  # covers mitochondrion/ia/ial

_ROMAN_VALUES = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10}


def _roman_to_int(s: str) -> int:
    """Convert a Roman numeral string to an integer."""
    s = s.upper()
    result = 0
    for i in range(len(s)):
        val = _ROMAN_VALUES[s[i]]
        if i + 1 < len(s) and val < _ROMAN_VALUES[s[i + 1]]:
            result -= val
        else:
            result += val
    return result


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def extract_new_name(header: str) -> str:
    """Derive a clean sequence name from a FASTA header line (without '>').

    Returns the first matching pattern name, or the original sequence ID
    (first token) if nothing matches.
    """
    # Mitochondrial check takes highest priority
    if _MITO_PATTERN.search(header):
        return f"{header.split()[0]}_mito"

    for pat in _PRIORITY_PATTERNS:
        m = pat.search(header)
        if m:
            return m.group(1)

    # Chromosome patterns only when no scaffold/contig keyword is present
    if not _SCAFFOLD_CONTIG_RE.search(header):
        m = _CHROMOSOME_PATTERN.search(header)
        if m:
            return f"chr_{m.group(1)}"

        m = _ROMAN_CHROMOSOME_PATTERN.search(header)
        if m:
            return f"chr_{_roman_to_int(m.group(1))}"

    # Fallback: original accession (first whitespace token)
    return header.split()[0]


def rename_fasta(in_path: Path, out_path: Path) -> dict[str, str]:
    """Read *in_path*, write renamed sequences to *out_path*.

    Returns a dict mapping old_name -> new_name for every sequence.
    """
    name_map: dict[str, str] = {}

    with open(in_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                header = line[1:].rstrip()
                old_name = header
                new_name = extract_new_name(header)
                name_map[old_name] = new_name
                fh_out.write(f">{new_name}\n")
            else:
                fh_out.write(line)

    return name_map


def rename_gff(in_path: Path, out_path: Path, name_map: dict[str, str]) -> int:
    """Rewrite *in_path* GFF replacing column-1 seqids using *name_map*.

    *name_map* maps full FASTA header -> new name (as produced by rename_fasta).
    Only column 1 (seqid) is touched; all other fields and formatting are
    preserved exactly as-is (comments and blank lines are passed through).

    Returns the number of data lines where the seqid was replaced.
    """
    # Build seqid lookup: first whitespace token of old header -> new name
    seqid_map = {old.split()[0]: new for old, new in name_map.items()}

    replaced = 0
    with open(in_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith("#") or not line.strip():
                fh_out.write(line)
                continue
            fields = line.split("\t")
            seqid = fields[0]
            if seqid in seqid_map:
                fields[0] = seqid_map[seqid]
                replaced += 1
            fh_out.write("\t".join(fields))
    return replaced


def write_name_map(name_map: dict[str, str], out_path: Path) -> None:
    """Write old->new name mapping as a two-column TSV."""
    with open(out_path, "w") as fh:
        fh.write("old_name\tnew_name\n")
        for old, new in name_map.items():
            fh.write(f"{old}\t{new}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rename FASTA sequence identifiers using standard naming rules.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input",  type=Path, help="Input FASTA file")
    parser.add_argument("output", type=Path, help="Output FASTA file with renamed headers")
    parser.add_argument(
        "--name-map",
        type=Path,
        default=None,
        metavar="FILE",
        help="TSV file for old->new name mapping (default: <output>.namemap.tsv)",
    )
    parser.add_argument(
        "--gff",
        type=Path,
        default=None,
        metavar="FILE",
        help="Input GFF3 file to rename alongside the FASTA (optional)",
    )
    parser.add_argument(
        "--gff-out",
        type=Path,
        default=None,
        metavar="FILE",
        help="Output GFF3 file with renamed seqids (default: <output stem>.renamed.gff3)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    name_map = rename_fasta(args.input, args.output)

    map_path = args.name_map if args.name_map else args.output.with_suffix(".namemap.tsv")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    write_name_map(name_map, map_path)

    print(f"Renamed {len(name_map)} sequence(s): {args.input} -> {args.output}")
    print(f"Name map written to {map_path}")

    if args.gff:
        if not args.gff.exists():
            sys.exit(f"Error: GFF file not found: {args.gff}")
        gff_out = args.gff_out if args.gff_out else args.output.with_suffix(".renamed.gff3")
        gff_out.parent.mkdir(parents=True, exist_ok=True)
        n_replaced = rename_gff(args.gff, gff_out, name_map)
        print(f"GFF seqids replaced: {n_replaced} lines -> {gff_out}")

    print("\nName mapping:")
    for old, new in name_map.items():
        print(f"  {old:30s} -> {new}")


if __name__ == "__main__":
    main()
