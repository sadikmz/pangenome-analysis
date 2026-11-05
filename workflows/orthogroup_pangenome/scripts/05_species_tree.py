#!/usr/bin/env python3
"""
05_species_tree.py
==================
Build a sequence-based species tree from single-copy orthologues, independent
of OrthoFinder's (slow) full gene-tree step.

Steps:
  1. From Orthogroups.txt + proteomes, find strict single-copy orthogroups
     (exactly one gene in every species).
  2. Sample N of them, extract protein sequences per orthogroup.
  3. Align each orthogroup with MUSCLE.
  4. Concatenate alignments into a supermatrix (species x sites).
  5. Infer the species tree with FastTree (LG model) and midpoint-root it.

Usage:
    python 05_species_tree.py [n_orthogroups]
Outputs:
    ../results/tables/species_tree.nwk
    ../results/figures/fig6_species_tree.png
    ../data/single_copy_supermatrix.fasta
"""
import os
import sys
import glob
import subprocess
import random
import dendropy
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
TAB = os.path.join(HERE, "..", "results", "tables")
FIG = os.path.join(HERE, "..", "results", "figures")
WORK = os.path.join(DATA, "_sctree")
os.makedirs(WORK, exist_ok=True)

N_OG = int(sys.argv[1]) if len(sys.argv) > 1 else 60
random.seed(42)


def load_proteomes(proteome_dir):
    """Return {species: {acc: seq}} and {acc: species}."""
    seqs, prot2sp = {}, {}
    for faa in sorted(glob.glob(os.path.join(proteome_dir, "*.faa"))):
        sp = os.path.splitext(os.path.basename(faa))[0]
        seqs[sp] = {}
        acc, buf = None, []
        with open(faa) as fh:
            for line in fh:
                if line.startswith(">"):
                    if acc:
                        seqs[sp][acc] = "".join(buf)
                    acc = line[1:].split()[0]
                    prot2sp[acc] = sp
                    buf = []
                else:
                    buf.append(line.strip())
            if acc:
                seqs[sp][acc] = "".join(buf)
    return seqs, prot2sp


def find_single_copy(orthogroups_txt, prot2sp, species):
    """Orthogroups with exactly one gene per species (strict single-copy)."""
    scos = {}
    with open(orthogroups_txt) as fh:
        for line in fh:
            og, _, members = line.partition(":")
            genes = members.split()
            by_sp = {}
            for g in genes:
                sp = prot2sp.get(g)
                if sp:
                    by_sp.setdefault(sp, []).append(g)
            if len(by_sp) == len(species) and all(len(v) == 1 for v in by_sp.values()):
                scos[og] = {sp: by_sp[sp][0] for sp in species}
    return scos


def muscle_align(in_fa, out_fa):
    subprocess.run(["muscle", "-align", in_fa, "-output", out_fa],
                   check=True, capture_output=True)


def read_fasta(path):
    d, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    d[acc] = "".join(buf)
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if acc:
            d[acc] = "".join(buf)
    return d


def main():
    proteome_dir = os.path.join(DATA, "proteomes")
    seqs, prot2sp = load_proteomes(proteome_dir)
    species = sorted(seqs)
    scos = find_single_copy(os.path.join(DATA, "Orthogroups.txt"), prot2sp, species)
    print(f"Found {len(scos)} strict single-copy orthogroups; sampling {N_OG}.")

    chosen = random.sample(list(scos), min(N_OG, len(scos)))
    concat = {sp: [] for sp in species}
    used = 0
    for og in chosen:
        mapping = scos[og]
        raw = os.path.join(WORK, f"{og}.fa")
        aln = os.path.join(WORK, f"{og}.aln")
        with open(raw, "w") as out:
            for sp in species:
                out.write(f">{sp}\n{seqs[sp][mapping[sp]]}\n")
        try:
            muscle_align(raw, aln)
        except subprocess.CalledProcessError:
            continue
        a = read_fasta(aln)
        if set(a) != set(species):
            continue
        L = len(next(iter(a.values())))
        for sp in species:
            concat[sp].append(a[sp])
        used += 1

    # Write supermatrix
    superfa = os.path.join(DATA, "single_copy_supermatrix.fasta")
    with open(superfa, "w") as out:
        for sp in species:
            out.write(f">{sp}\n{''.join(concat[sp])}\n")
    total_len = len("".join(concat[species[0]]))
    print(f"Concatenated {used} orthogroups -> supermatrix {total_len} aa x {len(species)} species")

    # FastTree (LG model for proteins)
    tree_nwk = os.path.join(TAB, "species_tree.nwk")
    with open(tree_nwk, "w") as out:
        subprocess.run(["fasttree", "-lg", superfa], check=True, stdout=out,
                       stderr=subprocess.DEVNULL)

    # Midpoint root + render
    tree = dendropy.Tree.get(path=tree_nwk, schema="newick")
    # midpoint rooting
    pdm = tree.phylogenetic_distance_matrix()
    tree.reroot_at_midpoint(update_bipartitions=False)
    tree.ladderize()
    print("\nSpecies tree (concatenated single-copy orthologues, LG+FastTree):")
    print(tree.as_ascii_plot(plot_metric="length"))

    leaves = list(tree.leaf_node_iter())
    ypos = {l: i for i, l in enumerate(leaves)}

    def y_of(node):
        if node.is_leaf():
            return ypos[node]
        ys = [y_of(c) for c in node.child_nodes()]
        return sum(ys) / len(ys)

    fig, ax = plt.subplots(figsize=(6.5, 3.8))

    def draw(node, x0):
        y = y_of(node)
        for child in node.child_nodes():
            el = child.edge.length or 0.0
            x1 = x0 + el
            yc = y_of(child)
            ax.plot([x0, x0], [y, yc], color="black", lw=1.3)
            ax.plot([x0, x1], [yc, yc], color="black", lw=1.3)
            # bootstrap-like support label from FastTree (node.label)
            if not child.is_leaf() and child.label:
                ax.text(x1, yc + 0.15, str(child.label), fontsize=7, color="#666")
            draw(child, x1)
        if node.is_leaf():
            ax.text(x0 + max(0.001, x0 * 0.02), y, node.taxon.label.replace("_", " "),
                    va="center", fontsize=10, fontstyle="italic")

    draw(tree.seed_node, 0.0)
    ax.set_yticks([])
    ax.set_xlabel("substitutions / site")
    ax.set_title(f"Aspergillus species tree\n({used} single-copy orthologues, "
                 f"{total_len} aa, LG+FastTree, midpoint-rooted)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlim(-0.005, ax.get_xlim()[1] * 1.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_species_tree.png"), dpi=200)
    plt.close(fig)
    print("\nWrote species_tree.nwk and fig6_species_tree.png")


if __name__ == "__main__":
    main()
