#!/usr/bin/env python3
"""
04_phylogenomics.py
===================
Comparative-genomics phylogenetics for the Aspergillus panel.

Two complementary trees:
  (A) GENE-CONTENT tree: neighbour-joining on a Jaccard distance computed from
      the orthogroup presence/absence matrix. This captures relationships in
      terms of shared gene repertoire (a pangenome-level signal) and needs only
      the presence/absence table produced by 03_pangenome_analysis.py.
  (B) OrthoFinder SPECIES tree: the STAG/STRIDE tree inferred from single-copy
      orthologue alignments (sequence-level signal). Rendered if present.

Also writes a pairwise shared-orthogroup / Jaccard-similarity matrix, a useful
comparative summary in its own right.

Inputs:
    ../results/tables/orthogroup_presence_absence.csv   (from script 03)
    ../data/SpeciesTree_rooted.txt   (optional; copy from OrthoFinder output)
Outputs:
    ../results/tables/pairwise_jaccard.csv
    ../results/tables/gene_content_tree.nwk
    ../results/figures/fig5_gene_content_tree.png
    ../results/figures/fig6_species_tree.png            (if OrthoFinder tree present)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import dendropy

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "..", "results", "tables")
FIG = os.path.join(HERE, "..", "results", "figures")
DATA = os.path.join(HERE, "..", "data")


def gene_content_nj(presence):
    """Neighbour-joining tree from Jaccard distances between species' gene sets."""
    import dendropy
    from dendropy.calculate import treecompare  # noqa
    species = list(presence.columns)
    X = presence.T.values  # species x orthogroups (0/1)
    # Jaccard distance between each pair of species
    dist = squareform(pdist(X, metric="jaccard"))
    dm_df = pd.DataFrame(dist, index=species, columns=species)
    dm_df.to_csv(os.path.join(TAB, "pairwise_jaccard.csv"))

    # Build a dendropy distance matrix and run neighbour-joining
    pdm = dendropy.PhylogeneticDistanceMatrix()
    taxa = dendropy.TaxonNamespace()
    taxon_objs = {s: taxa.new_taxon(label=s) for s in species}
    mapping = {}
    for i, si in enumerate(species):
        mapping[taxon_objs[si]] = {}
        for j, sj in enumerate(species):
            if i != j:
                mapping[taxon_objs[si]][taxon_objs[sj]] = float(dist[i, j])
    pdm.compile_from_dict(mapping, taxa)
    tree = pdm.nj_tree()
    return tree, dm_df


def ascii_and_plot(tree, title, out_png):
    """Render a dendropy tree: ASCII to stdout + a simple matplotlib cladogram."""
    print(f"\n{title}\n" + tree.as_ascii_plot(plot_metric="length"))
    # Minimal matplotlib rendering via a recursive layout.
    tree.ladderize()
    leaves = [l for l in tree.leaf_node_iter()]
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
            ax.plot([x0, x0], [y, yc], color="black", lw=1.2)      # vertical
            ax.plot([x0, x1], [yc, yc], color="black", lw=1.2)     # horizontal
            draw(child, x1)
        if node.is_leaf():
            ax.text(x0 + 0.002, y, node.taxon.label.replace("_", " "),
                    va="center", fontsize=10, fontstyle="italic")

    seed = tree.seed_node
    draw(seed, 0.0)
    ax.set_yticks([])
    ax.set_xlabel("branch length")
    ax.set_title(title)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlim(-0.005, ax.get_xlim()[1] * 1.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    import dendropy
    presence = pd.read_csv(os.path.join(TAB, "orthogroup_presence_absence.csv"),
                           index_col=0)

    # (A) Gene-content NJ tree
    tree, dm_df = gene_content_nj(presence)
    tree.write(path=os.path.join(TAB, "gene_content_tree.nwk"), schema="newick")
    print("Pairwise Jaccard distance matrix:")
    print(dm_df.round(3).to_string())
    ascii_and_plot(tree, "Gene-content neighbour-joining tree (Jaccard)",
                   os.path.join(FIG, "fig5_gene_content_tree.png"))

    # (B) OrthoFinder species tree, if available
    sp_tree_path = os.path.join(DATA, "SpeciesTree_rooted.txt")
    if os.path.exists(sp_tree_path):
        st = dendropy.Tree.get(path=sp_tree_path, schema="newick")
        st.write(path=os.path.join(TAB, "species_tree.nwk"), schema="newick")
        ascii_and_plot(st, "OrthoFinder species tree (single-copy orthologues)",
                       os.path.join(FIG, "fig6_species_tree.png"))
    else:
        print("\n[info] OrthoFinder species tree not found at", sp_tree_path)
        print("       (copy Species_Tree/SpeciesTree_rooted.txt here when ready)")


if __name__ == "__main__":
    main()
