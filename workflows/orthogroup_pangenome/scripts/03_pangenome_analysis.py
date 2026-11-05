#!/usr/bin/env python3
"""
03_pangenome_analysis.py
========================
Pan-gene / pangenome and comparative-genomics analysis of the Aspergillus
orthogroups produced by OrthoFinder.

What it computes
----------------
1. A complete orthogroup x species GENE-COUNT matrix, reconstructed from
   OrthoFinder's Orthogroups.txt (multi-gene orthogroups) plus every
   unassigned gene treated as a species-specific singleton orthogroup.
   This makes the matrix a true pangenome: every gene falls in exactly one OG.

2. Pangenome partitioning (presence/absence across N species):
      core        : present in ALL N species
      soft-core   : present in >= N-1 species
      accessory   : present in 2 .. N-1 species (a.k.a. dispensable/shell)
      species-specific : present in exactly 1 species (cloud)

3. Copy-number view: single-copy core vs multi-copy core (gene families that
   have expanded), and per-species gene-family expansion counts.

4. Pangenome / core-genome accumulation curves over all species orderings,
   with a Heaps'-law fit (n_genes = kappa * N^gamma) to test_out whether the
   pangenome is 'open' (gamma > 0) or 'closed'.

5. Per-species contribution table (genes, orthogroups, unique/strain-specific).

Inputs (all resolved relative to this script):
    <data-dir>/Orthogroups.txt        OrthoFinder MCL orthogroups
    <data-dir>/proteomes/*.faa        proteomes (for protein->species mapping)
Outputs:
    <run-dir>/tables/*.csv        all matrices and summary tables
    <run-dir>/figures/*.png       all figures
"""
import argparse
import os
import glob
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="Pangenome partition tables and figures from OrthoFinder orthogroups.")
parser.add_argument("--data-dir", required=True, help="Folder with Orthogroups.txt and proteomes/*.faa")
parser.add_argument("--run-dir", required=True, help="Run folder: output/<workflow>/<dataset>_<YYYYMMDD>")
ARGS = parser.parse_args()

DATA = ARGS.data_dir
TAB = os.path.join(ARGS.run_dir, "tables")
FIG = os.path.join(ARGS.run_dir, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ----------------------------------------------------------------------------
# Step 1. Map every protein accession -> species from the proteome FASTAs.
# ----------------------------------------------------------------------------
def build_protein_species_map(proteome_dir):
    prot2sp, sp_sizes = {}, {}
    for faa in sorted(glob.glob(os.path.join(proteome_dir, "*.faa"))):
        sp = os.path.splitext(os.path.basename(faa))[0]
        n = 0
        with open(faa) as fh:
            for line in fh:
                if line.startswith(">"):
                    acc = line[1:].split()[0]
                    prot2sp[acc] = sp
                    n += 1
        sp_sizes[sp] = n
    return prot2sp, sp_sizes


# ----------------------------------------------------------------------------
# Step 2. Build the full gene-count matrix (orthogroups x species).
#   - multi-gene orthogroups come from Orthogroups.txt
#   - every protein NOT assigned to such an OG becomes its own singleton OG
# ----------------------------------------------------------------------------
def build_gene_count_matrix(orthogroups_txt, prot2sp, species):
    rows, assigned = {}, set()
    with open(orthogroups_txt) as fh:
        for line in fh:
            og, _, members = line.partition(":")
            genes = members.split()
            counts = {sp: 0 for sp in species}
            for g in genes:
                sp = prot2sp.get(g)
                if sp is not None:
                    counts[sp] += 1
                    assigned.add(g)
            rows[og] = counts
    # Unassigned proteins -> one singleton orthogroup each
    singleton_i = 0
    for g, sp in prot2sp.items():
        if g not in assigned:
            counts = {s: 0 for s in species}
            counts[sp] = 1
            rows[f"SINGLETON{singleton_i:06d}"] = counts
            singleton_i += 1
    mat = pd.DataFrame.from_dict(rows, orient="index")[species]
    mat.index.name = "Orthogroup"
    return mat


# ----------------------------------------------------------------------------
# Step 3. Partition orthogroups into pangenome categories.
# ----------------------------------------------------------------------------
def classify(mat):
    N = mat.shape[1]
    presence = (mat > 0)
    nsp = presence.sum(axis=1)                 # species per orthogroup
    total_copies = mat.sum(axis=1)
    cat = pd.Series(index=mat.index, dtype=object)
    cat[nsp == N] = "core"
    cat[nsp == N - 1] = "softcore"
    cat[(nsp >= 2) & (nsp <= N - 2)] = "accessory"
    cat[nsp == 1] = "species_specific"
    # single- vs multi-copy for the core set
    single_copy_core = (nsp == N) & (mat.max(axis=1) == 1) & (mat.min(axis=1) == 1)
    out = pd.DataFrame({
        "n_species": nsp,
        "total_genes": total_copies,
        "category": cat,
        "single_copy_core": single_copy_core,
    })
    return out


# ----------------------------------------------------------------------------
# Step 4. Accumulation curves + Heaps' law fit.
#   For many random species orderings, record cumulative pan (union of OGs seen)
#   and core (OGs present in ALL sampled species) as species are added.
# ----------------------------------------------------------------------------
def accumulation_curves(presence, n_perm=500, seed=0):
    rng = np.random.default_rng(seed)
    species = list(presence.columns)
    N = len(species)
    P = presence.values.astype(int)
    col = {sp: i for i, sp in enumerate(species)}
    pan = np.zeros((n_perm, N), dtype=int)
    core = np.zeros((n_perm, N), dtype=int)
    for p in range(n_perm):
        order = rng.permutation(species)
        idx = [col[s] for s in order]
        seen_any = np.zeros(P.shape[0], dtype=bool)
        seen_all = np.ones(P.shape[0], dtype=bool)
        for k, ci in enumerate(idx):
            present = P[:, ci] > 0
            seen_any |= present
            seen_all &= present
            pan[p, k] = seen_any.sum()
            core[p, k] = seen_all.sum()
    return pan, core, N


def heaps_law(pan_mean):
    # Fit pan(N) = kappa * N^gamma  ->  log pan = log kappa + gamma log N
    N = np.arange(1, len(pan_mean) + 1)
    coef = np.polyfit(np.log(N), np.log(pan_mean), 1)
    gamma, log_kappa = coef[0], coef[1]
    return gamma, np.exp(log_kappa)


def main():
    proteome_dir = os.path.join(DATA, "proteomes")
    prot2sp, sp_sizes = build_protein_species_map(proteome_dir)
    species = sorted(sp_sizes)
    print("Species and proteome sizes:")
    for s in species:
        print(f"  {s:12s} {sp_sizes[s]:6d} proteins")

    mat = build_gene_count_matrix(os.path.join(DATA, "Orthogroups.txt"),
                                  prot2sp, species)
    mat.to_csv(os.path.join(TAB, "orthogroup_gene_counts.csv"))
    presence = (mat > 0).astype(int)
    presence.to_csv(os.path.join(TAB, "orthogroup_presence_absence.csv"))

    cls = classify(mat)
    cls.to_csv(os.path.join(TAB, "orthogroup_classification.csv"))

    # --- Pangenome partition summary ---
    summary = cls["category"].value_counts().reindex(
        ["core", "softcore", "accessory", "species_specific"]).fillna(0).astype(int)
    n_scc = int(cls["single_copy_core"].sum())
    summary_out = summary.to_frame("n_orthogroups")
    summary_out.loc["single_copy_core (subset of core)"] = n_scc
    summary_out.loc["TOTAL pangenome orthogroups"] = len(mat)
    summary_out.to_csv(os.path.join(TAB, "pangenome_summary.csv"))
    print("\nPangenome partition (orthogroups):")
    print(summary_out.to_string())

    # --- Per-species contribution ---
    per_sp = pd.DataFrame(index=species)
    per_sp["proteins"] = [sp_sizes[s] for s in species]
    per_sp["orthogroups_present"] = presence.sum(axis=0)
    # species-specific OGs contributed by each species
    ss = cls["category"] == "species_specific"
    per_sp["species_specific_OGs"] = [
        int(((presence[s] > 0) & ss).sum()) for s in species]
    # core OGs where this species has >1 copy (lineage-specific expansion)
    core_mask = cls["category"] == "core"
    per_sp["expanded_core_families"] = [
        int(((mat[s] > 1) & core_mask).sum()) for s in species]
    per_sp.to_csv(os.path.join(TAB, "per_species_contribution.csv"))
    print("\nPer-species contribution:")
    print(per_sp.to_string())

    # --- Accumulation curves ---
    pan, core, N = accumulation_curves(presence, n_perm=500)
    xs = np.arange(1, N + 1)
    pan_mean, pan_sd = pan.mean(0), pan.std(0)
    core_mean, core_sd = core.mean(0), core.std(0)
    gamma, kappa = heaps_law(pan_mean)
    curve_df = pd.DataFrame({
        "n_genomes": xs,
        "pan_mean": pan_mean, "pan_sd": pan_sd,
        "core_mean": core_mean, "core_sd": core_sd,
    })
    curve_df.to_csv(os.path.join(TAB, "accumulation_curves.csv"), index=False)
    print(f"\nHeaps' law fit: pan(N) = {kappa:.0f} * N^{gamma:.3f}  "
          f"({'OPEN' if gamma > 0.1 else 'CLOSED'} pangenome)")

    # =========================== FIGURES ===========================
    sns.set_style("whitegrid")

    # Fig 1: pangenome partition bar
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    cats = ["core", "softcore", "accessory", "species_specific"]
    vals = [summary[c] for c in cats]
    colors = ["#2166ac", "#67a9cf", "#f4a582", "#b2182b"]
    bars = ax.bar(cats, vals, color=colors)
    ax.bar_label(bars, padding=2)
    ax.set_ylabel("Number of orthogroups")
    ax.set_title(f"Aspergillus pangenome partition ({N} species)")
    ax.set_xticklabels(["Core\n(all %d)" % N, "Soft-core\n(%d)" % (N - 1),
                        "Accessory\n(2–%d)" % (N - 2), "Species-\nspecific (1)"])
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_pangenome_partition.png"), dpi=200)
    plt.close(fig)

    # Fig 2: accumulation curves
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(xs, pan_mean, yerr=pan_sd, marker="o", capsize=3,
                color="#b2182b", label="Pangenome (total)")
    ax.errorbar(xs, core_mean, yerr=core_sd, marker="s", capsize=3,
                color="#2166ac", label="Core genome")
    xf = np.linspace(1, N, 100)
    ax.plot(xf, kappa * xf ** gamma, "--", color="#b2182b", alpha=0.6,
            label=f"Heaps fit γ={gamma:.2f}")
    ax.set_xlabel("Number of genomes sampled")
    ax.set_ylabel("Number of orthogroups")
    ax.set_title("Pangenome & core-genome accumulation")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_accumulation_curves.png"), dpi=200)
    plt.close(fig)

    # Fig 3: presence/absence clustered heatmap (subsample accessory for clarity)
    acc_idx = cls.index[cls["category"].isin(["accessory", "softcore"])]
    show = presence.loc[acc_idx]
    if len(show) > 0:
        cg = sns.clustermap(show.T, cmap=["#f0f0f0", "#2166ac"], figsize=(9, 3.5),
                            col_cluster=True, row_cluster=True, cbar_pos=None,
                            xticklabels=False, linewidths=0)
        cg.ax_heatmap.set_xlabel(f"{len(show)} accessory + soft-core orthogroups")
        cg.ax_heatmap.set_ylabel("")
        cg.fig.suptitle("Accessory-genome presence/absence", y=1.02)
        cg.savefig(os.path.join(FIG, "fig3_accessory_heatmap.png"), dpi=200)
        plt.close(cg.fig)

    # Fig 4: per-species stacked composition (core vs accessory vs specific)
    comp = pd.DataFrame(index=species)
    for c in cats:
        m = cls["category"] == c
        comp[c] = [int(((presence[s] > 0) & m).sum()) for s in species]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    bottom = np.zeros(len(species))
    for c, col in zip(cats, colors):
        ax.bar(species, comp[c], bottom=bottom, label=c, color=col)
        bottom += comp[c].values
    ax.set_ylabel("Orthogroups present")
    ax.set_title("Per-species orthogroup composition")
    ax.legend(title="Category", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_per_species_composition.png"), dpi=200)
    plt.close(fig)
    comp.to_csv(os.path.join(TAB, "per_species_category_composition.csv"))

    print("\nDone. Tables ->", TAB)
    print("Figures ->", FIG)


if __name__ == "__main__":
    main()
