#!/usr/bin/env python3
"""
Phase 3.3: Class I vs Class II oscillatory pathway expression by cell type.

Hypothesis: Immune/communicating cells preferentially express Class II (symmetric)
pathways, while epithelial/autonomous cells preferentially express Class I (asymmetric).

Data: Human Protein Atlas v25, single cell type nCPM (154 cell types).
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# Gene classes from Paper 1
CLASS_I = {  # Autonomous/detection, fast spikes, A < 0.35
    "RELA": "NF-κB",
    "NFKBIA": "NF-κB",
    "TP53": "p53",
    "MDM2": "p53",
    "MAPK1": "ERK",
    "DUSP1": "ERK",
    "PLCG1": "Calcium",
    "ATP2A2": "Calcium",
    "MTOR": "mTOR",
    "TSC1": "mTOR",
    "CDK2": "Cell Cycle",
    "CDKN1A": "Cell Cycle",
}

CLASS_II = {  # Communication/intercellular, symmetric, A ≥ 0.45
    "CTNNB1": "Wnt",
    "AXIN2": "Wnt",
    "NOTCH1": "Notch",
    "FBXW7": "Notch",
    "CLOCK": "Circadian",
    "ARNTL": "Circadian",
    "PER2": "Circadian",
    "JAK1": "JAK-STAT",
    "STAT3": "JAK-STAT",
    "SOCS3": "JAK-STAT",
}

CLASS_III = {  # Shuttling, intermediate
    "YAP1": "Hippo",
    "LATS1": "Hippo",
    "SMAD3": "TGF-β",
    "SMAD7": "TGF-β",
    "NFE2L2": "NRF2",
    "KEAP1": "NRF2",
    "NFATC1": "NFAT",
    "RCAN1": "NFAT",
}

# Cell type categories for our hypothesis
IMMUNE_KEYWORDS = [
    "blood and immune cells",
]
EPITHELIAL_KEYWORDS = [
    "specialized epithelial cells",
    "glandular epithelial cells",
]


def load_expression():
    """Load HPA single cell type expression (long format)."""
    print("Loading HPA expression data...")
    expr = defaultdict(dict)  # gene_name -> {cell_type: nCPM}
    all_genes = set()
    all_celltypes = set()

    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row["Gene name"]
            ct = row["Cell type"]
            ncpm = float(row["nCPM"])
            expr[gene][ct] = ncpm
            all_genes.add(gene)
            all_celltypes.add(ct)

    print(f"  {len(all_genes)} genes, {len(all_celltypes)} cell types")
    return expr, all_celltypes


def load_cell_type_classes():
    """Load HPA cell type classification."""
    ct_class = {}
    ct_group = {}
    with open(DATA_DIR / "rna_single_cell_type_cell_types.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ct_class[row["Cell type"]] = row["Cell type class"]
            ct_group[row["Cell type"]] = row["Cell type group"]
    return ct_class, ct_group


def compute_class_means(expr, gene_class, cell_types):
    """Compute mean nCPM of a gene class per cell type."""
    means = {}
    for ct in cell_types:
        values = []
        for gene in gene_class:
            if gene in expr and ct in expr[gene]:
                values.append(expr[gene][ct])
        if values:
            means[ct] = np.mean(values)
    return means


def main():
    expr, all_celltypes = load_expression()
    ct_class, ct_group = load_cell_type_classes()

    # Check gene coverage
    print(f"\nGene coverage:")
    for label, genes in [("Class I", CLASS_I), ("Class II", CLASS_II), ("Class III", CLASS_III)]:
        found = [g for g in genes if g in expr]
        missing = [g for g in genes if g not in expr]
        print(f"  {label}: {len(found)}/{len(genes)} found", end="")
        if missing:
            print(f" (missing: {', '.join(missing)})", end="")
        print()

    # Compute per-cell-type means
    mean_I = compute_class_means(expr, CLASS_I, all_celltypes)
    mean_II = compute_class_means(expr, CLASS_II, all_celltypes)
    mean_III = compute_class_means(expr, CLASS_III, all_celltypes)

    # Compute ratio: Class II / Class I
    ratios = {}
    for ct in all_celltypes:
        if ct in mean_I and ct in mean_II and mean_I[ct] > 0:
            ratios[ct] = mean_II[ct] / mean_I[ct]

    # Classify cell types
    immune_cts = []
    epithelial_cts = []
    other_cts = []

    for ct in ratios:
        cls = ct_class.get(ct, "unknown")
        if cls in IMMUNE_KEYWORDS:
            immune_cts.append(ct)
        elif cls in EPITHELIAL_KEYWORDS:
            epithelial_cts.append(ct)
        else:
            other_cts.append(ct)

    print(f"\nCell type classification:")
    print(f"  Immune: {len(immune_cts)} cell types")
    print(f"  Epithelial: {len(epithelial_cts)} cell types")
    print(f"  Other: {len(other_cts)} cell types")

    # ===================================================================
    # MAIN TEST: Mann-Whitney U
    # ===================================================================
    print(f"\n{'='*70}")
    print("MAIN TEST: Class II/Class I ratio — Immune vs Epithelial")
    print("=" * 70)

    immune_ratios = [ratios[ct] for ct in immune_cts]
    epithelial_ratios = [ratios[ct] for ct in epithelial_cts]

    print(f"\n  Immune cells (n={len(immune_ratios)}):")
    print(f"    Mean ratio: {np.mean(immune_ratios):.3f}")
    print(f"    Median ratio: {np.median(immune_ratios):.3f}")
    for ct in sorted(immune_cts, key=lambda x: -ratios[x]):
        grp = ct_group.get(ct, "")
        print(f"      {ct:35s} ({grp:20s})  ratio = {ratios[ct]:.3f}")

    print(f"\n  Epithelial cells (n={len(epithelial_ratios)}):")
    print(f"    Mean ratio: {np.mean(epithelial_ratios):.3f}")
    print(f"    Median ratio: {np.median(epithelial_ratios):.3f}")
    for ct in sorted(epithelial_cts, key=lambda x: -ratios[x]):
        grp = ct_group.get(ct, "")
        print(f"      {ct:35s} ({grp:20s})  ratio = {ratios[ct]:.3f}")

    # Mann-Whitney U (one-sided: immune > epithelial)
    U, p_two = stats.mannwhitneyu(immune_ratios, epithelial_ratios, alternative='two-sided')
    _, p_greater = stats.mannwhitneyu(immune_ratios, epithelial_ratios, alternative='greater')

    effect_size = U / (len(immune_ratios) * len(epithelial_ratios))  # rank-biserial

    print(f"\n  Mann-Whitney U = {U:.0f}")
    print(f"  p (two-sided) = {p_two:.4f}")
    print(f"  p (one-sided, immune > epithelial) = {p_greater:.4f}")
    print(f"  Rank-biserial r = {2*effect_size - 1:.3f}")

    if p_greater < 0.05:
        print(f"\n  → SIGNIFICANT: Immune cells have higher Class II/Class I ratio")
    elif p_greater < 0.10:
        print(f"\n  → SUGGESTIVE: Trend towards higher ratio in immune cells")
    else:
        print(f"\n  → NOT SIGNIFICANT")

    # ===================================================================
    # Per-gene breakdown
    # ===================================================================
    print(f"\n{'='*70}")
    print("PER-GENE: Immune vs Epithelial mean nCPM")
    print("=" * 70)

    all_target_genes = {**CLASS_I, **CLASS_II, **CLASS_III}
    print(f"\n  {'Gene':12s} {'Class':8s} {'Pathway':12s} {'Immune':>10s} {'Epithelial':>10s} {'Ratio I/E':>10s}")
    print("  " + "-" * 65)

    for gene in sorted(all_target_genes.keys()):
        if gene not in expr:
            continue
        cls_label = "I" if gene in CLASS_I else "II" if gene in CLASS_II else "III"
        pw = all_target_genes[gene]

        imm_vals = [expr[gene].get(ct, 0) for ct in immune_cts]
        epi_vals = [expr[gene].get(ct, 0) for ct in epithelial_cts]

        imm_mean = np.mean(imm_vals)
        epi_mean = np.mean(epi_vals)
        ie_ratio = imm_mean / epi_mean if epi_mean > 0 else float('inf')

        print(f"  {gene:12s} {cls_label:8s} {pw:12s} {imm_mean:10.1f} {epi_mean:10.1f} {ie_ratio:10.2f}")

    # ===================================================================
    # All cell types ranked by ratio
    # ===================================================================
    print(f"\n{'='*70}")
    print("ALL CELL TYPES: Ranked by Class II / Class I ratio")
    print("=" * 70)

    print(f"\n  {'Rank':>4s} {'Cell type':35s} {'Class':25s} {'Ratio':>8s} {'ClassI':>8s} {'ClassII':>8s}")
    print("  " + "-" * 90)
    for rank, (ct, ratio) in enumerate(sorted(ratios.items(), key=lambda x: -x[1]), 1):
        cls = ct_class.get(ct, "unknown")
        marker = ""
        if cls in IMMUNE_KEYWORDS:
            marker = " ◆IMM"
        elif cls in EPITHELIAL_KEYWORDS:
            marker = " ◇EPI"
        print(f"  {rank:4d} {ct:35s} {cls:25s} {ratio:8.3f} {mean_I.get(ct,0):8.1f} {mean_II.get(ct,0):8.1f}{marker}")

    # ===================================================================
    # Bootstrap CI for difference
    # ===================================================================
    print(f"\n{'='*70}")
    print("Bootstrap 95% CI for mean ratio difference (Immune − Epithelial)")
    print("=" * 70)

    n_boot = 10000
    rng = np.random.default_rng(42)
    diffs = []
    for _ in range(n_boot):
        imm_sample = rng.choice(immune_ratios, size=len(immune_ratios), replace=True)
        epi_sample = rng.choice(epithelial_ratios, size=len(epithelial_ratios), replace=True)
        diffs.append(np.mean(imm_sample) - np.mean(epi_sample))

    ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])
    print(f"  Mean difference: {np.mean(immune_ratios) - np.mean(epithelial_ratios):.4f}")
    print(f"  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    if ci_lo > 0:
        print(f"  → CI excludes 0: difference is significant")
    else:
        print(f"  → CI includes 0: difference not significant")

    # ===================================================================
    # SECONDARY: Class III as intermediate?
    # ===================================================================
    print(f"\n{'='*70}")
    print("SECONDARY: Where does Class III fall?")
    print("=" * 70)

    ratio_III_I = {}
    for ct in all_celltypes:
        if ct in mean_I and ct in mean_III and mean_I[ct] > 0:
            ratio_III_I[ct] = mean_III[ct] / mean_I[ct]

    imm_III = [ratio_III_I[ct] for ct in immune_cts if ct in ratio_III_I]
    epi_III = [ratio_III_I[ct] for ct in epithelial_cts if ct in ratio_III_I]

    print(f"  Class III/Class I ratio:")
    print(f"    Immune: mean = {np.mean(imm_III):.3f}, median = {np.median(imm_III):.3f}")
    print(f"    Epithelial: mean = {np.mean(epi_III):.3f}, median = {np.median(epi_III):.3f}")
    U3, p3 = stats.mannwhitneyu(imm_III, epi_III, alternative='greater')
    print(f"    Mann-Whitney p (immune > epi) = {p3:.4f}")

    # Save results
    results = {
        "n_immune": len(immune_cts),
        "n_epithelial": len(epithelial_cts),
        "immune_mean_ratio": float(np.mean(immune_ratios)),
        "epithelial_mean_ratio": float(np.mean(epithelial_ratios)),
        "mannwhitney_U": float(U),
        "p_twosided": float(p_two),
        "p_onesided": float(p_greater),
        "bootstrap_ci": [float(ci_lo), float(ci_hi)],
    }

    import json
    with open(DATA_DIR / "celltype_expression_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved results to {DATA_DIR / 'celltype_expression_results.json'}")


if __name__ == "__main__":
    main()
