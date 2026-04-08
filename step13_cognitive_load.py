#!/usr/bin/env python3
"""
Paper 3 exploration: Cellular "cognitive load" — number of active
oscillatory feedback channels per cell type.

14 oscillatory modules from Paper 2. For each cell type (HPA 154):
- Compute mean expression of each module's genes
- Threshold: module "active" if mean nCPM > median across all cell types
- Cognitive load = number of active modules
- Mutual information between module activities

Predictions:
- Stem cells: low cognitive load (quiescent)
- Differentiated: higher (sensing specific niche)
- Immune cells: high (active sensing)
- Neurons: high (complex environment)
"""

import csv
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# 14 oscillatory modules — core genes per pathway
MODULES = {
    "NF-κB":      {"RELA", "NFKBIA", "TNFAIP3", "NFKB1"},
    "p53":        {"TP53", "MDM2", "CDKN1A"},
    "ERK/MAPK":   {"MAPK1", "MAPK3", "DUSP1", "DUSP6"},
    "Wnt":        {"CTNNB1", "AXIN2", "APC", "GSK3B"},
    "Notch":      {"NOTCH1", "FBXW7", "HES1", "DLL1"},
    "Circadian":  {"CLOCK", "PER2", "CRY1", "NR1D1"},
    "Hippo":      {"YAP1", "LATS1", "LATS2", "STK3"},
    "mTOR":       {"MTOR", "TSC1", "TSC2", "RPS6KB1"},
    "JAK-STAT":   {"JAK1", "JAK2", "STAT3", "SOCS3", "SOCS1"},
    "TGF-β":      {"TGFBR1", "SMAD2", "SMAD3", "SMAD7"},
    "NRF2":       {"NFE2L2", "KEAP1", "HMOX1"},
    "Calcium":    {"PLCG1", "ATP2A2", "ITPR1", "CALM1"},
    "Cell Cycle":  {"CDK2", "CDKN1A", "CCNE1", "RB1", "E2F1"},
    "PI3K/PTEN":  {"PIK3CA", "PTEN", "AKT1"},
}

# Waveform asymmetry class from Paper 1
MODULE_CLASS = {
    "NF-κB": "I", "p53": "I", "ERK/MAPK": "I", "Calcium": "I",
    "mTOR": "I", "Cell Cycle": "I",
    "Wnt": "II", "Notch": "II", "Circadian": "II", "JAK-STAT": "II",
    "Hippo": "III", "TGF-β": "III", "NRF2": "III", "PI3K/PTEN": "III",
}

# Module evolutionary age (Mya) from Phase 1.3
MODULE_AGE = {
    "ERK/MAPK": 1500, "Calcium": 1105, "Cell Cycle": 797,
    "NF-κB": 797, "Wnt": 797, "Circadian": 797, "Hippo": 797,
    "mTOR": 1105, "p53": 684, "TGF-β": 684,
    "Notch": 435, "JAK-STAT": 435, "NRF2": 1500, "PI3K/PTEN": 797,
}


def load_expression():
    """Load HPA nCPM per gene per cell type."""
    print("Loading HPA expression data...")
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row["Gene name"]
            ct = row["Cell type"]
            ncpm = float(row["nCPM"])
            expr[gene][ct] = ncpm
    cell_types = sorted(set(ct for gene_data in expr.values() for ct in gene_data))
    print(f"  {len(expr)} genes, {len(cell_types)} cell types")
    return expr, cell_types


def load_cell_classes():
    ct_class = {}
    ct_group = {}
    with open(DATA_DIR / "rna_single_cell_type_cell_types.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ct_class[row["Cell type"]] = row["Cell type class"]
            ct_group[row["Cell type"]] = row["Cell type group"]
    return ct_class, ct_group


def compute_module_activity(expr, cell_types):
    """Compute mean nCPM per module per cell type → activity matrix."""
    module_names = sorted(MODULES.keys())
    n_mod = len(module_names)
    n_ct = len(cell_types)

    # Matrix: modules × cell_types
    activity = np.zeros((n_mod, n_ct))
    gene_coverage = {}

    for i, mod in enumerate(module_names):
        genes = MODULES[mod]
        found = [g for g in genes if g in expr]
        gene_coverage[mod] = (len(found), len(genes))

        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in found]
            activity[i, j] = np.mean(vals) if vals else 0

    return activity, module_names, gene_coverage


def compute_cognitive_load(activity, module_names, cell_types):
    """
    For each cell type: count modules with above-median expression.
    Cognitive load = number of "active" modules.
    """
    n_mod, n_ct = activity.shape

    # Threshold: per module, median across cell types
    thresholds = np.median(activity, axis=1)

    # Binary: active or not
    active = activity > thresholds[:, np.newaxis]

    # Cognitive load per cell type
    cog_load = active.sum(axis=0)

    # Also compute continuous version: sum of z-scored activities
    z_activity = np.zeros_like(activity)
    for i in range(n_mod):
        mu = np.mean(activity[i, :])
        sd = np.std(activity[i, :])
        if sd > 0:
            z_activity[i, :] = (activity[i, :] - mu) / sd

    continuous_load = z_activity.sum(axis=0)

    return cog_load, active, continuous_load, thresholds


def compute_module_correlation(activity, module_names):
    """Correlation matrix between modules across cell types."""
    n_mod = len(module_names)
    corr = np.corrcoef(activity)
    return corr


def mutual_information_approx(x, y, n_bins=10):
    """Approximate MI using binned histogram."""
    c_xy, _, _ = np.histogram2d(x, y, bins=n_bins)
    c_xy = c_xy / c_xy.sum()
    c_x = c_xy.sum(axis=1)
    c_y = c_xy.sum(axis=0)

    mi = 0
    for i in range(n_bins):
        for j in range(n_bins):
            if c_xy[i, j] > 0 and c_x[i] > 0 and c_y[j] > 0:
                mi += c_xy[i, j] * np.log2(c_xy[i, j] / (c_x[i] * c_y[j]))
    return mi


def main():
    expr, cell_types = load_expression()
    ct_class, ct_group = load_cell_classes()

    # Module activity
    activity, module_names, gene_coverage = compute_module_activity(expr, cell_types)

    print(f"\nModule gene coverage:")
    for mod in module_names:
        found, total = gene_coverage[mod]
        print(f"  {mod:15s}: {found}/{total} genes found")

    # Cognitive load
    cog_load, active_binary, continuous_load, thresholds = compute_cognitive_load(
        activity, module_names, cell_types)

    # ===================================================================
    # RESULTS
    # ===================================================================
    print(f"\n{'='*80}")
    print("COGNITIVE LOAD: Number of active oscillatory modules per cell type")
    print("=" * 80)

    # Rank by cognitive load
    ranked = sorted(range(len(cell_types)),
                    key=lambda j: (-cog_load[j], -continuous_load[j]))

    print(f"\n{'Rank':>4} {'Cell type':35s} {'Class':25s} {'Load':>5} {'Cont':>7} {'Active modules'}")
    print("-" * 120)

    for rank, j in enumerate(ranked, 1):
        ct = cell_types[j]
        cls = ct_class.get(ct, "unknown")
        active_mods = [module_names[i] for i in range(len(module_names)) if active_binary[i, j]]
        mod_str = ", ".join(active_mods[:6])
        if len(active_mods) > 6:
            mod_str += f" +{len(active_mods)-6}"
        print(f"{rank:4d} {ct:35s} {cls:25s} {cog_load[j]:5.0f} {continuous_load[j]:7.1f} {mod_str}")

    # ===================================================================
    # CELL CLASS COMPARISON
    # ===================================================================
    print(f"\n{'='*80}")
    print("COGNITIVE LOAD BY CELL CLASS")
    print("=" * 80)

    class_loads = defaultdict(list)
    for j, ct in enumerate(cell_types):
        cls = ct_class.get(ct, "unknown")
        class_loads[cls].append(cog_load[j])

    print(f"\n{'Cell class':30s} {'N':>4} {'Mean':>6} {'Median':>7} {'SD':>6}")
    print("-" * 60)
    for cls in sorted(class_loads.keys(), key=lambda c: -np.mean(class_loads[c])):
        vals = class_loads[cls]
        print(f"{cls:30s} {len(vals):4d} {np.mean(vals):6.1f} {np.median(vals):7.1f} {np.std(vals):6.1f}")

    # Statistical tests
    immune = [cog_load[j] for j, ct in enumerate(cell_types)
              if ct_class.get(ct) == "blood and immune cells"]
    epithelial = [cog_load[j] for j, ct in enumerate(cell_types)
                  if ct_class.get(ct) in ("specialized epithelial cells", "glandular epithelial cells")]
    neuronal = [cog_load[j] for j, ct in enumerate(cell_types)
                if ct_class.get(ct) == "neuronal cells"]
    stem = [cog_load[j] for j, ct in enumerate(cell_types)
            if ct_class.get(ct) == "stem and proliferating cells"]
    muscle = [cog_load[j] for j, ct in enumerate(cell_types)
              if ct_class.get(ct) == "muscle cells"]

    print(f"\n--- Pairwise comparisons (Mann-Whitney) ---")
    comparisons = [
        ("Immune", immune, "Epithelial", epithelial),
        ("Neuronal", neuronal, "Epithelial", epithelial),
        ("Immune", immune, "Stem", stem),
        ("Immune", immune, "Muscle", muscle),
        ("Neuronal", neuronal, "Muscle", muscle),
    ]
    for name1, vals1, name2, vals2 in comparisons:
        if vals1 and vals2:
            U, p = stats.mannwhitneyu(vals1, vals2, alternative='two-sided')
            d = np.mean(vals1) - np.mean(vals2)
            print(f"  {name1:10s} vs {name2:10s}: Δ = {d:+.1f}, U = {U:.0f}, p = {p:.4f}")

    # Kruskal-Wallis across all classes
    all_groups = [v for v in class_loads.values() if len(v) >= 3]
    if len(all_groups) >= 3:
        H, p_kw = stats.kruskal(*all_groups)
        print(f"\n  Kruskal-Wallis across all classes: H = {H:.1f}, p = {p_kw:.4f}")

    # ===================================================================
    # MODULE CO-ACTIVATION PATTERNS
    # ===================================================================
    print(f"\n{'='*80}")
    print("MODULE CORRELATION MATRIX (Pearson r across 154 cell types)")
    print("=" * 80)

    corr = compute_module_correlation(activity, module_names)

    # Print top positive and negative correlations
    pairs = []
    for i in range(len(module_names)):
        for j in range(i+1, len(module_names)):
            pairs.append((module_names[i], module_names[j], corr[i, j]))

    pairs.sort(key=lambda x: -x[2])
    print(f"\nTop 10 positive correlations (co-active modules):")
    for m1, m2, r in pairs[:10]:
        c1, c2 = MODULE_CLASS[m1], MODULE_CLASS[m2]
        print(f"  {m1:15s} ({c1}) × {m2:15s} ({c2}): r = {r:+.3f}")

    print(f"\nTop 10 negative correlations (anti-correlated modules):")
    for m1, m2, r in pairs[-10:]:
        c1, c2 = MODULE_CLASS[m1], MODULE_CLASS[m2]
        print(f"  {m1:15s} ({c1}) × {m2:15s} ({c2}): r = {r:+.3f}")

    # Mean within-class vs between-class correlation
    within = []
    between = []
    for i in range(len(module_names)):
        for j in range(i+1, len(module_names)):
            ci = MODULE_CLASS[module_names[i]]
            cj = MODULE_CLASS[module_names[j]]
            if ci == cj:
                within.append(corr[i, j])
            else:
                between.append(corr[i, j])

    print(f"\n  Within-class mean r: {np.mean(within):.3f} (n={len(within)})")
    print(f"  Between-class mean r: {np.mean(between):.3f} (n={len(between)})")
    if within and between:
        U_c, p_c = stats.mannwhitneyu(within, between, alternative='greater')
        print(f"  Mann-Whitney (within > between): p = {p_c:.4f}")

    # ===================================================================
    # EVOLUTIONARY AGE vs EXPRESSION BREADTH
    # ===================================================================
    print(f"\n{'='*80}")
    print("MODULE AGE vs EXPRESSION BREADTH")
    print("=" * 80)

    breadths = []
    ages = []
    for i, mod in enumerate(module_names):
        # Breadth = fraction of cell types where module is active
        breadth = active_binary[i, :].sum() / len(cell_types)
        breadths.append(breadth)
        ages.append(MODULE_AGE.get(mod, 0))
        print(f"  {mod:15s}: age = {MODULE_AGE.get(mod, 0):5d} Mya, "
              f"breadth = {breadth:.2f} ({active_binary[i,:].sum():.0f}/{len(cell_types)} cell types)")

    rho_ab, p_ab = stats.spearmanr(ages, breadths)
    print(f"\n  Spearman ρ(age, breadth) = {rho_ab:.3f}, p = {p_ab:.4f}")
    if p_ab < 0.05:
        print(f"  → SIGNIFICANT: older modules expressed in more cell types")

    # ===================================================================
    # WHICH MODULES DISTINGUISH CELL CLASSES?
    # ===================================================================
    print(f"\n{'='*80}")
    print("DISCRIMINATIVE MODULES: Which modules best separate cell classes?")
    print("=" * 80)

    immune_idx = [j for j, ct in enumerate(cell_types)
                  if ct_class.get(ct) == "blood and immune cells"]
    epi_idx = [j for j, ct in enumerate(cell_types)
               if ct_class.get(ct) in ("specialized epithelial cells", "glandular epithelial cells")]

    print(f"\n  {'Module':15s} {'Class':>5s} {'Immune':>8s} {'Epithelial':>10s} {'Ratio':>7s} {'p':>8s}")
    print("  " + "-" * 60)

    for i, mod in enumerate(module_names):
        imm_act = activity[i, immune_idx]
        epi_act = activity[i, epi_idx]
        ratio = np.mean(imm_act) / np.mean(epi_act) if np.mean(epi_act) > 0 else float('inf')
        _, p_mod = stats.mannwhitneyu(imm_act, epi_act, alternative='two-sided')
        cls = MODULE_CLASS[mod]
        sig = "*" if p_mod < 0.05 else " "
        print(f"  {mod:15s} {cls:>5s} {np.mean(imm_act):8.1f} {np.mean(epi_act):10.1f} "
              f"{ratio:7.2f} {p_mod:8.4f} {sig}")

    # Save
    save_results(cell_types, cog_load, continuous_load, active_binary, module_names,
                 ct_class, corr, ages, breadths)
    print("\nDone!")


def save_results(cell_types, cog_load, continuous_load, active_binary, module_names,
                 ct_class, corr, ages, breadths):
    with open(DATA_DIR / "cognitive_load.csv", "w", newline="") as f:
        w = csv.writer(f)
        header = ["cell_type", "cell_class", "cognitive_load", "continuous_load"] + module_names
        w.writerow(header)
        for j, ct in enumerate(cell_types):
            row = [ct, ct_class.get(ct, ""), int(cog_load[j]), f"{continuous_load[j]:.2f}"]
            row += [int(active_binary[i, j]) for i in range(len(module_names))]
            w.writerow(row)

    # Correlation matrix
    corr_dict = {}
    for i in range(len(module_names)):
        for j in range(len(module_names)):
            corr_dict[f"{module_names[i]}|{module_names[j]}"] = round(float(corr[i, j]), 3)

    summary = {
        "n_modules": len(module_names),
        "n_cell_types": len(cell_types),
        "age_breadth_rho": float(stats.spearmanr(ages, breadths)[0]),
        "age_breadth_p": float(stats.spearmanr(ages, breadths)[1]),
    }
    with open(DATA_DIR / "cognitive_load_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: cognitive_load.csv, cognitive_load_summary.json")


if __name__ == "__main__":
    main()
