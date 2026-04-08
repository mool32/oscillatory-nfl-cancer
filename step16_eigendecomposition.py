#!/usr/bin/env python3
"""
Phase B.1/B.2: Module interaction matrix and eigendecomposition.

Build 14×14 expression co-variation matrix across cell types.
Eigendecompose to find effective dimensionality of cellular perception.
Compare with random gene sets (null model).

Pre-registered predictions:
  v1: general signaling intensity (all positive)
  v2: danger (NF-κB,JAK-STAT,ERK) vs tissue (Hippo,Wnt,Circadian)
  v3: metabolic (mTOR) vs inflammatory
  3 eigenvalues explain >80% variance → perception is 3-dimensional
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json
import random

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")
random.seed(42)
np.random.seed(42)

MODULES_CORE = {
    "NF-κB": {"RELA", "NFKBIA", "TNFAIP3", "NFKB1"},
    "p53": {"TP53", "MDM2", "CDKN1A"},
    "ERK/MAPK": {"MAPK1", "MAPK3", "DUSP1", "DUSP6"},
    "Wnt": {"CTNNB1", "AXIN2", "APC", "GSK3B"},
    "Notch": {"NOTCH1", "FBXW7", "HES1", "DLL1"},
    "Circadian": {"CLOCK", "PER2", "CRY1", "NR1D1"},
    "Hippo": {"YAP1", "LATS1", "LATS2", "STK3"},
    "mTOR": {"MTOR", "TSC1", "TSC2", "RPS6KB1"},
    "JAK-STAT": {"JAK1", "JAK2", "STAT3", "SOCS3", "SOCS1"},
    "TGF-β": {"TGFBR1", "SMAD2", "SMAD3", "SMAD7"},
    "NRF2": {"NFE2L2", "KEAP1", "HMOX1"},
    "Calcium": {"PLCG1", "ATP2A2", "ITPR1", "CALM1"},
    "Cell Cycle": {"CDK2", "CDKN1A", "CCNE1", "RB1", "E2F1"},
    "PI3K/PTEN": {"PIK3CA", "PTEN", "AKT1"},
}

MODULE_ORDER = ["NF-κB", "ERK/MAPK", "JAK-STAT", "p53", "Wnt", "Notch",
                "Hippo", "TGF-β", "mTOR", "Calcium", "Cell Cycle",
                "Circadian", "NRF2", "PI3K/PTEN"]


def load_expression():
    print("Loading HPA expression data...")
    expr = defaultdict(dict)
    genes_set = set()
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row["Gene name"]
            expr[gene][row["Cell type"]] = float(row["nCPM"])
            genes_set.add(gene)
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types, sorted(genes_set)


def module_activity_matrix(expr, cell_types):
    """Build modules × cell_types matrix."""
    N = len(MODULE_ORDER)
    M = len(cell_types)
    mat = np.zeros((N, M))
    for i, mod in enumerate(MODULE_ORDER):
        genes = [g for g in MODULES_CORE[mod] if g in expr]
        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in genes]
            mat[i, j] = np.mean(vals) if vals else 0
    return mat


def random_gene_set_matrix(expr, cell_types, all_genes, n_sets=14, set_size=4, seed=42):
    """Build random gene sets × cell_types matrix."""
    rng = np.random.RandomState(seed)
    available = [g for g in all_genes if len(expr[g]) > 50]  # expressed in many cell types
    N = n_sets
    M = len(cell_types)
    mat = np.zeros((N, M))
    for i in range(N):
        genes = rng.choice(available, size=set_size, replace=False)
        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in genes]
            mat[i, j] = np.mean(vals) if vals else 0
    return mat


def eigendecompose(mat, labels, name=""):
    """Compute correlation matrix and eigendecompose."""
    # Correlation matrix (modules × modules, computed across cell types)
    corr = np.corrcoef(mat)  # N×N correlation matrix

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(corr)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Variance explained
    total = np.sum(np.maximum(eigenvalues, 0))
    var_explained = np.maximum(eigenvalues, 0) / total * 100
    cumvar = np.cumsum(var_explained)

    return corr, eigenvalues, eigenvectors, var_explained, cumvar


def main():
    expr, cell_types, all_genes = load_expression()

    # ===================================================================
    # PART 1: Module correlation matrix
    # ===================================================================
    print("=" * 80)
    print("PART 1: Module co-variation matrix (14 modules × 154 cell types)")
    print("=" * 80)

    mat = module_activity_matrix(expr, cell_types)
    corr, eigenvalues, eigenvectors, var_explained, cumvar = eigendecompose(mat, MODULE_ORDER)

    # Print correlation matrix
    print(f"\nCorrelation matrix (Pearson, across {len(cell_types)} cell types):")
    header = "            " + "".join(f"{m[:7]:>8s}" for m in MODULE_ORDER)
    print(header)
    for i, mod in enumerate(MODULE_ORDER):
        row = f"  {mod:10s}" + "".join(f"{corr[i,j]:+8.2f}" for j in range(len(MODULE_ORDER)))
        print(row)

    # ===================================================================
    # PART 2: Eigendecomposition
    # ===================================================================
    print(f"\n{'='*80}")
    print("PART 2: Eigendecomposition — effective dimensionality")
    print("=" * 80)

    print(f"\n  {'PC':>3s} {'Eigenvalue':>11s} {'Var %':>7s} {'Cum %':>7s} {'> Kaiser?':>10s}")
    print("  " + "-" * 42)
    for i in range(len(eigenvalues)):
        kaiser = "✓" if eigenvalues[i] > 1.0 else ""
        print(f"  {i+1:3d} {eigenvalues[i]:11.3f} {var_explained[i]:7.1f} {cumvar[i]:7.1f} {kaiser:>10s}")

    n_kaiser = sum(1 for ev in eigenvalues if ev > 1.0)
    n_90pct = next(i+1 for i, cv in enumerate(cumvar) if cv >= 90)
    n_80pct = next(i+1 for i, cv in enumerate(cumvar) if cv >= 80)

    print(f"\n  Kaiser criterion (eigenvalue > 1): {n_kaiser} components")
    print(f"  Components for 80% variance: {n_80pct}")
    print(f"  Components for 90% variance: {n_90pct}")

    # ===================================================================
    # PART 3: Eigenvector interpretation
    # ===================================================================
    print(f"\n{'='*80}")
    print("PART 3: Eigenvector loadings — what does each dimension mean?")
    print("=" * 80)

    for pc in range(min(5, len(eigenvalues))):
        if eigenvalues[pc] < 0.5:
            break
        vec = eigenvectors[:, pc]
        print(f"\n  PC{pc+1} (λ={eigenvalues[pc]:.2f}, {var_explained[pc]:.1f}% variance):")

        # Sort by absolute loading
        sorted_idx = np.argsort(np.abs(vec))[::-1]
        for j in sorted_idx:
            bar = "█" * int(abs(vec[j]) * 20)
            sign = "+" if vec[j] > 0 else "−"
            print(f"    {sign} {MODULE_ORDER[j]:12s} {vec[j]:+.3f}  {bar}")

        # Interpret
        pos = [MODULE_ORDER[j] for j in range(len(vec)) if vec[j] > 0.15]
        neg = [MODULE_ORDER[j] for j in range(len(vec)) if vec[j] < -0.15]
        if pos and neg:
            print(f"    → Axis: [{', '.join(pos)}] vs [{', '.join(neg)}]")
        elif pos:
            print(f"    → All positive: general {'|'.join(pos[:3])} intensity")

    # ===================================================================
    # PART 4: Pre-registered prediction check
    # ===================================================================
    print(f"\n{'='*80}")
    print("PART 4: Pre-registered prediction check")
    print("=" * 80)

    # Prediction 1: v1 all positive
    v1 = eigenvectors[:, 0]
    all_pos_v1 = all(v > -0.05 for v in v1)
    print(f"\n  Prediction 1: PC1 = general signaling intensity (all positive)")
    print(f"    Min loading: {min(v1):.3f}")
    print(f"    Result: {'✓ CONFIRMED' if all_pos_v1 else '✗ NOT confirmed'}")

    # Prediction 2: v2 = danger vs tissue
    v2 = eigenvectors[:, 1]
    danger_mods = {"NF-κB", "JAK-STAT", "ERK/MAPK"}
    tissue_mods = {"Hippo", "Wnt", "Circadian"}
    danger_idx = [MODULE_ORDER.index(m) for m in danger_mods]
    tissue_idx = [MODULE_ORDER.index(m) for m in tissue_mods]

    # Check if danger and tissue have OPPOSITE signs
    danger_signs = [np.sign(v2[i]) for i in danger_idx]
    tissue_signs = [np.sign(v2[i]) for i in tissue_idx]

    # They should have opposite signs (don't know which is + or -)
    danger_mean = np.mean([v2[i] for i in danger_idx])
    tissue_mean = np.mean([v2[i] for i in tissue_idx])
    opposite = (danger_mean * tissue_mean) < 0

    print(f"\n  Prediction 2: PC2 = danger vs tissue maintenance")
    print(f"    Danger modules (NF-κB, JAK-STAT, ERK): mean loading = {danger_mean:.3f}")
    print(f"    Tissue modules (Hippo, Wnt, Circadian): mean loading = {tissue_mean:.3f}")
    print(f"    Result: {'✓ CONFIRMED — opposite signs' if opposite else '✗ NOT confirmed'}")

    # Prediction 3: mTOR anti-correlates with danger
    v3 = eigenvectors[:, 2] if len(eigenvectors[0]) > 2 else None
    if v3 is not None:
        mtor_idx = MODULE_ORDER.index("mTOR")
        mtor_v3 = v3[mtor_idx]
        danger_v3 = np.mean([v3[i] for i in danger_idx])
        opposite_v3 = (mtor_v3 * danger_v3) < 0
        print(f"\n  Prediction 3: mTOR vs danger in some PC")
        print(f"    mTOR loading on PC3: {mtor_v3:.3f}")
        print(f"    Danger mean on PC3: {danger_v3:.3f}")
        print(f"    Result: {'✓ Opposite' if opposite_v3 else '✗ Same sign'}")

        # Also check PC2 for mTOR
        mtor_v2 = v2[mtor_idx]
        print(f"    mTOR on PC2: {mtor_v2:.3f} (danger on PC2: {danger_mean:.3f})")

    # Prediction 4: 3 PCs explain >80%
    print(f"\n  Prediction 4: 3 PCs explain >80% variance")
    print(f"    PC1-3 cumulative: {cumvar[2]:.1f}%")
    print(f"    Result: {'✓ CONFIRMED' if cumvar[2] > 80 else '✗ NOT confirmed — need ' + str(n_80pct) + ' PCs'}")

    # ===================================================================
    # PART 5: Null model — random gene sets
    # ===================================================================
    print(f"\n{'='*80}")
    print("PART 5: Null model — 14 random gene sets")
    print("=" * 80)

    N_PERM = 100
    rand_kaiser_counts = []
    rand_80pct_counts = []
    rand_top3_var = []

    for perm in range(N_PERM):
        rmat = random_gene_set_matrix(expr, cell_types, all_genes, seed=perm)
        _, rev, _, rvar, rcum = eigendecompose(rmat, [f"rand_{i}" for i in range(14)])
        rand_kaiser_counts.append(sum(1 for ev in rev if ev > 1.0))
        rand_80pct_counts.append(next(i+1 for i, cv in enumerate(rcum) if cv >= 80))
        rand_top3_var.append(rcum[2])

    print(f"\n  Real modules:  Kaiser={n_kaiser}, PCs for 80%={n_80pct}, Top3 var={cumvar[2]:.1f}%")
    print(f"  Random (N={N_PERM}):")
    print(f"    Kaiser:    mean={np.mean(rand_kaiser_counts):.1f} ± {np.std(rand_kaiser_counts):.1f}")
    print(f"    PCs for 80%: mean={np.mean(rand_80pct_counts):.1f} ± {np.std(rand_80pct_counts):.1f}")
    print(f"    Top3 var:  mean={np.mean(rand_top3_var):.1f}% ± {np.std(rand_top3_var):.1f}%")

    # Z-scores
    z_kaiser = (n_kaiser - np.mean(rand_kaiser_counts)) / (np.std(rand_kaiser_counts) + 0.01)
    z_top3 = (cumvar[2] - np.mean(rand_top3_var)) / (np.std(rand_top3_var) + 0.01)

    print(f"\n  Z-score (Kaiser): {z_kaiser:.2f}")
    print(f"  Z-score (Top3 var): {z_top3:.2f}")

    if cumvar[2] > np.mean(rand_top3_var) + 2*np.std(rand_top3_var):
        print(f"  ✓ Real modules MORE structured than random (top3 var > random + 2σ)")
    elif cumvar[2] < np.mean(rand_top3_var) - 2*np.std(rand_top3_var):
        print(f"  Real modules LESS structured than random — unexpected")
    else:
        print(f"  ≈ Real modules similar to random — structure may be trivial")

    if n_kaiser < np.mean(rand_kaiser_counts) - 2*np.std(rand_kaiser_counts):
        print(f"  ✓ Fewer significant PCs than random → more concentrated structure")

    # ===================================================================
    # PART 6: Module clustering in PC space
    # ===================================================================
    print(f"\n{'='*80}")
    print("PART 6: Module positions in PC2 × PC3 space")
    print("=" * 80)

    print(f"\n  {'Module':15s} {'PC1':>7s} {'PC2':>7s} {'PC3':>7s}")
    print("  " + "-" * 40)
    for i, mod in enumerate(MODULE_ORDER):
        print(f"  {mod:15s} {eigenvectors[i,0]:+7.3f} {eigenvectors[i,1]:+7.3f} {eigenvectors[i,2]:+7.3f}")

    # Identify clusters using PC2-PC3
    # Danger cluster: high PC2 (or low, depending on sign)
    # Tissue cluster: opposite PC2
    # Metabolic: extreme PC3

    print(f"\n  Cluster assignment (by PC2/PC3 quadrant):")
    for i, mod in enumerate(MODULE_ORDER):
        pc2, pc3 = eigenvectors[i, 1], eigenvectors[i, 2]
        if pc2 > 0.15 and pc3 > 0:
            q = "Q1 (PC2+, PC3+)"
        elif pc2 > 0.15 and pc3 <= 0:
            q = "Q2 (PC2+, PC3−)"
        elif pc2 <= 0.15 and pc2 >= -0.15:
            q = "Center"
        elif pc2 < -0.15 and pc3 > 0:
            q = "Q3 (PC2−, PC3+)"
        else:
            q = "Q4 (PC2−, PC3−)"
        print(f"    {mod:15s} → {q}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("EIGENDECOMPOSITION SUMMARY")
    print("=" * 80)

    results = {
        "n_modules": 14,
        "n_cell_types": len(cell_types),
        "kaiser_n_components": n_kaiser,
        "pcs_for_80pct": n_80pct,
        "pcs_for_90pct": n_90pct,
        "top3_cumvar": round(float(cumvar[2]), 1),
        "eigenvalues": [round(float(ev), 3) for ev in eigenvalues],
        "var_explained": [round(float(ve), 1) for ve in var_explained],
        "null_kaiser_mean": round(float(np.mean(rand_kaiser_counts)), 1),
        "null_top3_mean": round(float(np.mean(rand_top3_var)), 1),
        "null_top3_std": round(float(np.std(rand_top3_var)), 1),
        "predictions": {
            "v1_all_positive": bool(all_pos_v1),
            "v2_danger_vs_tissue": bool(opposite),
            "top3_gt_80pct": bool(cumvar[2] > 80),
        },
        "correlation_matrix": [[round(float(corr[i,j]), 3) for j in range(14)] for i in range(14)],
        "module_order": MODULE_ORDER,
    }

    with open(DATA_DIR / "eigendecomposition.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n  Kaiser: {n_kaiser} significant PCs")
    print(f"  Top 3 PCs: {cumvar[2]:.1f}% variance")
    print(f"  Null model (random): {np.mean(rand_top3_var):.1f}% ± {np.std(rand_top3_var):.1f}%")
    print(f"\n  Predictions confirmed: {sum(results['predictions'].values())}/3")
    print(f"\nSaved to eigendecomposition.json")


if __name__ == "__main__":
    main()
