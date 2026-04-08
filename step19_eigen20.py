#!/usr/bin/env python3
"""
#56: Expanded 20×20 eigendecomposition with all modules.

Same method as step16 but with 6 new modules added.
Key questions:
  - Do the same four axes emerge?
  - Where do AMPK, SREBP, Rho/ROCK, Autophagy load?
  - Is structure more or less concentrated?
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

ALL_MODULES = {
    # Original 14
    "NF-κB": {"RELA", "NFKBIA", "TNFAIP3", "NFKB1"},
    "ERK/MAPK": {"MAPK1", "MAPK3", "DUSP1", "DUSP6"},
    "JAK-STAT": {"JAK1", "JAK2", "STAT3", "SOCS3", "SOCS1"},
    "p53": {"TP53", "MDM2", "CDKN1A"},
    "Wnt": {"CTNNB1", "AXIN2", "APC", "GSK3B"},
    "Notch": {"NOTCH1", "FBXW7", "HES1", "DLL1"},
    "Hippo": {"YAP1", "LATS1", "LATS2", "STK3"},
    "TGF-β": {"TGFBR1", "SMAD2", "SMAD3", "SMAD7"},
    "mTOR": {"MTOR", "TSC1", "TSC2", "RPS6KB1"},
    "Calcium": {"PLCG1", "ATP2A2", "ITPR1", "CALM1"},
    "Cell Cycle": {"CDK2", "CDKN1A", "CCNE1", "RB1", "E2F1"},
    "Circadian": {"CLOCK", "PER2", "CRY1", "NR1D1"},
    "NRF2": {"NFE2L2", "KEAP1", "HMOX1"},
    "PI3K/PTEN": {"PIK3CA", "PTEN", "AKT1"},
    # New 6
    "AMPK": {"PRKAA1", "PRKAA2", "STK11", "ACACB"},
    "SREBP": {"SREBF1", "SREBF2", "INSIG1", "SCAP"},
    "ATR/CHK1": {"ATR", "CHEK1", "WEE1", "CDC25A"},
    "Rho/ROCK": {"RHOA", "ROCK1", "MKL1", "CFL1"},
    "PPAR/LXR": {"PPARA", "PPARG", "NR1H3", "NCOR1"},
    "Autophagy": {"TFEB", "BECN1", "ULK1", "ATG7"},
}

MODULE_ORDER = [
    # Original 14 in same order
    "NF-κB", "ERK/MAPK", "JAK-STAT", "p53", "Wnt", "Notch",
    "Hippo", "TGF-β", "mTOR", "Calcium", "Cell Cycle",
    "Circadian", "NRF2", "PI3K/PTEN",
    # New 6
    "AMPK", "SREBP", "ATR/CHK1", "Rho/ROCK", "PPAR/LXR", "Autophagy",
]


def load_expression():
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types


def main():
    print("Loading expression data...")
    expr, cell_types = load_expression()
    N = len(MODULE_ORDER)
    M = len(cell_types)

    # Build activity matrix
    mat = np.zeros((N, M))
    for i, mod in enumerate(MODULE_ORDER):
        genes = [g for g in ALL_MODULES[mod] if g in expr]
        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in genes]
            mat[i, j] = np.mean(vals) if vals else 0

    # Correlation matrix
    corr = np.corrcoef(mat)
    eigenvalues, eigenvectors = np.linalg.eigh(corr)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    total = np.sum(np.maximum(eigenvalues, 0))
    var_explained = np.maximum(eigenvalues, 0) / total * 100
    cumvar = np.cumsum(var_explained)

    # ===================================================================
    # EIGENVALUES
    # ===================================================================
    print(f"\n{'='*80}")
    print(f"20×20 EIGENDECOMPOSITION ({N} modules × {M} cell types)")
    print("=" * 80)

    n_kaiser = sum(1 for ev in eigenvalues if ev > 1.0)
    n_80 = next(i+1 for i, cv in enumerate(cumvar) if cv >= 80)

    print(f"\n  {'PC':>3s} {'Eigenvalue':>11s} {'Var %':>7s} {'Cum %':>7s} {'Kaiser':>7s}")
    print("  " + "-" * 40)
    for i in range(N):
        k = "✓" if eigenvalues[i] > 1.0 else ""
        print(f"  {i+1:3d} {eigenvalues[i]:11.3f} {var_explained[i]:7.1f} {cumvar[i]:7.1f} {k:>7s}")

    print(f"\n  Kaiser: {n_kaiser} components (was 3 with 14 modules)")
    print(f"  PCs for 80%: {n_80} (was 7)")
    print(f"  Top 3 cumvar: {cumvar[2]:.1f}% (was 56.7%)")

    # ===================================================================
    # EIGENVECTORS — top 5 PCs
    # ===================================================================
    print(f"\n{'='*80}")
    print("EIGENVECTOR LOADINGS")
    print("=" * 80)

    for pc in range(min(5, N)):
        if eigenvalues[pc] < 0.5:
            break
        vec = eigenvectors[:, pc]
        print(f"\n  PC{pc+1} (λ={eigenvalues[pc]:.2f}, {var_explained[pc]:.1f}%):")

        sorted_idx = np.argsort(np.abs(vec))[::-1]
        for j in sorted_idx:
            bar = "█" * int(abs(vec[j]) * 20)
            sign = "+" if vec[j] > 0 else "−"
            new = " ★" if MODULE_ORDER[j] in ["AMPK", "SREBP", "ATR/CHK1", "Rho/ROCK", "PPAR/LXR", "Autophagy"] else ""
            print(f"    {sign} {MODULE_ORDER[j]:12s} {vec[j]:+.3f}  {bar}{new}")

    # ===================================================================
    # STABILITY CHECK: do original axes survive?
    # ===================================================================
    print(f"\n{'='*80}")
    print("STABILITY CHECK: Do original 14-module axes survive?")
    print("=" * 80)

    # Check PC1: is it still general signaling (all negative except mTOR)?
    v1 = eigenvectors[:, 0]
    orig_14_idx = list(range(14))
    mtor_idx = MODULE_ORDER.index("mTOR")
    non_mtor_neg = sum(1 for i in orig_14_idx if i != mtor_idx and v1[i] < 0)
    print(f"\n  PC1 stability:")
    print(f"    Original 14 (excl mTOR) with same sign: {non_mtor_neg}/13")
    print(f"    mTOR loading: {v1[mtor_idx]:+.3f}")

    # Check PC3 (danger vs tissue): NF-κB/ERK/JAK-STAT vs Hippo/TGF-β/Circadian
    # Find which PC best matches this pattern
    danger_idx = [MODULE_ORDER.index(m) for m in ["NF-κB", "ERK/MAPK", "JAK-STAT"]]
    tissue_idx = [MODULE_ORDER.index(m) for m in ["Hippo", "TGF-β", "Circadian"]]

    best_pc = -1
    best_sep = 0
    for pc in range(min(6, N)):
        vec = eigenvectors[:, pc]
        danger_mean = np.mean([vec[i] for i in danger_idx])
        tissue_mean = np.mean([vec[i] for i in tissue_idx])
        sep = danger_mean - tissue_mean
        if abs(sep) > abs(best_sep):
            best_sep = sep
            best_pc = pc

    vec_best = eigenvectors[:, best_pc]
    print(f"\n  Danger vs tissue axis: best match = PC{best_pc+1}")
    print(f"    Danger mean: {np.mean([vec_best[i] for i in danger_idx]):+.3f}")
    print(f"    Tissue mean: {np.mean([vec_best[i] for i in tissue_idx]):+.3f}")
    print(f"    Separation: {abs(best_sep):.3f}")

    # ===================================================================
    # WHERE DO NEW MODULES LOAD?
    # ===================================================================
    print(f"\n{'='*80}")
    print("NEW MODULE POSITIONS IN PERCEPTION SPACE")
    print("=" * 80)

    new_mods = ["AMPK", "SREBP", "ATR/CHK1", "Rho/ROCK", "PPAR/LXR", "Autophagy"]
    print(f"\n  {'Module':15s} {'PC1':>7s} {'PC2':>7s} {'PC3':>7s} {'PC4':>7s}  Interpretation")
    print("  " + "-" * 75)

    for mod in new_mods:
        i = MODULE_ORDER.index(mod)
        pc1, pc2, pc3, pc4 = eigenvectors[i, 0], eigenvectors[i, 1], eigenvectors[i, 2], eigenvectors[i, 3]

        # Interpretation based on loading pattern
        interp = []
        if abs(pc1) > 0.2:
            interp.append(f"{'signaling' if pc1 < 0 else 'metabolic'}-aligned (PC1)")
        if abs(pc2) > 0.2:
            interp.append(f"{'growth' if pc2 > 0 else 'checkpoint'} (PC2)")
        if abs(pc3) > 0.2:
            interp.append(f"{'danger' if pc3 > 0 else 'tissue'} (PC3)")
        if not interp:
            interp.append("weak loading on main axes")

        print(f"  {mod:15s} {pc1:+7.3f} {pc2:+7.3f} {pc3:+7.3f} {pc4:+7.3f}  {'; '.join(interp)}")

    # ===================================================================
    # CORRELATION of new modules with original
    # ===================================================================
    print(f"\n{'='*80}")
    print("NEW MODULE CORRELATIONS WITH ORIGINAL 14")
    print("=" * 80)

    for mod in new_mods:
        i = MODULE_ORDER.index(mod)
        top_corrs = []
        for j in range(14):
            top_corrs.append((MODULE_ORDER[j], corr[i, j]))
        top_corrs.sort(key=lambda x: -abs(x[1]))
        top3 = top_corrs[:3]
        print(f"  {mod:15s} → {', '.join(f'{n}({r:+.2f})' for n, r in top3)}")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "n_modules": N,
        "n_cell_types": M,
        "kaiser": n_kaiser,
        "pcs_80": n_80,
        "top3_cumvar": round(float(cumvar[2]), 1),
        "eigenvalues": [round(float(ev), 3) for ev in eigenvalues],
        "var_explained": [round(float(ve), 1) for ve in var_explained],
        "module_order": MODULE_ORDER,
        "danger_tissue_best_pc": best_pc + 1,
    }
    with open(DATA_DIR / "eigen20.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)
    print(f"  14-module: Kaiser=3, top3={56.7}%, 80% at PC7")
    print(f"  20-module: Kaiser={n_kaiser}, top3={cumvar[2]:.1f}%, 80% at PC{n_80}")
    print(f"\nSaved to eigen20.json")


if __name__ == "__main__":
    main()
