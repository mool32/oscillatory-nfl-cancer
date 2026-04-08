#!/usr/bin/env python3
"""
Test A: Maintenance cost per eigenmode → cancer vulnerability.

For each of 6 Kaiser eigenmodes:
  - Compute "maintenance cost" = mean expression of genes heavily loaded on that mode
  - Use HPA tissue data as proxy for protein production cost

Prediction: PC5 (where CGC genes cluster, p=0.0001) should have
            HIGHEST maintenance cost among 6 eigenmodes.

Falsification: if PC5 has average/low cost → energetic principle doesn't explain
               cancer gene clustering on this axis.
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

ALL_MODULES = {
    "NF-κB": {"RELA", "RELB", "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "TNFAIP3",
              "TRAF2", "TRAF6", "IKBKB", "IKBKG", "CHUK", "MAP3K7", "TAB1", "TAB2"},
    "ERK/MAPK": {"MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BRAF", "RAF1", "ARAF",
                 "HRAS", "KRAS", "NRAS", "SOS1", "GRB2", "DUSP1", "DUSP6", "SPRY2"},
    "JAK-STAT": {"JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT3", "STAT5A", "STAT5B",
                 "SOCS1", "SOCS3", "CISH", "PIAS1"},
    "p53": {"TP53", "MDM2", "MDM4", "CDKN1A", "BAX", "BBC3", "PMAIP1",
            "ATM", "ATR", "CHEK1", "CHEK2"},
    "Wnt": {"CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B", "DVL1",
            "TCF7L2", "LEF1", "RNF43", "ZNRF3"},
    "Notch": {"NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "FBXW7", "HES1", "HEY1",
              "MAML1", "RBPJ", "DLL1", "DLL4", "JAG1", "JAG2"},
    "Hippo": {"YAP1", "WWTR1", "LATS1", "LATS2", "STK3", "STK4", "SAV1",
              "MOB1A", "NF2", "TEAD1", "TEAD4"},
    "TGF-β": {"TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7",
              "SMURF1", "SMURF2", "BMPR1A", "BMPR2", "ACVR1"},
    "mTOR": {"MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "RHEB", "RPS6KB1",
             "EIF4EBP1", "DEPTOR", "MLST8"},
    "Calcium": {"PLCG1", "PLCG2", "ITPR1", "ITPR2", "ATP2A2", "CALM1",
                "NFATC1", "NFATC2", "CAMK2A", "CAMK2B"},
    "Cell Cycle": {"CDK2", "CDK4", "CDK6", "CCND1", "CCNE1", "CCNA2", "CCNB1",
                   "RB1", "E2F1", "CDKN1A", "CDKN2A", "CDKN1B", "CDC25A"},
    "Circadian": {"CLOCK", "ARNTL", "PER1", "PER2", "CRY1", "CRY2", "NR1D1",
                  "NR1D2", "CSNK1D", "CSNK1E", "FBXL3"},
    "NRF2": {"NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1",
             "SOD2", "CAT", "GPX1"},
    "PI3K/PTEN": {"PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "AKT1", "AKT2",
                  "PDK1", "INPP4B"},
    "AMPK": {"PRKAA1", "PRKAA2", "PRKAB1", "PRKAG1", "STK11", "ACACB",
             "PPARGC1A", "FOXO3", "TSC2", "CREB1"},
    "SREBP": {"SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2", "HMGCR",
              "FASN", "SCD", "ACLY"},
    "ATR/CHK1": {"ATR", "CHEK1", "ATRIP", "TOPBP1", "WEE1", "CDC25A",
                 "RAD17", "RAD9A", "HUS1"},
    "Rho/ROCK": {"RHOA", "RHOB", "RHOC", "ROCK1", "ROCK2", "MKL1",
                 "LIMK1", "CFL1", "ARHGAP1", "ARHGAP5"},
    "PPAR/LXR": {"PPARA", "PPARG", "PPARD", "NR1H3", "NR1H2", "RXRA",
                 "NCOR1", "NCOR2", "NCOA1", "NCOA2"},
    "Autophagy": {"ULK1", "ULK2", "BECN1", "ATG5", "ATG7", "ATG12",
                  "TFEB", "SQSTM1", "MAP1LC3B"},
}

MODULE_ORDER = [
    "NF-κB", "ERK/MAPK", "JAK-STAT", "p53", "Wnt", "Notch",
    "Hippo", "TGF-β", "mTOR", "Calcium", "Cell Cycle",
    "Circadian", "NRF2", "PI3K/PTEN",
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


def load_cgc():
    cgc = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                cgc.add(g)
    return cgc


def main():
    expr, cell_types = load_expression()
    cgc = load_cgc()

    # Recompute eigendecomposition
    N = len(MODULE_ORDER)
    M = len(cell_types)
    mat = np.zeros((N, M))
    for i, mod in enumerate(MODULE_ORDER):
        genes = [g for g in ALL_MODULES[mod] if g in expr]
        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in genes]
            mat[i, j] = np.mean(vals) if vals else 0

    corr = np.corrcoef(mat)
    eigenvalues, eigenvectors = np.linalg.eigh(corr)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    N_PCS = 6

    # ===================================================================
    # Maintenance cost per eigenmode
    # ===================================================================
    print("=" * 80)
    print("TEST A: MAINTENANCE COST PER EIGENMODE")
    print("=" * 80)

    # Method 1: Weighted module expression
    # For each PC: identify modules with |loading| > 0.2 (substantially loaded)
    # Maintenance cost = mean expression of those modules' genes across all cell types

    # First: compute mean expression per module across all cell types
    module_mean_expr = {}
    module_total_expr = {}  # total = mean × n_genes (proxy for total protein)
    for mod in MODULE_ORDER:
        genes = [g for g in ALL_MODULES[mod] if g in expr]
        if not genes:
            module_mean_expr[mod] = 0
            module_total_expr[mod] = 0
            continue
        vals = []
        for g in genes:
            gene_mean = np.mean([expr[g].get(ct, 0) for ct in cell_types])
            vals.append(gene_mean)
        module_mean_expr[mod] = np.mean(vals)
        module_total_expr[mod] = np.sum(vals)

    print(f"\n  Module expression levels (mean nCPM across cell types):")
    print(f"  {'Module':15s} {'Mean nCPM':>10s} {'Total nCPM':>11s} {'N genes':>7s}")
    print("  " + "-" * 47)
    for mod in sorted(MODULE_ORDER, key=lambda m: -module_mean_expr[m]):
        n = len([g for g in ALL_MODULES[mod] if g in expr])
        print(f"  {mod:15s} {module_mean_expr[mod]:10.1f} {module_total_expr[mod]:11.1f} {n:>7d}")

    # Method 2: Eigenmode maintenance cost
    # For each PC: sum of |loading_i| × module_expression_i
    # This weights each module's expression by how much it contributes to the mode

    print(f"\n{'='*80}")
    print("EIGENMODE MAINTENANCE COST")
    print("=" * 80)

    pc_costs = []
    pc_cgc_density = []

    # Also compute CGC density per eigenmode
    # Genes loaded on each PC: weight = sum over modules of |loading| × (gene in module)
    all_module_genes = set()
    for genes in ALL_MODULES.values():
        all_module_genes.update(genes)

    print(f"\n  {'PC':>4s} {'λ':>6s} {'Var%':>6s} {'Maint. Cost':>12s} {'CGC density':>12s} {'Top modules'}")
    print("  " + "-" * 80)

    for pc in range(N_PCS):
        vec = eigenvectors[:, pc]

        # Maintenance cost: weighted sum of module expression
        cost = 0
        for i, mod in enumerate(MODULE_ORDER):
            cost += abs(vec[i]) * module_total_expr[mod]
        pc_costs.append(cost)

        # CGC density weighted by loading
        # For each gene: compute its "loading" on this PC based on module membership
        n_cgc_weighted = 0
        n_total_weighted = 0
        for gene in all_module_genes:
            # Gene's loading = max |module loading| across its modules
            gene_loading = 0
            for i, mod in enumerate(MODULE_ORDER):
                if gene in ALL_MODULES[mod]:
                    gene_loading = max(gene_loading, abs(vec[i]))

            if gene_loading > 0.15:  # substantially loaded on this PC
                n_total_weighted += 1
                if gene in cgc:
                    n_cgc_weighted += 1

        cgc_dens = n_cgc_weighted / n_total_weighted if n_total_weighted > 0 else 0
        pc_cgc_density.append(cgc_dens)

        # Top loaded modules
        top_mods = sorted(range(N), key=lambda i: -abs(vec[i]))[:3]
        top_str = ", ".join(f"{MODULE_ORDER[i]}({vec[i]:+.2f})" for i in top_mods)

        print(f"  PC{pc+1:1d} {eigenvalues[pc]:6.2f} {eigenvalues[pc]/sum(eigenvalues[:N])*100:5.1f}% "
              f"{cost:12.0f} {cgc_dens:11.1%}    {top_str}")

    # ===================================================================
    # KEY TEST: Does PC5 have highest maintenance cost?
    # ===================================================================
    print(f"\n{'='*80}")
    print("PRE-REGISTERED TEST: Does PC5 (CGC cluster) have highest maintenance cost?")
    print("=" * 80)

    pc5_rank = sorted(range(N_PCS), key=lambda i: -pc_costs[i]).index(4) + 1
    pc5_cost = pc_costs[4]
    max_cost_pc = np.argmax(pc_costs) + 1
    max_cost = max(pc_costs)

    print(f"\n  PC5 maintenance cost: {pc5_cost:.0f}")
    print(f"  PC5 rank among 6 PCs: {pc5_rank}/6")
    print(f"  Highest cost: PC{max_cost_pc} ({max_cost:.0f})")

    if pc5_rank == 1:
        print(f"\n  ✓ PREDICTION CONFIRMED: PC5 has highest maintenance cost")
    else:
        print(f"\n  ✗ Prediction failed: PC{max_cost_pc} has highest cost, not PC5")

    # ===================================================================
    # Correlation: maintenance cost vs CGC density per eigenmode
    # ===================================================================
    print(f"\n{'='*80}")
    print("CORRELATION: Maintenance cost ↔ CGC density per eigenmode")
    print("=" * 80)

    costs = np.array(pc_costs)
    densities = np.array(pc_cgc_density)

    rho, p = stats.spearmanr(costs, densities)
    r, p_r = stats.pearsonr(costs, densities)

    print(f"\n  Spearman ρ(cost, CGC density): {rho:.3f}, p = {p:.4f}")
    print(f"  Pearson r(cost, CGC density): {r:.3f}, p = {p_r:.4f}")

    print(f"\n  Per-PC summary:")
    for pc in range(N_PCS):
        print(f"    PC{pc+1}: cost={pc_costs[pc]:8.0f}, CGC density={pc_cgc_density[pc]:.1%}")

    # ===================================================================
    # Alternative: per-gene analysis
    # ===================================================================
    print(f"\n{'='*80}")
    print("PER-GENE: Expression level of CGC vs non-CGC genes per eigenmode")
    print("=" * 80)

    for pc in range(N_PCS):
        vec = eigenvectors[:, pc]

        cgc_expr_vals = []
        noncgc_expr_vals = []

        for gene in all_module_genes:
            gene_loading = 0
            for i, mod in enumerate(MODULE_ORDER):
                if gene in ALL_MODULES[mod]:
                    gene_loading = max(gene_loading, abs(vec[i]))

            if gene_loading < 0.15:
                continue

            # Mean expression across cell types
            if gene in expr:
                gene_mean = np.mean([expr[gene].get(ct, 0) for ct in cell_types])
            else:
                continue

            if gene in cgc:
                cgc_expr_vals.append(gene_mean)
            else:
                noncgc_expr_vals.append(gene_mean)

        if cgc_expr_vals and noncgc_expr_vals:
            t, p = stats.mannwhitneyu(cgc_expr_vals, noncgc_expr_vals, alternative='greater')
            print(f"  PC{pc+1}: CGC expr={np.median(cgc_expr_vals):.1f}, "
                  f"non-CGC={np.median(noncgc_expr_vals):.1f}, "
                  f"U-test p={p:.4f} (CGC > non-CGC?)")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "pc5_cost": round(float(pc5_cost), 0),
        "pc5_rank": pc5_rank,
        "max_cost_pc": int(max_cost_pc),
        "costs": [round(float(c), 0) for c in pc_costs],
        "cgc_densities": [round(float(d), 3) for d in pc_cgc_density],
        "cost_cgc_rho": round(float(rho), 3),
        "cost_cgc_p": round(float(p), 4),
    }

    with open(DATA_DIR / "maintenance_cost.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to maintenance_cost.json")


if __name__ == "__main__":
    main()
