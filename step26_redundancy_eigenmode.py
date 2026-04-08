#!/usr/bin/env python3
"""
Test A': Redundancy per eigenmode → cancer vulnerability.

Replaces failed Test A (maintenance cost).
Hypothesis: PC5 (where CGC clusters) has LOWEST redundancy among 6 eigenmodes.

Redundancy measures per module:
  1. Paralog count: number of paralogous gene pairs within module
  2. Module redundancy score (from step21): assigned redundancy per module
  3. Mean degree in inter-module network: how connected/backed-up

For each eigenmode: weight module redundancy by |loading| on that PC.
Low redundancy eigenmode → high CGC density.
"""

import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json
import csv

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

MODULE_ORDER = [
    "NF-κB", "ERK/MAPK", "JAK-STAT", "p53", "Wnt", "Notch",
    "Hippo", "TGF-β", "mTOR", "Calcium", "Cell Cycle",
    "Circadian", "NRF2", "PI3K/PTEN",
    "AMPK", "SREBP", "ATR/CHK1", "Rho/ROCK", "PPAR/LXR", "Autophagy",
]

# Module redundancy scores (from step21 analysis)
# = number of paralogous gene families that provide backup
MODULE_REDUNDANCY = {
    "NF-κB": 4,       # RELA/RELB/NFKB1/NFKB2
    "ERK/MAPK": 3,    # ERK1/2, RAF1/BRAF/ARAF, RAS×3
    "JAK-STAT": 4,    # JAK1/2/3/TYK2, STAT1/3/5A/5B/6
    "p53": 1,         # TP53 only
    "Wnt": 2,         # AXIN1/2, DVL1/2/3, but APC irreplaceable
    "Notch": 4,       # NOTCH1/2/3/4
    "Hippo": 3,       # YAP/TAZ, LATS1/2, MST1/2
    "TGF-β": 2,       # SMAD2/3, but SMAD4 irreplaceable
    "mTOR": 1,        # MTOR singular, TSC1/2 obligate
    "Calcium": 4,     # ITPR1/2/3, PLCG1/2, NFATC×4, CaMKII×4
    "Cell Cycle": 2,  # CDK4/6, CyclinD1/2/3, but RB1 irreplaceable
    "Circadian": 3,   # PER1/2, CRY1/2, NR1D1/2
    "NRF2": 1,        # NFE2L2 and KEAP1 both singular
    "PI3K/PTEN": 2,   # PIK3CA/CB, AKT1/2/3, but PTEN irreplaceable
    "AMPK": 2,        # AMPKα1/α2, but STK11 irreplaceable
    "SREBP": 2,       # SREBF1/2, INSIG1/2
    "ATR/CHK1": 1,    # ATR singular
    "Rho/ROCK": 3,    # RhoA/B/C, ROCK1/2, LIMK1/2
    "PPAR/LXR": 3,    # PPARα/γ/δ, LXRα/β, NCOR1/2
    "Autophagy": 2,   # ULK1/2, but BECN1/ATG5 singular
}

# Number of paralogous pairs per module (more granular)
MODULE_PARALOG_PAIRS = {
    "NF-κB": 6,       # RELA-RELB, NFKB1-NFKB2, NFKBIA-NFKBIB, TRAF2-TRAF6, IKBKB-CHUK, TAB1-TAB2
    "ERK/MAPK": 8,    # MAPK1-MAPK3, MAP2K1-MAP2K2, BRAF-RAF1-ARAF(3), HRAS-KRAS-NRAS(3), DUSP1-DUSP6
    "JAK-STAT": 9,    # JAK1-JAK2-JAK3-TYK2(6), STAT3-STAT5A-STAT5B(3), SOCS1-SOCS3
    "p53": 1,         # MDM2-MDM4 only; TP53 no paralog
    "Wnt": 2,         # AXIN1-AXIN2, TCF7L2-LEF1
    "Notch": 7,       # NOTCH1-4(6), DLL1-DLL4, JAG1-JAG2, HES1-HEY1
    "Hippo": 4,       # YAP1-WWTR1, LATS1-LATS2, STK3-STK4, TEAD1-TEAD4
    "TGF-β": 3,       # SMAD2-SMAD3, SMURF1-SMURF2, BMPR1A-BMPR2
    "mTOR": 1,        # RPTOR-RICTOR (different complexes, not true backup)
    "Calcium": 7,     # PLCG1-PLCG2, ITPR1-ITPR2(1), NFATC1-NFATC2, CAMK2A-CAMK2B, etc
    "Cell Cycle": 5,  # CDK4-CDK6, CCND1-CCNE1-CCNA2-CCNB1(partial), CDKN1A-CDKN2A-CDKN1B
    "Circadian": 4,   # PER1-PER2, CRY1-CRY2, NR1D1-NR1D2, CSNK1D-CSNK1E
    "NRF2": 1,        # GCLC-GCLM only; NFE2L2/KEAP1 singular
    "PI3K/PTEN": 3,   # PIK3CA-PIK3CB, PIK3R1 singular, AKT1-AKT2
    "AMPK": 2,        # PRKAA1-PRKAA2, PRKAB1-PRKAG1 (subunits, partial)
    "SREBP": 3,       # SREBF1-SREBF2, INSIG1-INSIG2, FASN-SCD(partial)
    "ATR/CHK1": 1,    # RAD9A-HUS1 (partial); ATR singular
    "Rho/ROCK": 5,    # RHOA-RHOB-RHOC(3), ROCK1-ROCK2, LIMK1-LIMK2(partial)
    "PPAR/LXR": 5,    # PPARA-PPARG-PPARD(3), NR1H3-NR1H2, NCOR1-NCOR2, NCOA1-NCOA2
    "Autophagy": 2,   # ULK1-ULK2, ATG5-ATG12(conjugation pair, not paralog)
}

# Has irreplaceable bottleneck? (gene with no paralog that entire module depends on)
MODULE_HAS_BOTTLENECK = {
    "NF-κB": False,    # Multiple subunits compensate
    "ERK/MAPK": False, # Parallel cascades
    "JAK-STAT": False, # Multiple JAKs and STATs
    "p53": True,       # TP53 irreplaceable
    "Wnt": True,       # APC, CTNNB1 irreplaceable
    "Notch": False,    # 4 paralogs
    "Hippo": False,    # Redundant pairs
    "TGF-β": True,     # SMAD4 irreplaceable
    "mTOR": True,      # MTOR singular
    "Calcium": False,  # Extensive paralogy
    "Cell Cycle": True, # RB1 irreplaceable
    "Circadian": False, # Redundant pairs
    "NRF2": True,      # NFE2L2 and KEAP1 singular
    "PI3K/PTEN": True, # PTEN irreplaceable
    "AMPK": True,      # STK11 irreplaceable
    "SREBP": False,    # SREBF1/2 partially redundant
    "ATR/CHK1": True,  # ATR singular
    "Rho/ROCK": False, # Triple paralogs
    "PPAR/LXR": False, # Triple paralogs
    "Autophagy": True, # BECN1 singular
}

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


def load_cgc():
    cgc = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                cgc.add(g)
    return cgc


def load_expression():
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types


def main():
    cgc = load_cgc()
    expr, cell_types = load_expression()

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
    # Compute redundancy metrics per eigenmode
    # ===================================================================
    print("=" * 80)
    print("TEST A': REDUNDANCY PER EIGENMODE → CANCER VULNERABILITY")
    print("=" * 80)

    # Three redundancy measures
    redun_scores = np.array([MODULE_REDUNDANCY[m] for m in MODULE_ORDER], dtype=float)
    paralog_scores = np.array([MODULE_PARALOG_PAIRS[m] for m in MODULE_ORDER], dtype=float)
    bottleneck_scores = np.array([1.0 if MODULE_HAS_BOTTLENECK[m] else 0.0 for m in MODULE_ORDER])

    # Compute per-eigenmode weighted redundancy
    print(f"\n  Module redundancy scores:")
    print(f"  {'Module':15s} {'Redun':>6s} {'Paralogs':>8s} {'Bottleneck':>10s}")
    print("  " + "-" * 45)
    for i, mod in enumerate(MODULE_ORDER):
        bn = "✓" if MODULE_HAS_BOTTLENECK[mod] else ""
        print(f"  {mod:15s} {MODULE_REDUNDANCY[mod]:>6d} {MODULE_PARALOG_PAIRS[mod]:>8d} {bn:>10s}")

    # Weighted redundancy per eigenmode
    pc_redun = []
    pc_paralogs = []
    pc_bottleneck_frac = []
    pc_cgc_density = []

    all_module_genes = set()
    for genes in ALL_MODULES.values():
        all_module_genes.update(genes)

    print(f"\n{'='*80}")
    print("EIGENMODE REDUNDANCY vs CGC DENSITY")
    print("=" * 80)

    print(f"\n  {'PC':>4s} {'λ':>6s} {'Wt.Redun':>9s} {'Wt.Paral':>9s} {'Btlnck%':>8s} {'CGC%':>6s}")
    print("  " + "-" * 50)

    for pc in range(N_PCS):
        vec = eigenvectors[:, pc]
        weights = np.abs(vec)
        weights = weights / weights.sum()  # normalize

        wr = np.dot(weights, redun_scores)
        wp = np.dot(weights, paralog_scores)
        wb = np.dot(weights, bottleneck_scores)

        pc_redun.append(wr)
        pc_paralogs.append(wp)
        pc_bottleneck_frac.append(wb)

        # CGC density for genes loaded on this PC
        n_cgc = 0
        n_total = 0
        for gene in all_module_genes:
            gene_loading = max(abs(vec[MODULE_ORDER.index(m)]) for m in MODULE_ORDER if gene in ALL_MODULES[m])
            if gene_loading > 0.15:
                n_total += 1
                if gene in cgc:
                    n_cgc += 1
        cgc_dens = n_cgc / n_total if n_total > 0 else 0
        pc_cgc_density.append(cgc_dens)

        print(f"  PC{pc+1:1d} {eigenvalues[pc]:6.2f} {wr:9.2f} {wp:9.1f} {wb*100:7.0f}% {cgc_dens*100:5.1f}%")

    # ===================================================================
    # CORRELATIONS
    # ===================================================================
    print(f"\n{'='*80}")
    print("CORRELATIONS: Redundancy metrics ↔ CGC density (N=6 eigenmodes)")
    print("=" * 80)

    cgc_arr = np.array(pc_cgc_density)

    metrics = {
        "Weighted redundancy": np.array(pc_redun),
        "Weighted paralogs": np.array(pc_paralogs),
        "Bottleneck fraction": np.array(pc_bottleneck_frac),
    }

    for name, arr in metrics.items():
        rho, p = stats.spearmanr(arr, cgc_arr)
        r, p_r = stats.pearsonr(arr, cgc_arr)
        print(f"\n  {name}:")
        print(f"    Spearman ρ = {rho:+.3f}, p = {p:.4f}")
        print(f"    Pearson  r = {r:+.3f}, p = {p_r:.4f}")
        direction = "Low redundancy → more cancer" if rho < 0 else "High redundancy → more cancer"
        print(f"    Direction: {direction}")

    # ===================================================================
    # KEY TEST: PC5 redundancy rank
    # ===================================================================
    print(f"\n{'='*80}")
    print("PRE-REGISTERED: Does PC5 have LOWEST redundancy?")
    print("=" * 80)

    pc5_redun_rank = sorted(range(N_PCS), key=lambda i: pc_redun[i]).index(4) + 1
    pc5_para_rank = sorted(range(N_PCS), key=lambda i: pc_paralogs[i]).index(4) + 1
    pc5_btlnck_rank = sorted(range(N_PCS), key=lambda i: -pc_bottleneck_frac[i]).index(4) + 1

    print(f"\n  PC5 weighted redundancy: {pc_redun[4]:.2f} (rank {pc5_redun_rank}/6, 1=lowest)")
    print(f"  PC5 weighted paralogs: {pc_paralogs[4]:.1f} (rank {pc5_para_rank}/6, 1=lowest)")
    print(f"  PC5 bottleneck fraction: {pc_bottleneck_frac[4]*100:.0f}% (rank {pc5_btlnck_rank}/6, 1=highest)")

    lowest_redun_pc = np.argmin(pc_redun) + 1
    highest_btlnck_pc = np.argmax(pc_bottleneck_frac) + 1

    print(f"\n  Lowest redundancy eigenmode: PC{lowest_redun_pc} ({pc_redun[lowest_redun_pc-1]:.2f})")
    print(f"  Highest bottleneck eigenmode: PC{highest_btlnck_pc} ({pc_bottleneck_frac[highest_btlnck_pc-1]*100:.0f}%)")

    # ===================================================================
    # MODULE-LEVEL: redundancy vs CGC per module (not eigenmode)
    # ===================================================================
    print(f"\n{'='*80}")
    print("MODULE-LEVEL: Redundancy → CGC (across 20 modules)")
    print("=" * 80)

    mod_cgc_frac = []
    for mod in MODULE_ORDER:
        genes = ALL_MODULES[mod]
        n_cgc = len(genes & cgc)
        mod_cgc_frac.append(n_cgc / len(genes) if genes else 0)

    mod_cgc_arr = np.array(mod_cgc_frac)

    rho_r, p_r = stats.spearmanr(redun_scores, mod_cgc_arr)
    rho_p, p_p = stats.spearmanr(paralog_scores, mod_cgc_arr)
    rho_b, p_b = stats.spearmanr(bottleneck_scores, mod_cgc_arr)

    print(f"\n  Redundancy score ↔ CGC: ρ = {rho_r:+.3f}, p = {p_r:.4f}")
    print(f"  Paralog pairs ↔ CGC: ρ = {rho_p:+.3f}, p = {p_p:.4f}")
    print(f"  Has bottleneck ↔ CGC: ρ = {rho_b:+.3f}, p = {p_b:.4f}")

    # Binary split: bottleneck vs no bottleneck
    bn_cgc = [mod_cgc_frac[i] for i, m in enumerate(MODULE_ORDER) if MODULE_HAS_BOTTLENECK[m]]
    no_bn_cgc = [mod_cgc_frac[i] for i, m in enumerate(MODULE_ORDER) if not MODULE_HAS_BOTTLENECK[m]]

    u, p_u = stats.mannwhitneyu(bn_cgc, no_bn_cgc, alternative='greater')
    print(f"\n  Bottleneck modules (N={len(bn_cgc)}): mean CGC = {np.mean(bn_cgc)*100:.1f}%")
    print(f"  No bottleneck (N={len(no_bn_cgc)}): mean CGC = {np.mean(no_bn_cgc)*100:.1f}%")
    print(f"  Mann-Whitney (bottleneck > no): U={u:.0f}, p = {p_u:.4f}")

    # ===================================================================
    # COMBINED: IA × bottleneck
    # ===================================================================
    print(f"\n{'='*80}")
    print("COMBINED: Irreversible Authority × Bottleneck")
    print("=" * 80)

    # From step22
    IA = {
        "p53": 3, "Wnt": 2, "TGF-β": 2, "JAK-STAT": 2,
        "PI3K/PTEN": 1, "Cell Cycle": 1, "ERK/MAPK": 1, "PPAR/LXR": 1,
        "Hippo": 1, "Notch": 1, "NF-κB": 1,
        "AMPK": 0, "mTOR": 0, "ATR/CHK1": 0, "Calcium": 0,
        "NRF2": 0, "Rho/ROCK": 0, "Autophagy": 0, "Circadian": 0, "SREBP": 0,
    }

    ia_arr = np.array([IA[m] for m in MODULE_ORDER], dtype=float)
    combined = ia_arr + bottleneck_scores  # IA + bottleneck (additive)

    rho_ia, p_ia = stats.spearmanr(ia_arr, mod_cgc_arr)
    rho_comb, p_comb = stats.spearmanr(combined, mod_cgc_arr)
    rho_btn, p_btn = stats.spearmanr(bottleneck_scores, mod_cgc_arr)

    print(f"\n  IA alone → CGC: ρ = {rho_ia:+.3f}, p = {p_ia:.4f}")
    print(f"  Bottleneck alone → CGC: ρ = {rho_btn:+.3f}, p = {p_btn:.4f}")
    print(f"  IA + Bottleneck → CGC: ρ = {rho_comb:+.3f}, p = {p_comb:.4f}")

    if abs(rho_comb) > abs(rho_ia):
        print(f"\n  ✓ Combined (IA + Bottleneck) outperforms IA alone")
    else:
        print(f"\n  IA alone is sufficient; bottleneck doesn't add")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "pc5_redundancy_rank": pc5_redun_rank,
        "pc5_bottleneck_rank": pc5_btlnck_rank,
        "eigenmode_redun_cgc_rho": round(float(stats.spearmanr(pc_redun, cgc_arr)[0]), 3),
        "eigenmode_redun_cgc_p": round(float(stats.spearmanr(pc_redun, cgc_arr)[1]), 4),
        "module_redun_cgc_rho": round(float(rho_r), 3),
        "module_bottleneck_cgc_rho": round(float(rho_b), 3),
        "bottleneck_mean_cgc": round(float(np.mean(bn_cgc)), 3),
        "no_bottleneck_mean_cgc": round(float(np.mean(no_bn_cgc)), 3),
        "bottleneck_p": round(float(p_u), 4),
        "ia_rho": round(float(rho_ia), 3),
        "ia_plus_bottleneck_rho": round(float(rho_comb), 3),
    }

    with open(DATA_DIR / "redundancy_eigenmode.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to redundancy_eigenmode.json")


if __name__ == "__main__":
    main()
