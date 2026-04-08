#!/usr/bin/env python3
"""
Test C: Mutation order in cancer progression → eigenspace trajectory.

Pre-registered prediction:
  Order of mutations follows DESCENDING IA+Bottleneck score.
  First mutation hits module with highest IA+BN, then next highest, etc.

Well-established progression sequences:
  Colorectal: APC → KRAS → TP53 → SMAD4
  Pancreatic: KRAS → CDKN2A → TP53 → SMAD4
  Endometrial: PTEN → PIK3CA → CTNNB1 → TP53
  Lung adeno: KRAS → TP53 → STK11 → KEAP1
  Melanoma: BRAF → CDKN2A → TP53 → PTEN

For each: map gene → module → IA+BN score → check if descending.
Also: project into eigenspace → visualize trajectory.
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

# IA + Bottleneck score per module (from step22 + step26)
IA_BN = {
    "p53": 4,         # IA=3 + BN=1
    "Wnt": 3,         # IA=2 + BN=1
    "TGF-β": 3,       # IA=2 + BN=1
    "JAK-STAT": 2,    # IA=2 + BN=0
    "PI3K/PTEN": 2,   # IA=1 + BN=1
    "Cell Cycle": 2,  # IA=1 + BN=1
    "mTOR": 1,        # IA=0 + BN=1
    "NRF2": 1,        # IA=0 + BN=1
    "AMPK": 1,        # IA=0 + BN=1
    "ATR/CHK1": 1,    # IA=0 + BN=1
    "ERK/MAPK": 1,    # IA=1 + BN=0
    "PPAR/LXR": 1,    # IA=1 + BN=0 (adipogenesis irrev)
    "Hippo": 1,       # IA=1 + BN=0
    "Notch": 1,       # IA=1 + BN=0
    "NF-κB": 1,       # IA=1 + BN=0
    "Autophagy": 1,   # IA=0 + BN=1 (BECN1)
    "Calcium": 0,
    "Circadian": 0,
    "Rho/ROCK": 0,
    "SREBP": 0,
}


# Cancer progression sequences (well-established from literature)
PROGRESSIONS = {
    "Colorectal (Vogelstein)": [
        ("APC", "Wnt", "initiating — loss of Wnt brake"),
        ("KRAS", "ERK/MAPK", "early progression — constitutive growth"),
        ("TP53", "p53", "late — loss of apoptosis/senescence"),
        ("SMAD4", "TGF-β", "late — loss of growth inhibition"),
    ],
    "Pancreatic (PDAC)": [
        ("KRAS", "ERK/MAPK", "initiating — G12D in >90%"),
        ("CDKN2A", "Cell Cycle", "early — p16 loss"),
        ("TP53", "p53", "intermediate"),
        ("SMAD4", "TGF-β", "late — in ~55% of cases"),
    ],
    "Endometrial (Type I)": [
        ("PTEN", "PI3K/PTEN", "initiating — PTEN loss"),
        ("PIK3CA", "PI3K/PTEN", "early — PI3K activation"),
        ("CTNNB1", "Wnt", "intermediate — β-catenin"),
        ("TP53", "p53", "late — rare in Type I, more in Type II"),
    ],
    "Lung adenocarcinoma": [
        ("KRAS", "ERK/MAPK", "initiating driver"),
        ("TP53", "p53", "early co-mutation"),
        ("STK11", "AMPK", "progression — LKB1 loss"),
        ("KEAP1", "NRF2", "late — metabolic adaptation"),
    ],
    "Melanoma": [
        ("BRAF", "ERK/MAPK", "initiating — V600E"),
        ("CDKN2A", "Cell Cycle", "early — p16 loss"),
        ("TP53", "p53", "variable timing"),
        ("PTEN", "PI3K/PTEN", "late — PI3K/AKT activation"),
    ],
    "CML/Hematological": [
        ("JAK2", "JAK-STAT", "initiating — V617F"),
        ("TP53", "p53", "blast crisis / transformation"),
        ("NOTCH1", "Notch", "T-ALL progression"),
    ],
}


def gene_to_module(gene):
    """Find which module(s) a gene belongs to."""
    modules = []
    for mod, genes in ALL_MODULES.items():
        if gene in genes:
            modules.append(mod)
    return modules


def main():
    print("=" * 80)
    print("TEST C: MUTATION ORDER vs IA+BOTTLENECK SCORE")
    print("=" * 80)

    print(f"\nPrediction: first mutations hit modules with HIGHEST IA+BN score.")
    print(f"Perfect ordering = descending IA+BN with each successive mutation.\n")

    all_concordances = []
    all_sequences = []

    for cancer_type, sequence in PROGRESSIONS.items():
        print(f"\n{'─'*60}")
        print(f"  {cancer_type}")
        print(f"{'─'*60}")

        print(f"  {'Step':>5s} {'Gene':>8s} {'Module':>12s} {'IA+BN':>6s} {'Timing'}")
        print("  " + "-" * 55)

        scores = []
        for i, (gene, module, timing) in enumerate(sequence):
            score = IA_BN.get(module, 0)
            scores.append(score)
            print(f"  {i+1:>5d} {gene:>8s} {module:>12s} {score:>6d}   {timing}")

        # Check if scores are monotonically descending
        n_concordant = 0
        n_discordant = 0
        n_pairs = 0
        for i in range(len(scores)):
            for j in range(i+1, len(scores)):
                n_pairs += 1
                if scores[i] > scores[j]:
                    n_concordant += 1
                elif scores[i] < scores[j]:
                    n_discordant += 1
                # ties don't count

        kendall_tau = (n_concordant - n_discordant) / n_pairs if n_pairs > 0 else 0

        print(f"\n  Score sequence: {scores}")
        print(f"  Concordant pairs: {n_concordant}/{n_pairs}")
        print(f"  Kendall τ (descending): {kendall_tau:+.2f}")

        all_concordances.append(kendall_tau)
        all_sequences.append({
            "cancer": cancer_type,
            "genes": [s[0] for s in sequence],
            "modules": [s[1] for s in sequence],
            "scores": scores,
            "kendall_tau": kendall_tau,
        })

    # ===================================================================
    # AGGREGATE across all cancer types
    # ===================================================================
    print(f"\n{'='*80}")
    print("AGGREGATE: IA+BN score vs mutation order position")
    print("=" * 80)

    # Pool all mutation positions
    all_positions = []  # (normalized_position, IA_BN_score)
    for seq in all_sequences:
        n = len(seq["scores"])
        for i, score in enumerate(seq["scores"]):
            norm_pos = i / (n - 1) if n > 1 else 0  # 0=first, 1=last
            all_positions.append((norm_pos, score))

    positions = np.array([p[0] for p in all_positions])
    ia_scores = np.array([p[1] for p in all_positions])

    rho, p = stats.spearmanr(positions, ia_scores)
    print(f"\n  All mutations pooled (N={len(all_positions)}):")
    print(f"  Spearman ρ(position, IA+BN): {rho:+.3f}, p = {p:.4f}")
    print(f"  (Negative ρ = high score mutated first = prediction confirmed)")

    # Mean Kendall tau across cancer types
    mean_tau = np.mean(all_concordances)
    t_tau, p_tau = stats.ttest_1samp(all_concordances, 0)
    print(f"\n  Mean Kendall τ across {len(all_concordances)} cancer types: {mean_tau:+.3f}")
    print(f"  One-sample t-test (τ > 0 = descending): t = {t_tau:.2f}, p = {p_tau:.4f}")

    # ===================================================================
    # Per mutation position: what's the mean IA+BN?
    # ===================================================================
    print(f"\n{'='*80}")
    print("MEAN IA+BN BY MUTATION POSITION")
    print("=" * 80)

    for pos in range(1, 5):
        scores_at_pos = [seq["scores"][pos-1] for seq in all_sequences if len(seq["scores"]) >= pos]
        if scores_at_pos:
            print(f"  Position {pos} (mutation #{pos}): mean IA+BN = {np.mean(scores_at_pos):.2f} "
                  f"(N={len(scores_at_pos)}, range {min(scores_at_pos)}-{max(scores_at_pos)})")

    # ===================================================================
    # Alternative test: initiating mutation has highest score?
    # ===================================================================
    print(f"\n{'='*80}")
    print("INITIATING MUTATION = HIGHEST IA+BN?")
    print("=" * 80)

    n_first_max = 0
    n_total = 0
    for seq in all_sequences:
        if len(seq["scores"]) < 2:
            continue
        n_total += 1
        if seq["scores"][0] >= max(seq["scores"]):
            n_first_max += 1
            print(f"  ✓ {seq['cancer']}: first ({seq['genes'][0]}, IA+BN={seq['scores'][0]}) ≥ rest")
        else:
            max_idx = seq["scores"].index(max(seq["scores"]))
            print(f"  ✗ {seq['cancer']}: first ({seq['genes'][0]}, IA+BN={seq['scores'][0]}) "
                  f"< {seq['genes'][max_idx]} (IA+BN={seq['scores'][max_idx]})")

    print(f"\n  First mutation is highest-score: {n_first_max}/{n_total} ({100*n_first_max/n_total:.0f}%)")

    # Binomial test (chance = 1/N_modules ≈ 5%)
    # Under null: first mutation hits any module randomly
    # P(first = highest among 4 sequential) ≈ 0.25
    p_binom = stats.binomtest(n_first_max, n_total, 0.25, alternative='greater')
    print(f"  Binomial test (vs chance 25%): p = {p_binom.pvalue:.4f}")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "n_cancer_types": len(all_sequences),
        "n_total_mutations": len(all_positions),
        "pooled_rho": round(float(rho), 3),
        "pooled_p": round(float(p), 4),
        "mean_kendall_tau": round(float(mean_tau), 3),
        "tau_ttest_p": round(float(p_tau), 4),
        "first_is_highest": n_first_max,
        "first_is_highest_n": n_total,
        "sequences": [{
            "cancer": s["cancer"],
            "genes": s["genes"],
            "modules": s["modules"],
            "scores": s["scores"],
            "tau": round(s["kendall_tau"], 2),
        } for s in all_sequences],
    }

    with open(DATA_DIR / "mutation_trajectories.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  Pooled position vs IA+BN: ρ = {rho:+.3f}, p = {p:.4f}")
    print(f"  Mean Kendall τ: {mean_tau:+.3f}, p = {p_tau:.4f}")
    print(f"  First mutation = highest score: {n_first_max}/{n_total}")
    if rho < 0 and p < 0.05:
        print(f"\n  ✓ CONFIRMED: High IA+BN modules mutated first in cancer progression")
    elif rho < 0 and p < 0.10:
        print(f"\n  ~ SUGGESTIVE: Trend toward high-score-first, not significant")
    else:
        print(f"\n  ✗ NOT CONFIRMED: Mutation order does not follow IA+BN")

    print(f"\nSaved to mutation_trajectories.json")


if __name__ == "__main__":
    main()
