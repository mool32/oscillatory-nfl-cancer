#!/usr/bin/env python3
"""
Cancer Vulnerability Score (CVS) — Two-axis metric.

CVS = S × P = |RDS| × Paper1_polarity_sign

S = Structural vulnerability = |RDS| (number of irreplaceable feedback loops)
P = Phenotypic polarity = Paper 1 manual classification with growth-inhibitory inversion
    +1 = mutation → proliferative advantage (oncogene direction)
    -1 = mutation → growth arrest/death advantage (TSG direction)

Tests:
A. sign(CVS) vs actual CGC role (accuracy)
B. |CVS| vs COSMIC mutation frequency (Spearman)
C. |CVS| vs drug target status (Mann-Whitney)
D. ROC AUC comparison: |RDS| vs Paper1 binary vs |CVS|
E. Blind test on 12 NFL-only genes
"""

import csv
import json
import numpy as np
from scipy import stats
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# ---------------------------------------------------------------
# Paper 1 polarity sign: +1 = oncogene direction, -1 = TSG direction
# This is the MANUAL classification from Paper 1 including
# growth-inhibitory pathway inversion.
#
# Rise arm + growth-promoting pathway → GoF = proliferation → +1 (oncogene)
# Recovery arm + growth-promoting pathway → LoF = no brake → +1 actually wait...
#
# Let me be precise:
# For ONCOGENE prediction: gene whose mutation PROMOTES cancer
#   - Rise arm in growth-promoting: GoF → constitutive activation → +1
#   - Recovery arm in growth-inhibitory: LoF → removes growth arrest → +1
# For TSG prediction: gene whose mutation ALLOWS cancer
#   - Recovery arm in growth-promoting: LoF → removes brake → -1
#   - Rise arm in growth-inhibitory: LoF → removes growth suppression → -1
#
# Paper1_sign encodes: if this gene is mutated (GoF for rise, LoF for recovery),
# does it promote cancer (+1) or suppress cancer (-1)?
# ---------------------------------------------------------------

PAPER1_SIGN = {
    # NF-κB (growth-promoting) — rise = OG, recovery = TSG
    "RELA": +1, "NFKBIA": -1, "TNFAIP3": -1, "TRAF6": +1,

    # p53 (growth-INHIBITORY) — rise = TSG (activates growth arrest), recovery = OG (removes arrest)
    # TP53: recovery arm of growth cycle, but RISE arm of growth-inhibitory checkpoint
    # TP53 LoF → no checkpoint → cancer. So TP53 sign = -1 (TSG)
    # MDM2: recovery arm of p53 (degrades p53). MDM2 GoF → p53 removed → cancer. MDM2 sign = +1 (OG)
    "TP53": -1, "MDM2": +1,
    "PMAIP1": -1, "SIVA1": -1, "BID": -1, "BBC3": -1, "TP73": -1,
    "BCL2": +1,  # BCL2: anti-apoptotic. GoF → blocks apoptosis → cancer

    # ERK/MAPK (growth-promoting)
    "MAPK1": +1, "MAPK14": +1, "DUSP1": -1,

    # JAK-STAT (growth-promoting)
    "JAK1": +1, "JAK2": +1, "TYK2": +1,
    "STAT1": +1, "STAT2": +1, "STAT3": +1, "STAT5A": +1, "STAT5B": +1,
    "SOCS1": -1, "SOCS3": -1, "SOCS4": -1, "SOCS5": -1, "SOCS6": -1, "SOCS7": -1,
    "CISH": -1, "CSF1R": +1, "LEPR": +1, "PRLR": +1,

    # Wnt (growth-promoting)
    "CTNNB1": +1, "CDH1": -1, "AXIN1": -1, "AXIN2": -1,

    # Notch (context-dependent but mostly oncogenic in hematologic)
    "NOTCH1": +1, "FBXW7": -1,

    # mTOR/PI3K (growth-promoting)
    "MTOR": +1, "IRS1": +1, "IRS2": +1,
    "PIK3CA": +1, "PIK3R1": -1, "PIK3R3": +1,
    "AKT1": +1, "AKT2": +1, "AKT3": +1, "PTK2": +1,
    "PTEN": -1, "TSC1": -1,

    # TGF-β (growth-INHIBITORY in epithelial context)
    # SMAD2/3/4: rise arm of growth inhibition → LoF = lose growth arrest → TSG
    "TGFBR1": -1, "SMAD2": -1, "SMAD3": -1, "SMAD4": -1,
    "SMAD1": -1, "SMAD5": -1, "SMAD7": +1,  # SMAD7 recovery of GI → removes arrest → OG direction

    # BMP (growth-inhibitory context)
    "BMPR1A": -1, "BMPR1B": -1, "BMPR2": -1,

    # NRF2 (growth-promoting: cytoprotective)
    "NFE2L2": +1, "KEAP1": -1,

    # Hedgehog (growth-promoting)
    "GLI1": +1, "GLI2": +1, "GLI3": -1, "SUFU": -1,

    # Hippo (growth-INHIBITORY: LATS/MST suppress YAP)
    # YAP1: effector, GoF → proliferation → OG
    # LATS1: rise arm of growth inhibition → LoF = YAP active → TSG
    "YAP1": +1, "LATS1": -1, "CSNK1E": -1,

    # Cell cycle (growth-promoting for CDK/cyclin, inhibitory for CKIs/Rb)
    "CDK2": +1, "CCNE1": +1, "E2F1": +1, "E2F2": +1, "E2F3": +1, "TFDP1": +1,
    "CDKN1A": -1, "CDKN1B": -1, "CDKN2A": -1,
    "RB1": -1, "RBL1": -1, "RBL2": -1,

    # FOXO (growth-INHIBITORY: FOXO activates cell cycle arrest)
    "FOXO1": -1, "FOXO3": -1, "FOXO4": -1, "FOXO6": -1,

    # Calcium (signaling, context-dependent)
    "PLCG1": +1, "ATP2A2": -1,

    # NFAT
    "NFATC1": +1, "RCAN1": -1,
}

# CGC classification
CGC_ROLE = {
    "AKT1": "oncogene", "AKT2": "oncogene", "AKT3": "oncogene",
    "AXIN1": "TSG", "BCL2": "oncogene", "BMPR1A": "TSG",
    "CCNE1": "oncogene", "CDH1": "TSG", "CDK2": "oncogene",
    "CDKN1A": "TSG", "CDKN1B": "TSG", "CDKN2A": "TSG",
    "CSF1R": "oncogene", "CTNNB1": "oncogene",
    "E2F1": "oncogene",
    "FOXO1": "both", "FOXO3": "TSG", "FOXO4": "oncogene",
    "GLI1": "oncogene", "GLI2": "oncogene", "GLI3": "TSG",
    "JAK1": "oncogene", "JAK2": "oncogene",
    "MAPK1": "oncogene", "MDM2": "oncogene", "MTOR": "oncogene",
    "PIK3CA": "oncogene", "PIK3R1": "TSG",
    "PTEN": "TSG", "RB1": "TSG", "RELA": "oncogene",
    "SMAD2": "TSG", "SMAD3": "TSG", "SMAD4": "TSG",
    "SOCS1": "TSG", "STAT3": "oncogene", "STAT5B": "oncogene",
    "TGFBR1": "TSG", "TNFAIP3": "TSG", "TP53": "TSG",
}

# Approximate COSMIC mutation frequencies (% tumors with mutation)
# From COSMIC top mutated genes across all cancer types
COSMIC_FREQ = {
    "TP53": 28.6, "PIK3CA": 10.2, "PTEN": 5.8, "RB1": 4.8,
    "CDKN2A": 6.2, "AKT1": 1.8, "CTNNB1": 3.5, "MTOR": 2.1,
    "SMAD4": 4.1, "CDK2": 0.3, "STAT3": 1.2, "JAK2": 1.5,
    "MAPK1": 0.8, "MDM2": 3.2, "CCNE1": 2.8, "E2F1": 0.4,
    "CDH1": 3.0, "AXIN1": 1.5, "FBXW7": 4.5, "SOCS1": 1.0,
    "GLI1": 0.5, "JAK1": 0.8, "AKT2": 0.5, "BCL2": 1.8,
    "BMPR1A": 0.6, "TGFBR1": 0.8, "SMAD2": 1.8, "SMAD3": 0.7,
    "FOXO1": 0.6, "RELA": 0.3, "PIK3R1": 1.5,
    "CDKN1A": 0.8, "CDKN1B": 1.2, "STAT5B": 0.4,
    "GLI2": 0.6, "GLI3": 0.4, "CSF1R": 0.5,
    "TNFAIP3": 1.0, "RBL1": 0.1, "RBL2": 0.1,
    "FOXO3": 0.3, "FOXO4": 0.1, "PMAIP1": 0.1,
    "NFKBIA": 0.5, "NFE2L2": 2.5, "KEAP1": 2.0,
    "YAP1": 0.3, "LATS1": 0.5, "NOTCH1": 3.5,
}

# Drug target status
HAS_DRUG = {
    "MTOR": True, "PIK3CA": True, "AKT1": True, "CDK2": True,
    "JAK1": True, "JAK2": True, "STAT3": True, "MAPK1": True,
    "NOTCH1": True, "BCL2": True, "PTEN": False, "RB1": False,
    "TP53": False, "SMAD2": False, "SMAD3": False, "SMAD4": False,
    "CDKN2A": False, "CDKN1A": False, "CDKN1B": False,
    "CDH1": False, "AXIN1": False, "SOCS1": False,
    "TGFBR1": True, "BMPR1A": False,
    "GLI1": True, "GLI2": True, "NFE2L2": False, "KEAP1": False,
    "CCNE1": False, "E2F1": False, "MDM2": True,
    "CTNNB1": False, "CSF1R": True,
}

# Paper 1 genes (for overlap detection)
PAPER1_GENES = {
    "RELA", "RELB", "NFKB1", "NFKB2", "REL", "IKBKB", "CHUK",
    "NFKBIA", "NFKBIB", "TNFAIP3", "CYLD", "OTULIN",
    "TP53", "ATM", "ATR", "CHEK1", "CHEK2", "MDM2", "MDM4", "PPM1D", "USP7",
    "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2",
    "MAPK1", "MAPK3", "DUSP1", "DUSP2", "DUSP3", "DUSP4", "DUSP5", "DUSP6",
    "CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B",
    "NOTCH1", "NOTCH2", "FBXW7",
    "CLOCK", "ARNTL", "PER1", "PER2", "PER3", "CRY1", "CRY2",
    "YAP1", "WWTR1", "LATS1", "LATS2", "STK3", "STK4", "NF2",
    "MTOR", "RPTOR", "RHEB", "AKT1", "AKT2", "AKT3",
    "PIK3CA", "PIK3CB", "PIK3R1", "TSC1", "TSC2", "PTEN",
    "JAK1", "JAK2", "JAK3", "TYK2", "STAT3", "STAT5A", "STAT5B",
    "SOCS1", "SOCS3",
    "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD6", "SMAD7",
    "NFE2L2", "KEAP1",
    "GLI1", "GLI2", "PTCH1", "SUFU",
    "ERN1", "XBP1", "HSPA5",
    "PLCG1", "PLCG2", "ATP2A1", "ATP2A2",
}


def load_rds():
    with open(DATA_DIR / "rds_scores.csv") as f:
        reader = csv.DictReader(f)
        rds = {}
        for row in reader:
            rds[row["gene"]] = {
                "abs_RDS": int(row["abs_RDS"]),
                "RDS": int(row["RDS"]),
                "n_loops": int(row["n_loops"]),
                "n_irreplaceable": int(row["n_irreplaceable"]),
            }
    return rds


def load_cgc():
    genes = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g: genes.add(g)
    return genes


def compute_cvs(rds):
    """CVS = |RDS| × Paper1_sign."""
    cvs = {}
    for gene, r in rds.items():
        s = r["abs_RDS"]  # structural vulnerability
        p = PAPER1_SIGN.get(gene, 0)  # phenotypic polarity
        cvs[gene] = {
            "S": s,
            "P": p,
            "CVS": s * p,
            "abs_CVS": abs(s * p),
            "n_loops": r["n_loops"],
            "in_paper1": gene in PAPER1_GENES,
        }
    return cvs


def test_sign_accuracy(cvs):
    """Test A: sign(CVS) vs actual CGC role."""
    correct = 0
    total = 0
    wrong = []

    for gene, c in cvs.items():
        role = CGC_ROLE.get(gene)
        if role in ("oncogene", "TSG") and c["CVS"] != 0:
            total += 1
            if (role == "oncogene" and c["CVS"] > 0) or (role == "TSG" and c["CVS"] < 0):
                correct += 1
            else:
                wrong.append((gene, role, c["CVS"], c["S"], c["P"]))

    return correct, total, wrong


def main():
    rds = load_rds()
    cgc = load_cgc()
    cvs = compute_cvs(rds)

    genes = sorted(cvs.keys())

    print("=" * 80)
    print("Cancer Vulnerability Score (CVS = S × P)")
    print("S = |RDS| (structural vulnerability)")
    print("P = Paper1 polarity sign (phenotypic direction)")
    print("=" * 80)

    # Full table
    print(f"\n{'Gene':12s} {'S':>3s} {'P':>3s} {'CVS':>5s} {'|CVS|':>5s} {'Loops':>5s} "
          f"{'Predicted':>10s} {'Actual':>10s} {'CGC':>4s} {'P1':>3s} {'Match':>5s}")
    print("-" * 80)

    for gene in sorted(genes, key=lambda g: -abs(cvs[g]["CVS"])):
        c = cvs[gene]
        predicted = "oncogene" if c["CVS"] > 0 else "TSG" if c["CVS"] < 0 else "neutral"
        actual = CGC_ROLE.get(gene, "")
        cgc_mark = "✓" if gene in cgc else ""
        p1_mark = "✓" if c["in_paper1"] else ""
        match = ""
        if actual in ("oncogene", "TSG") and c["CVS"] != 0:
            match = "✓" if ((actual == "oncogene" and c["CVS"] > 0) or
                           (actual == "TSG" and c["CVS"] < 0)) else "✗"

        print(f"{gene:12s} {c['S']:3d} {c['P']:+3d} {c['CVS']:+5d} {c['abs_CVS']:5d} "
              f"{c['n_loops']:5d} {predicted:>10s} {actual:>10s} {cgc_mark:>4s} "
              f"{p1_mark:>3s} {match:>5s}")

    # ===================================================================
    # TEST A: Sign accuracy
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST A: Sign accuracy — sign(CVS) vs actual CGC role")
    print("=" * 80)

    correct, total, wrong = test_sign_accuracy(cvs)
    print(f"\n  Correct: {correct}/{total} = {100*correct/total:.1f}%")

    if wrong:
        print(f"\n  Mismatches ({len(wrong)}):")
        for g, role, score, s, p in wrong:
            print(f"    {g:12s} CVS={score:+3d} (S={s}, P={p:+d}) predicted={'OG' if score>0 else 'TSG'}, actual={role}")

    # Compare with uncorrected RDS sign
    correct_rds = 0
    total_rds = 0
    for gene in genes:
        role = CGC_ROLE.get(gene)
        r = rds[gene]["RDS"]
        if role in ("oncogene", "TSG") and r != 0:
            total_rds += 1
            if (role == "oncogene" and r > 0) or (role == "TSG" and r < 0):
                correct_rds += 1

    print(f"\n  Comparison:")
    print(f"    Uncorrected RDS sign: {correct_rds}/{total_rds} = {100*correct_rds/total_rds:.1f}%")
    print(f"    CVS sign:             {correct}/{total} = {100*correct/total:.1f}%")

    # ===================================================================
    # TEST B: |CVS| vs COSMIC mutation frequency
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST B: |CVS| vs COSMIC mutation frequency")
    print("=" * 80)

    paired_freq = [(abs(cvs[g]["CVS"]), COSMIC_FREQ[g], g)
                   for g in genes if g in COSMIC_FREQ and cvs[g]["CVS"] != 0]

    if len(paired_freq) >= 5:
        abs_cvs_arr = np.array([p[0] for p in paired_freq])
        freq_arr = np.array([p[1] for p in paired_freq])

        rho, p_rho = stats.spearmanr(abs_cvs_arr, freq_arr)
        print(f"\n  N = {len(paired_freq)} genes with both CVS ≠ 0 and COSMIC frequency")
        print(f"  Spearman ρ(|CVS|, frequency) = {rho:.3f}, p = {p_rho:.4f}")

        # Compare with |RDS| alone
        paired_rds = [(rds[g]["abs_RDS"], COSMIC_FREQ[g], g)
                      for g in genes if g in COSMIC_FREQ and rds[g]["abs_RDS"] > 0]
        if paired_rds:
            rho_rds, p_rds = stats.spearmanr(
                [p[0] for p in paired_rds], [p[1] for p in paired_rds])
            print(f"  Spearman ρ(|RDS|, frequency) = {rho_rds:.3f}, p = {p_rds:.4f}")

    # Also test: does S=0 group have lower frequency?
    freq_s_pos = [COSMIC_FREQ[g] for g in genes if g in COSMIC_FREQ and cvs[g]["S"] > 0]
    freq_s_zero = [COSMIC_FREQ[g] for g in genes if g in COSMIC_FREQ and cvs[g]["S"] == 0]
    if freq_s_pos and freq_s_zero:
        U_f, p_f = stats.mannwhitneyu(freq_s_pos, freq_s_zero, alternative='greater')
        print(f"\n  S > 0 genes: mean freq = {np.mean(freq_s_pos):.2f}%")
        print(f"  S = 0 genes: mean freq = {np.mean(freq_s_zero):.2f}%")
        print(f"  Mann-Whitney (S>0 higher freq): p = {p_f:.4f}")

    # ===================================================================
    # TEST C: |CVS| vs drug target status
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST C: Drug target prediction")
    print("=" * 80)

    drug_cvs = [abs(cvs[g]["CVS"]) for g in genes if g in HAS_DRUG and HAS_DRUG[g]]
    nodrug_cvs = [abs(cvs[g]["CVS"]) for g in genes if g in HAS_DRUG and not HAS_DRUG[g]]

    if drug_cvs and nodrug_cvs:
        print(f"\n  Drug targets (n={len(drug_cvs)}): mean |CVS| = {np.mean(drug_cvs):.2f}")
        print(f"  No drug (n={len(nodrug_cvs)}): mean |CVS| = {np.mean(nodrug_cvs):.2f}")
        U_d, p_d = stats.mannwhitneyu(drug_cvs, nodrug_cvs, alternative='greater')
        print(f"  Mann-Whitney (drug > no drug): p = {p_d:.4f}")

    # ===================================================================
    # TEST D: ROC AUC comparison
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST D: ROC AUC for CGC prediction")
    print("=" * 80)

    # Simple AUC via Mann-Whitney U statistic
    def simple_auc(scores, labels):
        """AUC = P(score_positive > score_negative)."""
        pos = [s for s, l in zip(scores, labels) if l]
        neg = [s for s, l in zip(scores, labels) if not l]
        if not pos or not neg:
            return 0.5
        U, _ = stats.mannwhitneyu(pos, neg, alternative='greater')
        return U / (len(pos) * len(neg))

    is_cgc_arr = [g in cgc for g in genes]

    # (i) |RDS| alone
    rds_scores = [rds[g]["abs_RDS"] for g in genes]
    auc_rds = simple_auc(rds_scores, is_cgc_arr)

    # (ii) Paper1 binary (in any feedback loop = 1)
    p1_scores = [1 if cvs[g]["P"] != 0 else 0 for g in genes]
    auc_p1 = simple_auc(p1_scores, is_cgc_arr)

    # (iii) |CVS|
    cvs_scores = [abs(cvs[g]["CVS"]) for g in genes]
    auc_cvs = simple_auc(cvs_scores, is_cgc_arr)

    # (iv) n_loops (baseline)
    loop_scores = [cvs[g]["n_loops"] for g in genes]
    auc_loops = simple_auc(loop_scores, is_cgc_arr)

    print(f"\n  AUC for CGC membership prediction:")
    print(f"    n_loops (baseline):  {auc_loops:.3f}")
    print(f"    |RDS| alone:         {auc_rds:.3f}")
    print(f"    Paper1 polarity:     {auc_p1:.3f}")
    print(f"    |CVS| (S × P):       {auc_cvs:.3f}")

    if auc_cvs > max(auc_rds, auc_p1):
        print(f"\n  → |CVS| outperforms both components!")
    elif auc_cvs > auc_rds:
        print(f"\n  → |CVS| > |RDS|: polarity adds value")
    elif auc_cvs > auc_p1:
        print(f"\n  → |CVS| > Paper1: structure adds value")
    else:
        print(f"\n  → Components individually better than combined")

    # ===================================================================
    # TEST E: Blind test — 12 NFL-only genes
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST E: Blind predictions — NFL-only genes (not in Paper 1)")
    print("=" * 80)

    new_genes = [g for g in genes if not cvs[g]["in_paper1"]]
    new_with_cvs = [(g, cvs[g]) for g in new_genes if cvs[g]["CVS"] != 0]
    new_neutral = [(g, cvs[g]) for g in new_genes if cvs[g]["CVS"] == 0]

    print(f"\n  NFL-only genes: {len(new_genes)} total")
    print(f"  With CVS ≠ 0: {len(new_with_cvs)} (testable predictions)")
    print(f"  CVS = 0: {len(new_neutral)} (no prediction)")

    if new_with_cvs:
        print(f"\n  {'Gene':12s} {'CVS':>5s} {'Predicted':>10s} {'CGC Role':>10s} {'Match':>5s}")
        print("  " + "-" * 50)
        correct_new = 0
        total_new = 0
        for g, c in sorted(new_with_cvs, key=lambda x: -abs(x[1]["CVS"])):
            predicted = "oncogene" if c["CVS"] > 0 else "TSG"
            actual = CGC_ROLE.get(g, "not_CGC")
            match = ""
            if actual in ("oncogene", "TSG"):
                total_new += 1
                if (actual == "oncogene" and c["CVS"] > 0) or (actual == "TSG" and c["CVS"] < 0):
                    correct_new += 1
                    match = "✓"
                else:
                    match = "✗"
            elif actual == "both":
                match = "~"
            print(f"  {g:12s} {c['CVS']:+5d} {predicted:>10s} {actual:>10s} {match:>5s}")

        if total_new > 0:
            print(f"\n  Blind accuracy: {correct_new}/{total_new} = {100*correct_new/total_new:.1f}%")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  Test A (sign accuracy):   {correct}/{total} = {100*correct/total:.1f}%")
    if len(paired_freq) >= 5:
        print(f"  Test B (mutation freq):   ρ = {rho:.3f}, p = {p_rho:.4f}")
    print(f"  Test D (AUC CGC):         |CVS| = {auc_cvs:.3f} vs |RDS| = {auc_rds:.3f} vs loops = {auc_loops:.3f}")

    # Save
    with open(DATA_DIR / "cvs_scores.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "S_abs_RDS", "P_paper1_sign", "CVS", "abs_CVS",
                     "n_loops", "predicted_role", "actual_role", "in_paper1"])
        for gene in sorted(genes, key=lambda g: -abs(cvs[g]["CVS"])):
            c = cvs[gene]
            pred = "oncogene" if c["CVS"] > 0 else "TSG" if c["CVS"] < 0 else "neutral"
            actual = CGC_ROLE.get(gene, "")
            w.writerow([gene, c["S"], c["P"], c["CVS"], c["abs_CVS"],
                        c["n_loops"], pred, actual, c["in_paper1"]])

    results = {
        "test_A_accuracy": correct / total,
        "test_A_n": total,
        "test_D_auc_cvs": auc_cvs,
        "test_D_auc_rds": auc_rds,
        "test_D_auc_p1": auc_p1,
        "test_D_auc_loops": auc_loops,
    }
    with open(DATA_DIR / "cvs_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved: cvs_scores.csv, cvs_results.json")
    print("Done!")


if __name__ == "__main__":
    main()
