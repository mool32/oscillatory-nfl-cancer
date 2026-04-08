#!/usr/bin/env python3
"""
Polarity-Corrected RDS.

Corrected_RDS(G) = Σ_l Polarity(l) × Sign(G,l) × R(G,l)

Where Polarity(l) = +1 for growth-promoting pathways,
                    −1 for growth-inhibitory pathways.

Growth-inhibitory: p53, TGF-β, BMP (tumor suppressive signaling)
  → rise arm = growth INHIBITION → TSG
  → recovery arm = growth PROMOTION → oncogene

Growth-promoting: everything else (NF-κB, ERK, Wnt, mTOR, JAK-STAT, etc.)
  → rise arm = growth PROMOTION → oncogene
  → recovery arm = growth INHIBITION → TSG
"""

import csv
import json
import numpy as np
from scipy import stats
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# Load RDS details from step10
def load_rds_and_motifs():
    """Load unique motifs and rebuild gene-loop mapping."""
    motifs = []
    with open(DATA_DIR / "unique_nfl_motifs.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes = row["genes_hgnc"].split("|")
            motifs.append({
                "genes": genes,
                "gene_set": frozenset(genes),
                "length": int(row["length"]),
            })
    return motifs

# Position sign: +1 = rise, -1 = recovery
POSITION_SIGN = {
    "RELA": +1, "NFKB1": +1, "RELB": +1,
    "NFKBIA": -1, "TNFAIP3": -1, "TRAF6": +1,
    "TP53": -1, "MDM2": +1, "MDM4": +1, "ATM": -1, "CHEK2": -1,
    "MAPK1": +1, "MAPK3": +1, "MAPK14": +1,
    "DUSP1": -1,
    "SOCS1": -1, "SOCS3": -1, "SOCS4": -1, "SOCS5": -1, "SOCS6": -1, "SOCS7": -1,
    "CISH": -1,
    "CTNNB1": +1, "CDH1": -1, "AXIN1": -1, "AXIN2": -1,
    "NOTCH1": +1, "FBXW7": -1,
    "CLOCK": +1, "ARNTL": +1, "PER2": -1, "CRY1": -1,
    "YAP1": +1, "LATS1": -1, "LATS2": -1, "STK3": -1,
    "MTOR": +1, "TSC1": -1, "TSC2": -1, "IRS1": +1, "IRS2": +1,
    "PIK3CA": +1, "PIK3R1": -1, "PIK3R3": +1, "PTEN": -1,
    "AKT1": +1, "AKT2": +1, "AKT3": +1, "PTK2": +1,
    "JAK1": +1, "JAK2": +1, "TYK2": +1,
    "STAT1": +1, "STAT2": +1, "STAT3": +1, "STAT5A": +1, "STAT5B": +1,
    "LEPR": +1, "CSF1R": +1, "PRLR": +1,
    "TGFBR1": +1, "SMAD2": +1, "SMAD3": +1, "SMAD4": +1,
    "SMAD1": +1, "SMAD5": +1, "SMAD7": -1,
    "BMPR1A": +1, "BMPR1B": +1, "BMPR2": +1,
    "NFE2L2": +1, "KEAP1": -1,
    "GLI1": +1, "GLI2": +1, "GLI3": -1, "SUFU": -1,
    "ERN1": +1, "HSPA5": -1,
    "PLCG1": +1, "ATP2A2": -1,
    "CDK2": +1, "CDKN1A": -1, "CDKN1B": -1, "CDKN2A": -1,
    "CCNE1": +1, "E2F1": +1, "E2F2": +1, "E2F3": +1, "TFDP1": +1,
    "RB1": -1, "RBL1": -1, "RBL2": -1,
    "FOXO1": -1, "FOXO3": -1, "FOXO4": -1, "FOXO6": -1,
    "BCL2": +1, "PMAIP1": -1, "SIVA1": -1, "BID": -1, "BBC3": -1, "TP73": -1,
    "NFATC1": +1, "RCAN1": -1,
    "IL6": +1, "NFKBIZ": -1, "MYC": +1,
    "HIF1A": +1, "VHL": -1, "CSNK1E": -1,
}

# Pathway polarity: which loops are growth-inhibitory?
# Growth-inhibitory pathways: activation of the pathway SUPPRESSES growth
# → rise arm genes are TSGs, recovery arm genes are oncogenes
GROWTH_INHIBITORY_GENES = {
    # p53 pathway: p53 activation = growth arrest/apoptosis
    "TP53", "ATM", "CHEK2",  # rise arm (activators of growth inhibition)
    "MDM2", "MDM4",           # recovery arm (turn off growth inhibition)
    # p53-apoptosis: pro-apoptotic arm = growth inhibitory
    "PMAIP1", "SIVA1", "BID", "BBC3", "TP73",  # rise (pro-apoptotic)
    "BCL2",                    # recovery (anti-apoptotic)
    # TGF-β/SMAD: TGF-β = growth inhibitory in epithelial cells
    "TGFBR1", "SMAD2", "SMAD3", "SMAD4", "SMAD1", "SMAD5",  # rise
    "SMAD7",                   # recovery
    # BMP (in context of SMAD-mediated growth arrest)
    "BMPR1A", "BMPR1B", "BMPR2",  # rise
    # Rb/E2F: Rb = growth inhibitory (blocks E2F)
    # But E2F/CDK2/CCNE1 = growth promoting. Rb is recovery arm OF growth-promoting cycle.
    # So Rb loop is growth-PROMOTING (cell cycle), not inhibitory.
    # FOXO: growth-inhibitory (FOXO activates cell cycle arrest genes)
    "FOXO1", "FOXO3", "FOXO4", "FOXO6",  # rise (growth inhibition)
    # CDK2 in context of FOXO loop: CDK2 phosphorylates FOXO → inactivates → recovery
}

# For each gene, determine if it's in a growth-inhibitory context
# Polarity = -1 if gene is in growth-inhibitory pathway
def get_loop_polarity(loop_genes):
    """
    Determine if a loop is growth-inhibitory.
    If ANY gene in the loop is in GROWTH_INHIBITORY_GENES → polarity = -1.
    Otherwise → polarity = +1.
    """
    for g in loop_genes:
        if g in GROWTH_INHIBITORY_GENES:
            return -1
    return +1


def load_irreplaceability():
    """Load irreplaceability from step10 RDS results."""
    with open(DATA_DIR / "rds_full.json") as f:
        data = json.load(f)
    return data["rds"]


# Full CGC
def load_cgc():
    genes = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g: genes.add(g)
    return genes

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


def compute_corrected_rds(motifs, rds_old):
    """Compute polarity-corrected RDS."""
    gene_loops = defaultdict(list)
    for m in motifs:
        for g in m["genes"]:
            gene_loops[g].append(m)

    rds_new = {}

    for gene in sorted(gene_loops.keys()):
        loops = gene_loops[gene]
        sign = POSITION_SIGN.get(gene, 0)
        n_irreplaceable = rds_old.get(gene, {}).get("n_irreplaceable", 0)

        # For each loop, compute polarity-corrected contribution
        corrected_rds = 0
        uncorrected_rds = 0
        loop_details = []

        for i, loop in enumerate(loops):
            # Check irreplaceability (reuse from step10)
            # Simple heuristic: if n_irreplaceable > 0, distribute
            # For proper computation, we'd need per-loop irreplaceability
            # For now: use the ratio from step10
            n_loops = len(loops)
            # Approximate: if gene is irreplaceable in k of n loops,
            # assume the first k loops are the irreplaceable ones
            # (This is a simplification — proper would need per-loop check)
            is_irreplaceable = i < n_irreplaceable

            r = 1 if is_irreplaceable else 0
            polarity = get_loop_polarity(loop["genes"])

            uncorrected = sign * r
            corrected = polarity * sign * r

            uncorrected_rds += uncorrected
            corrected_rds += corrected

            loop_details.append({
                "genes": loop["genes"],
                "polarity": polarity,
                "irreplaceable": is_irreplaceable,
                "uncorrected": uncorrected,
                "corrected": corrected,
            })

        # Predicted role from corrected RDS
        if corrected_rds > 0:
            predicted = "oncogene"
        elif corrected_rds < 0:
            predicted = "TSG"
        elif corrected_rds == 0 and n_irreplaceable > 0:
            predicted = "both"
        else:
            predicted = "neutral"

        actual = CGC_ROLE.get(gene, "")

        rds_new[gene] = {
            "RDS_uncorrected": uncorrected_rds,
            "RDS_corrected": corrected_rds,
            "abs_corrected": abs(corrected_rds),
            "n_loops": len(loops),
            "n_irreplaceable": n_irreplaceable,
            "sign": sign,
            "predicted_role": predicted,
            "actual_role": actual,
            "loops": loop_details,
        }

    return rds_new


def test_all_predictions(rds, cgc):
    """Test all predictions with corrected RDS."""
    genes = sorted(rds.keys())
    corrected = np.array([rds[g]["RDS_corrected"] for g in genes])
    uncorrected = np.array([rds[g]["RDS_uncorrected"] for g in genes])
    abs_c = np.array([abs(c) for c in corrected])
    is_cgc = np.array([g in cgc for g in genes])

    print(f"\n{'='*70}")
    print("POLARITY-CORRECTED RDS — ALL PREDICTIONS")
    print("=" * 70)

    # ===================================================================
    # Prediction 1: |corrected RDS| > 0 → CGC
    # ===================================================================
    print(f"\n--- Prediction 1: |RDS_corrected| > 0 → CGC ---")
    pos = abs_c > 0
    zero = abs_c == 0
    cgc_pos = sum(is_cgc[pos])
    cgc_zero = sum(is_cgc[zero])
    n_pos = sum(pos)
    n_zero = sum(zero)

    print(f"  |RDS| > 0: {cgc_pos}/{n_pos} CGC ({100*cgc_pos/n_pos:.1f}%)")
    print(f"  |RDS| = 0: {cgc_zero}/{n_zero} CGC ({100*cgc_zero/n_zero:.1f}%)")

    table1 = [[cgc_pos, n_pos - cgc_pos], [cgc_zero, n_zero - cgc_zero]]
    OR1, p1 = stats.fisher_exact(table1)
    print(f"  Fisher: OR = {OR1:.2f}, p = {p1:.4f}")
    rpb, p_rpb = stats.pointbiserialr(is_cgc.astype(int), abs_c)
    print(f"  Point-biserial r = {rpb:.3f}, p = {p_rpb:.4f}")

    # ===================================================================
    # Prediction 2: Sign accuracy (THE MAIN TEST)
    # ===================================================================
    print(f"\n{'='*70}")
    print("--- Prediction 2: Sign accuracy (CORRECTED vs UNCORRECTED) ---")
    print("=" * 70)

    # Uncorrected
    correct_unc = 0
    total_unc = 0
    for g in genes:
        role = CGC_ROLE.get(g, "")
        if role in ("oncogene", "TSG"):
            total_unc += 1
            u = rds[g]["RDS_uncorrected"]
            if (role == "oncogene" and u > 0) or (role == "TSG" and u < 0):
                correct_unc += 1

    # Corrected
    correct_cor = 0
    total_cor = 0
    mismatches = []
    for g in genes:
        role = CGC_ROLE.get(g, "")
        if role in ("oncogene", "TSG"):
            total_cor += 1
            c = rds[g]["RDS_corrected"]
            if (role == "oncogene" and c > 0) or (role == "TSG" and c < 0):
                correct_cor += 1
            elif c != 0:  # wrong sign (exclude RDS=0 from wrong count)
                mismatches.append((g, role, c, rds[g]["RDS_uncorrected"]))

    print(f"\n  UNCORRECTED: {correct_unc}/{total_unc} = {100*correct_unc/total_unc:.1f}% accuracy")
    print(f"  CORRECTED:   {correct_cor}/{total_cor} = {100*correct_cor/total_cor:.1f}% accuracy")
    print(f"  Improvement: +{correct_cor - correct_unc} genes corrected")

    # Detailed: which genes flipped?
    print(f"\n  Genes where polarity correction changed prediction:")
    for g in genes:
        u = rds[g]["RDS_uncorrected"]
        c = rds[g]["RDS_corrected"]
        if u != c and (u != 0 or c != 0):
            role = CGC_ROLE.get(g, "n/a")
            u_pred = "OG" if u > 0 else "TSG" if u < 0 else "—"
            c_pred = "OG" if c > 0 else "TSG" if c < 0 else "—"
            correct_mark = "✓" if ((role == "oncogene" and c > 0) or (role == "TSG" and c < 0)) else "✗" if c != 0 and role in ("oncogene", "TSG") else "?"
            print(f"    {g:12s} uncorr={u:+3d}({u_pred:3s}) → corr={c:+3d}({c_pred:3s})  actual={role:10s} {correct_mark}")

    # Remaining mismatches
    if mismatches:
        print(f"\n  Remaining mismatches after correction:")
        for g, role, c, u in mismatches:
            print(f"    {g:12s} RDS_corr={c:+3d} predicted={'OG' if c>0 else 'TSG'}, actual={role}")

    # ===================================================================
    # Mann-Whitney: oncogenes vs TSGs on corrected RDS
    # ===================================================================
    og_rds = [rds[g]["RDS_corrected"] for g in genes if CGC_ROLE.get(g) == "oncogene"]
    tsg_rds = [rds[g]["RDS_corrected"] for g in genes if CGC_ROLE.get(g) == "TSG"]
    both_rds = [rds[g]["RDS_corrected"] for g in genes if CGC_ROLE.get(g) == "both"]

    print(f"\n--- Distribution ---")
    print(f"  Oncogenes (n={len(og_rds)}): mean = {np.mean(og_rds):+.2f}, median = {np.median(og_rds):+.1f}")
    print(f"  TSGs (n={len(tsg_rds)}): mean = {np.mean(tsg_rds):+.2f}, median = {np.median(tsg_rds):+.1f}")
    if both_rds:
        print(f"  Both (n={len(both_rds)}): mean = {np.mean(both_rds):+.2f}")

    U, p_u = stats.mannwhitneyu(og_rds, tsg_rds, alternative='greater')
    print(f"  Mann-Whitney (OG > TSG): U = {U:.0f}, p = {p_u:.6f}")

    # ===================================================================
    # Comparison: corrected vs uncorrected separation
    # ===================================================================
    print(f"\n--- Effect size comparison ---")
    og_unc = [rds[g]["RDS_uncorrected"] for g in genes if CGC_ROLE.get(g) == "oncogene"]
    tsg_unc = [rds[g]["RDS_uncorrected"] for g in genes if CGC_ROLE.get(g) == "TSG"]

    d_unc = (np.mean(og_unc) - np.mean(tsg_unc))
    d_cor = (np.mean(og_rds) - np.mean(tsg_rds))
    print(f"  Uncorrected: OG−TSG mean diff = {d_unc:.2f}")
    print(f"  Corrected:   OG−TSG mean diff = {d_cor:.2f}")
    print(f"  Improvement ratio: {d_cor/d_unc:.2f}× wider separation" if d_unc != 0 else "")

    # ===================================================================
    # Prediction 5: |RDS| ≈ 0 for "both" genes
    # ===================================================================
    if both_rds:
        print(f"\n--- Prediction 5: 'Both' genes ---")
        print(f"  Both: |RDS| = {[abs(r) for r in both_rds]}")
        print(f"  OG mean |RDS| = {np.mean([abs(r) for r in og_rds]):.2f}")
        print(f"  TSG mean |RDS| = {np.mean([abs(r) for r in tsg_rds]):.2f}")

    return {
        "accuracy_uncorrected": correct_unc / total_unc,
        "accuracy_corrected": correct_cor / total_cor,
        "n_genes_tested": total_cor,
        "mannwhitney_p": float(p_u),
        "effect_size_uncorrected": d_unc,
        "effect_size_corrected": d_cor,
    }


def main():
    print("=" * 70)
    print("Polarity-Corrected Reversibility Dependency Score")
    print("=" * 70)

    motifs = load_rds_and_motifs()
    rds_old = load_irreplaceability()
    cgc = load_cgc()

    rds = compute_corrected_rds(motifs, rds_old)

    # Print full table
    print(f"\n{'Gene':15s} {'Uncorr':>7s} {'Corr':>6s} {'|C|':>4s} {'Loops':>6s} {'Irr':>4s} "
          f"{'Predicted':>10s} {'Actual':>10s} {'CGC':>4s} {'Match':>6s}")
    print("-" * 85)

    for gene in sorted(rds.keys(), key=lambda g: -abs(rds[g]["RDS_corrected"])):
        r = rds[gene]
        cgc_mark = "✓" if gene in cgc else ""
        match = ""
        if r["actual_role"] in ("oncogene", "TSG") and r["RDS_corrected"] != 0:
            if (r["actual_role"] == "oncogene" and r["RDS_corrected"] > 0) or \
               (r["actual_role"] == "TSG" and r["RDS_corrected"] < 0):
                match = "✓"
            else:
                match = "✗"
        elif r["actual_role"] == "both" and r["RDS_corrected"] == 0:
            match = "✓"

        print(f"{gene:15s} {r['RDS_uncorrected']:+7d} {r['RDS_corrected']:+6d} {r['abs_corrected']:4d} "
              f"{r['n_loops']:6d} {r['n_irreplaceable']:4d} "
              f"{r['predicted_role']:>10s} {r['actual_role']:>10s} {cgc_mark:>4s} {match:>6s}")

    # Test predictions
    results = test_all_predictions(rds, cgc)

    # Save
    with open(DATA_DIR / "rds_corrected.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "RDS_uncorrected", "RDS_corrected", "abs_corrected",
                     "n_loops", "n_irreplaceable", "predicted_role", "actual_role"])
        for gene in sorted(rds.keys(), key=lambda g: -abs(rds[g]["RDS_corrected"])):
            r = rds[gene]
            w.writerow([gene, r["RDS_uncorrected"], r["RDS_corrected"],
                        r["abs_corrected"], r["n_loops"], r["n_irreplaceable"],
                        r["predicted_role"], r["actual_role"]])

    with open(DATA_DIR / "rds_corrected_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved: rds_corrected.csv, rds_corrected_results.json")
    print("Done!")


if __name__ == "__main__":
    main()
