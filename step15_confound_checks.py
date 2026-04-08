#!/usr/bin/env python3
"""
Confound checks for module-receptor coupling and cognitive load results.

Check 1: Gene list overlap — are module genes and receptor genes the same?
Check 2: Partial correlation — module↔receptor controlling for total expression
Check 3: Partial correlation — cognitive load↔receptor diversity controlling for total genes
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# Import gene lists from step14
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

MODULE_RECEPTORS = {
    "NF-κB": {"TNFRSF1A", "TNFRSF1B", "TLR4", "TLR2", "TLR7", "TLR9",
              "IL1R1", "IRAK1", "NFKB1", "MYD88", "RIPK1"},
    "ERK/MAPK": {"EGFR", "FGFR1", "FGFR2", "PDGFRA", "PDGFRB", "KIT",
                 "MET", "ERBB2", "ERBB3", "IGF1R", "KRAS", "GRB2", "SOS1"},
    "JAK-STAT": {"IL6R", "IL6ST", "IL2RA", "IL2RB", "IL2RG", "IL7R",
                 "IL10RA", "IL12RB1", "IFNAR1", "IFNAR2", "IFNGR1",
                 "IFNGR2", "EPOR", "CSF2RA", "LEPR"},
    "Wnt": {"FZD1", "FZD2", "FZD4", "FZD5", "FZD7", "LRP5", "LRP6", "ROR2", "RNF43"},
    "Notch": {"NOTCH1", "NOTCH2", "NOTCH3", "DLL1", "DLL4", "JAG1", "JAG2"},
    "TGF-β": {"TGFBR1", "TGFBR2", "ACVR1", "ACVR1B", "ACVR2A",
              "BMPR1A", "BMPR1B", "BMPR2"},
    "Hippo": {"CDH1", "CDH2", "ITGA6", "ITGB1", "ITGB4", "NF2",
              "FAT4", "DCHS1", "AMOT", "GPR87"},
    "mTOR": {"SLC1A5", "SLC7A5", "SLC3A2", "INSR", "IGF1R",
             "SESN2", "CASTOR1", "FLCN"},
    "Calcium": {"ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "TRPC1",
                "TRPV4", "ORAI1", "STIM1", "P2RX7"},
    "p53": {"ATM", "ATR", "CHEK1", "CHEK2", "PARP1", "H2AFX"},
    "Cell Cycle": {"CCND1", "CCND2", "CCND3", "CDK4", "CDK6", "MYC"},
    "Circadian": {"OPN4", "MTNR1A", "MTNR1B", "PER1", "NR1D1", "ARNTL"},
    "NRF2": {"NFE2L2", "NQO1", "GCLC", "HMOX1", "TXNRD1", "GPX2"},
    "PI3K/PTEN": {"EGFR", "ERBB2", "IGF1R", "INSR", "PDGFRA", "ITGB1"},
}


def partial_corr(x, y, z):
    """Partial Pearson correlation of x,y controlling for z."""
    # Residuals of x ~ z and y ~ z
    slope_xz, intercept_xz, _, _, _ = stats.linregress(z, x)
    slope_yz, intercept_yz, _, _, _ = stats.linregress(z, y)
    res_x = x - (slope_xz * z + intercept_xz)
    res_y = y - (slope_yz * z + intercept_yz)
    r, p = stats.pearsonr(res_x, res_y)
    return r, p


def load_expression():
    print("Loading HPA expression data...")
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types


def load_cognitive_load():
    cog = {}
    with open(DATA_DIR / "cognitive_load.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cog[row["cell_type"]] = {
                "load": int(row["cognitive_load"]),
                "continuous": float(row["continuous_load"]),
            }
    return cog


def main():
    expr, cell_types = load_expression()
    cog = load_cognitive_load()

    # ===================================================================
    # CHECK 1: Gene list overlap
    # ===================================================================
    print("=" * 80)
    print("CHECK 1: Gene overlap between module core genes and receptor genes")
    print("=" * 80)

    any_overlap = False
    for mod in sorted(MODULES_CORE.keys()):
        if mod not in MODULE_RECEPTORS:
            continue
        core = MODULES_CORE[mod]
        recs_dict = MODULE_RECEPTORS[mod]
        recs = set(recs_dict.keys()) if isinstance(recs_dict, dict) else recs_dict
        overlap = core & recs
        if overlap:
            any_overlap = True
            print(f"\n  ⚠️  {mod}: OVERLAP detected!")
            print(f"      Core genes: {sorted(core)}")
            print(f"      Receptor genes: {sorted(recs)}")
            print(f"      Shared: {sorted(overlap)}")

    if not any_overlap:
        print("\n  ✓ No gene overlap between module core and receptor gene lists.")

    # ===================================================================
    # CHECK 1b: Recompute NF-κB WITHOUT overlapping genes
    # ===================================================================
    print(f"\n{'='*80}")
    print("CHECK 1b: Recompute with overlapping genes removed")
    print("=" * 80)

    for mod in sorted(MODULES_CORE.keys()):
        if mod not in MODULE_RECEPTORS:
            continue
        core = MODULES_CORE[mod]
        recs_raw = MODULE_RECEPTORS[mod]
        recs = set(recs_raw.keys()) if isinstance(recs_raw, dict) else recs_raw
        overlap = core & recs

        if not overlap:
            continue

        # Compute with and without overlap
        core_clean = core - overlap
        recs_clean = recs - overlap

        core_found = [g for g in core if g in expr]
        core_clean_found = [g for g in core_clean if g in expr]
        recs_found = [g for g in recs if g in expr]
        recs_clean_found = [g for g in recs_clean if g in expr]

        # Original
        mod_scores = np.array([np.mean([expr[g].get(ct, 0) for g in core_found])
                               for ct in cell_types])
        rec_scores = np.array([np.mean([expr[g].get(ct, 0) for g in recs_found])
                               for ct in cell_types])
        r_orig, p_orig = stats.pearsonr(mod_scores, rec_scores)

        # Clean (no overlap)
        if core_clean_found and recs_clean_found:
            mod_clean = np.array([np.mean([expr[g].get(ct, 0) for g in core_clean_found])
                                  for ct in cell_types])
            rec_clean = np.array([np.mean([expr[g].get(ct, 0) for g in recs_clean_found])
                                  for ct in cell_types])
            r_clean, p_clean = stats.pearsonr(mod_clean, rec_clean)
            print(f"\n  {mod}:")
            print(f"    Overlapping genes: {sorted(overlap)}")
            print(f"    Original:        r = {r_orig:.3f}, p = {p_orig:.2e}")
            print(f"    Without overlap: r = {r_clean:.3f}, p = {p_clean:.2e}")
            print(f"    Core genes used: {sorted(core_clean_found)}")
            print(f"    Receptor genes used: {sorted(recs_clean_found)}")

    # ===================================================================
    # CHECK 2: Partial correlation controlling for total expression
    # ===================================================================
    print(f"\n{'='*80}")
    print("CHECK 2: Module↔receptor PARTIAL correlation (controlling for total expression)")
    print("=" * 80)

    # Compute total expression per cell type (mean across ALL genes)
    print("\nComputing total expression per cell type...")
    all_genes = sorted(expr.keys())
    total_expr = np.array([np.mean([expr[g].get(ct, 0) for g in all_genes[:2000]])  # sample 2000 for speed
                           for ct in cell_types])

    # Also compute total detected genes per cell type
    n_detected = np.array([sum(1 for g in all_genes[:5000] if expr[g].get(ct, 0) > 1.0)
                           for ct in cell_types])

    print(f"  Total expression range: {total_expr.min():.1f} - {total_expr.max():.1f}")
    print(f"  Detected genes range: {n_detected.min()} - {n_detected.max()}")

    print(f"\n  {'Module':15s} {'Raw r':>8s} {'Partial r':>10s} {'Part. p':>10s} {'Δr':>8s} {'Survives?':>10s}")
    print("  " + "-" * 65)

    partial_results = []
    for mod in sorted(MODULES_CORE.keys()):
        if mod not in MODULE_RECEPTORS:
            continue
        core = MODULES_CORE[mod]
        recs_raw = MODULE_RECEPTORS[mod]
        recs = set(recs_raw.keys()) if isinstance(recs_raw, dict) else recs_raw

        # Use clean (non-overlapping) gene lists
        core_clean = [g for g in (core - recs) if g in expr] or [g for g in core if g in expr]
        recs_clean = [g for g in (recs - core) if g in expr] or [g for g in recs if g in expr]

        mod_scores = np.array([np.mean([expr[g].get(ct, 0) for g in core_clean])
                               for ct in cell_types])
        rec_scores = np.array([np.mean([expr[g].get(ct, 0) for g in recs_clean])
                               for ct in cell_types])

        r_raw, _ = stats.pearsonr(mod_scores, rec_scores)
        r_part, p_part = partial_corr(mod_scores, rec_scores, total_expr)

        survives = "✓" if (p_part < 0.05 and r_part > 0.1) else "✗"
        delta = r_part - r_raw

        partial_results.append(r_part)
        print(f"  {mod:15s} {r_raw:+8.3f} {r_part:+10.3f} {p_part:10.2e} {delta:+8.3f} {survives:>10s}")

    mean_partial = np.mean(partial_results)
    t_part, p_tpart = stats.ttest_1samp(partial_results, 0)
    print(f"\n  Mean partial r: {mean_partial:.3f}")
    print(f"  One-sample t-test (mean partial r > 0): t = {t_part:.2f}, p = {p_tpart:.4f}")

    # ===================================================================
    # CHECK 3: Cognitive load ↔ receptor diversity, controlling for total genes
    # ===================================================================
    print(f"\n{'='*80}")
    print("CHECK 3: Cognitive load ↔ receptor diversity (controlling for total detected genes)")
    print("=" * 80)

    # Receptor diversity from step14
    from step14_receptor_module_coupling import ALL_RECEPTORS
    all_recs_found = [g for g in ALL_RECEPTORS if g in expr]
    diversity = np.array([sum(1 for g in all_recs_found if expr[g].get(ct, 0) > 1.0)
                          for ct in cell_types])

    valid = [j for j, ct in enumerate(cell_types) if ct in cog]
    load_arr = np.array([cog[cell_types[j]]["load"] for j in valid]).astype(float)
    div_arr = diversity[valid].astype(float)
    total_arr = n_detected[valid].astype(float)

    r_raw, p_raw = stats.pearsonr(load_arr, div_arr)
    r_part, p_part = partial_corr(load_arr, div_arr, total_arr)

    print(f"\n  Raw: r(load, diversity) = {r_raw:.3f}, p = {p_raw:.2e}")
    print(f"  Partial: r(load, diversity | total_genes) = {r_part:.3f}, p = {p_part:.2e}")

    # Also: receptor diversity / total genes = FRACTION of receptors
    frac = div_arr / (total_arr + 1)
    rho_frac, p_frac = stats.spearmanr(load_arr, frac)
    print(f"  Fraction: ρ(load, receptor_frac) = {rho_frac:.3f}, p = {p_frac:.2e}")

    if r_part > 0.1 and p_part < 0.05:
        print(f"\n  ✓ Cognitive load ↔ receptor diversity SURVIVES after controlling for total expression")
    else:
        print(f"\n  ✗ Relationship driven by total expression level (confound)")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("CONFOUND CHECK SUMMARY")
    print("=" * 80)

    if any_overlap:
        print("\n  ⚠️  Gene overlap found — corrected in Check 1b")
    else:
        print("\n  ✓ No gene list overlap")

    print(f"  Module-receptor partial r (mean): {mean_partial:.3f}, t-test p = {p_tpart:.4f}")
    if p_tpart < 0.05:
        print(f"  ✓ Module-receptor coupling SURVIVES after controlling for total expression")
    else:
        print(f"  ✗ Module-receptor coupling is confounded by total expression")

    print(f"  Cognitive load ↔ diversity partial r: {r_part:.3f}, p = {p_part:.2e}")

    # Save
    results = {
        "gene_overlap_found": any_overlap,
        "partial_r_mean": round(float(mean_partial), 3),
        "partial_r_ttest_p": float(p_tpart),
        "load_diversity_raw_r": round(float(r_raw), 3),
        "load_diversity_partial_r": round(float(r_part), 3),
        "load_diversity_partial_p": float(p_part),
        "load_diversity_fraction_rho": round(float(rho_frac), 3),
    }
    with open(DATA_DIR / "confound_checks.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to confound_checks.json")


if __name__ == "__main__":
    main()
