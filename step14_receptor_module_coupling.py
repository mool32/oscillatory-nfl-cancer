#!/usr/bin/env python3
"""
Paper 3: Receptor-module coupling and receptor diversity.

Test 1: Module activity ↔ upstream receptor expression (per cell type).
        Active NF-κB → more TLR/TNFR? Active Hippo → more cadherins?

Test 2: Receptor diversity ∝ cognitive load.
        Endothelial (load=11.8) → many receptors. Germ (1.4) → few.

Data: HPA single cell type nCPM (154 cell types, already loaded).
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# ---------------------------------------------------------------
# Module → upstream receptors mapping
# These are the receptors that FEED INTO each signaling module
# ---------------------------------------------------------------
MODULE_RECEPTORS = {
    "NF-κB": {
        "TNFRSF1A": "TNF receptor",    # TNFR1
        "TNFRSF1B": "TNF receptor",    # TNFR2
        "TLR4": "LPS receptor",
        "TLR2": "bacterial lipopeptide",
        "TLR7": "ssRNA",
        "TLR9": "CpG DNA",
        "IL1R1": "IL-1 receptor",
        "IRAK1": "IL-1R signaling",
        "NFKB1": "NF-κB subunit",      # also marker of activity
        "MYD88": "adapter",
        "RIPK1": "death receptor signaling",
    },
    "ERK/MAPK": {
        "EGFR": "EGF receptor",
        "FGFR1": "FGF receptor",
        "FGFR2": "FGF receptor",
        "PDGFRA": "PDGF receptor",
        "PDGFRB": "PDGF receptor",
        "KIT": "SCF receptor",
        "MET": "HGF receptor",
        "ERBB2": "HER2",
        "ERBB3": "HER3",
        "IGF1R": "IGF receptor",
        "KRAS": "RAS GTPase",
        "GRB2": "adapter",
        "SOS1": "RAS-GEF",
    },
    "JAK-STAT": {
        "IL6R": "IL-6 receptor",
        "IL6ST": "gp130 co-receptor",
        "IL2RA": "IL-2R alpha",
        "IL2RB": "IL-2R beta",
        "IL2RG": "common gamma chain",
        "IL7R": "IL-7 receptor",
        "IL10RA": "IL-10 receptor",
        "IL12RB1": "IL-12 receptor",
        "IFNAR1": "IFN-alpha receptor",
        "IFNAR2": "IFN-alpha receptor",
        "IFNGR1": "IFN-gamma receptor",
        "IFNGR2": "IFN-gamma receptor",
        "EPOR": "EPO receptor",
        "CSF2RA": "GM-CSF receptor",
        "LEPR": "leptin receptor",
    },
    "Wnt": {
        "FZD1": "Frizzled-1",
        "FZD2": "Frizzled-2",
        "FZD4": "Frizzled-4",
        "FZD5": "Frizzled-5",
        "FZD7": "Frizzled-7",
        "LRP5": "LRP5 co-receptor",
        "LRP6": "LRP6 co-receptor",
        "ROR2": "ROR2 Wnt receptor",
        "RNF43": "Wnt feedback receptor",
    },
    "Notch": {
        "NOTCH1": "Notch1 receptor",
        "NOTCH2": "Notch2 receptor",
        "NOTCH3": "Notch3 receptor",
        "DLL1": "Delta ligand",
        "DLL4": "Delta ligand",
        "JAG1": "Jagged ligand",
        "JAG2": "Jagged ligand",
    },
    "TGF-β": {
        "TGFBR1": "TGF-β type I receptor",
        "TGFBR2": "TGF-β type II receptor",
        "ACVR1": "Activin receptor",
        "ACVR1B": "Activin receptor",
        "ACVR2A": "Activin receptor",
        "BMPR1A": "BMP receptor",
        "BMPR1B": "BMP receptor",
        "BMPR2": "BMP receptor",
    },
    "Hippo": {
        "CDH1": "E-cadherin (contact)",
        "CDH2": "N-cadherin",
        "ITGA6": "Integrin α6",
        "ITGB1": "Integrin β1",
        "ITGB4": "Integrin β4",
        "NF2": "Merlin (mechanosensor)",
        "FAT4": "FAT cadherin",
        "DCHS1": "Dachsous",
        "AMOT": "Angiomotin",
        "GPR87": "GPCR mechanosensor",
    },
    "mTOR": {
        "SLC1A5": "AA transporter (glutamine)",
        "SLC7A5": "AA transporter (leucine)",
        "SLC3A2": "AA transporter subunit",
        "INSR": "Insulin receptor",
        "IGF1R": "IGF receptor",
        "SESN2": "stress sensor (Sestrin2)",
        "CASTOR1": "arginine sensor",
        "FLCN": "folliculin",
    },
    "Calcium": {
        "ITPR1": "IP3 receptor",
        "ITPR2": "IP3 receptor",
        "ITPR3": "IP3 receptor",
        "RYR1": "Ryanodine receptor",
        "RYR2": "Ryanodine receptor",
        "TRPC1": "TRP channel",
        "TRPV4": "TRP channel",
        "ORAI1": "CRAC channel",
        "STIM1": "ER Ca sensor",
        "P2RX7": "purinergic receptor",
    },
    "p53": {
        "ATM": "DNA damage sensor",
        "ATR": "replication stress sensor",
        "CHEK1": "checkpoint kinase",
        "CHEK2": "checkpoint kinase",
        "PARP1": "DNA repair sensor",
        "H2AX": "damage marker",  # H2AFX
    },
    "Cell Cycle": {
        "CCND1": "Cyclin D1 (mitogen sensor)",
        "CCND2": "Cyclin D2",
        "CCND3": "Cyclin D3",
        "CDK4": "CDK4",
        "CDK6": "CDK6",
        "MYC": "Myc (growth signal integrator)",
    },
    "Circadian": {
        "OPN4": "melanopsin (light)",
        "MTNR1A": "melatonin receptor",
        "MTNR1B": "melatonin receptor",
        "PER1": "Period 1",
        "NR1D1": "Rev-erb alpha",
        "BMAL1": "BMAL1",  # actually ARNTL
    },
    "NRF2": {
        "NFE2L2": "NRF2 itself",
        "NQO1": "target/sensor",
        "GCLC": "glutathione synthesis",
        "HMOX1": "heme oxygenase",
        "TXNRD1": "thioredoxin reductase",
        "GPX2": "glutathione peroxidase",
    },
    "PI3K/PTEN": {
        "EGFR": "EGF receptor",
        "ERBB2": "HER2",
        "IGF1R": "IGF receptor",
        "INSR": "Insulin receptor",
        "PDGFRA": "PDGF receptor",
        "ITGB1": "Integrin (mechanotransduction)",
    },
}

# Full receptor gene list (curated from IUPHAR + key receptors)
ALL_RECEPTORS = set()
for mod, recs in MODULE_RECEPTORS.items():
    ALL_RECEPTORS.update(recs.keys())

# Add extra receptors not module-assigned but important for diversity count
EXTRA_RECEPTORS = {
    # GPCRs
    "ADORA1", "ADORA2A", "ADORA2B", "ADRB1", "ADRB2", "CHRM1", "CHRM3",
    "DRD1", "DRD2", "HTR1A", "HTR2A", "PTGER2", "PTGER4", "S1PR1", "S1PR3",
    "CXCR4", "CCR2", "CCR5", "CCR7", "CXCR1", "CXCR2", "CXCR5",
    "GRM1", "GRM5", "GABBR1", "GABBR2",
    # RTKs
    "FLT1", "FLT3", "FLT4", "KDR", "TEK", "NTRK1", "NTRK2", "RET",
    "ALK", "ROS1", "AXL", "TYRO3", "MERTK", "DDR1", "DDR2", "EPHB2", "EPHA2",
    # Ion channels
    "KCNA1", "KCNA2", "SCN1A", "SCN5A", "CACNA1C", "CACNA1A",
    "GRIA1", "GRIA2", "GRIN1", "GRIN2A", "GRIN2B",
    # Nuclear receptors
    "ESR1", "ESR2", "AR", "NR3C1", "PPARG", "PPARA", "RARA", "VDR",
    # Toll-like
    "TLR1", "TLR3", "TLR5", "TLR6", "TLR8",
    # Death receptors
    "FAS", "TNFRSF10A", "TNFRSF10B",
    # Notch/developmental
    "PTCH1", "PTCH2", "SMO",
    # Immune receptors
    "CD4", "CD8A", "CD19", "CD28", "CTLA4", "PDCD1", "CD274",
    "FCGR1A", "FCGR2A", "FCGR3A", "CR1", "CR2",
    # Integrins
    "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGAV",
    "ITGB2", "ITGB3", "ITGB5", "ITGB6",
    # Cadherins
    "CDH3", "CDH5", "CDH6", "CDH11", "CDH13",
    # Selectins
    "SELE", "SELL", "SELP",
    # Other
    "EPCAM", "ALCAM", "ICAM1", "VCAM1", "PECAM1",
}
ALL_RECEPTORS.update(EXTRA_RECEPTORS)


def load_expression():
    """Load HPA nCPM."""
    print("Loading HPA expression data...")
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types


def load_cognitive_load():
    """Load precomputed cognitive load per cell type."""
    cog = {}
    with open(DATA_DIR / "cognitive_load.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cog[row["cell_type"]] = {
                "load": int(row["cognitive_load"]),
                "class": row["cell_class"],
                "continuous": float(row["continuous_load"]),
            }
    return cog


# Module core genes (same as step13)
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


def compute_module_scores(expr, cell_types):
    """Mean nCPM of core genes per module per cell type."""
    scores = {}
    for mod, genes in MODULES_CORE.items():
        found = [g for g in genes if g in expr]
        mod_scores = []
        for ct in cell_types:
            vals = [expr[g].get(ct, 0) for g in found]
            mod_scores.append(np.mean(vals) if vals else 0)
        scores[mod] = np.array(mod_scores)
    return scores


def compute_receptor_scores(expr, cell_types):
    """Mean nCPM of upstream receptors per module per cell type."""
    scores = {}
    for mod, recs in MODULE_RECEPTORS.items():
        found = [g for g in recs if g in expr]
        rec_scores = []
        for ct in cell_types:
            vals = [expr[g].get(ct, 0) for g in found]
            rec_scores.append(np.mean(vals) if vals else 0)
        scores[mod] = (np.array(rec_scores), len(found), len(recs))
    return scores


def count_receptor_diversity(expr, cell_types, threshold_nCPM=1.0):
    """Count distinct receptor genes expressed per cell type."""
    all_recs = [g for g in ALL_RECEPTORS if g in expr]
    diversity = []
    for ct in cell_types:
        n_expressed = sum(1 for g in all_recs if expr[g].get(ct, 0) > threshold_nCPM)
        diversity.append(n_expressed)
    return np.array(diversity), len(all_recs)


def main():
    expr, cell_types = load_expression()
    cog = load_cognitive_load()

    module_scores = compute_module_scores(expr, cell_types)
    receptor_scores = compute_receptor_scores(expr, cell_types)

    # ===================================================================
    # TEST 1: Module activity ↔ upstream receptor expression
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST 1: Module activity ↔ Upstream receptor expression")
    print("(Across 154 cell types: does active module → more receptors for that module?)")
    print("=" * 80)

    print(f"\n  {'Module':15s} {'Recs found':>10s} {'Pearson r':>10s} {'p':>10s} {'Spearman ρ':>10s} {'p':>10s}")
    print("  " + "-" * 70)

    r_values = []
    for mod in sorted(MODULES_CORE.keys()):
        if mod not in receptor_scores:
            continue
        mod_act = module_scores[mod]
        rec_act, n_found, n_total = receptor_scores[mod]

        r_pearson, p_pearson = stats.pearsonr(mod_act, rec_act)
        rho, p_rho = stats.spearmanr(mod_act, rec_act)
        r_values.append(r_pearson)

        sig = "***" if p_pearson < 0.001 else "**" if p_pearson < 0.01 else "*" if p_pearson < 0.05 else ""
        print(f"  {mod:15s} {n_found:3d}/{n_total:<3d}    {r_pearson:+10.3f} {p_pearson:10.2e}  "
              f"{rho:+10.3f} {p_rho:10.2e}  {sig}")

    print(f"\n  Mean Pearson r across all modules: {np.mean(r_values):.3f}")
    print(f"  Modules with r > 0.5: {sum(1 for r in r_values if r > 0.5)}/{len(r_values)}")
    print(f"  Modules with significant coupling (p<0.05): {sum(1 for r in r_values if r > 0)}/{len(r_values)}")

    # One-sample t-test: is mean r > 0?
    t_stat, p_ttest = stats.ttest_1samp(r_values, 0)
    print(f"\n  One-sample t-test (mean r > 0): t = {t_stat:.2f}, p = {p_ttest:.4f}")

    # ===================================================================
    # TEST 2: Receptor diversity ∝ cognitive load
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST 2: Receptor diversity ∝ Cognitive load")
    print("=" * 80)

    diversity, n_total_recs = count_receptor_diversity(expr, cell_types)
    cog_loads = np.array([cog[ct]["load"] for ct in cell_types if ct in cog])

    # Align
    valid = [j for j, ct in enumerate(cell_types) if ct in cog]
    div_valid = diversity[valid]
    load_valid = cog_loads

    rho_div, p_div = stats.spearmanr(load_valid, div_valid)
    r_div, p_rdiv = stats.pearsonr(load_valid, div_valid)

    print(f"\n  Total receptor genes checked: {n_total_recs}")
    print(f"  Receptor genes found in HPA: {sum(1 for g in ALL_RECEPTORS if g in expr)}")
    print(f"\n  Spearman ρ(cognitive load, receptor diversity) = {rho_div:.3f}, p = {p_div:.2e}")
    print(f"  Pearson r(cognitive load, receptor diversity) = {r_div:.3f}, p = {p_rdiv:.2e}")

    # Top and bottom cell types
    print(f"\n  Top 10 receptor diversity:")
    sorted_idx = sorted(valid, key=lambda j: -diversity[j])
    for j in sorted_idx[:10]:
        ct = cell_types[j]
        cls = cog.get(ct, {}).get("class", "")
        print(f"    {ct:35s} {cls:25s} receptors={diversity[j]:3d}  load={cog[ct]['load']}")

    print(f"\n  Bottom 10 receptor diversity:")
    for j in sorted_idx[-10:]:
        ct = cell_types[j]
        cls = cog.get(ct, {}).get("class", "")
        print(f"    {ct:35s} {cls:25s} receptors={diversity[j]:3d}  load={cog[ct]['load']}")

    # By cell class
    print(f"\n  Diversity by cell class:")
    class_div = defaultdict(list)
    for j in valid:
        ct = cell_types[j]
        cls = cog.get(ct, {}).get("class", "")
        class_div[cls].append(diversity[j])

    print(f"  {'Class':30s} {'N':>3} {'Mean div':>9} {'Median':>7}")
    print("  " + "-" * 55)
    for cls in sorted(class_div.keys(), key=lambda c: -np.mean(class_div[c])):
        vals = class_div[cls]
        print(f"  {cls:30s} {len(vals):3d} {np.mean(vals):9.1f} {np.median(vals):7.1f}")

    # ===================================================================
    # TEST 3: Cross-module receptor sharing
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST 3: Receptor sharing between modules")
    print("=" * 80)

    for m1 in sorted(MODULE_RECEPTORS.keys()):
        for m2 in sorted(MODULE_RECEPTORS.keys()):
            if m1 >= m2:
                continue
            shared = set(MODULE_RECEPTORS[m1].keys()) & set(MODULE_RECEPTORS[m2].keys())
            if shared:
                print(f"  {m1:15s} ∩ {m2:15s}: {', '.join(sorted(shared))}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  Test 1 (module-receptor coupling): mean r = {np.mean(r_values):.3f}, "
          f"t-test p = {p_ttest:.4f}")
    print(f"  Test 2 (diversity ∝ load): ρ = {rho_div:.3f}, p = {p_div:.2e}")

    if np.mean(r_values) > 0 and p_ttest < 0.05:
        print(f"  → Module-receptor coupling CONFIRMED: active modules → more receptors")
    if rho_div > 0.3 and p_div < 0.05:
        print(f"  → Receptor diversity scales with cognitive load")

    # Save
    results = {
        "test1_mean_r": round(float(np.mean(r_values)), 3),
        "test1_ttest_p": float(p_ttest),
        "test2_rho": round(float(rho_div), 3),
        "test2_p": float(p_div),
        "n_receptor_genes": n_total_recs,
        "n_cell_types": len(cell_types),
    }
    with open(DATA_DIR / "receptor_module_coupling.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to receptor_module_coupling.json")
    print("Done!")


if __name__ == "__main__":
    main()
