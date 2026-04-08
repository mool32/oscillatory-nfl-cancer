#!/usr/bin/env python3
"""
step29: Differentiation trajectory in the oscillatory module eigenspace.

Question: As cells differentiate (E9.5→E13.5), do they move systematically
through the 6D eigenspace? Do lineages heading toward "cancer-vulnerable"
cell types show distinct trajectory shapes?

Data available:
  - moca_group_entropy.csv: 185 rows, per Main_Cluster×trajectory×day
    Fields: E_intra_mean, E_inter (oscillatory entropy)
  - moca_per_cell.csv: 300K cells with trajectory, day, cell_type, E_intra
  - HPA eigenspace: 20-module eigenvectors (step19 eigen20.json)

Analysis layers:
  1. E_intra trajectory per lineage: does complexity increase/decrease with development?
  2. Lineage → adult cell type mapping: where do mature versions land in eigenspace?
  3. Cross-lineage comparison: neural vs hematopoietic vs mesenchymal signatures
  4. Cancer relevance: do cancer-relevant cell types (epithelial, hematopoietic)
     show specific developmental patterns?

Key biological question:
  "Is cancer vulnerability acquired during differentiation (ascending trajectory)
   or intrinsic to adult cell type architecture?"

Connection to K5 (mutation trajectory):
  Cancer = ascending IA+BN (break accelerators first, then checkpoints).
  Does development itself build up IA+BN? If yes → cancer exploits the
  developmental ascending strategy in reverse.
"""

import csv
import json
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")
MOCA_DIR = Path("/Users/teo/Desktop/research/oscilatory/results/embryogenesis/analysis")

# ── Trajectory → adult HPA cell type mapping ──────────────────────────────
# Based on MOCA trajectory labels and HPA cell type names
TRAJECTORY_TO_HPA = {
    "Neural":             ["brain excitatory neurons", "brain inhibitory neurons",
                           "astrocytes", "oligodendrocytes", "microglia"],
    "Epithelial":         ["colonocytes", "enterocytes", "hepatocytes",
                           "respiratory ciliated cells", "alveolar cells type 1"],
    "Hematopoietic":      ["b-cells", "t-cells", "monocytes", "neutrophils",
                           "macrophages", "nk-cells", "erythrocytes"],
    "Mesenchyme":         ["fibroblasts", "smooth muscle cells", "pericytes",
                           "fibro-adipogenic progenitors"],
    "Endothelial":        ["vascular endothelial cells", "lymphatic endothelial cells"],
    "Intermediate_mesoderm": ["podocytes", "proximal tubule cells",
                               "loop of henle epithelial cells"],
    "Hepatic":            ["hepatocytes", "cholangiocytes", "kupffer cells"],
    "Neural_crest":       ["melanocytes", "schwann cells", "adrenal medulla cells"],
    "Notochord":          [],  # no direct adult equivalent
    "Eye":                ["retinal ganglion cells", "rod photoreceptor cells",
                           "cone photoreceptor cells", "retinal pigment epithelial cells"],
    "Lens":               ["ocular epithelial cells"],
    "Kidney":             ["podocytes", "proximal tubule cells",
                           "renal collecting duct principal cells"],
}

# Cancer relevance from CGC density per module (from paper results)
# Which trajectories lead to cancer-relevant adult cell types?
CANCER_RELEVANCE = {
    "Hematopoietic":      "HIGH",    # leukemias, lymphomas — major cancer burden
    "Epithelial":         "HIGH",    # carcinomas — largest cancer burden
    "Hepatic":            "HIGH",    # HCC, cholangiocarcinoma
    "Neural":             "MEDIUM",  # glioma, neuroblastoma
    "Neural_crest":       "MEDIUM",  # melanoma, neuroblastoma
    "Mesenchyme":         "MEDIUM",  # sarcoma
    "Kidney":             "MEDIUM",  # RCC
    "Intermediate_mesoderm": "LOW",
    "Endothelial":        "LOW",
    "Eye":                "LOW",
    "Notochord":          "LOW",
    "Lens":               "LOW",
}

# IRreversible Authority scores (from paper step22, IA metric)
# These are the module IA scores that we'll use to weight the eigenspace position
MODULE_IA = {
    "NF-κB": 2, "ERK/MAPK": 2, "JAK-STAT": 2, "p53": 3, "Wnt": 3,
    "Notch": 3, "Hippo": 2, "TGF-β": 2, "mTOR": 1, "Calcium": 1,
    "Cell Cycle": 2, "Circadian": 0, "NRF2": 0, "PI3K/PTEN": 2,
    "AMPK": 1, "SREBP": 0, "ATR/CHK1": 2, "Rho/ROCK": 1,
    "PPAR/LXR": 0, "Autophagy": 1,
}

MODULE_ORDER = [
    "NF-κB","ERK/MAPK","JAK-STAT","p53","Wnt","Notch",
    "Hippo","TGF-β","mTOR","Calcium","Cell Cycle",
    "Circadian","NRF2","PI3K/PTEN",
    "AMPK","SREBP","ATR/CHK1","Rho/ROCK","PPAR/LXR","Autophagy",
]

ALL_MODULES = {
    "NF-κB":      {"RELA","RELB","NFKB1","NFKB2","NFKBIA","NFKBIB","TNFAIP3",
                   "TRAF2","TRAF6","IKBKB","IKBKG","CHUK","MAP3K7","TAB1","TAB2"},
    "ERK/MAPK":   {"MAPK1","MAPK3","MAP2K1","MAP2K2","BRAF","RAF1","ARAF",
                   "HRAS","KRAS","NRAS","SOS1","GRB2","DUSP1","DUSP6","SPRY2"},
    "JAK-STAT":   {"JAK1","JAK2","JAK3","TYK2","STAT1","STAT3","STAT5A","STAT5B",
                   "SOCS1","SOCS3","CISH","PIAS1"},
    "p53":        {"TP53","MDM2","MDM4","CDKN1A","BAX","BBC3","PMAIP1",
                   "ATM","ATR","CHEK1","CHEK2"},
    "Wnt":        {"CTNNB1","APC","AXIN1","AXIN2","GSK3B","DVL1",
                   "TCF7L2","LEF1","RNF43","ZNRF3"},
    "Notch":      {"NOTCH1","NOTCH2","NOTCH3","NOTCH4","FBXW7","HES1","HEY1",
                   "MAML1","RBPJ","DLL1","DLL4","JAG1","JAG2"},
    "Hippo":      {"YAP1","WWTR1","LATS1","LATS2","STK3","STK4","SAV1",
                   "MOB1A","NF2","TEAD1","TEAD4"},
    "TGF-β":      {"TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                   "SMURF1","SMURF2","BMPR1A","BMPR2","ACVR1"},
    "mTOR":       {"MTOR","RPTOR","RICTOR","TSC1","TSC2","RHEB","RPS6KB1",
                   "EIF4EBP1","DEPTOR","MLST8"},
    "Calcium":    {"PLCG1","PLCG2","ITPR1","ITPR2","ATP2A2","CALM1",
                   "NFATC1","NFATC2","CAMK2A","CAMK2B"},
    "Cell Cycle": {"CDK2","CDK4","CDK6","CCND1","CCNE1","CCNA2","CCNB1",
                   "RB1","E2F1","CDKN1A","CDKN2A","CDKN1B","CDC25A"},
    "Circadian":  {"CLOCK","ARNTL","PER1","PER2","CRY1","CRY2","NR1D1",
                   "NR1D2","CSNK1D","CSNK1E","FBXL3"},
    "NRF2":       {"NFE2L2","KEAP1","HMOX1","NQO1","GCLC","GCLM","TXNRD1",
                   "SOD2","CAT","GPX1"},
    "PI3K/PTEN":  {"PIK3CA","PIK3CB","PIK3R1","PTEN","AKT1","AKT2",
                   "PDK1","INPP4B"},
    "AMPK":       {"PRKAA1","PRKAA2","PRKAB1","PRKAG1","STK11","ACACB",
                   "PPARGC1A","FOXO3","TSC2","CREB1"},
    "SREBP":      {"SREBF1","SREBF2","SCAP","INSIG1","INSIG2","HMGCR",
                   "FASN","SCD","ACLY"},
    "ATR/CHK1":   {"ATR","CHEK1","ATRIP","TOPBP1","WEE1","CDC25A",
                   "RAD17","RAD9A","HUS1"},
    "Rho/ROCK":   {"RHOA","RHOB","RHOC","ROCK1","ROCK2","MKL1",
                   "LIMK1","CFL1","ARHGAP1","ARHGAP5"},
    "PPAR/LXR":   {"PPARA","PPARG","PPARD","NR1H3","NR1H2","RXRA",
                   "NCOR1","NCOR2","NCOA1","NCOA2"},
    "Autophagy":  {"ULK1","ULK2","BECN1","ATG5","ATG7","ATG12",
                   "TFEB","SQSTM1","MAP1LC3B"},
}


def load_entropy_data():
    """Load moca_group_entropy.csv."""
    data = []
    with open(MOCA_DIR / "moca_group_entropy.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                "cluster": float(row["Main_Cluster"]),
                "cell_type": row["cell_type"],
                "trajectory": row["trajectory"],
                "day": float(row["day"]),
                "n_cells": int(row["n_cells"]),
                "E_intra": float(row["E_intra_mean"]),
                "E_inter": float(row["E_inter"]),
            })
    return data


def load_hpa_expression():
    """Load HPA expression for adult module scoring."""
    from collections import defaultdict
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    return expr


def compute_module_activity(expr, cell_types):
    """Compute mean activity per module for a set of cell types."""
    scores = {}
    for mod in MODULE_ORDER:
        genes = [g for g in ALL_MODULES[mod] if g in expr]
        if not genes:
            scores[mod] = 0.0
            continue
        vals = []
        for ct in cell_types:
            for g in genes:
                vals.append(expr[g].get(ct, 0))
        scores[mod] = np.mean(vals) if vals else 0.0
    return scores


def compute_ia_weighted_score(module_scores):
    """Compute IA-weighted mean module activity — proxy for vulnerability."""
    total_ia = sum(MODULE_IA.values())
    weighted = sum(module_scores.get(m, 0) * MODULE_IA[m] for m in MODULE_ORDER)
    return weighted / total_ia if total_ia > 0 else 0.0


def project_to_eigenspace(module_scores_dict, eigenvectors, module_order):
    """
    Project module activity scores onto eigenspace axes.
    Returns coordinates on each PC.
    """
    vec = np.array([module_scores_dict.get(m, 0.0) for m in module_order])
    # Normalize
    vec = (vec - np.mean(vec)) / (np.std(vec) + 1e-12)
    coords = eigenvectors.T @ vec
    return coords


def main():
    print("Loading MOCA entropy trajectories...")
    entropy_data = load_entropy_data()

    # ===================================================================
    # LAYER 1: E_intra trajectory per lineage
    # ===================================================================
    print("\n" + "="*80)
    print("LAYER 1: E_INTRA TRAJECTORY PER DEVELOPMENTAL LINEAGE")
    print("="*80)

    # Group by trajectory × day
    traj_day = defaultdict(lambda: defaultdict(list))
    for row in entropy_data:
        traj_day[row["trajectory"]][row["day"]].append(row["E_intra"])

    traj_stats = {}
    for traj in sorted(traj_day.keys()):
        days_sorted = sorted(traj_day[traj].keys())
        means = [np.mean(traj_day[traj][d]) for d in days_sorted]
        sems = [np.std(traj_day[traj][d]) / np.sqrt(len(traj_day[traj][d]))
                for d in days_sorted]
        n_days = len(days_sorted)

        if n_days >= 3:
            rho, p = stats.spearmanr(days_sorted, means)
        else:
            rho, p = float('nan'), float('nan')

        direction = "↑ increasing" if rho > 0.3 else ("↓ decreasing" if rho < -0.3 else "→ stable")
        cancer_rel = CANCER_RELEVANCE.get(traj, "UNKNOWN")

        print(f"\n  {traj:25s} [cancer: {cancer_rel}]")
        print(f"    Days: " + "  ".join(f"E{d:.1f}:{np.mean(traj_day[traj][d]):.3f}" for d in days_sorted))
        print(f"    ρ(day,E_intra)={rho:+.3f} p={p:.4f} → {direction}")

        traj_stats[traj] = {
            "days": days_sorted,
            "E_intra_means": means,
            "E_intra_sems": sems,
            "rho_day": rho if not np.isnan(rho) else None,
            "p_day": p if not np.isnan(p) else None,
            "direction": direction.split()[1] if direction else "stable",
            "cancer_relevance": cancer_rel,
        }

    # ===================================================================
    # LAYER 2: E_inter trajectory (between-cell diversity)
    # ===================================================================
    print("\n" + "="*80)
    print("LAYER 2: E_INTER TRAJECTORY (BETWEEN-CELL DIVERSITY)")
    print("="*80)

    traj_inter = defaultdict(lambda: defaultdict(list))
    for row in entropy_data:
        traj_inter[row["trajectory"]][row["day"]].append(row["E_inter"])

    print(f"\n  {'Trajectory':25s} {'ρ(day,E_inter)':>15s} {'p':>8s} {'Direction'}")
    print("  " + "-" * 60)
    for traj in sorted(traj_inter.keys()):
        days_sorted = sorted(traj_inter[traj].keys())
        means = [np.mean(traj_inter[traj][d]) for d in days_sorted]
        if len(days_sorted) >= 3:
            rho, p = stats.spearmanr(days_sorted, means)
        else:
            rho, p = float('nan'), float('nan')
        direction = "↑" if rho > 0.3 else ("↓" if rho < -0.3 else "→")
        print(f"  {traj:25s} {rho:+15.3f} {p:8.4f} {direction}")
        if traj in traj_stats:
            traj_stats[traj]["E_inter_rho"] = float(rho) if not np.isnan(rho) else None

    # ===================================================================
    # LAYER 3: Does ΔE_intra and ΔE_inter anti-correlate per trajectory?
    # (Replication of aging result in development context)
    # ===================================================================
    print("\n" + "="*80)
    print("LAYER 3: ρ(ΔE_intra, ΔE_inter) — REPLICATION IN DEVELOPMENT")
    print("="*80)

    delta_intra_list = []
    delta_inter_list = []

    for traj in sorted(traj_day.keys()):
        days_intra = sorted(traj_day[traj].keys())
        days_inter = sorted(traj_inter[traj].keys())
        if len(days_intra) < 3 or len(days_inter) < 3:
            continue
        means_intra = [np.mean(traj_day[traj][d]) for d in days_intra]
        means_inter = [np.mean(traj_inter[traj][d]) for d in days_inter]
        if len(days_intra) >= 3:
            rho_i, _ = stats.spearmanr(days_intra, means_intra)
            rho_e, _ = stats.spearmanr(days_inter, means_inter)
            delta_intra_list.append(rho_i)
            delta_inter_list.append(rho_e)

    if len(delta_intra_list) >= 4:
        rho_rep, p_rep = stats.spearmanr(delta_intra_list, delta_inter_list)
        print(f"\n  ρ(ΔE_intra_slope, ΔE_inter_slope) = {rho_rep:+.3f}, p = {p_rep:.4f}")
        print(f"  N trajectories = {len(delta_intra_list)}")
        if rho_rep < -0.3:
            print(f"  ✓ Anti-correlation confirmed in development (consistent with aging result)")
        elif rho_rep > 0.3:
            print(f"  ✗ Positive correlation — development differs from aging")
        else:
            print(f"  → No clear direction — development may not replicate aging pattern")
    else:
        print(f"\n  Insufficient data for cross-trajectory correlation")
        rho_rep, p_rep = None, None

    # ===================================================================
    # LAYER 4: Adult eigenspace position per trajectory
    # ===================================================================
    print("\n" + "="*80)
    print("LAYER 4: ADULT EIGENSPACE POSITION PER LINEAGE")
    print("="*80)

    # Load HPA expression and eigenvectors
    print("  Loading HPA expression data...")
    expr = load_hpa_expression()

    # Load eigendecomposition from step19
    eigen_file = DATA_DIR / "eigen20.json"
    if not eigen_file.exists():
        # Try alternative paths
        for fname in ["eigendecomposition.json", "eigenspace_cancer.json"]:
            if (DATA_DIR / fname).exists():
                eigen_file = DATA_DIR / fname
                break

    if eigen_file.exists():
        with open(eigen_file) as f:
            eigen_data = json.load(f)
        # Reconstruct eigenvectors
        # eigen20.json format depends on how step19 saved it
        if "eigenvectors" in eigen_data:
            vecs_global = np.array(eigen_data["eigenvectors"])
        elif "pcs" in eigen_data:
            vecs_global = np.array([pc["loadings"] for pc in eigen_data["pcs"]]).T
        else:
            vecs_global = None
        print(f"  Loaded eigenspace from {eigen_file.name}")
    else:
        # Recompute from HPA data
        print("  Recomputing eigenspace from HPA data...")
        all_cts = list(set(ct for gene_d in expr.values() for ct in gene_d))
        N = len(MODULE_ORDER)
        M = len(all_cts)
        mat = np.zeros((N, M))
        from collections import defaultdict as dd
        for i, mod in enumerate(MODULE_ORDER):
            genes = [g for g in ALL_MODULES[mod] if g in expr]
            for j, ct in enumerate(all_cts):
                vals = [expr[g].get(ct, 0) for g in genes]
                mat[i, j] = np.mean(vals) if vals else 0
        corr = np.corrcoef(mat)
        corr = np.nan_to_num(corr, nan=0.0)
        np.fill_diagonal(corr, 1.0)
        ev, vecs_global = np.linalg.eigh(corr)
        idx = np.argsort(ev)[::-1]
        vecs_global = vecs_global[:, idx]
        print(f"  Eigenspace recomputed (N={M} cell types)")

    # Compute module activity for each trajectory's adult cell types
    print(f"\n  {'Trajectory':25s} {'Cancer':>8s} {'IA-score':>10s} {'PC1':>7s} {'PC2':>7s} {'PC3':>7s}")
    print("  " + "-" * 70)

    traj_eigen_data = {}
    for traj in sorted(TRAJECTORY_TO_HPA.keys()):
        hpa_cts = TRAJECTORY_TO_HPA[traj]
        # Filter to HPA cell types that exist in our data
        available_cts = [ct for ct in hpa_cts if any(ct in d for d in expr.values())]
        # Actually check if gene has this ct
        available_cts = []
        for ct in hpa_cts:
            # Check if this ct appears in any gene's expression
            found = any(ct in gene_dict for gene_dict in list(expr.values())[:100])
            if found:
                available_cts.append(ct)

        if not available_cts:
            print(f"  {traj:25s} {'':>8s} {'N/A (no HPA match)':>10s}")
            continue

        mod_scores = compute_module_activity(expr, available_cts)
        ia_score = compute_ia_weighted_score(mod_scores)
        cancer_rel = CANCER_RELEVANCE.get(traj, "UNKNOWN")

        if vecs_global is not None:
            coords = project_to_eigenspace(mod_scores, vecs_global, MODULE_ORDER)
            pc1, pc2, pc3 = coords[0], coords[1], coords[2]
            print(f"  {traj:25s} {cancer_rel:>8s} {ia_score:10.2f} {pc1:+7.3f} {pc2:+7.3f} {pc3:+7.3f}")
        else:
            pc1, pc2, pc3 = None, None, None
            print(f"  {traj:25s} {cancer_rel:>8s} {ia_score:10.2f}")

        traj_eigen_data[traj] = {
            "hpa_cell_types": available_cts,
            "ia_score": round(float(ia_score), 3),
            "pc1": round(float(pc1), 3) if pc1 is not None else None,
            "pc2": round(float(pc2), 3) if pc2 is not None else None,
            "pc3": round(float(pc3), 3) if pc3 is not None else None,
            "cancer_relevance": cancer_rel,
        }

    # ===================================================================
    # LAYER 5: Do cancer-relevant trajectories have distinct E_intra dynamics?
    # ===================================================================
    print("\n" + "="*80)
    print("LAYER 5: CANCER RELEVANCE × E_INTRA TRAJECTORY")
    print("="*80)

    high_cancer_rhos = []
    low_cancer_rhos = []

    for traj, st in traj_stats.items():
        if st["rho_day"] is None:
            continue
        rho = st["rho_day"]
        if st["cancer_relevance"] == "HIGH":
            high_cancer_rhos.append(rho)
        elif st["cancer_relevance"] == "LOW":
            low_cancer_rhos.append(rho)

    if high_cancer_rhos and low_cancer_rhos:
        print(f"\n  HIGH cancer trajectories (N={len(high_cancer_rhos)}):")
        print(f"    ρ(day, E_intra): {high_cancer_rhos}")
        print(f"    Mean: {np.mean(high_cancer_rhos):+.3f}")
        print(f"\n  LOW cancer trajectories (N={len(low_cancer_rhos)}):")
        print(f"    ρ(day, E_intra): {low_cancer_rhos}")
        print(f"    Mean: {np.mean(low_cancer_rhos):+.3f}")

        if len(high_cancer_rhos) >= 2 and len(low_cancer_rhos) >= 2:
            u, p_mwu = stats.mannwhitneyu(high_cancer_rhos, low_cancer_rhos,
                                           alternative='two-sided')
            print(f"\n  Mann-Whitney U test: U={u:.0f}, p={p_mwu:.4f}")
            if p_mwu < 0.1:
                direction = "higher" if np.mean(high_cancer_rhos) > np.mean(low_cancer_rhos) else "lower"
                print(f"  → Cancer-relevant trajectories show {direction} E_intra slope during development")

    # ===================================================================
    # CONNECTION TO K5 (mutation order)
    # ===================================================================
    print("\n" + "="*80)
    print("CONNECTION TO K5: DEVELOPMENT AS ASCENDING IA+BN?")
    print("="*80)

    # Load IA+BN scores from step26 if available
    ia_bn_file = DATA_DIR / "vulnerability_metric.json"
    if ia_bn_file.exists():
        with open(ia_bn_file) as f:
            vuln_data = json.load(f)
        print("\n  Module IA+BN scores available:")
        if "scores" in vuln_data:
            scores = vuln_data["scores"]
            sorted_modules = sorted(scores.items(), key=lambda x: x[1], reverse=True)
            for mod, score in sorted_modules[:5]:
                print(f"    {mod}: {score:.2f}")

    print("\n  Hypothesis: If development = ascending IA+BN (mirror of K5),")
    print("  then E_intra should INCREASE during differentiation as cells")
    print("  gain more complex (high-IA) regulatory architectures.")
    print()

    increasing_trjs = [t for t, s in traj_stats.items()
                       if s.get("rho_day") and s["rho_day"] > 0.3]
    decreasing_trjs = [t for t, s in traj_stats.items()
                       if s.get("rho_day") and s["rho_day"] < -0.3]
    stable_trjs = [t for t, s in traj_stats.items()
                   if s.get("rho_day") and abs(s["rho_day"]) <= 0.3]

    print(f"  ↑ Increasing E_intra ({len(increasing_trjs)}): {increasing_trjs}")
    print(f"  ↓ Decreasing E_intra ({len(decreasing_trjs)}): {decreasing_trjs}")
    print(f"  → Stable E_intra ({len(stable_trjs)}): {stable_trjs}")
    print()

    if len(increasing_trjs) > len(decreasing_trjs):
        print("  ✓ CONSISTENT with developmental ascending hypothesis:")
        print("    Most lineages INCREASE transcriptional complexity during differentiation")
        print("    → Development builds up the regulatory architecture that cancer later exploits")
    elif len(decreasing_trjs) > len(increasing_trjs):
        print("  ✗ INCONSISTENT: Most lineages DECREASE complexity (pruning/specialization dominates)")
        print("    → Cancer vulnerability is NOT a simple mirror of developmental trajectory")
    else:
        print("  → Mixed: lineage-specific patterns, no universal developmental direction")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "n_trajectories": len(traj_stats),
        "n_stages": 5,
        "stage_range": "E9.5-E13.5",
        "trajectory_stats": {
            t: {
                "rho_day_eintra": s.get("rho_day"),
                "p_day_eintra": s.get("p_day"),
                "direction": s.get("direction"),
                "cancer_relevance": s.get("cancer_relevance"),
                "E_intra_means": [round(v, 4) for v in s.get("E_intra_means", [])],
            }
            for t, s in traj_stats.items()
        },
        "eigenspace_positions": traj_eigen_data,
        "dev_aging_replication": {
            "rho": float(rho_rep) if rho_rep is not None else None,
            "p": float(p_rep) if p_rep is not None else None,
            "n": len(delta_intra_list),
        },
        "increasing_trajectories": increasing_trjs,
        "decreasing_trajectories": decreasing_trjs,
        "cancer_high_mean_rho": round(float(np.mean(high_cancer_rhos)), 3) if high_cancer_rhos else None,
        "cancer_low_mean_rho": round(float(np.mean(low_cancer_rhos)), 3) if low_cancer_rhos else None,
    }

    with open(DATA_DIR / "moca_differentiation.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to moca_differentiation.json")


if __name__ == "__main__":
    main()
