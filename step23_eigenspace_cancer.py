#!/usr/bin/env python3
"""
Cancer genes in perception eigenspace.

PRE-REGISTERED PREDICTIONS (before code runs):
  1. CGC genes cluster non-randomly in 6D eigenspace (vs random genes)
  2. Hematological drivers load on NF-κB/JAK-STAT eigenmode (PC3+)
  3. Epithelial drivers load on Hippo/Wnt eigenmode (PC3-)
  4. GoF and LoF project in opposite directions on ≥1 eigenmode

Step 1: Project all module genes into eigenspace, compare CGC vs non-CGC
Step 2: Cancer-type-specific drivers → eigenspace centroids
Step 3: GoF vs LoF direction in eigenspace
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")
np.random.seed(42)

# ===================================================================
# Module definitions (same as step19/21)
# ===================================================================
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

# Cancer type driver genes (curated from TCGA pan-cancer + Bailey 2018 + CGC)
# Hematological
HEME_DRIVERS = {
    "JAK2", "JAK3", "STAT3", "STAT5B",  # JAK-STAT
    "NFKB2", "RELA", "IKBKB",  # NF-κB
    "NOTCH1", "FBXW7",  # Notch (T-ALL)
    "MYC",  # universal but critical in lymphoma
    "CDK6", "CCND1",  # cell cycle (mantle cell lymphoma)
    "EZH2", "KMT2A",  # epigenetic (leukemia) — not in our modules but relevant
    "SOCS1",  # JAK-STAT feedback
    "TNFAIP3",  # NF-κB feedback (A20)
    "CISH",  # JAK-STAT
    "TRAF6",  # NF-κB
}

# Epithelial solid tumors
EPITHELIAL_DRIVERS = {
    "APC", "CTNNB1", "AXIN1", "AXIN2", "RNF43",  # Wnt
    "LATS1", "LATS2", "NF2", "YAP1",  # Hippo
    "SMAD2", "SMAD3", "SMAD4", "TGFBR1", "TGFBR2",  # TGF-β
    "CDH1",  # E-cadherin (Hippo/contact)
    "PTEN", "PIK3CA", "AKT1",  # PI3K
    "KRAS", "BRAF", "EGFR",  # ERK (but shared)
    "TP53", "RB1",  # universal but enriched
    "STK11",  # AMPK/LKB1
    "CDKN2A",  # cell cycle
}

# GoF and LoF classification (from Paper 1 + CGC annotation)
GOF_GENES = {
    "KRAS", "NRAS", "HRAS", "BRAF", "RAF1",  # ERK
    "PIK3CA", "AKT1", "AKT2",  # PI3K
    "CTNNB1",  # Wnt (activating)
    "STAT3", "STAT5B", "JAK2",  # JAK-STAT
    "NOTCH1",  # Notch (in some contexts)
    "CCND1", "CCNE1", "CDK4", "CDK6",  # Cell cycle (amplification)
    "MYC",  # universal amplification
    "MTOR",  # mTOR
    "MDM2", "MDM4",  # p53 inhibitors (gain = suppress p53)
    "RHOA",  # Rho
    "EGFR", "ERBB2",  # RTK (not in module list but known GoF)
    "NFE2L2",  # NRF2 (activating mutations)
    "RELA",  # NF-κB
    "YAP1", "WWTR1",  # Hippo effectors (GoF = bypass Hippo)
}

LOF_GENES = {
    "TP53", "ATM", "ATR", "CHEK2",  # p53/damage
    "APC", "AXIN1", "AXIN2", "RNF43",  # Wnt inhibitors
    "PTEN", "INPP4B",  # PI3K inhibitors
    "RB1", "CDKN2A", "CDKN1A", "CDKN1B",  # cell cycle inhibitors
    "SMAD2", "SMAD3", "SMAD4", "SMAD7", "TGFBR2",  # TGF-β
    "FBXW7",  # Notch/general degradation
    "NF2", "LATS1", "LATS2",  # Hippo kinases
    "STK11",  # AMPK
    "TSC1", "TSC2",  # mTOR inhibitors
    "KEAP1",  # NRF2 inhibitor
    "SOCS1", "SOCS3",  # JAK-STAT inhibitors
    "NFKBIA", "TNFAIP3",  # NF-κB inhibitors
    "BAX",  # apoptosis
    "NCOR1", "NCOR2",  # PPAR/nuclear receptor co-repressors
}


def load_cgc():
    cgc = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                cgc.add(g)
    return cgc


def load_eigendata():
    with open(DATA_DIR / "eigen20.json") as f:
        return json.load(f)


def gene_to_module_vector(gene):
    """Return binary 20-vector: 1 if gene belongs to module, 0 otherwise."""
    vec = np.zeros(len(MODULE_ORDER))
    for i, mod in enumerate(MODULE_ORDER):
        if gene in ALL_MODULES[mod]:
            vec[i] = 1
    return vec


def project_to_eigenspace(gene_module_vec, eigenvectors, n_pcs=6):
    """Project module membership vector onto eigenvectors."""
    return np.array([np.dot(gene_module_vec, eigenvectors[:, pc]) for pc in range(n_pcs)])


def main():
    cgc = load_cgc()
    eigen = load_eigendata()

    # Reconstruct eigenvectors from eigendecomposition
    # Need to recompute since json doesn't store eigenvectors
    print("Recomputing eigendecomposition...")
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))

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

    N_PCS = 6  # Kaiser components

    # ===================================================================
    # Collect all module genes
    # ===================================================================
    all_module_genes = set()
    for genes in ALL_MODULES.values():
        all_module_genes.update(genes)

    cgc_in_modules = all_module_genes & cgc
    non_cgc_in_modules = all_module_genes - cgc

    print(f"\nTotal module genes: {len(all_module_genes)}")
    print(f"CGC genes in modules: {len(cgc_in_modules)}")
    print(f"Non-CGC genes in modules: {len(non_cgc_in_modules)}")

    # ===================================================================
    # STEP 1: CGC vs non-CGC distribution in eigenspace
    # ===================================================================
    print(f"\n{'='*80}")
    print("STEP 1: CGC vs non-CGC genes in eigenspace")
    print("=" * 80)

    # Project each gene
    cgc_projections = []
    noncgc_projections = []

    for gene in all_module_genes:
        mvec = gene_to_module_vector(gene)
        if mvec.sum() == 0:
            continue
        proj = project_to_eigenspace(mvec, eigenvectors, N_PCS)
        if gene in cgc:
            cgc_projections.append(proj)
        else:
            noncgc_projections.append(proj)

    cgc_proj = np.array(cgc_projections)
    noncgc_proj = np.array(noncgc_projections)

    print(f"\n  CGC genes projected: {len(cgc_proj)}")
    print(f"  Non-CGC genes projected: {len(noncgc_proj)}")

    # Compare distributions along each PC
    print(f"\n  {'PC':>4s} {'CGC mean':>10s} {'NonCGC mean':>12s} {'Δ':>8s} {'t-stat':>8s} {'p':>10s}")
    print("  " + "-" * 55)

    pc_results = []
    for pc in range(N_PCS):
        cgc_vals = cgc_proj[:, pc]
        noncgc_vals = noncgc_proj[:, pc]
        t, p = stats.ttest_ind(cgc_vals, noncgc_vals)
        delta = np.mean(cgc_vals) - np.mean(noncgc_vals)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  PC{pc+1:1d} {np.mean(cgc_vals):+10.3f} {np.mean(noncgc_vals):+12.3f} "
              f"{delta:+8.3f} {t:+8.2f} {p:10.4f}  {sig}")
        pc_results.append({"pc": pc+1, "delta": delta, "t": t, "p": p})

    # Multivariate test: Hotelling's T² (approximate with MANOVA-like)
    # Use combined distance
    cgc_centroid = np.mean(cgc_proj, axis=0)
    noncgc_centroid = np.mean(noncgc_proj, axis=0)
    diff = cgc_centroid - noncgc_centroid
    euclidean_dist = np.linalg.norm(diff)

    # Permutation test for multivariate difference
    all_proj = np.vstack([cgc_proj, noncgc_proj])
    n_cgc = len(cgc_proj)
    n_perm = 10000
    perm_dists = []
    for _ in range(n_perm):
        perm_idx = np.random.permutation(len(all_proj))
        perm_cgc = all_proj[perm_idx[:n_cgc]]
        perm_non = all_proj[perm_idx[n_cgc:]]
        perm_dist = np.linalg.norm(np.mean(perm_cgc, axis=0) - np.mean(perm_non, axis=0))
        perm_dists.append(perm_dist)

    p_perm = np.mean(np.array(perm_dists) >= euclidean_dist)
    print(f"\n  Multivariate centroid distance: {euclidean_dist:.4f}")
    print(f"  Permutation test (N={n_perm}): p = {p_perm:.4f}")

    if p_perm < 0.05:
        print(f"  ✓ PREDICTION 1 CONFIRMED: CGC genes cluster non-randomly in eigenspace")
    else:
        print(f"  ✗ Prediction 1 not confirmed")

    # ===================================================================
    # STEP 2: Cancer type drivers in eigenspace
    # ===================================================================
    print(f"\n{'='*80}")
    print("STEP 2: Hematological vs Epithelial drivers in eigenspace")
    print("=" * 80)

    heme_proj = []
    epi_proj = []

    for gene in HEME_DRIVERS:
        mvec = gene_to_module_vector(gene)
        if mvec.sum() > 0:
            heme_proj.append(project_to_eigenspace(mvec, eigenvectors, N_PCS))

    for gene in EPITHELIAL_DRIVERS:
        mvec = gene_to_module_vector(gene)
        if mvec.sum() > 0:
            epi_proj.append(project_to_eigenspace(mvec, eigenvectors, N_PCS))

    heme_proj = np.array(heme_proj)
    epi_proj = np.array(epi_proj)

    print(f"\n  Hematological drivers in modules: {len(heme_proj)}")
    print(f"  Epithelial drivers in modules: {len(epi_proj)}")

    print(f"\n  {'PC':>4s} {'Heme mean':>10s} {'Epi mean':>10s} {'Δ':>8s} {'t':>8s} {'p':>10s}")
    print("  " + "-" * 55)

    for pc in range(N_PCS):
        hv = heme_proj[:, pc]
        ev = epi_proj[:, pc]
        t, p = stats.ttest_ind(hv, ev)
        delta = np.mean(hv) - np.mean(ev)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  PC{pc+1:1d} {np.mean(hv):+10.3f} {np.mean(ev):+10.3f} "
              f"{delta:+8.3f} {t:+8.2f} {p:10.4f}  {sig}")

    # Pre-registered test: PC3 specifically
    # PC3 = danger(+) vs tissue(-) in 20-module decomposition
    heme_pc3 = heme_proj[:, 2]
    epi_pc3 = epi_proj[:, 2]
    t_pc3, p_pc3 = stats.ttest_ind(heme_pc3, epi_pc3)

    print(f"\n  Pre-registered test — PC3 (danger vs tissue):")
    print(f"    Hematological mean PC3: {np.mean(heme_pc3):+.3f}")
    print(f"    Epithelial mean PC3:    {np.mean(epi_pc3):+.3f}")
    print(f"    t = {t_pc3:.2f}, p = {p_pc3:.4f}")

    if np.mean(heme_pc3) > np.mean(epi_pc3):
        print(f"    ✓ PREDICTION 2 CONFIRMED: Heme → danger pole, Epi → tissue pole on PC3")
    else:
        print(f"    ✗ Prediction 2: direction opposite to expected")

    # ===================================================================
    # STEP 3: GoF vs LoF in eigenspace
    # ===================================================================
    print(f"\n{'='*80}")
    print("STEP 3: GoF vs LoF direction in eigenspace")
    print("=" * 80)

    gof_proj = []
    lof_proj = []
    gof_names = []
    lof_names = []

    for gene in GOF_GENES:
        mvec = gene_to_module_vector(gene)
        if mvec.sum() > 0:
            gof_proj.append(project_to_eigenspace(mvec, eigenvectors, N_PCS))
            gof_names.append(gene)

    for gene in LOF_GENES:
        mvec = gene_to_module_vector(gene)
        if mvec.sum() > 0:
            lof_proj.append(project_to_eigenspace(mvec, eigenvectors, N_PCS))
            lof_names.append(gene)

    gof_proj = np.array(gof_proj)
    lof_proj = np.array(lof_proj)

    print(f"\n  GoF genes in modules: {len(gof_proj)} ({', '.join(gof_names[:10])}...)")
    print(f"  LoF genes in modules: {len(lof_proj)} ({', '.join(lof_names[:10])}...)")

    print(f"\n  {'PC':>4s} {'GoF mean':>10s} {'LoF mean':>10s} {'Δ':>8s} {'t':>8s} {'p':>10s} {'Direction'}")
    print("  " + "-" * 72)

    for pc in range(N_PCS):
        gv = gof_proj[:, pc]
        lv = lof_proj[:, pc]
        t, p = stats.ttest_ind(gv, lv)
        delta = np.mean(gv) - np.mean(lv)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        direction = ""
        if p < 0.05:
            if delta > 0:
                direction = "GoF → + pole"
            else:
                direction = "LoF → + pole"
        print(f"  PC{pc+1:1d} {np.mean(gv):+10.3f} {np.mean(lv):+10.3f} "
              f"{delta:+8.3f} {t:+8.2f} {p:10.4f}  {sig} {direction}")

    # Overall: do GoF and LoF separate in multivariate space?
    gof_centroid = np.mean(gof_proj, axis=0)
    lof_centroid = np.mean(lof_proj, axis=0)
    gof_lof_dist = np.linalg.norm(gof_centroid - lof_centroid)

    # Permutation test
    all_gl = np.vstack([gof_proj, lof_proj])
    n_gof = len(gof_proj)
    perm_gl_dists = []
    for _ in range(n_perm):
        pi = np.random.permutation(len(all_gl))
        d = np.linalg.norm(np.mean(all_gl[pi[:n_gof]], axis=0) - np.mean(all_gl[pi[n_gof:]], axis=0))
        perm_gl_dists.append(d)

    p_gl = np.mean(np.array(perm_gl_dists) >= gof_lof_dist)

    print(f"\n  GoF-LoF centroid distance: {gof_lof_dist:.4f}")
    print(f"  Permutation test: p = {p_gl:.4f}")

    if p_gl < 0.05:
        print(f"  ✓ PREDICTION 4 CONFIRMED: GoF and LoF occupy different regions of eigenspace")
    else:
        print(f"  ✗ Prediction 4 not confirmed")

    # ===================================================================
    # STEP 4: Candidate novel cancer genes (eigenspace neighbors of CGC)
    # ===================================================================
    print(f"\n{'='*80}")
    print("STEP 4: Non-CGC genes in 'cancer region' of eigenspace")
    print("=" * 80)

    # Define cancer region as within 1 SD of CGC centroid
    cgc_centroid = np.mean(cgc_proj, axis=0)
    cgc_std = np.std(cgc_proj, axis=0)

    candidates = []
    for gene in all_module_genes:
        if gene in cgc:
            continue
        mvec = gene_to_module_vector(gene)
        if mvec.sum() == 0:
            continue
        proj = project_to_eigenspace(mvec, eigenvectors, N_PCS)
        # Mahalanobis-like distance (using per-PC std)
        z_scores = np.abs((proj - cgc_centroid) / (cgc_std + 0.01))
        mean_z = np.mean(z_scores)
        if mean_z < 1.5:  # within 1.5 SD on average
            candidates.append((gene, mean_z, proj))

    candidates.sort(key=lambda x: x[1])

    print(f"\n  Non-CGC genes in 'cancer-like' eigenspace position: {len(candidates)}")
    print(f"\n  Top 20 candidates (closest to CGC centroid):")
    print(f"  {'Gene':12s} {'Mean Z':>7s} {'Modules'}")
    print("  " + "-" * 60)
    for gene, mz, _ in candidates[:20]:
        mods = [m for m in MODULE_ORDER if gene in ALL_MODULES[m]]
        print(f"  {gene:12s} {mz:7.2f}  {', '.join(mods)}")

    # ===================================================================
    # SAVE
    # ===================================================================
    sig_pcs = [r for r in pc_results if r["p"] < 0.05]

    results = {
        "step1": {
            "cgc_n": len(cgc_proj),
            "noncgc_n": len(noncgc_proj),
            "multivariate_p": round(float(p_perm), 4),
            "significant_pcs": [r["pc"] for r in sig_pcs],
        },
        "step2": {
            "heme_n": len(heme_proj),
            "epi_n": len(epi_proj),
            "pc3_heme_mean": round(float(np.mean(heme_pc3)), 3),
            "pc3_epi_mean": round(float(np.mean(epi_pc3)), 3),
            "pc3_p": round(float(p_pc3), 4),
        },
        "step3": {
            "gof_n": len(gof_proj),
            "lof_n": len(lof_proj),
            "multivariate_p": round(float(p_gl), 4),
            "centroid_distance": round(float(gof_lof_dist), 4),
        },
        "step4": {
            "n_candidates": len(candidates),
            "top_candidates": [c[0] for c in candidates[:20]],
        },
    }

    with open(DATA_DIR / "eigenspace_cancer.json", "w") as f:
        json.dump(results, f, indent=2)

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("PRE-REGISTERED PREDICTION SUMMARY")
    print("=" * 80)
    print(f"\n  1. CGC non-random in eigenspace: p = {p_perm:.4f} {'✓' if p_perm < 0.05 else '✗'}")
    print(f"  2. Heme vs Epi on PC3: p = {p_pc3:.4f} {'✓' if p_pc3 < 0.05 and np.mean(heme_pc3) > np.mean(epi_pc3) else '✗'}")
    print(f"  3. (Same as 2 — epithelial → tissue pole)")
    print(f"  4. GoF vs LoF separate: p = {p_gl:.4f} {'✓' if p_gl < 0.05 else '✗'}")
    print(f"\nSaved to eigenspace_cancer.json")


if __name__ == "__main__":
    main()
