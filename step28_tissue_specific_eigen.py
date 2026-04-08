#!/usr/bin/env python3
"""
step28: Tissue-specific eigendecomposition of 20 oscillatory modules.

Question: Are the 6 Kaiser eigenmodes universal across tissue classes,
or are they tissue-specific?

Design:
- HPA data: 154 cell types × 20 modules
- Split cell types by HPA class (blood/immune, epithelial, neuronal, etc.)
- Compute eigendecomposition separately per class
- Compare: are PC axes similar? Do the same module pairs co-vary?
- Use Procrustes / axis correlation to quantify similarity

Biological meaning:
  - Universal axes → cancer vulnerability is a universal architectural feature
  - Tissue-specific axes → vulnerability landscape is remodeled per tissue
  - Intermediate → some axes universal (core), some tissue-specific (context)
"""

import csv
import json
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

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

MODULE_ORDER = [
    "NF-κB","ERK/MAPK","JAK-STAT","p53","Wnt","Notch",
    "Hippo","TGF-β","mTOR","Calcium","Cell Cycle",
    "Circadian","NRF2","PI3K/PTEN",
    "AMPK","SREBP","ATR/CHK1","Rho/ROCK","PPAR/LXR","Autophagy",
]

# Tissue class groupings — merge small classes for stability
TISSUE_CLASSES = {
    "blood_immune":    ["blood and immune cells"],
    "epithelial":      ["glandular epithelial cells", "specialized epithelial cells",
                        "squamous epithelial cells", "ciliated cells"],
    "neuronal":        ["neuronal cells", "glial cells"],
    "mesenchymal":     ["mesenchymal cells", "muscle cells"],
    "endocrine":       ["endocrine cells"],
    "stem_germ":       ["stem and proliferating cells", "germ cells"],
    "vascular":        ["endothelial and mural cells", "trophoblast cells", "pigment cells"],
}


def load_cell_type_classes():
    """Map cell type → tissue class."""
    ct_to_hpa_class = {}
    with open(DATA_DIR / "rna_single_cell_type_cell_types.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ct_to_hpa_class[row["Cell type"]] = row["Cell type class"]
    # Invert tissue classes
    hpa_class_to_tissue = {}
    for tissue, hpa_classes in TISSUE_CLASSES.items():
        for hc in hpa_classes:
            hpa_class_to_tissue[hc] = tissue
    ct_to_tissue = {}
    for ct, hc in ct_to_hpa_class.items():
        tissue = hpa_class_to_tissue.get(hc, "other")
        ct_to_tissue[ct] = tissue
    return ct_to_tissue


def load_expression():
    """Load HPA nCPM data."""
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    return expr


def build_module_matrix(cell_types, expr):
    """Build 20 × N_celltypes module expression matrix."""
    N = len(MODULE_ORDER)
    M = len(cell_types)
    mat = np.zeros((N, M))
    for i, mod in enumerate(MODULE_ORDER):
        genes = [g for g in ALL_MODULES[mod] if g in expr]
        for j, ct in enumerate(cell_types):
            vals = [expr[g].get(ct, 0) for g in genes]
            mat[i, j] = np.mean(vals) if vals else 0.0
    return mat


def eigen_decompose(mat):
    """Correlation matrix → eigenvectors (sorted by descending eigenvalue)."""
    if mat.shape[1] < 4:
        return None, None, None
    corr = np.corrcoef(mat)
    # Handle NaN (constant rows)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    eigenvalues, eigenvectors = np.linalg.eigh(corr)
    idx = np.argsort(eigenvalues)[::-1]
    return eigenvalues[idx], eigenvectors[:, idx], corr


def kaiser_n(eigenvalues):
    """Number of Kaiser components (eigenvalue > 1)."""
    return int(np.sum(eigenvalues > 1))


def axis_similarity(v1, v2):
    """Cosine similarity between two eigenvectors (absolute, sign-invariant)."""
    return abs(float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)))


def top_modules(vec, n=4):
    """Return top n modules by |loading| with sign."""
    top = np.argsort(np.abs(vec))[::-1][:n]
    return [(MODULE_ORDER[i], float(vec[i])) for i in top]


def main():
    print("Loading HPA expression data...")
    expr = load_expression()
    ct_to_tissue = load_cell_type_classes()

    # Collect cell types per tissue class
    tissue_cts = defaultdict(list)
    for ct, tissue in ct_to_tissue.items():
        tissue_cts[tissue].append(ct)

    print(f"\nTissue class cell type counts:")
    for tissue, cts in sorted(tissue_cts.items()):
        print(f"  {tissue:20s}: {len(cts):3d} cell types")

    # Global eigendecomposition (all 154 cell types — reference)
    print("\n" + "="*80)
    print("GLOBAL REFERENCE (all 154 cell types)")
    print("="*80)
    all_cts = sorted(ct_to_tissue.keys())
    mat_global = build_module_matrix(all_cts, expr)
    ev_global, vecs_global, corr_global = eigen_decompose(mat_global)
    n_kaiser_global = kaiser_n(ev_global)
    var_explained = ev_global / ev_global.sum() * 100

    print(f"\n  Kaiser components: {n_kaiser_global}")
    print(f"  Variance explained (PC1-6): {var_explained[:6]}")
    print(f"\n  Global PC axes:")
    for pc in range(min(6, n_kaiser_global)):
        tm = top_modules(vecs_global[:, pc], n=3)
        print(f"    PC{pc+1} ({var_explained[pc]:.1f}%): "
              + ", ".join(f"{m}({v:+.2f})" for m, v in tm))

    # Per-tissue eigendecomposition
    tissue_results = {}
    print("\n" + "="*80)
    print("TISSUE-SPECIFIC EIGENDECOMPOSITION")
    print("="*80)

    for tissue in sorted(tissue_cts.keys()):
        cts = tissue_cts[tissue]
        if len(cts) < 4:
            print(f"\n  {tissue}: SKIPPED (N={len(cts)} < 4)")
            continue

        mat = build_module_matrix(cts, expr)
        ev, vecs, corr = eigen_decompose(mat)
        if ev is None:
            continue

        n_kaiser = kaiser_n(ev)
        var = ev / ev.sum() * 100

        print(f"\n  {'─'*60}")
        print(f"  {tissue.upper()} (N={len(cts)} cell types, Kaiser={n_kaiser})")
        print(f"  {'─'*60}")
        print(f"  Variance: {var[:6]}")

        # Top axes
        tissue_pcs = []
        for pc in range(min(6, max(n_kaiser, 3))):
            tm = top_modules(vecs[:, pc], n=3)
            sim_global = axis_similarity(vecs[:, pc], vecs_global[:, pc])
            print(f"    PC{pc+1} ({var[pc]:.1f}%): "
                  + ", ".join(f"{m}({v:+.2f})" for m, v in tm)
                  + f"  | sim_global={sim_global:.3f}")
            tissue_pcs.append({
                "var": round(float(var[pc]), 1),
                "top_modules": [(m, round(v, 3)) for m, v in tm],
                "sim_global": round(sim_global, 3),
            })

        tissue_results[tissue] = {
            "n_cell_types": len(cts),
            "n_kaiser": n_kaiser,
            "variance_pct": [round(float(v), 1) for v in var[:6]],
            "pcs": tissue_pcs,
        }

    # ===================================================================
    # Cross-tissue PC1 similarity matrix
    # ===================================================================
    print("\n" + "="*80)
    print("PC1 SIMILARITY ACROSS TISSUES (axis conservation)")
    print("="*80)

    tissues_with_results = [t for t in sorted(tissue_cts.keys())
                             if t in tissue_results]
    vecs_per_tissue = {}
    for tissue in tissues_with_results:
        cts = tissue_cts[tissue]
        mat = build_module_matrix(cts, expr)
        ev, vecs, _ = eigen_decompose(mat)
        if vecs is not None:
            vecs_per_tissue[tissue] = vecs

    print(f"\n  PC1 cosine similarity matrix:")
    header = f"  {'':20s}" + "".join(f" {t[:10]:>10s}" for t in vecs_per_tissue)
    print(header)
    for t1, v1 in vecs_per_tissue.items():
        row = f"  {t1[:20]:20s}"
        for t2, v2 in vecs_per_tissue.items():
            sim = axis_similarity(v1[:, 0], v2[:, 0])
            row += f" {sim:10.3f}"
        print(row)

    # PC1 vs global similarity summary
    print(f"\n  PC1 similarity to global reference:")
    for tissue, vecs in vecs_per_tissue.items():
        sim = axis_similarity(vecs[:, 0], vecs_global[:, 0])
        stability = "★★★" if sim > 0.85 else ("★★" if sim > 0.70 else "★")
        print(f"    {tissue:20s}: {sim:.3f} {stability}")

    # ===================================================================
    # Module co-variation stability: do specific pairs always co-vary?
    # ===================================================================
    print("\n" + "="*80)
    print("MODULE PAIR CO-VARIATION STABILITY")
    print("="*80)

    # Key pairs from global analysis
    key_pairs = [
        ("ERK/MAPK", "Wnt"),
        ("ERK/MAPK", "Notch"),
        ("p53", "Cell Cycle"),
        ("mTOR", "AMPK"),
        ("JAK-STAT", "NF-κB"),
        ("Hippo", "TGF-β"),
        ("SREBP", "PPAR/LXR"),
        ("Circadian", "ATR/CHK1"),
    ]

    print(f"\n  {'Module pair':35s} {'Global':>7s}", end="")
    for t in vecs_per_tissue:
        print(f" {t[:8]:>8s}", end="")
    print()
    print("  " + "-" * (35 + 7 + 8*len(vecs_per_tissue) + 5))

    # Compute global correlations first
    global_corrs = {}
    for (m1, m2) in key_pairs:
        i1 = MODULE_ORDER.index(m1)
        i2 = MODULE_ORDER.index(m2)
        global_corrs[(m1, m2)] = corr_global[i1, i2]

    pair_data = {}
    for (m1, m2) in key_pairs:
        i1 = MODULE_ORDER.index(m1)
        i2 = MODULE_ORDER.index(m2)
        tissue_cors = {}
        for tissue, cts in tissue_cts.items():
            if tissue not in vecs_per_tissue:
                continue
            mat = build_module_matrix(cts, expr)
            if mat.shape[1] < 4:
                continue
            corr_t = np.corrcoef(mat)
            tissue_cors[tissue] = corr_t[i1, i2]
        pair_data[(m1, m2)] = tissue_cors

    for (m1, m2) in key_pairs:
        pair_str = f"{m1} ~ {m2}"
        global_r = global_corrs[(m1, m2)]
        print(f"  {pair_str:35s} {global_r:+7.3f}", end="")
        for tissue in vecs_per_tissue:
            r = pair_data.get((m1, m2), {}).get(tissue, float('nan'))
            print(f" {r:+8.3f}", end="")
        print()

    # ===================================================================
    # SUMMARY CONCLUSION
    # ===================================================================
    print("\n" + "="*80)
    print("CONCLUSION: UNIVERSAL vs TISSUE-SPECIFIC")
    print("="*80)

    # Average PC1 similarity across tissues
    pc1_sims = [axis_similarity(v[:, 0], vecs_global[:, 0])
                for v in vecs_per_tissue.values()]
    mean_pc1_sim = np.mean(pc1_sims) if pc1_sims else 0

    print(f"\n  Mean PC1 similarity to global: {mean_pc1_sim:.3f}")
    print(f"  Range: {min(pc1_sims):.3f} – {max(pc1_sims):.3f}")

    if mean_pc1_sim > 0.80:
        verdict = "HIGHLY UNIVERSAL — core oscillatory axes conserved across tissues"
    elif mean_pc1_sim > 0.60:
        verdict = "PARTIALLY UNIVERSAL — axis conserved in some tissues, remodeled in others"
    else:
        verdict = "TISSUE-SPECIFIC — oscillatory architecture is tissue-dependent"

    print(f"\n  Verdict: {verdict}")

    # ===================================================================
    # SAVE
    # ===================================================================
    results = {
        "global": {
            "n_cell_types": len(all_cts),
            "n_kaiser": n_kaiser_global,
            "variance_pct": [round(float(v), 1) for v in var_explained[:6]],
        },
        "tissues": tissue_results,
        "pc1_mean_similarity": round(float(mean_pc1_sim), 3),
        "pc1_sim_range": [round(float(min(pc1_sims)), 3), round(float(max(pc1_sims)), 3)],
        "verdict": verdict,
        "pair_correlations": {
            f"{m1}_{m2}": {
                "global": round(float(global_corrs[(m1, m2)]), 3),
                "by_tissue": {t: round(float(r), 3) for t, r in pair_data.get((m1, m2), {}).items()},
            }
            for (m1, m2) in key_pairs
        },
    }

    with open(DATA_DIR / "tissue_specific_eigen.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to tissue_specific_eigen.json")


if __name__ == "__main__":
    main()
