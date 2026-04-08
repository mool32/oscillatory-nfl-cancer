#!/usr/bin/env python3
"""
CORRECTED CGC Enrichment Analysis — using full COSMIC CGC (740 genes).
Replaces step06 results which used incomplete (~200) CGC set.

Test 1: NFL genes vs whole genome (hypergeometric)
Test 2: NFL genes vs all KEGG signaling genes (hypergeometric)
Test 3: NFL genes vs non-NFL KEGG signaling genes (Fisher exact) ← CRITICAL
Test 4: Overlap with Paper 1
"""

import csv
import json
import numpy as np
from scipy import stats
from pathlib import Path
from collections import Counter

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

N_GENOME = 20000  # protein-coding genes


def load_cgc():
    """Load full COSMIC CGC (740 genes)."""
    genes = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    print(f"Loaded {len(genes)} CGC genes")
    return genes


def load_nfl_genes():
    """Load NFL motif genes."""
    genes = set()
    with open(DATA_DIR / "unique_nfl_motifs.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for g in row["genes_hgnc"].split("|"):
                genes.add(g)
    return genes


def load_kegg_background():
    """Load all KEGG signaling genes."""
    genes = set()
    with open(DATA_DIR / "kegg_all_signaling_genes.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    return genes


# Paper 1 genes (157 across 14 pathways)
PAPER1_GENES = {
    "RELA", "RELB", "NFKB1", "NFKB2", "REL", "IKBKB", "CHUK",
    "NFKBIA", "NFKBIB", "TNFAIP3", "CYLD", "OTULIN",
    "TP53", "ATM", "ATR", "CHEK1", "CHEK2", "MDM2", "MDM4", "PPM1D", "USP7",
    "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2",
    "MAPK1", "MAPK3", "DUSP1", "DUSP2", "DUSP3", "DUSP4", "DUSP5", "DUSP6",
    "SPRY1", "SPRY2", "NF1",
    "CTNNB1", "WNT1", "WNT3A", "LRP5", "LRP6", "DVL1", "DVL2",
    "APC", "AXIN1", "AXIN2", "GSK3B", "RNF43", "ZNRF3",
    "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "JAG1",
    "FBXW7", "NUMB", "NUMBL", "ITCH", "SEL1L", "MAML1",
    "CLOCK", "ARNTL", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
    "CSNK1D", "CSNK1E",
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
    "LATS1", "LATS2", "STK3", "STK4", "SAV1", "NF2",
    "MTOR", "RPTOR", "RHEB", "AKT1", "AKT2", "AKT3",
    "PIK3CA", "PIK3CB", "PIK3R1", "TSC1", "TSC2", "PTEN", "STK11", "DEPDC5",
    "JAK1", "JAK2", "JAK3", "TYK2", "STAT3", "STAT5A", "STAT5B",
    "SOCS1", "SOCS3", "PTPN1", "PTPN2", "PTPN6", "PIAS1", "SOCS2",
    "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD6", "SMAD7",
    "SMURF1", "SMURF2", "SKI", "SKIL",
    "NFE2L2", "MAF", "MAFG", "KEAP1", "CUL3", "RBX1", "BACH1",
    "SHH", "IHH", "SMO", "GLI1", "GLI2", "PTCH1", "PTCH2", "SUFU", "GPR161",
    "ERN1", "XBP1", "ATF6", "HSPA5", "EDEM1", "DNAJB9",
    "ITPR1", "ITPR2", "ITPR3", "PLCG1", "PLCG2", "ATP2A1", "ATP2A2", "CALM1",
}


def hypergeometric_test(k, n, K, N):
    """P(X >= k) in hypergeometric."""
    p = stats.hypergeom.sf(k - 1, N, K, n)
    expected = n * K / N
    fold = k / expected if expected > 0 else float('inf')
    return p, fold, expected


def main():
    cgc = load_cgc()
    nfl_genes = load_nfl_genes()
    kegg_genes = load_kegg_background()

    nfl_cgc = nfl_genes & cgc
    kegg_cgc = kegg_genes & cgc
    non_nfl_kegg = kegg_genes - nfl_genes
    non_nfl_cgc = non_nfl_kegg & cgc

    print("=" * 70)
    print("CORRECTED CGC ENRICHMENT — Full COSMIC CGC (740 genes)")
    print("=" * 70)
    print(f"\nNFL motif genes: {len(nfl_genes)}")
    print(f"NFL ∩ CGC: {len(nfl_cgc)} ({100*len(nfl_cgc)/len(nfl_genes):.1f}%)")
    print(f"NFL CGC genes: {', '.join(sorted(nfl_cgc))}")
    print(f"\nKEGG signaling genes: {len(kegg_genes)}")
    print(f"KEGG ∩ CGC: {len(kegg_cgc)} ({100*len(kegg_cgc)/len(kegg_genes):.1f}%)")
    print(f"\nNon-NFL KEGG genes: {len(non_nfl_kegg)}")
    print(f"Non-NFL ∩ CGC: {len(non_nfl_cgc)} ({100*len(non_nfl_cgc)/len(non_nfl_kegg):.1f}%)")

    # ===================================================================
    # TEST 1: NFL vs genome
    # ===================================================================
    print(f"\n{'='*70}")
    print("TEST 1: NFL genes vs Whole Genome")
    print("=" * 70)
    k1 = len(nfl_cgc)
    n1 = len(nfl_genes)
    K1 = len(cgc)  # 740
    N1 = N_GENOME
    p1, fold1, exp1 = hypergeometric_test(k1, n1, K1, N1)
    print(f"  {k1} CGC in {n1} NFL genes")
    print(f"  Background: {K1} CGC in {N1} genome = {100*K1/N1:.1f}%")
    print(f"  Expected: {exp1:.1f}")
    print(f"  Fold enrichment: {fold1:.1f}×")
    print(f"  p = {p1:.2e}")

    # ===================================================================
    # TEST 2: NFL vs KEGG signaling
    # ===================================================================
    print(f"\n{'='*70}")
    print("TEST 2: NFL genes vs All KEGG Signaling Genes")
    print("=" * 70)
    k2 = len(nfl_cgc)
    n2 = len(nfl_genes)
    K2 = len(kegg_cgc)
    N2 = len(kegg_genes)
    p2, fold2, exp2 = hypergeometric_test(k2, n2, K2, N2)
    print(f"  {k2} CGC in {n2} NFL genes")
    print(f"  Background: {K2} CGC in {N2} KEGG signaling = {100*K2/N2:.1f}%")
    print(f"  Expected: {exp2:.1f}")
    print(f"  Fold enrichment: {fold2:.1f}×")
    print(f"  p = {p2:.2e}")

    # ===================================================================
    # TEST 3: NFL vs non-NFL in KEGG (Fisher exact) ← CRITICAL
    # ===================================================================
    print(f"\n{'='*70}")
    print("TEST 3: NFL vs Non-NFL KEGG Signaling (Fisher Exact) ← CRITICAL")
    print("=" * 70)

    nfl_non_cgc = nfl_genes - cgc
    non_nfl_non_cgc = non_nfl_kegg - cgc

    table3 = [
        [len(nfl_cgc), len(nfl_non_cgc)],
        [len(non_nfl_cgc), len(non_nfl_non_cgc)]
    ]
    OR3, p3 = stats.fisher_exact(table3)

    print(f"  NFL: {len(nfl_cgc)} CGC / {len(nfl_non_cgc)} non-CGC = {100*len(nfl_cgc)/len(nfl_genes):.1f}% CGC")
    print(f"  Non-NFL KEGG: {len(non_nfl_cgc)} CGC / {len(non_nfl_non_cgc)} non-CGC = {100*len(non_nfl_cgc)/len(non_nfl_kegg):.1f}% CGC")
    print(f"\n  2×2 table:")
    print(f"              CGC    non-CGC")
    print(f"    NFL:      {table3[0][0]:5d}    {table3[0][1]:5d}")
    print(f"    Non-NFL:  {table3[1][0]:5d}    {table3[1][1]:5d}")
    print(f"\n  Fisher exact OR = {OR3:.2f}")
    print(f"  p = {p3:.2e}")

    # Chi-square for comparison
    chi2, p_chi, _, _ = stats.chi2_contingency(table3)
    print(f"  Chi-square: χ² = {chi2:.1f}, p = {p_chi:.2e}")

    if p3 < 0.05:
        print(f"\n  → SIGNIFICANT: NFL membership specifically enriches for cancer")
    else:
        print(f"\n  → NOT significant: enrichment driven by pathway membership")

    # ===================================================================
    # TEST 4: Paper 1 overlap
    # ===================================================================
    print(f"\n{'='*70}")
    print("TEST 4: Paper 1 Overlap")
    print("=" * 70)

    paper1_cgc = PAPER1_GENES & cgc

    both_cgc = nfl_cgc & paper1_cgc
    nfl_only_cgc = nfl_cgc - paper1_cgc
    paper1_only_cgc = paper1_cgc - nfl_cgc

    print(f"\n  Paper 1: {len(PAPER1_GENES)} genes, {len(paper1_cgc)} CGC ({100*len(paper1_cgc)/len(PAPER1_GENES):.1f}%)")
    print(f"  Paper 2: {len(nfl_genes)} genes, {len(nfl_cgc)} CGC ({100*len(nfl_cgc)/len(nfl_genes):.1f}%)")

    p1_fold = (len(paper1_cgc) / len(PAPER1_GENES)) / (len(cgc) / N_GENOME)
    p2_fold = (len(nfl_cgc) / len(nfl_genes)) / (len(cgc) / N_GENOME)

    print(f"\n  Paper 1 vs genome: {p1_fold:.1f}×")
    print(f"  Paper 2 vs genome: {p2_fold:.1f}×")

    print(f"\n  CGC Venn:")
    print(f"    Shared: {len(both_cgc)}")
    print(f"    NFL-only: {len(nfl_only_cgc)} — {', '.join(sorted(nfl_only_cgc))}")
    print(f"    Paper1-only: {len(paper1_only_cgc)}")

    print(f"\n  Gene overlap (all): {len(nfl_genes & PAPER1_GENES)}/{len(nfl_genes)} NFL in Paper1 ({100*len(nfl_genes & PAPER1_GENES)/len(nfl_genes):.0f}%)")
    print(f"  New NFL genes: {len(nfl_genes - PAPER1_GENES)}")

    # ===================================================================
    # DECISION: Does Paper 2 stand alone?
    # ===================================================================
    print(f"\n{'='*70}")
    print("DECISION SUMMARY")
    print("=" * 70)
    print(f"\n  Test 1 (vs genome):      {fold1:.1f}× enrichment, p = {p1:.2e}")
    print(f"  Test 2 (vs KEGG):        {fold2:.1f}× enrichment, p = {p2:.2e}")
    print(f"  Test 3 (NFL vs non-NFL): OR = {OR3:.1f}, p = {p3:.2e}")
    print(f"  Paper 1 convergence:     {p1_fold:.1f}× vs {p2_fold:.1f}×")

    if OR3 > 5 and p3 < 0.001:
        print(f"\n  → STRONG: NFL membership specifically predicts cancer (OR={OR3:.1f})")
        print(f"    Paper 2 stands alone.")
    elif OR3 > 2 and p3 < 0.05:
        print(f"\n  → MODERATE: NFL enrichment beyond pathway membership (OR={OR3:.1f})")
        print(f"    Paper 2 viable but weaker standalone case.")
    else:
        print(f"\n  → WEAK: NFL enrichment ≈ pathway membership (OR={OR3:.1f})")
        print(f"    Consider merging key findings into Paper 1.")

    # Save
    results = {
        "cgc_version": "COSMIC_CGC_v103_740genes",
        "nfl_genes": len(nfl_genes),
        "nfl_cgc": len(nfl_cgc),
        "nfl_cgc_pct": round(100 * len(nfl_cgc) / len(nfl_genes), 1),
        "kegg_genes": len(kegg_genes),
        "kegg_cgc": len(kegg_cgc),
        "kegg_cgc_pct": round(100 * len(kegg_cgc) / len(kegg_genes), 1),
        "non_nfl_cgc": len(non_nfl_cgc),
        "non_nfl_cgc_pct": round(100 * len(non_nfl_cgc) / len(non_nfl_kegg), 1),
        "test1_fold": round(fold1, 1),
        "test1_p": float(p1),
        "test2_fold": round(fold2, 1),
        "test2_p": float(p2),
        "test3_OR": round(float(OR3), 2),
        "test3_p": float(p3),
        "test3_table": table3,
        "paper1_fold": round(p1_fold, 1),
        "paper2_fold": round(p2_fold, 1),
        "venn_both": len(both_cgc),
        "venn_nfl_only": len(nfl_only_cgc),
        "venn_paper1_only": len(paper1_only_cgc),
    }
    with open(DATA_DIR / "cgc_enrichment_CORRECTED.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {DATA_DIR / 'cgc_enrichment_CORRECTED.json'}")


if __name__ == "__main__":
    main()
