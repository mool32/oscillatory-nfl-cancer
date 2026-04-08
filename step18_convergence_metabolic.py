#!/usr/bin/env python3
"""
Phase C.2: Input convergence → cancer vulnerability (test #44)
Phase E: Metabolic regulatory vs enzymatic cancer enrichment (test #45)

C.2: For each module — count upstream receptors from KEGG-derived mapping.
     Test: more receptors → more cancer genes?

E: Classify metabolic genes as regulatory (with feedback) vs enzymatic (without).
   Test: regulatory metabolic genes more CGC-enriched than enzymatic?
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")


def load_cgc():
    cgc = set()
    # Try plain text file (one gene per line)
    cgc_path = DATA_DIR / "cosmic_cgc_full.txt"
    if cgc_path.exists():
        with open(cgc_path) as f:
            for line in f:
                gene = line.strip()
                if gene:
                    cgc.add(gene)
    # Fallback to CSV
    if not cgc:
        cgc_path = DATA_DIR / "cancer_gene_census.csv"
        if cgc_path.exists():
            with open(cgc_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    name_field = "Gene Symbol" if "Gene Symbol" in row else list(row.keys())[0]
                    cgc.add(row[name_field].strip())
    print(f"  Loaded {len(cgc)} CGC genes")
    return cgc


# ===================================================================
# MODULE DATA: all 20 modules with receptor counts and gene lists
# ===================================================================

# Core + extended gene lists per module (for CGC density calculation)
MODULE_GENES = {
    "NF-κB": {"RELA", "RELB", "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "TNFAIP3",
              "TRAF2", "TRAF6", "IKBKB", "IKBKG", "CHUK", "MAP3K7", "TAB1", "TAB2"},
    "ERK/MAPK": {"MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BRAF", "RAF1", "ARAF",
                 "HRAS", "KRAS", "NRAS", "SOS1", "GRB2", "DUSP1", "DUSP6", "SPRY2"},
    "JAK-STAT": {"JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT3", "STAT5A", "STAT5B",
                 "SOCS1", "SOCS3", "CISH", "PIAS1"},
    "p53": {"TP53", "MDM2", "MDM4", "CDKN1A", "BAX", "BBC3", "PMAIP1",
            "ATM", "ATR", "CHEK1", "CHEK2"},
    "Wnt": {"CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B", "CK1", "DVL1",
            "TCF7L2", "LEF1", "WNT3A", "RNF43", "ZNRF3"},
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
    # New modules
    "AMPK": {"PRKAA1", "PRKAA2", "PRKAB1", "PRKAG1", "STK11", "ACACB",
             "PPARGC1A", "NUAK1"},
    "SREBP": {"SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2", "HMGCR",
              "FASN", "SCD", "ACLY"},
    "ATR/CHK1": {"ATR", "ATRIP", "CHEK1", "TOPBP1", "WEE1", "CDC25A",
                 "CLSPN", "CLASPIN", "RAD17"},
    "Rho/ROCK": {"RHOA", "RHOB", "RHOC", "ROCK1", "ROCK2", "MKL1",
                 "LIMK1", "CFL1", "ARHGAP1"},
    "PPAR/LXR": {"PPARA", "PPARG", "PPARD", "NR1H3", "NR1H2", "RXRA",
                 "NCOR1", "NCOR2"},
    "Autophagy": {"TFEB", "ULK1", "ULK2", "BECN1", "ATG5", "ATG7", "ATG12",
                  "SQSTM1", "MAP1LC3B"},
}

# Upstream receptor count per module (curated from KEGG + literature)
MODULE_RECEPTOR_COUNT = {
    # Original 14
    "NF-κB": {"receptors": ["TNFRSF1A", "TNFRSF1B", "TLR1", "TLR2", "TLR3", "TLR4",
                             "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "IL1R1", "IL1R2",
                             "IL18R1", "NOD1", "NOD2", "RIPK1", "RIPK2", "CD40",
                             "LTBR", "RANK"], "count": 21},
    "ERK/MAPK": {"receptors": ["EGFR", "ERBB2", "ERBB3", "FGFR1", "FGFR2", "FGFR3",
                                "PDGFRA", "PDGFRB", "KIT", "MET", "RET", "ALK",
                                "NTRK1", "NTRK2", "IGF1R", "INSR", "FLT3"], "count": 17},
    "JAK-STAT": {"receptors": ["IL6R", "IL6ST", "IL2RA", "IL2RB", "IL2RG", "IL7R",
                                "IL10RA", "IL12RB1", "IL12RB2", "IL21R", "IL23R",
                                "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
                                "EPOR", "TPOR", "CSF2RA", "CSF3R", "LEPR",
                                "GHR", "PRLR"], "count": 22},
    "p53": {"receptors": ["ATM", "ATR", "PARP1", "H2AFX", "DNAPK"], "count": 5},
    "Wnt": {"receptors": ["FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7",
                           "FZD8", "FZD9", "FZD10", "LRP5", "LRP6", "ROR1", "ROR2",
                           "RYK"], "count": 15},
    "Notch": {"receptors": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"], "count": 4},
    "Hippo": {"receptors": ["CDH1", "CDH2", "CDH3", "NF2", "FAT4", "DCHS1",
                             "AMOT", "GPR87", "PIEZO1"], "count": 9},
    "TGF-β": {"receptors": ["TGFBR1", "TGFBR2", "BMPR1A", "BMPR1B", "BMPR2",
                              "ACVR1", "ACVR1B", "ACVR2A", "ACVR2B"], "count": 9},
    "mTOR": {"receptors": ["SLC1A5", "SLC7A5", "SLC3A2", "SLC38A2", "INSR",
                            "IGF1R", "SESN2", "CASTOR1", "FLCN", "LAMTOR1"], "count": 10},
    "Calcium": {"receptors": ["ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "RYR3",
                               "TRPC1", "TRPV4", "TRPM7", "ORAI1", "P2RX7",
                               "P2RX4", "CACNA1C"], "count": 13},
    "Cell Cycle": {"receptors": ["CCND1", "CDK4", "CDK6", "MYC"], "count": 4},
    "Circadian": {"receptors": ["OPN4", "MTNR1A", "MTNR1B"], "count": 3},
    "NRF2": {"receptors": ["KEAP1"], "count": 1},  # Direct sensor, few receptors
    "PI3K/PTEN": {"receptors": ["EGFR", "ERBB2", "IGF1R", "INSR", "PDGFRA",
                                 "PDGFRB", "ITGB1", "ITGA5", "KIT", "MET",
                                 "FGFR1"], "count": 11},
    # New modules
    "AMPK": {"receptors": ["ADIPOR1", "ADIPOR2", "INSR", "SLC2A1", "SLC2A4",
                            "CAMKK2"], "count": 6},
    "SREBP": {"receptors": ["LDLR", "SCARB1", "NPC1", "NPC1L1", "INSR"], "count": 5},
    "ATR/CHK1": {"receptors": ["RPA1", "RPA2", "RAD17", "RFC2"], "count": 4},
    "Rho/ROCK": {"receptors": ["ITGB1", "ITGA5", "ITGAV", "CDH1", "CDH2",
                                "LPAR1", "S1PR1", "PIEZO1", "PIEZO2"], "count": 9},
    "PPAR/LXR": {"receptors": ["CD36", "FABP1", "FABP4", "SLC27A1",
                                "SCARB1"], "count": 5},
    "Autophagy": {"receptors": ["LAMP1", "LAMP2", "MCOLN1", "P2RX7",
                                 "SQSTM1"], "count": 5},
}


# ===================================================================
# METABOLIC GENES: regulatory vs enzymatic
# ===================================================================

METABOLIC_REGULATORY = {
    # Genes that SENSE metabolic state and REGULATE via feedback
    "PRKAA1", "PRKAA2",  # AMPK catalytic
    "STK11",  # LKB1 → AMPK
    "MTOR", "RPTOR", "TSC1", "TSC2",  # mTOR pathway
    "HIF1A", "EPAS1",  # Hypoxia sensors
    "SIRT1", "SIRT3", "SIRT6",  # Sirtuins (NAD+ sensors)
    "PPARGC1A",  # PGC-1α (metabolic regulator)
    "PPARA", "PPARG", "PPARD",  # PPARs
    "SREBF1", "SREBF2",  # SREBP
    "NR1H3", "NR1H2",  # LXR
    "MLXIPL",  # ChREBP (glucose sensor/TF)
    "NFE2L2",  # NRF2 (redox sensor)
    "KEAP1",  # NRF2 inhibitor
    "PRKAB1", "PRKAG1",  # AMPK regulatory
    "DEPTOR",  # mTOR inhibitor
    "INSIG1", "INSIG2",  # SREBP regulators
    "NCOR1",  # Nuclear co-repressor
    "FOXO1", "FOXO3",  # FOXO TFs (metabolic integration)
    "TIGAR",  # p53-regulated glycolysis
    "PKM",  # Pyruvate kinase M (splice-regulated)
    "TP53",  # p53 regulates metabolism
}

METABOLIC_ENZYMATIC = {
    # Glycolysis enzymes (no feedback regulation of their own)
    "HK1", "HK2",  # Hexokinase
    "GPI",  # Glucose-6-phosphate isomerase
    "PFKL", "PFKM",  # Phosphofructokinase
    "ALDOA", "ALDOB",  # Aldolase
    "GAPDH",  # GAPDH
    "PGK1",  # Phosphoglycerate kinase
    "ENO1", "ENO2",  # Enolase
    "LDHA", "LDHB",  # Lactate dehydrogenase
    # TCA cycle
    "CS",  # Citrate synthase
    "ACO1", "ACO2",  # Aconitase
    "IDH1", "IDH2", "IDH3A",  # Isocitrate DH
    "OGDH",  # α-ketoglutarate DH
    "SUCLG1",  # Succinyl-CoA ligase
    "SDHA", "SDHB", "SDHC", "SDHD",  # Succinate DH
    "FH",  # Fumarase
    "MDH1", "MDH2",  # Malate DH
    # Oxidative phosphorylation
    "NDUFS1", "NDUFV1",  # Complex I
    "UQCRC1", "UQCRC2",  # Complex III
    "COX4I1", "COX5A",  # Complex IV
    "ATP5F1A", "ATP5F1B",  # ATP synthase
    # Fatty acid oxidation
    "CPT1A", "CPT2",  # Carnitine palmitoyltransferase
    "ACADM", "ACADL",  # Acyl-CoA dehydrogenase
    "HADHA", "HADHB",  # Trifunctional enzyme
    # Fatty acid synthesis
    "FASN",  # Fatty acid synthase
    "ACLY",  # ATP citrate lyase
    "ACACA",  # ACC1
    # Amino acid metabolism
    "GOT1", "GOT2",  # Aspartate aminotransferase
    "GLS", "GLS2",  # Glutaminase
    "GLUD1",  # Glutamate DH
}


def main():
    cgc = load_cgc()

    # ===================================================================
    # TEST C.2: Input convergence → cancer vulnerability
    # ===================================================================
    print("=" * 80)
    print("TEST C.2: Input convergence (receptor count) → CGC density per module")
    print("=" * 80)

    modules_data = []
    for mod, info in MODULE_RECEPTOR_COUNT.items():
        genes = MODULE_GENES.get(mod, set())
        n_recs = info["count"]
        n_genes = len(genes)
        n_cgc = len(genes & cgc)
        cgc_frac = n_cgc / n_genes if n_genes > 0 else 0
        modules_data.append({
            "module": mod,
            "n_receptors": n_recs,
            "n_genes": n_genes,
            "n_cgc": n_cgc,
            "cgc_fraction": cgc_frac,
            "cgc_members": sorted(genes & cgc),
        })

    # Sort by receptor count
    modules_data.sort(key=lambda x: -x["n_receptors"])

    print(f"\n  {'Module':15s} {'Receptors':>9s} {'Genes':>5s} {'CGC':>4s} {'CGC%':>6s} {'CGC members'}")
    print("  " + "-" * 85)
    for d in modules_data:
        print(f"  {d['module']:15s} {d['n_receptors']:>9d} {d['n_genes']:>5d} {d['n_cgc']:>4d} "
              f"{d['cgc_fraction']*100:>5.0f}% {', '.join(d['cgc_members'][:5])}")

    # Correlation: receptor count vs CGC fraction
    recs = np.array([d["n_receptors"] for d in modules_data])
    cgc_frac = np.array([d["cgc_fraction"] for d in modules_data])
    n_genes = np.array([d["n_genes"] for d in modules_data])

    rho, p_rho = stats.spearmanr(recs, cgc_frac)
    r, p_r = stats.pearsonr(recs, cgc_frac)

    print(f"\n  Spearman ρ(receptors, CGC fraction): {rho:.3f}, p = {p_rho:.4f}")
    print(f"  Pearson r(receptors, CGC fraction): {r:.3f}, p = {p_r:.4f}")

    # Partial correlation controlling for module size
    from step15_confound_checks import partial_corr
    r_part, p_part = partial_corr(recs.astype(float), cgc_frac, n_genes.astype(float))
    print(f"  Partial r (controlling for module size): {r_part:.3f}, p = {p_part:.4f}")

    # Split: high convergence (>10 receptors) vs low (<= 10)
    high_conv = [d for d in modules_data if d["n_receptors"] > 10]
    low_conv = [d for d in modules_data if d["n_receptors"] <= 10]
    high_cgc = [d["cgc_fraction"] for d in high_conv]
    low_cgc = [d["cgc_fraction"] for d in low_conv]

    u_stat, p_mann = stats.mannwhitneyu(high_cgc, low_cgc, alternative='greater')
    print(f"\n  High convergence (>10 receptors): N={len(high_conv)}, mean CGC={np.mean(high_cgc)*100:.1f}%")
    print(f"  Low convergence (≤10 receptors): N={len(low_conv)}, mean CGC={np.mean(low_cgc)*100:.1f}%")
    print(f"  Mann-Whitney U (high > low): U={u_stat:.0f}, p={p_mann:.4f}")

    # ===================================================================
    # TEST E: Metabolic regulatory vs enzymatic
    # ===================================================================
    print(f"\n{'='*80}")
    print("TEST E: Metabolic regulatory (feedback) vs enzymatic (no feedback) → CGC")
    print("=" * 80)

    reg_in_cgc = METABOLIC_REGULATORY & cgc
    reg_not_cgc = METABOLIC_REGULATORY - cgc
    enz_in_cgc = METABOLIC_ENZYMATIC & cgc
    enz_not_cgc = METABOLIC_ENZYMATIC - cgc

    print(f"\n  Regulatory metabolic genes: {len(METABOLIC_REGULATORY)}")
    print(f"    In CGC: {len(reg_in_cgc)} ({100*len(reg_in_cgc)/len(METABOLIC_REGULATORY):.1f}%)")
    print(f"    CGC members: {sorted(reg_in_cgc)}")
    print(f"\n  Enzymatic metabolic genes: {len(METABOLIC_ENZYMATIC)}")
    print(f"    In CGC: {len(enz_in_cgc)} ({100*len(enz_in_cgc)/len(METABOLIC_ENZYMATIC):.1f}%)")
    print(f"    CGC members: {sorted(enz_in_cgc)}")

    # Fisher exact test
    table = [[len(reg_in_cgc), len(reg_not_cgc)],
             [len(enz_in_cgc), len(enz_not_cgc)]]
    odds, p_fisher = stats.fisher_exact(table, alternative='greater')
    print(f"\n  Fisher exact test (regulatory > enzymatic):")
    print(f"    Contingency: {table}")
    print(f"    Odds ratio: {odds:.2f}")
    print(f"    p = {p_fisher:.4f}")

    # Baseline CGC rate
    # ~700 CGC genes out of ~20,000 protein-coding = ~3.5%
    baseline_rate = 700 / 20000
    print(f"\n  Baseline CGC rate: ~{baseline_rate*100:.1f}%")
    print(f"  Regulatory: {100*len(reg_in_cgc)/len(METABOLIC_REGULATORY):.1f}% (enrichment: {len(reg_in_cgc)/len(METABOLIC_REGULATORY)/baseline_rate:.1f}×)")
    print(f"  Enzymatic:  {100*len(enz_in_cgc)/len(METABOLIC_ENZYMATIC):.1f}% (enrichment: {len(enz_in_cgc)/len(METABOLIC_ENZYMATIC)/baseline_rate:.1f}×)")

    # Binomial test for each
    p_reg = stats.binomtest(len(reg_in_cgc), len(METABOLIC_REGULATORY), baseline_rate, alternative='greater')
    p_enz = stats.binomtest(len(enz_in_cgc), len(METABOLIC_ENZYMATIC), baseline_rate, alternative='greater')
    print(f"\n  Binomial test (vs baseline):")
    print(f"    Regulatory: p = {p_reg.pvalue:.4f}")
    print(f"    Enzymatic:  p = {p_enz.pvalue:.4f}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)

    results = {
        "convergence": {
            "spearman_rho": round(float(rho), 3),
            "spearman_p": round(float(p_rho), 4),
            "partial_r": round(float(r_part), 3),
            "partial_p": round(float(p_part), 4),
            "high_conv_mean_cgc": round(float(np.mean(high_cgc)), 3),
            "low_conv_mean_cgc": round(float(np.mean(low_cgc)), 3),
            "mann_whitney_p": round(float(p_mann), 4),
        },
        "metabolic": {
            "regulatory_n": len(METABOLIC_REGULATORY),
            "regulatory_cgc": len(reg_in_cgc),
            "enzymatic_n": len(METABOLIC_ENZYMATIC),
            "enzymatic_cgc": len(enz_in_cgc),
            "fisher_odds": round(float(odds), 2),
            "fisher_p": round(float(p_fisher), 4),
        },
    }

    with open(DATA_DIR / "convergence_metabolic.json", "w") as f:
        json.dump(results, f, indent=2)

    if p_rho < 0.05:
        print(f"\n  ✓ Convergence → cancer: more receptors → more cancer genes (ρ={rho:.3f})")
    else:
        print(f"\n  ✗ Convergence → cancer: no significant correlation (ρ={rho:.3f})")

    if p_fisher < 0.05:
        print(f"  ✓ Regulatory metabolic > enzymatic in CGC (OR={odds:.1f})")
    else:
        print(f"  ✗ Regulatory vs enzymatic: no significant difference (OR={odds:.1f})")

    print(f"\nSaved to convergence_metabolic.json")


if __name__ == "__main__":
    main()
