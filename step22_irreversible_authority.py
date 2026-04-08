#!/usr/bin/env python3
"""
Binary Irreversible Authority (IA) vs raw authority vs vulnerability.

For each module: count ONLY clearly irreversible decisions:
  - Apoptosis (cell death)
  - Senescence (permanent arrest)
  - Terminal differentiation (permanent exit from stem state)
  - Ferroptosis / necroptosis (cell death variants)

No subjective 0.1-1.0 scoring. Binary: irreversible or not.
Compare: IA vs raw authority vs convergence → CGC density.
"""

import numpy as np
from scipy import stats
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")


def load_cgc():
    cgc = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            gene = line.strip()
            if gene:
                cgc.add(gene)
    return cgc


# Module data: genes, all decisions, irreversible decisions marked
MODULES = {
    "p53": {
        "genes": {"TP53", "MDM2", "MDM4", "CDKN1A", "BAX", "BBC3", "PMAIP1",
                  "ATM", "ATR", "CHEK1", "CHEK2"},
        "all_decisions": [
            ("Apoptosis (BAX/BAK)", True),
            ("Senescence (sustained p21+RB)", True),
            ("Ferroptosis (SLC7A11 suppression)", True),
            ("G1 cell cycle arrest", False),  # reversible — cell can re-enter
            ("G2 cell cycle arrest", False),
            ("DNA repair activation", False),
            ("Metabolic reprogramming", False),
            ("Autophagy induction", False),
            ("Anti-angiogenesis", False),
        ],
        "authority": 9,
        "n_receptors": 5,
    },
    "Cell Cycle": {
        "genes": {"CDK2", "CDK4", "CDK6", "CCND1", "CCNE1", "CCNA2", "CCNB1",
                  "RB1", "E2F1", "CDKN1A", "CDKN2A", "CDKN1B", "CDC25A"},
        "all_decisions": [
            ("G1/S commitment (restriction point)", True),  # past restriction point = committed
            ("S phase progression", False),
            ("G2/M transition", False),
            ("Mitotic exit", False),
        ],
        "authority": 4,
        "n_receptors": 4,
    },
    "Wnt": {
        "genes": {"CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B", "DVL1",
                  "TCF7L2", "LEF1", "RNF43", "ZNRF3"},
        "all_decisions": [
            ("Stem cell self-renewal vs differentiation", True),  # differentiation largely irreversible
            ("Cell fate specification", True),  # lineage commitment irreversible
            ("Tissue patterning", False),
        ],
        "authority": 3,
        "n_receptors": 15,
    },
    "TGF-β": {
        "genes": {"TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7",
                  "SMURF1", "SMURF2", "BMPR1A", "BMPR2", "ACVR1"},
        "all_decisions": [
            ("EMT (epithelial-mesenchymal)", True),  # involves chromatin remodeling, slow to reverse
            ("Growth arrest in epithelial", False),  # reversible
            ("Fibrosis/ECM production", False),
            ("Immune suppression (Treg)", True),  # Treg commitment largely stable
        ],
        "authority": 4,
        "n_receptors": 9,
    },
    "PI3K/PTEN": {
        "genes": {"PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "AKT1", "AKT2",
                  "PDK1", "INPP4B"},
        "all_decisions": [
            ("Cell survival (block apoptosis)", True),  # preventing death = irreversible consequence
            ("Growth factor response", False),
            ("Glucose uptake", False),
        ],
        "authority": 3,
        "n_receptors": 11,
    },
    "NF-κB": {
        "genes": {"RELA", "RELB", "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "TNFAIP3",
                  "TRAF2", "TRAF6", "IKBKB", "IKBKG", "CHUK", "MAP3K7", "TAB1", "TAB2"},
        "all_decisions": [
            ("Anti-apoptotic program (survival)", True),  # blocking death = irreversible consequence
            ("Inflammatory cytokine production", False),
            ("Immune cell activation", False),
            ("Acute phase response", False),
        ],
        "authority": 4,
        "n_receptors": 21,
    },
    "ERK/MAPK": {
        "genes": {"MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BRAF", "RAF1", "ARAF",
                  "HRAS", "KRAS", "NRAS", "SOS1", "GRB2", "DUSP1", "DUSP6", "SPRY2"},
        "all_decisions": [
            ("Proliferation commitment", False),  # each cycle is new decision
            ("Differentiation (sustained ERK)", True),  # terminal differentiation
            ("Migration", False),
        ],
        "authority": 3,
        "n_receptors": 17,
    },
    "JAK-STAT": {
        "genes": {"JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT3", "STAT5A", "STAT5B",
                  "SOCS1", "SOCS3", "CISH", "PIAS1"},
        "all_decisions": [
            ("Immune differentiation (Th lineage)", True),  # T cell commitment irreversible
            ("Hematopoietic commitment", True),  # lineage choice irreversible
            ("Antiviral IFN response", False),
        ],
        "authority": 3,
        "n_receptors": 22,
    },
    "Notch": {
        "genes": {"NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "FBXW7", "HES1", "HEY1",
                  "MAML1", "RBPJ", "DLL1", "DLL4", "JAG1", "JAG2"},
        "all_decisions": [
            ("Lateral inhibition (binary fate)", True),  # cell fate choice irreversible
            ("Boundary formation", False),
        ],
        "authority": 2,
        "n_receptors": 4,
    },
    "Hippo": {
        "genes": {"YAP1", "WWTR1", "LATS1", "LATS2", "STK3", "STK4", "SAV1",
                  "MOB1A", "NF2", "TEAD1", "TEAD4"},
        "all_decisions": [
            ("Organ size control (terminal growth arrest)", True),
            ("Contact inhibition", False),
        ],
        "authority": 2,
        "n_receptors": 9,
    },
    "mTOR": {
        "genes": {"MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "RHEB", "RPS6KB1",
                  "EIF4EBP1", "DEPTOR", "MLST8"},
        "all_decisions": [
            ("Protein synthesis activation", False),
            ("Lipid synthesis", False),
            ("Autophagy inhibition", False),
        ],
        "authority": 3,
        "n_receptors": 10,
    },
    "Calcium": {
        "genes": {"PLCG1", "PLCG2", "ITPR1", "ITPR2", "ATP2A2", "CALM1",
                  "NFATC1", "NFATC2", "CAMK2A", "CAMK2B"},
        "all_decisions": [
            ("NFAT transcription (immune)", False),  # transient
            ("Muscle contraction", False),
        ],
        "authority": 2,
        "n_receptors": 13,
    },
    "Circadian": {
        "genes": {"CLOCK", "ARNTL", "PER1", "PER2", "CRY1", "CRY2", "NR1D1",
                  "NR1D2", "CSNK1D", "CSNK1E", "FBXL3"},
        "all_decisions": [
            ("Temporal gating of division", False),
            ("Metabolic rhythm", False),
        ],
        "authority": 2,
        "n_receptors": 3,
    },
    "NRF2": {
        "genes": {"NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1",
                  "SOD2", "CAT", "GPX1"},
        "all_decisions": [
            ("Antioxidant response", False),
            ("Xenobiotic metabolism", False),
        ],
        "authority": 2,
        "n_receptors": 1,
    },
    # New modules
    "AMPK": {
        "genes": {"PRKAA1", "PRKAA2", "PRKAB1", "PRKAG1", "STK11", "ACACB",
                  "PPARGC1A", "FOXO3", "TSC2", "CREB1"},
        "all_decisions": [
            ("Energy homeostasis", False),
            ("mTOR inhibition", False),
            ("Autophagy activation", False),
        ],
        "authority": 3,
        "n_receptors": 6,
    },
    "SREBP": {
        "genes": {"SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2", "HMGCR",
                  "FASN", "SCD", "ACLY"},
        "all_decisions": [
            ("Lipogenesis program", False),
            ("Cholesterol homeostasis", False),
        ],
        "authority": 2,
        "n_receptors": 5,
    },
    "ATR/CHK1": {
        "genes": {"ATR", "CHEK1", "ATRIP", "TOPBP1", "WEE1", "CDC25A",
                  "RAD17", "RAD9A", "HUS1"},
        "all_decisions": [
            ("S-phase checkpoint", False),  # transient arrest
            ("Fork stabilization", False),
        ],
        "authority": 2,
        "n_receptors": 4,
    },
    "Rho/ROCK": {
        "genes": {"RHOA", "RHOB", "RHOC", "ROCK1", "ROCK2", "MKL1",
                  "LIMK1", "CFL1", "ARHGAP1", "ARHGAP5"},
        "all_decisions": [
            ("Cytoskeletal remodeling", False),
            ("MRTF/SRF transcription", False),
        ],
        "authority": 2,
        "n_receptors": 9,
    },
    "PPAR/LXR": {
        "genes": {"PPARA", "PPARG", "PPARD", "NR1H3", "NR1H2", "RXRA",
                  "NCOR1", "NCOR2", "NCOA1", "NCOA2"},
        "all_decisions": [
            ("Fatty acid oxidation program", False),
            ("Adipogenesis", True),  # adipocyte differentiation is terminal
            ("Cholesterol efflux", False),
        ],
        "authority": 3,
        "n_receptors": 5,
    },
    "Autophagy": {
        "genes": {"ULK1", "ULK2", "BECN1", "ATG5", "ATG7", "ATG12",
                  "TFEB", "SQSTM1", "MAP1LC3B"},
        "all_decisions": [
            ("Bulk macroautophagy", False),
            ("Selective autophagy", False),
        ],
        "authority": 2,
        "n_receptors": 5,
    },
}


def main():
    cgc = load_cgc()

    # Compute metrics
    data = []
    for name, info in MODULES.items():
        genes = info["genes"]
        n_cgc = len(genes & cgc)
        cgc_frac = n_cgc / len(genes)

        irrev_decisions = [d for d, is_irrev in info["all_decisions"] if is_irrev]
        n_irrev = len(irrev_decisions)
        n_total = len(info["all_decisions"])

        data.append({
            "module": name,
            "n_genes": len(genes),
            "n_cgc": n_cgc,
            "cgc_frac": cgc_frac,
            "authority": info["authority"],
            "irrev_authority": n_irrev,
            "n_decisions": n_total,
            "irrev_fraction": n_irrev / n_total if n_total > 0 else 0,
            "n_receptors": info["n_receptors"],
            "irrev_decisions": irrev_decisions,
        })

    # Sort by irreversible authority
    data.sort(key=lambda x: (-x["irrev_authority"], -x["cgc_frac"]))

    # ===================================================================
    # TABLE
    # ===================================================================
    print("=" * 80)
    print("IRREVERSIBLE AUTHORITY (IA) vs RAW AUTHORITY vs CONVERGENCE → CGC")
    print("=" * 80)

    print(f"\n  {'Module':15s} {'IA':>3s} {'Auth':>5s} {'Recv':>5s} {'CGC%':>6s} {'Irreversible decisions'}")
    print("  " + "-" * 80)
    for d in data:
        irrev_str = "; ".join(d["irrev_decisions"][:3]) if d["irrev_decisions"] else "none"
        print(f"  {d['module']:15s} {d['irrev_authority']:>3d} {d['authority']:>5d} "
              f"{d['n_receptors']:>5d} {d['cgc_frac']*100:>5.0f}%  {irrev_str}")

    # ===================================================================
    # CORRELATIONS
    # ===================================================================
    print(f"\n{'='*80}")
    print("CORRELATION COMPARISON")
    print("=" * 80)

    cgc_arr = np.array([d["cgc_frac"] for d in data])
    ia_arr = np.array([d["irrev_authority"] for d in data], dtype=float)
    auth_arr = np.array([d["authority"] for d in data], dtype=float)
    conv_arr = np.array([d["n_receptors"] for d in data], dtype=float)

    metrics = {
        "Irreversible Authority (IA)": ia_arr,
        "Raw Authority": auth_arr,
        "Convergence (receptors)": conv_arr,
        "IA + Convergence (sum)": ia_arr + conv_arr / 10,  # scaled
    }

    print(f"\n  {'Metric':35s} {'Spearman ρ':>11s} {'p':>10s} {'Pearson r':>10s} {'R²':>6s}")
    print("  " + "-" * 75)

    best_rho = 0
    best_name = ""
    for name, arr in metrics.items():
        rho, p_rho = stats.spearmanr(arr, cgc_arr)
        r, p_r = stats.pearsonr(arr, cgc_arr)
        sig = "***" if p_rho < 0.001 else "**" if p_rho < 0.01 else "*" if p_rho < 0.05 else ""
        print(f"  {name:35s} {rho:+11.3f} {p_rho:10.4f} {r:+10.3f} {r**2:6.3f}  {sig}")
        if abs(rho) > abs(best_rho):
            best_rho = rho
            best_name = name

    # ===================================================================
    # IA = 0 vs IA > 0 split
    # ===================================================================
    print(f"\n{'='*80}")
    print("BINARY SPLIT: modules with ANY irreversible decision vs NONE")
    print("=" * 80)

    has_irrev = [d for d in data if d["irrev_authority"] > 0]
    no_irrev = [d for d in data if d["irrev_authority"] == 0]

    cgc_has = [d["cgc_frac"] for d in has_irrev]
    cgc_no = [d["cgc_frac"] for d in no_irrev]

    u, p_u = stats.mannwhitneyu(cgc_has, cgc_no, alternative='greater')

    print(f"\n  With irreversible decisions (N={len(has_irrev)}):")
    print(f"    Mean CGC: {np.mean(cgc_has)*100:.1f}%")
    print(f"    Modules: {', '.join(d['module'] for d in has_irrev)}")
    print(f"\n  Without irreversible decisions (N={len(no_irrev)}):")
    print(f"    Mean CGC: {np.mean(cgc_no)*100:.1f}%")
    print(f"    Modules: {', '.join(d['module'] for d in no_irrev)}")
    print(f"\n  Mann-Whitney U (irrev > no irrev): U={u:.0f}, p = {p_u:.4f}")

    # ===================================================================
    # DOSE-RESPONSE: CGC by IA level
    # ===================================================================
    print(f"\n{'='*80}")
    print("DOSE-RESPONSE: CGC% by number of irreversible decisions")
    print("=" * 80)

    for ia_level in sorted(set(d["irrev_authority"] for d in data)):
        mods = [d for d in data if d["irrev_authority"] == ia_level]
        mean_cgc = np.mean([d["cgc_frac"] for d in mods])
        names = ", ".join(d["module"] for d in mods)
        print(f"  IA={ia_level}: mean CGC = {mean_cgc*100:.1f}% (N={len(mods)}) — {names}")

    # ===================================================================
    # RESIDUALS from IA prediction
    # ===================================================================
    print(f"\n{'='*80}")
    print("RESIDUALS: What IA misses")
    print("=" * 80)

    slope, intercept, _, _, _ = stats.linregress(ia_arr, cgc_arr)
    predicted = slope * ia_arr + intercept
    residuals = cgc_arr - predicted

    outliers = [(data[i]["module"], residuals[i]) for i in range(len(data)) if abs(residuals[i]) > 0.15]
    if outliers:
        print(f"\n  Modules with |residual| > 15%:")
        for mod, resid in sorted(outliers, key=lambda x: -abs(x[1])):
            direction = "MORE" if resid > 0 else "LESS"
            print(f"    {mod:15s} residual = {resid*100:+.0f}% ({direction} than IA predicts)")
    else:
        print(f"\n  No large outliers — IA captures most variance")

    # ===================================================================
    # SAVE
    # ===================================================================
    rho_ia, p_ia = stats.spearmanr(ia_arr, cgc_arr)
    rho_auth, p_auth = stats.spearmanr(auth_arr, cgc_arr)

    results = {
        "n_modules": len(data),
        "ia_spearman_rho": round(float(rho_ia), 3),
        "ia_spearman_p": round(float(p_ia), 4),
        "auth_spearman_rho": round(float(rho_auth), 3),
        "auth_spearman_p": round(float(p_auth), 4),
        "binary_split_p": round(float(p_u), 4),
        "irrev_mean_cgc": round(float(np.mean(cgc_has)), 3),
        "rev_mean_cgc": round(float(np.mean(cgc_no)), 3),
        "best_predictor": best_name,
    }

    with open(DATA_DIR / "irreversible_authority.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*80}")
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  Best predictor: {best_name} (ρ = {best_rho:.3f})")
    print(f"  IA: ρ = {rho_ia:.3f}, p = {p_ia:.4f}")
    print(f"  Raw authority: ρ = {rho_auth:.3f}, p = {p_auth:.4f}")
    print(f"  Binary split: irrev {np.mean(cgc_has)*100:.0f}% vs rev {np.mean(cgc_no)*100:.0f}%, p = {p_u:.4f}")
    print(f"\nSaved to irreversible_authority.json")


if __name__ == "__main__":
    main()
