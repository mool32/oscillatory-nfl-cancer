#!/usr/bin/env python3
"""
Phase A.1: Screen new module candidates.

For each candidate, check 5 criteria:
1. Convergence: >3 different inputs
2. Negative feedback: documented feedback loop with core genes
3. TF output: transcriptional program
4. Signal transformation: input ≠ output format
5. Temporal competence: architecture supports cycling

Then compute all metrics (cognitive load contribution, CGC enrichment,
receptor coupling) for confirmed modules.

Candidates:
  AMPK, SREBP, ATR/CHK1, Rho/ROCK/MRTF, PPAR/LXR, Autophagy/TFEB
"""

import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import json

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# ===================================================================
# New module candidates with criteria assessment
# ===================================================================
# Based on established literature — each field references known biology

NEW_MODULES = {
    "AMPK": {
        "description": "Energy/metabolic stress sensor",
        "core_genes": {
            "activator": {"PRKAA1", "PRKAA2", "PRKAB1", "PRKAG1", "STK11"},  # AMPK subunits + LKB1
            "inhibitor": {"PPARGC1A", "ACACB"},  # PGC-1α (downstream target that feeds back), ACC
            "core_for_score": {"PRKAA1", "PRKAA2", "STK11", "ACACB"},
        },
        "criteria": {
            "convergence": True,  # AMP:ATP ratio, Ca2+, glucose deprivation, metformin, exercise
            "negative_feedback": True,  # AMPK→mTOR inhibition→reduced protein synthesis→ATP recovery→AMPK off
            "tf_output": True,  # CREB, FOXO3, PGC-1α transcriptional programs
            "signal_transform": True,  # AMP:ATP ratio (metabolic) → gene expression (transcriptional)
            "temporal_competence": True,  # AMPK shows pulsatile activation (Tsou et al 2011)
        },
        "upstream_receptors": {"SLC2A1", "SLC2A4", "ADIPOR1", "ADIPOR2", "INSR",
                               "CAMKK2", "STK11"},
        "n_criteria_met": 5,
    },

    "SREBP": {
        "description": "Lipid/cholesterol sensing and synthesis",
        "core_genes": {
            "activator": {"SREBF1", "SREBF2", "SCAP"},
            "inhibitor": {"INSIG1", "INSIG2", "HMGCR"},
            "core_for_score": {"SREBF1", "SREBF2", "INSIG1", "SCAP"},
        },
        "criteria": {
            "convergence": True,  # Cholesterol, oxysterols, fatty acids, insulin, ER stress
            "negative_feedback": True,  # SREBP→HMGCR→cholesterol→INSIG→sequesters SCAP→blocks SREBP
            "tf_output": True,  # SREBP is itself a TF (bHLH), drives LDLR, FASN, SCD, HMGCR
            "signal_transform": True,  # Cholesterol level (metabolic) → lipogenic gene program
            "temporal_competence": True,  # SREBP shows oscillatory processing in cholesterol-deprived cells
        },
        "upstream_receptors": {"LDLR", "SCARB1", "NPC1", "NPC1L1", "INSR", "SCAP"},
        "n_criteria_met": 5,
    },

    "ATR/CHK1": {
        "description": "Replication stress checkpoint",
        "core_genes": {
            "activator": {"ATR", "CHEK1", "ATRIP", "TOPBP1"},
            "inhibitor": {"WEE1", "CDC25A", "CLSPN"},  # CDC25A degradation = negative output
            "core_for_score": {"ATR", "CHEK1", "WEE1", "CDC25A"},
        },
        "criteria": {
            "convergence": True,  # ssDNA, stalled forks, nucleotide depletion, oncogene activation
            "negative_feedback": True,  # ATR→CHK1→CDC25A degradation→CDK inhibition→fork stabilization→ATR off
            "tf_output": False,  # Primarily post-translational (phosphorylation), not transcriptional
            "signal_transform": True,  # DNA structure → cell cycle arrest
            "temporal_competence": True,  # CHK1 shows pulsatile activation during S phase
        },
        "upstream_receptors": {"RPA1", "RPA2", "RAD17", "RFC2"},  # ssDNA sensors, not classic receptors
        "n_criteria_met": 4,  # Missing TF output
    },

    "Rho/ROCK/MRTF": {
        "description": "Mechanotransduction and cytoskeletal sensing",
        "core_genes": {
            "activator": {"RHOA", "ROCK1", "ROCK2", "MKL1"},  # MKL1 = MRTF-A
            "inhibitor": {"ARHGAP1", "LIMK1", "CFL1"},  # GAP, cofilin feedback
            "core_for_score": {"RHOA", "ROCK1", "MKL1", "CFL1"},
        },
        "criteria": {
            "convergence": True,  # ECM stiffness, cell-cell contact, GPCRs, integrins, growth factors
            "negative_feedback": True,  # Rho→ROCK→LIMK→cofilin-P→actin stabilization→reduced Rho activation
            "tf_output": True,  # MRTF-A/SRF axis drives transcriptional program (actin genes, etc)
            "signal_transform": True,  # Mechanical force → gene expression
            "temporal_competence": True,  # Rho oscillations documented in migrating cells (Machacek 2009)
        },
        "upstream_receptors": {"ITGB1", "ITGA5", "ITGAV", "CDH1", "CDH2",
                               "LPAR1", "S1PR1", "PIEZO1"},
        "n_criteria_met": 5,
    },

    "PPAR/LXR": {
        "description": "Lipid-activated nuclear receptor sensing",
        "core_genes": {
            "activator": {"PPARA", "PPARG", "NR1H3", "NR1H2"},  # NR1H3=LXRα, NR1H2=LXRβ
            "inhibitor": {"NCOR1", "NCOR2", "SIRT1"},  # Nuclear co-repressors
            "core_for_score": {"PPARA", "PPARG", "NR1H3", "NCOR1"},
        },
        "criteria": {
            "convergence": True,  # Fatty acids, eicosanoids, oxysterols, thiazolidinediones
            "negative_feedback": True,  # PPARγ→adipogenesis→lipid uptake→reduced free FA→reduced PPARγ activation
            "tf_output": True,  # PPARs ARE transcription factors (nuclear receptors)
            "signal_transform": True,  # Lipid metabolites → gene expression
            "temporal_competence": False,  # Slow (hours-days), no documented oscillation
        },
        "upstream_receptors": {"CD36", "FABP1", "FABP4", "SLC27A1", "SLC27A2",
                               "SCARB1"},
        "n_criteria_met": 4,  # Missing temporal competence
    },

    "Autophagy/TFEB": {
        "description": "Lysosomal sensing and autophagy regulation",
        "core_genes": {
            "activator": {"TFEB", "ATG5", "ATG7", "BECN1", "ULK1"},
            "inhibitor": {"MTOR", "RPTOR", "LAMP1"},  # mTOR inhibits TFEB; LAMP1 = lysosome
            "core_for_score": {"TFEB", "BECN1", "ULK1", "ATG7"},
        },
        "criteria": {
            "convergence": True,  # Starvation, ER stress, pathogen, damaged organelles, hypoxia
            "negative_feedback": True,  # TFEB→lysosome biogenesis→mTORC1 reactivation on lysosome→TFEB phosphorylation→nuclear export
            "tf_output": True,  # TFEB drives CLEAR gene network (>500 genes)
            "signal_transform": True,  # Lysosomal state → transcriptional program
            "temporal_competence": True,  # TFEB shows nuclear-cytoplasmic shuttling cycles
        },
        "upstream_receptors": {"LAMP1", "LAMP2", "ATP6V1A", "MCOLN1",
                               "P2RX7", "SQSTM1"},  # lysosomal sensors
        "n_criteria_met": 5,
    },
}


def load_expression():
    expr = defaultdict(dict)
    with open(DATA_DIR / "rna_single_cell_type.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            expr[row["Gene name"]][row["Cell type"]] = float(row["nCPM"])
    cell_types = sorted(set(ct for gd in expr.values() for ct in gd))
    return expr, cell_types


def load_cgc():
    cgc = set()
    cgc_path = DATA_DIR / "cancer_gene_census.csv"
    if cgc_path.exists():
        with open(cgc_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                name_field = "Gene Symbol" if "Gene Symbol" in row else list(row.keys())[0]
                cgc.add(row[name_field].strip())
    return cgc


def load_cognitive_load():
    cog = {}
    with open(DATA_DIR / "cognitive_load.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cog[row["cell_type"]] = int(row["cognitive_load"])
    return cog


def main():
    print("Loading data...")
    expr, cell_types = load_expression()
    cgc = load_cgc()
    cog = load_cognitive_load()

    # ===================================================================
    # SCREENING: 5 criteria per candidate
    # ===================================================================
    print(f"\n{'='*80}")
    print("NEW MODULE CANDIDATE SCREENING")
    print("=" * 80)

    print(f"\n  {'Module':20s} {'Conv':>5s} {'FB':>5s} {'TF':>5s} {'Trans':>5s} {'Temp':>5s} {'Score':>6s} {'Include?':>9s}")
    print("  " + "-" * 65)

    included = {}
    for name, info in NEW_MODULES.items():
        c = info["criteria"]
        score = sum(c.values())
        include = score >= 4
        sym = lambda v: "✓" if v else "✗"
        print(f"  {name:20s} {sym(c['convergence']):>5s} {sym(c['negative_feedback']):>5s} "
              f"{sym(c['tf_output']):>5s} {sym(c['signal_transform']):>5s} "
              f"{sym(c['temporal_competence']):>5s} {score:>4d}/5 {'→ INCLUDE' if include else '→ exclude'}")
        if include:
            included[name] = info

    print(f"\n  Included: {len(included)} new modules")

    # ===================================================================
    # METRICS for included modules
    # ===================================================================
    print(f"\n{'='*80}")
    print("METRICS FOR NEW MODULES")
    print("=" * 80)

    for name, info in included.items():
        print(f"\n  --- {name} ({info['description']}) ---")

        # Gene expression
        core = info["core_genes"]["core_for_score"]
        found = [g for g in core if g in expr]
        print(f"  Core genes: {sorted(core)} ({len(found)}/{len(core)} found)")

        # Activity per cell type (top/bottom)
        if found:
            scores = {}
            for ct in cell_types:
                vals = [expr[g].get(ct, 0) for g in found]
                scores[ct] = np.mean(vals)

            sorted_cts = sorted(scores.items(), key=lambda x: -x[1])
            print(f"  Top 5 expressing cell types:")
            for ct, sc in sorted_cts[:5]:
                print(f"    {ct:35s} {sc:8.1f} nCPM")
            print(f"  Bottom 3:")
            for ct, sc in sorted_cts[-3:]:
                print(f"    {ct:35s} {sc:8.1f} nCPM")

        # CGC enrichment
        all_genes = set()
        for cat in info["core_genes"].values():
            all_genes.update(cat)
        n_cgc = len(all_genes & cgc)
        print(f"  CGC genes: {n_cgc}/{len(all_genes)} ({100*n_cgc/len(all_genes):.0f}%)")
        if all_genes & cgc:
            print(f"    CGC members: {sorted(all_genes & cgc)}")

        # Receptor coupling
        recs = info.get("upstream_receptors", set())
        recs_found = [g for g in recs if g in expr]
        if found and recs_found:
            mod_scores = np.array([np.mean([expr[g].get(ct, 0) for g in found]) for ct in cell_types])
            rec_scores = np.array([np.mean([expr[g].get(ct, 0) for g in recs_found]) for ct in cell_types])

            # Remove overlap
            overlap = core & recs
            if overlap:
                print(f"  ⚠️  Gene overlap: {sorted(overlap)} — computing clean correlation")
                core_clean = [g for g in found if g not in recs]
                recs_clean = [g for g in recs_found if g not in core]
                if core_clean and recs_clean:
                    mod_clean = np.array([np.mean([expr[g].get(ct, 0) for g in core_clean]) for ct in cell_types])
                    rec_clean = np.array([np.mean([expr[g].get(ct, 0) for g in recs_clean]) for ct in cell_types])
                    r, p = stats.pearsonr(mod_clean, rec_clean)
                    print(f"  Receptor coupling (clean): r = {r:.3f}, p = {p:.2e}")
            else:
                r, p = stats.pearsonr(mod_scores, rec_scores)
                print(f"  Receptor coupling: r = {r:.3f}, p = {p:.2e}")

    # ===================================================================
    # COGNITIVE LOAD CONTRIBUTION
    # ===================================================================
    print(f"\n{'='*80}")
    print("COGNITIVE LOAD CONTRIBUTION OF NEW MODULES")
    print("=" * 80)

    ACTIVITY_THRESHOLD = 10.0  # nCPM threshold for "active"

    print(f"\n  {'Module':20s} {'N active CTs':>12s} {'Mean in active':>14s}")
    print("  " + "-" * 50)

    for name, info in included.items():
        core = info["core_genes"]["core_for_score"]
        found = [g for g in core if g in expr]
        if not found:
            continue

        n_active = 0
        active_scores = []
        for ct in cell_types:
            score = np.mean([expr[g].get(ct, 0) for g in found])
            if score > ACTIVITY_THRESHOLD:
                n_active += 1
                active_scores.append(score)

        mean_act = np.mean(active_scores) if active_scores else 0
        print(f"  {name:20s} {n_active:>8d}/{len(cell_types):<3d} {mean_act:>14.1f}")

    # ===================================================================
    # SAVE RESULTS
    # ===================================================================
    results = {
        "screened": len(NEW_MODULES),
        "included": len(included),
        "included_modules": list(included.keys()),
        "excluded_modules": [m for m in NEW_MODULES if m not in included],
        "criteria_summary": {
            name: {
                "score": info["n_criteria_met"],
                "missing": [k for k, v in info["criteria"].items() if not v]
            }
            for name, info in NEW_MODULES.items()
        },
    }

    with open(DATA_DIR / "new_modules_screening.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to new_modules_screening.json")
    print(f"\nTotal modules for expanded analysis: 14 + {len(included)} = {14 + len(included)}")


if __name__ == "__main__":
    main()
