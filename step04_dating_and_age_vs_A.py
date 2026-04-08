#!/usr/bin/env python3
"""
Phase 1.3 + 1.4 — Phylogenetic dating of feedback loops and Age vs A correlation.

For each loop: oldest species with BOTH components present → LCA divergence time.
Then: test correlation between loop age and waveform asymmetry A.

Pre-registered hypothesis:
  H1: Spearman ρ(age, A) > 0 (older loops = lower A = more asymmetric)
  H0: ρ = 0
  Threshold: p < 0.05
"""

import csv
import numpy as np
from scipy import stats
from pathlib import Path
import json

# ---------------------------------------------------------------------------
# TimeTree divergence times (Mya) — from TimeTree.org, human as reference
# Ordered from most distant to closest
# ---------------------------------------------------------------------------
DIVERGENCE_TIMES = {
    "D.discoideum":    1500,
    "A.thaliana":      1500,
    "S.pombe":         1105,
    "S.cerevisiae":    1105,
    "N.vectensis":      824,
    "C.elegans":        797,
    "D.melanogaster":   797,
    "S.purpuratus":     684,
    "C.intestinalis":   684,
    "D.rerio":          435,
    "X.tropicalis":     352,
    "G.gallus":         320,
    "C.l.familiaris":    96,
    "R.norvegicus":      90,
    "M.musculus":        90,
}

# Species ordered by divergence (most distant first)
SPECIES_ORDER = sorted(DIVERGENCE_TIMES.keys(), key=lambda s: -DIVERGENCE_TIMES[s])

# ---------------------------------------------------------------------------
# A values from Paper 1 (Phase 1 oscillator compilation)
# A = τ_rise / T — waveform asymmetry
# Only pathways with measured oscillatory A values
# ---------------------------------------------------------------------------
A_VALUES = {
    "NF-κB / IκBα":                          0.31,
    "p53 / Mdm2":                             0.32,
    "ERK / DUSP":                             0.28,
    "Wnt / APC-Axin":                         0.47,
    "Notch / FBXW7":                          0.47,  # Hes1 autorepressive ~0.47
    "Circadian (CLOCK-BMAL1 / PER-CRY)":     0.42,
    "Hippo / LATS":                           0.68,  # YAP shuttling
    "mTOR / TSC":                             0.35,  # estimated from pS6K dynamics
    "JAK-STAT / SOCS":                        0.50,
    "TGF-β / SMAD6-7":                        0.60,  # Smad shuttling
    "NRF2 / KEAP1":                           0.68,  # Nrf2 shuttling class
    "Calcium / SERCA":                        0.20,  # fast Ca2+ spikes
    "Cell Cycle CKI / CDK":                   0.35,  # CDK oscillation
    "NFAT / RCAN":                            0.68,  # NFAT shuttling class
}

# Feedback type (from Phase 1.1)
FEEDBACK_TYPE = {
    "NF-κB / IκBα": "SEQ",
    "p53 / Mdm2": "ENZ",
    "ERK / DUSP": "ENZ",
    "Wnt / APC-Axin": "SEQ",
    "Notch / FBXW7": "ENZ",
    "Circadian (CLOCK-BMAL1 / PER-CRY)": "SEQ",
    "Hippo / LATS": "ENZ",
    "mTOR / TSC": "ENZ",
    "JAK-STAT / SOCS": "ENZ",
    "TGF-β / SMAD6-7": "SEQ",
    "NRF2 / KEAP1": "SEQ",
    "Hedgehog / SUFU": "ENZ",
    "UPR / IRE1-BiP": "SEQ",
    "Calcium / SERCA": "ENZ",
    "Cell Cycle CKI / CDK": "SEQ",
    "Rb / E2F": "SEQ",
    "HIF / VHL": "SEQ",
    "NFAT / RCAN": "SEQ",
    "PI3K-AKT / PTEN": "ENZ",
    "Myc / FBXW7": "ENZ",
    "Id / bHLH": "SEQ",
    "IκBζ / IL-6": "SEQ",
}

# Cluster from Paper 1
CLUSTER = {
    "NF-κB / IκBα": "I (transcriptional NF)",
    "p53 / Mdm2": "I (transcriptional NF)",
    "ERK / DUSP": "I (transcriptional NF)",
    "Calcium / SERCA": "I (transcriptional NF)",
    "mTOR / TSC": "I (transcriptional NF)",
    "Cell Cycle CKI / CDK": "I (transcriptional NF)",
    "Wnt / APC-Axin": "II (autorepressive)",
    "Notch / FBXW7": "II (autorepressive)",
    "Circadian (CLOCK-BMAL1 / PER-CRY)": "II (autorepressive)",
    "JAK-STAT / SOCS": "II (autorepressive)",
    "Hippo / LATS": "III (shuttling)",
    "TGF-β / SMAD6-7": "III (shuttling)",
    "NRF2 / KEAP1": "III (shuttling)",
    "NFAT / RCAN": "III (shuttling)",
}


DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")


def load_loop_matrix():
    """Load loop ortholog matrix from Phase 1.2."""
    loops = []
    with open(DATA_DIR / "loop_ortholog_matrix.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            loops.append(row)
    return loops


def compute_loop_ages(loops):
    """For each loop, find the most distant species with both components."""
    results = []

    for loop in loops:
        pathway = loop["pathway_name"]

        oldest_species = None
        age_mya = 0

        # Check species from most distant to closest
        for sp in SPECIES_ORDER:
            if sp in loop and loop[sp] == "1":
                oldest_species = sp
                age_mya = DIVERGENCE_TIMES[sp]
                break  # found the most distant, stop

        # If no species found (shouldn't happen — at least mouse should be present)
        if oldest_species is None:
            # Check if any species has both
            for sp in SPECIES_ORDER:
                if sp in loop and loop[sp] == "1":
                    oldest_species = sp
                    age_mya = DIVERGENCE_TIMES[sp]
                    break

        # Confidence level
        n_species = int(loop.get("n_species_both", 0))
        if n_species >= 8:
            confidence = "high"
        elif n_species >= 5:
            confidence = "medium"
        else:
            confidence = "low"

        a_val = A_VALUES.get(pathway, None)
        fb_type = FEEDBACK_TYPE.get(pathway, "unknown")
        cluster = CLUSTER.get(pathway, "unclassified")

        results.append({
            "loop_id": loop["loop_id"],
            "pathway": pathway,
            "oldest_species": oldest_species,
            "age_mya": age_mya,
            "n_species_both": n_species,
            "confidence": confidence,
            "A": a_val,
            "feedback_type": fb_type,
            "cluster": cluster,
        })

    return results


def test_age_vs_A(results):
    """Pre-registered test: Spearman ρ(age, A)."""
    # Filter to loops with both age and A
    paired = [(r["age_mya"], r["A"], r["pathway"], r["cluster"])
              for r in results if r["A"] is not None and r["age_mya"] > 0]

    ages = np.array([p[0] for p in paired])
    As = np.array([p[1] for p in paired])
    names = [p[2] for p in paired]
    clusters = [p[3] for p in paired]

    print(f"\n{'='*60}")
    print("PHASE 1.4: Age vs A Correlation Test")
    print(f"{'='*60}")
    print(f"\nN = {len(paired)} loops with both age and A values")
    print(f"\nPre-registered hypothesis:")
    print(f"  H1: Spearman ρ(age, A) < 0 (older loops → lower A → more asymmetric)")
    print(f"  H0: ρ = 0")
    print(f"  α = 0.05, one-tailed")

    # Spearman correlation
    rho, p_two = stats.spearmanr(ages, As)
    p_one = p_two / 2 if rho < 0 else 1 - p_two / 2  # one-tailed for negative

    print(f"\n--- Main result ---")
    print(f"  Spearman ρ = {rho:.3f}")
    print(f"  p (two-tailed) = {p_two:.4f}")
    print(f"  p (one-tailed, ρ < 0) = {p_one:.4f}")

    # Pearson on log(age) vs A
    log_ages = np.log10(ages)
    r_pearson, p_pearson = stats.pearsonr(log_ages, As)
    print(f"\n  Pearson r(log10(age), A) = {r_pearson:.3f}, p = {p_pearson:.4f}")

    # Power analysis note
    # For N=14, can detect |ρ| > 0.55 at power=0.80
    print(f"\n  Power note: N={len(paired)}, detectable |ρ| > 0.53 at power=0.80")

    # Bootstrap CI
    n_boot = 10000
    rho_boot = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        idx = rng.choice(len(ages), size=len(ages), replace=True)
        r_b, _ = stats.spearmanr(ages[idx], As[idx])
        rho_boot.append(r_b)
    ci_lo, ci_hi = np.percentile(rho_boot, [2.5, 97.5])
    print(f"  Bootstrap 95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")

    # Leave-one-out sensitivity
    print(f"\n--- Leave-one-out sensitivity ---")
    for i in range(len(ages)):
        mask = np.ones(len(ages), dtype=bool)
        mask[i] = False
        r_loo, p_loo = stats.spearmanr(ages[mask], As[mask])
        print(f"  Without {names[i]:40s}: ρ = {r_loo:+.3f}, p = {p_loo:.4f}")

    # Test: age correlates with cluster membership?
    print(f"\n--- Age by cluster ---")
    cluster_ages = {}
    for a, c in zip(ages, clusters):
        cluster_ages.setdefault(c, []).append(a)
    for c, c_ages in sorted(cluster_ages.items()):
        print(f"  {c}: mean age = {np.mean(c_ages):.0f} Mya, "
              f"median = {np.median(c_ages):.0f}, n = {len(c_ages)}")

    # Kruskal-Wallis if ≥3 clusters
    if len(cluster_ages) >= 3:
        groups = [np.array(v) for v in cluster_ages.values()]
        H_stat, p_kw = stats.kruskal(*groups)
        print(f"  Kruskal-Wallis H = {H_stat:.2f}, p = {p_kw:.4f}")

    # Fisher exact: ancient (>700 Mya) = Class I vs young (<700 Mya) = Class III?
    print(f"\n--- Fisher exact: ancient vs Class membership ---")
    ancient_classI = sum(1 for a, c in zip(ages, clusters) if a > 700 and "I (" in c and "III" not in c)
    ancient_notI = sum(1 for a, c in zip(ages, clusters) if a > 700 and ("I (" not in c or "III" in c))
    young_classI = sum(1 for a, c in zip(ages, clusters) if a <= 700 and "I (" in c and "III" not in c)
    young_notI = sum(1 for a, c in zip(ages, clusters) if a <= 700 and ("I (" not in c or "III" in c))

    table = [[ancient_classI, ancient_notI], [young_classI, young_notI]]
    print(f"  2×2 table (ancient>700 × ClassI):")
    print(f"    Ancient ClassI: {ancient_classI}, Ancient other: {ancient_notI}")
    print(f"    Young ClassI: {young_classI}, Young other: {young_notI}")
    OR_fisher, p_fisher = stats.fisher_exact(table)
    print(f"  Fisher exact OR = {OR_fisher:.2f}, p = {p_fisher:.4f}")

    # Conclusion
    print(f"\n{'='*60}")
    if p_one < 0.05:
        print(f"CONCLUSION: H1 SUPPORTED — older loops are more asymmetric (lower A)")
        print(f"  ρ = {rho:.3f}, p_one-tailed = {p_one:.4f}")
    elif p_two < 0.10:
        print(f"CONCLUSION: SUGGESTIVE but not significant at α=0.05")
        print(f"  ρ = {rho:.3f}, p_one-tailed = {p_one:.4f}")
    else:
        print(f"CONCLUSION: H0 NOT REJECTED — no significant age-A correlation")
        print(f"  ρ = {rho:.3f}, p_one-tailed = {p_one:.4f}")
        print(f"  This could be due to low power (N={len(paired)})")
    print(f"{'='*60}")

    return {
        "rho": rho,
        "p_two": p_two,
        "p_one": p_one,
        "r_pearson": r_pearson,
        "p_pearson": p_pearson,
        "ci": (ci_lo, ci_hi),
        "N": len(paired),
    }


def save_dating_results(results, stats_result):
    """Save dating table and statistics."""
    # Dating table
    csv_path = DATA_DIR / "loop_ages.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["loop_id", "pathway", "oldest_species", "age_mya",
                         "n_species_both", "confidence", "A", "feedback_type", "cluster"])
        for r in results:
            writer.writerow([r["loop_id"], r["pathway"], r["oldest_species"],
                            r["age_mya"], r["n_species_both"], r["confidence"],
                            r["A"] if r["A"] is not None else "",
                            r["feedback_type"], r["cluster"]])
    print(f"\nSaved dating results to {csv_path}")

    # Stats summary
    stats_path = DATA_DIR / "age_vs_A_stats.json"
    with open(stats_path, "w") as f:
        json.dump({k: float(v) if isinstance(v, (np.floating, float)) else
                   [float(x) for x in v] if isinstance(v, tuple) else v
                   for k, v in stats_result.items()}, f, indent=2)
    print(f"Saved statistics to {stats_path}")


def main():
    print("=" * 60)
    print("Phase 1.3: Phylogenetic Dating of Feedback Loops")
    print("=" * 60)

    loops = load_loop_matrix()
    results = compute_loop_ages(loops)

    # Print dating table
    print(f"\n{'Pathway':<45} {'Oldest species':<20} {'Age (Mya)':>10} {'N sp':>5} {'A':>6} {'Type':>5}")
    print("-" * 95)
    for r in sorted(results, key=lambda x: -x["age_mya"]):
        a_str = f"{r['A']:.2f}" if r["A"] is not None else "  —"
        print(f"{r['pathway']:<45} {r['oldest_species'] or 'N/A':<20} {r['age_mya']:>10} "
              f"{r['n_species_both']:>5} {a_str:>6} {r['feedback_type']:>5}")

    # Phase 1.4: Age vs A
    stats_result = test_age_vs_A(results)

    # Save
    save_dating_results(results, stats_result)

    print("\nDone!")


if __name__ == "__main__":
    main()
