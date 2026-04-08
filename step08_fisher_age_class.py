#!/usr/bin/env python3
"""
Fisher exact test: Ancient vs Young loops × Class I vs Class II/III.

Binary test — more power than Spearman correlation.
Cutoffs tested: 600, 800, 1000 Mya for sensitivity.
"""

import csv
import numpy as np
from scipy import stats
from pathlib import Path

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")


def load_data():
    loops = []
    with open(DATA_DIR / "loop_ages.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a_val = row["A"]
            if not a_val:
                continue  # skip loops without A
            loops.append({
                "pathway": row["pathway"],
                "age": int(row["age_mya"]),
                "A": float(a_val),
                "cluster": row["cluster"],
                "fb_type": row["feedback_type"],
            })
    return loops


def classify(loops, age_cutoff, a_cutoff_low=0.35, a_cutoff_high=0.45):
    """Classify into 2×2: ancient/young × classI/classII+III."""
    # Class I: A < a_cutoff_low (fast rise, strong asymmetry)
    # Class II+III: A >= a_cutoff_high (slower rise, more symmetric)
    # Middle band (0.35-0.45): exclude for clean separation, or include

    results = {"ancient_I": [], "ancient_II": [], "young_I": [], "young_II": []}
    excluded_middle = []

    for l in loops:
        is_ancient = l["age"] > age_cutoff

        if l["A"] < a_cutoff_low:
            a_class = "I"
        elif l["A"] >= a_cutoff_high:
            a_class = "II"
        else:
            excluded_middle.append(l)
            continue

        if is_ancient and a_class == "I":
            results["ancient_I"].append(l)
        elif is_ancient and a_class == "II":
            results["ancient_II"].append(l)
        elif not is_ancient and a_class == "I":
            results["young_I"].append(l)
        else:
            results["young_II"].append(l)

    return results, excluded_middle


def run_fisher(results, age_cutoff, label=""):
    a_I = len(results["ancient_I"])
    a_II = len(results["ancient_II"])
    y_I = len(results["young_I"])
    y_II = len(results["young_II"])

    table = [[a_I, a_II], [y_I, y_II]]

    # Handle zero cells
    if min(a_I + a_II, y_I + y_II, a_I + y_I, a_II + y_II) == 0:
        OR, p = float('inf'), 1.0
    else:
        OR, p = stats.fisher_exact(table, alternative='greater')

    n_total = a_I + a_II + y_I + y_II

    print(f"\n--- Cutoff: {age_cutoff} Mya {label}---")
    print(f"  {'':20s} {'Class I (A<0.35)':>16s} {'Class II+ (A≥0.45)':>18s}")
    print(f"  {'Ancient (>'+str(age_cutoff)+')':20s} {a_I:>16d} {a_II:>18d}")
    print(f"  {'Young (≤'+str(age_cutoff)+')':20s} {y_I:>16d} {y_II:>18d}")
    print(f"  N = {n_total}, Fisher OR = {OR:.2f}, p (one-sided) = {p:.4f}")

    if a_I > 0 and a_II >= 0 and y_I >= 0 and y_II > 0:
        pct_ancient_I = 100 * a_I / (a_I + a_II) if (a_I + a_II) > 0 else 0
        pct_young_I = 100 * y_I / (y_I + y_II) if (y_I + y_II) > 0 else 0
        print(f"  Ancient: {pct_ancient_I:.0f}% Class I | Young: {pct_young_I:.0f}% Class I")

    # List the loops
    for key in ["ancient_I", "ancient_II", "young_I", "young_II"]:
        if results[key]:
            names = [l["pathway"] for l in results[key]]
            print(f"  {key}: {', '.join(names)}")

    return OR, p, table


def run_inclusive(loops, age_cutoff):
    """Include middle band: Class I = A ≤ 0.42, Class II = A > 0.42 (median split)."""
    As = sorted([l["A"] for l in loops])
    median_A = np.median(As)

    a_I, a_II, y_I, y_II = 0, 0, 0, 0
    for l in loops:
        is_ancient = l["age"] > age_cutoff
        is_low_A = l["A"] <= median_A

        if is_ancient and is_low_A: a_I += 1
        elif is_ancient: a_II += 1
        elif is_low_A: y_I += 1
        else: y_II += 1

    table = [[a_I, a_II], [y_I, y_II]]
    OR, p = stats.fisher_exact(table, alternative='greater')

    print(f"\n--- Median split (A ≤ {median_A:.2f} vs > {median_A:.2f}), age cutoff {age_cutoff} ---")
    print(f"  {'':20s} {'Low A':>10s} {'High A':>10s}")
    print(f"  {'Ancient':20s} {a_I:>10d} {a_II:>10d}")
    print(f"  {'Young':20s} {y_I:>10d} {y_II:>10d}")
    print(f"  N = {sum(sum(r) for r in table)}, Fisher OR = {OR:.2f}, p = {p:.4f}")

    return OR, p


def main():
    loops = load_data()
    print("=" * 60)
    print("Fisher Exact: Ancient vs Young × Class I vs Class II+III")
    print("=" * 60)
    print(f"\nLoops with A values: {len(loops)}")

    for l in sorted(loops, key=lambda x: -x["age"]):
        print(f"  {l['pathway']:45s} age={l['age']:5d}  A={l['A']:.2f}  {l['cluster']}")

    # Strict classification (exclude 0.35-0.45 middle band)
    print(f"\n{'='*60}")
    print("STRICT: Class I (A < 0.35) vs Class II+ (A ≥ 0.45)")
    print("Excludes middle band 0.35-0.45")
    print("=" * 60)

    all_results = {}
    for cutoff in [600, 800, 1000]:
        results, excluded = classify(loops, cutoff)
        OR, p, table = run_fisher(results, cutoff)
        all_results[cutoff] = {"OR": OR, "p": p, "table": table}

    if excluded:
        print(f"\n  Middle band excluded ({len(excluded)}): " +
              ", ".join(f"{l['pathway']} (A={l['A']})" for l in excluded))

    # Inclusive: median split
    print(f"\n{'='*60}")
    print("INCLUSIVE: Median split (no exclusion)")
    print("=" * 60)

    for cutoff in [600, 800, 1000]:
        run_inclusive(loops, cutoff)

    # Without NRF2 (sensitivity to outlier)
    print(f"\n{'='*60}")
    print("SENSITIVITY: Without NRF2/KEAP1 (convergent evolution outlier)")
    print("=" * 60)

    loops_no_nrf2 = [l for l in loops if "NRF2" not in l["pathway"]]
    for cutoff in [600, 800, 1000]:
        results, _ = classify(loops_no_nrf2, cutoff)
        run_fisher(results, cutoff, label="(no NRF2) ")

    # Additional: by feedback type
    print(f"\n{'='*60}")
    print("EXPLORATORY: Age × Feedback Type (SEQ vs ENZ)")
    print("=" * 60)

    for cutoff in [800]:
        seq_ancient = sum(1 for l in loops if l["age"] > cutoff and l["fb_type"] == "SEQ")
        enz_ancient = sum(1 for l in loops if l["age"] > cutoff and l["fb_type"] == "ENZ")
        seq_young = sum(1 for l in loops if l["age"] <= cutoff and l["fb_type"] == "SEQ")
        enz_young = sum(1 for l in loops if l["age"] <= cutoff and l["fb_type"] == "ENZ")

        table = [[enz_ancient, seq_ancient], [enz_young, seq_young]]
        OR, p = stats.fisher_exact(table, alternative='two-sided')

        print(f"\n  Cutoff {cutoff} Mya:")
        print(f"  {'':15s} {'ENZ':>6s} {'SEQ':>6s}")
        print(f"  {'Ancient':15s} {enz_ancient:>6d} {seq_ancient:>6d}")
        print(f"  {'Young':15s} {enz_young:>6d} {seq_young:>6d}")
        print(f"  OR = {OR:.2f}, p = {p:.4f}")

    print(f"\n{'='*60}")
    print("SUMMARY")
    print("=" * 60)
    for cutoff, r in all_results.items():
        sig = "***" if r["p"] < 0.01 else "**" if r["p"] < 0.05 else "*" if r["p"] < 0.1 else "ns"
        print(f"  {cutoff} Mya: OR = {r['OR']:.1f}, p = {r['p']:.4f} {sig}")


if __name__ == "__main__":
    main()
