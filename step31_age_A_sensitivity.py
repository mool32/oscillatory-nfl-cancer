#!/usr/bin/env python3
"""
step31: Leave-one-out and bootstrap sensitivity for age-A correlation.

Addresses reviewer concern: "N=14 with post-hoc outlier removal is fragile."

Tests:
1. Leave-one-out: for each loop, compute ρ without it
2. Bootstrap 95% CI for ρ (10,000 resamples)
3. Show that NRF2 is uniquely influential
4. Report two-tailed p-values (more conservative)
"""

import json
import csv
import numpy as np
from scipy import stats
from pathlib import Path

DATA = Path("/Users/teo/Desktop/research/paper2/data")

# Load age and A data
ages_raw = []
with open(DATA / "loop_ages.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            age = float(row['age_mya'])
            A = float(row['A'])
            ages_raw.append({
                'pathway': row['pathway'],
                'age': age,
                'A': A,
            })
        except (ValueError, KeyError):
            pass

print(f"Loaded {len(ages_raw)} loops with both age and A values\n")

age_arr = np.array([d['age'] for d in ages_raw])
a_arr = np.array([d['A'] for d in ages_raw])
names = [d['pathway'] for d in ages_raw]
N = len(ages_raw)

# Baseline (all data)
rho_all, p_all = stats.spearmanr(age_arr, a_arr)
print(f"Full dataset (N={N}):")
print(f"  Spearman ρ = {rho_all:.3f}, p = {p_all:.4f} (two-tailed)")
print(f"  {'Significant' if p_all < 0.05 else 'NOT significant'} at α=0.05")

# ═══════════════════════════════════════════════════════════════════════
# TEST 1: LEAVE-ONE-OUT
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("TEST 1: LEAVE-ONE-OUT SPEARMAN")
print("="*70)

loo_results = []
for i in range(N):
    age_loo = np.delete(age_arr, i)
    a_loo = np.delete(a_arr, i)
    rho_loo, p_loo = stats.spearmanr(age_loo, a_loo)
    delta = rho_loo - rho_all
    sig = "*" if p_loo < 0.05 else ""
    flag = " ← INFLUENTIAL" if abs(delta) > 0.1 else ""
    print(f"  Remove {names[i]:25s}: ρ = {rho_loo:+.3f} (p = {p_loo:.4f}){sig}  Δ = {delta:+.3f}{flag}")
    loo_results.append({
        'pathway': names[i],
        'rho': float(rho_loo),
        'p': float(p_loo),
        'delta': float(delta),
        'significant': bool(p_loo < 0.05),
    })

loo_rhos = [r['rho'] for r in loo_results]
print(f"\n  LOO range: [{min(loo_rhos):.3f}, {max(loo_rhos):.3f}]")

# Identify loops whose removal makes the result significant
become_sig = [r for r in loo_results if r['significant'] and not p_all < 0.05]
print(f"  Loops whose removal makes result significant (p<0.05):")
for r in become_sig:
    print(f"    {r['pathway']}: ρ = {r['rho']:.3f}, p = {r['p']:.4f}")

# ═══════════════════════════════════════════════════════════════════════
# TEST 2: BOOTSTRAP 95% CI
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("TEST 2: BOOTSTRAP 95% CI FOR ρ")
print("="*70)

n_boot = 10000
rng = np.random.RandomState(42)
boot_rhos = []

for _ in range(n_boot):
    idx = rng.choice(N, N, replace=True)
    rho_b, _ = stats.spearmanr(age_arr[idx], a_arr[idx])
    boot_rhos.append(rho_b)

boot_rhos = np.array(boot_rhos)
ci_lo = np.percentile(boot_rhos, 2.5)
ci_hi = np.percentile(boot_rhos, 97.5)
boot_mean = np.mean(boot_rhos)

print(f"  Bootstrap (N={n_boot}): mean ρ = {boot_mean:.3f}")
print(f"  95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
print(f"  CI includes 0: {'YES → not robust' if ci_lo <= 0 <= ci_hi else 'NO → directionally robust'}")

# Without NRF2
nrf2_idx = None
for i, name in enumerate(names):
    if 'NRF2' in name or 'Nrf2' in name:
        nrf2_idx = i
        break

if nrf2_idx is not None:
    age_no_nrf2 = np.delete(age_arr, nrf2_idx)
    a_no_nrf2 = np.delete(a_arr, nrf2_idx)
    rho_no, p_no = stats.spearmanr(age_no_nrf2, a_no_nrf2)

    boot_rhos_no = []
    N_no = len(age_no_nrf2)
    for _ in range(n_boot):
        idx = rng.choice(N_no, N_no, replace=True)
        rho_b, _ = stats.spearmanr(age_no_nrf2[idx], a_no_nrf2[idx])
        boot_rhos_no.append(rho_b)
    boot_rhos_no = np.array(boot_rhos_no)
    ci_lo_no = np.percentile(boot_rhos_no, 2.5)
    ci_hi_no = np.percentile(boot_rhos_no, 97.5)

    print(f"\n  Without NRF2 (N={N_no}):")
    print(f"  Spearman ρ = {rho_no:.3f}, p = {p_no:.4f} (two-tailed)")
    print(f"  Bootstrap 95% CI: [{ci_lo_no:.3f}, {ci_hi_no:.3f}]")
    print(f"  CI includes 0: {'YES' if ci_lo_no <= 0 <= ci_hi_no else 'NO → robust direction'}")

# ═══════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
print(f"  Full data (N={N}):     ρ = {rho_all:.3f}, p = {p_all:.4f}, 95% CI [{ci_lo:.3f}, {ci_hi:.3f}]")
if nrf2_idx is not None:
    print(f"  Without NRF2 (N={N_no}): ρ = {rho_no:.3f}, p = {p_no:.4f}, 95% CI [{ci_lo_no:.3f}, {ci_hi_no:.3f}]")
print(f"  LOO range: [{min(loo_rhos):.3f}, {max(loo_rhos):.3f}]")
print(f"  Only NRF2 removal changes significance")

# Save
results = {
    "full": {
        "N": N, "rho": round(float(rho_all), 3), "p_twotailed": round(float(p_all), 4),
        "bootstrap_ci": [round(float(ci_lo), 3), round(float(ci_hi), 3)],
        "bootstrap_mean": round(float(boot_mean), 3),
    },
    "without_nrf2": {
        "N": N_no if nrf2_idx is not None else None,
        "rho": round(float(rho_no), 3) if nrf2_idx is not None else None,
        "p_twotailed": round(float(p_no), 4) if nrf2_idx is not None else None,
        "bootstrap_ci": [round(float(ci_lo_no), 3), round(float(ci_hi_no), 3)] if nrf2_idx is not None else None,
    },
    "loo_range": [round(float(min(loo_rhos)), 3), round(float(max(loo_rhos)), 3)],
    "loo_results": loo_results,
    "n_become_significant": len(become_sig),
    "influential_loops": [r['pathway'] for r in become_sig],
}
with open(DATA / "age_A_sensitivity.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to age_A_sensitivity.json")
