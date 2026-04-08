#!/usr/bin/env python3
"""
Step 03 — Ortholog mapping for 22 feedback loop core components across 15 species.

Queries Ensembl REST API for orthologs of 43 unique human genes,
builds presence/absence matrices per gene and per loop.

Uses two Ensembl Compara databases:
  - Default (vertebrate) compara: vertebrate + C. intestinalis + S. cerevisiae
  - pan_homology compara: adds S. pombe, D. discoideum, A. thaliana, D. melanogaster
    (cross-kingdom comparisons, includes many2many)

Phase 1.2 of the signaling oscillator asymmetry project.
"""

import csv
import json
import os
import sys
import time
import requests
import pandas as pd
from pathlib import Path
from collections import defaultdict

# ── Paths ───────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
LOOP_FILE   = DATA / "loop_components.csv"
CACHE_FILE  = DATA / "ortholog_details.json"
GENE_OUT    = DATA / "gene_orthologs.csv"
LOOP_OUT    = DATA / "loop_ortholog_matrix.csv"

# ── Species of interest (ordered by divergence from human) ──────────────────
SPECIES = [
    "mus_musculus",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "gallus_gallus",
    "xenopus_tropicalis",
    "danio_rerio",
    "ciona_intestinalis",
    "strongylocentrotus_purpuratus",
    "drosophila_melanogaster",
    "caenorhabditis_elegans",
    "nematostella_vectensis",
    "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
    "arabidopsis_thaliana",
    "dictyostelium_discoideum",
]

# Divergence times (Mya) for annotation
DIVERGENCE_MYA = {
    "mus_musculus": 90,
    "rattus_norvegicus": 90,
    "canis_lupus_familiaris": 96,
    "gallus_gallus": 320,
    "xenopus_tropicalis": 352,
    "danio_rerio": 435,
    "ciona_intestinalis": 684,
    "strongylocentrotus_purpuratus": 684,
    "drosophila_melanogaster": 797,
    "caenorhabditis_elegans": 797,
    "nematostella_vectensis": 824,
    "saccharomyces_cerevisiae": 1105,
    "schizosaccharomyces_pombe": 1105,
    "arabidopsis_thaliana": 1500,
    "dictyostelium_discoideum": 1500,
}

# Pretty names for output
SPECIES_LABELS = {
    "mus_musculus": "M.musculus",
    "rattus_norvegicus": "R.norvegicus",
    "canis_lupus_familiaris": "C.l.familiaris",
    "gallus_gallus": "G.gallus",
    "xenopus_tropicalis": "X.tropicalis",
    "danio_rerio": "D.rerio",
    "ciona_intestinalis": "C.intestinalis",
    "strongylocentrotus_purpuratus": "S.purpuratus",
    "drosophila_melanogaster": "D.melanogaster",
    "caenorhabditis_elegans": "C.elegans",
    "nematostella_vectensis": "N.vectensis",
    "saccharomyces_cerevisiae": "S.cerevisiae",
    "schizosaccharomyces_pombe": "S.pombe",
    "arabidopsis_thaliana": "A.thaliana",
    "dictyostelium_discoideum": "D.discoideum",
}

ENSEMBL_REST = "https://rest.ensembl.org"

# Orthology types to include
ALLOWED_OTYPES = {"ortholog_one2one", "ortholog_one2many", "ortholog_many2many"}


def load_loops(path):
    """Load loop_components.csv and return list of loop dicts + set of unique gene IDs."""
    loops = []
    gene_ids = set()
    gene_names = {}  # ensembl_id -> gene_name
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            loops.append(row)
            for role in ("activator", "inhibitor"):
                eid = row[f"{role}_ensembl"].strip()
                gname = row[f"{role}_gene"].strip()
                gene_ids.add(eid)
                gene_names[eid] = gname
    return loops, sorted(gene_ids), gene_names


def load_cache(path):
    """Load cached API responses."""
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(cache, path):
    """Save cache to disk."""
    with open(path, "w") as f:
        json.dump(cache, f, indent=2)


def query_ensembl(ensembl_id, session, compara=None):
    """
    Query Ensembl REST for orthologs of a single gene.
    compara: None for default vertebrate, or 'pan_homology' for cross-kingdom.
    """
    url = f"{ENSEMBL_REST}/homology/id/human/{ensembl_id}"
    params = {}
    if compara:
        params["compara"] = compara
    headers = {"Accept": "application/json"}
    resp = session.get(url, params=params, headers=headers, timeout=60)
    if resp.status_code == 429:
        wait = float(resp.headers.get("Retry-After", 2))
        print(f"    Rate limited, waiting {wait}s...")
        time.sleep(wait)
        return query_ensembl(ensembl_id, session, compara)
    resp.raise_for_status()
    return resp.json()


def extract_orthologs(api_data):
    """
    Extract ortholog info from Ensembl API response.
    Includes one2one, one2many, and many2many orthology types.
    Returns dict: species -> list of {target_id, orthology_type}
    """
    orthologs = defaultdict(list)
    species_set = set(SPECIES)

    for entry in api_data.get("data", []):
        for hom in entry.get("homologies", []):
            otype = hom.get("type", "")
            if otype not in ALLOWED_OTYPES:
                continue
            target = hom.get("target", {})
            sp = target.get("species", "").lower()
            if sp in species_set:
                orthologs[sp].append({
                    "target_id": target.get("id", ""),
                    "orthology_type": otype,
                })
    return dict(orthologs)


def merge_orthologs(existing, new):
    """Merge new ortholog dict into existing, avoiding duplicates by target_id."""
    for sp, orth_list in new.items():
        if sp not in existing:
            existing[sp] = orth_list
        else:
            existing_ids = {o["target_id"] for o in existing[sp]}
            for o in orth_list:
                if o["target_id"] not in existing_ids:
                    existing[sp].append(o)
                    existing_ids.add(o["target_id"])
    return existing


def fetch_all_orthologs(gene_ids, gene_names, cache):
    """
    Fetch orthologs for all genes using two compara passes:
    1. Default (vertebrate) compara — covers vertebrates + C. intestinalis + S. cerevisiae
    2. pan_homology compara — covers cross-kingdom (fungi, plants, protists, invertebrates)
    """
    session = requests.Session()
    session.headers.update({"User-Agent": "ortholog-mapping/1.0 (research)"})

    # Determine which genes need querying
    # Cache key includes "_v2" to distinguish from old single-pass cache
    to_query = [g for g in gene_ids if g not in cache or "compara_version" not in cache.get(g, {})]
    cached_count = len(gene_ids) - len(to_query)
    if cached_count > 0:
        print(f"  {cached_count} genes already cached (v2), {len(to_query)} to query")

    if not to_query:
        print("  All genes cached!")
        return cache

    for i, gid in enumerate(to_query):
        gname = gene_names.get(gid, gid)
        print(f"  [{i+1}/{len(to_query)}] Querying {gname} ({gid})...")

        ortho = {}
        errors = []

        # Pass 1: Default vertebrate compara
        try:
            raw1 = query_ensembl(gid, session, compara=None)
            o1 = extract_orthologs(raw1)
            ortho = merge_orthologs(ortho, o1)
            print(f"    vertebrate: {sum(len(v) for v in o1.values())} orthologs in {len(o1)} species")
        except Exception as e:
            errors.append(f"vertebrate: {e}")
            print(f"    vertebrate ERROR: {e}")

        time.sleep(0.15)

        # Pass 2: Pan-compara (cross-kingdom)
        try:
            raw2 = query_ensembl(gid, session, compara="pan_homology")
            o2 = extract_orthologs(raw2)
            ortho = merge_orthologs(ortho, o2)
            new_sp = set(o2.keys()) - set(o1.keys()) if 'o1' in dir() else set(o2.keys())
            print(f"    pan_homology: {sum(len(v) for v in o2.values())} orthologs in {len(o2)} species"
                  f" (+{len(new_sp)} new species)")
        except Exception as e:
            errors.append(f"pan_homology: {e}")
            print(f"    pan_homology ERROR: {e}")

        cache[gid] = {
            "gene_name": gname,
            "orthologs": ortho,
            "compara_version": "v2_dual",
            "query_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        if errors:
            cache[gid]["errors"] = errors

        # Save cache after each gene (resumable)
        save_cache(cache, CACHE_FILE)

        # Rate limiting
        if i < len(to_query) - 1:
            time.sleep(0.15)

    return cache


def build_gene_matrix(gene_ids, gene_names, cache):
    """Build genes x species presence/absence matrix."""
    rows = []
    for gid in gene_ids:
        gname = gene_names.get(gid, gid)
        entry = cache.get(gid, {})
        ortho = entry.get("orthologs", {})
        row = {"ensembl_id": gid, "gene_name": gname}
        for sp in SPECIES:
            row[SPECIES_LABELS[sp]] = 1 if sp in ortho and len(ortho[sp]) > 0 else 0
        # Total species count
        row["n_species"] = sum(1 for sp in SPECIES if sp in ortho and len(ortho[sp]) > 0)
        rows.append(row)

    df = pd.DataFrame(rows)
    # Sort by gene name
    df = df.sort_values("gene_name").reset_index(drop=True)
    return df


def build_loop_matrix(loops, cache):
    """Build loops x species both-present matrix."""
    rows = []
    for loop in loops:
        act_id = loop["activator_ensembl"].strip()
        inh_id = loop["inhibitor_ensembl"].strip()
        act_ortho = cache.get(act_id, {}).get("orthologs", {})
        inh_ortho = cache.get(inh_id, {}).get("orthologs", {})

        row = {
            "loop_id": loop["loop_id"],
            "pathway_name": loop["pathway_name"],
            "activator": loop["activator_gene"],
            "inhibitor": loop["inhibitor_gene"],
        }
        for sp in SPECIES:
            act_present = sp in act_ortho and len(act_ortho[sp]) > 0
            inh_present = sp in inh_ortho and len(inh_ortho[sp]) > 0
            row[SPECIES_LABELS[sp]] = 1 if (act_present and inh_present) else 0
        # Total species count where both present
        row["n_species_both"] = sum(
            1 for sp in SPECIES
            if sp in act_ortho and len(act_ortho[sp]) > 0
            and sp in inh_ortho and len(inh_ortho[sp]) > 0
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    return df


def main():
    print("=" * 70)
    print("Step 03 — Ortholog Mapping for Feedback Loop Components")
    print("  Dual compara: vertebrate + pan_homology (cross-kingdom)")
    print("  Orthology types: one2one, one2many, many2many")
    print("=" * 70)

    # 1. Load loops
    print("\n1. Loading loop components...")
    loops, gene_ids, gene_names = load_loops(LOOP_FILE)
    print(f"   {len(loops)} loops, {len(gene_ids)} unique genes")

    # 2. Load cache
    print("\n2. Loading cache...")
    cache = load_cache(CACHE_FILE)
    v2_count = sum(1 for g in cache.values() if g.get("compara_version") == "v2_dual")
    print(f"   {len(cache)} genes in cache ({v2_count} with dual compara)")

    # 3. Fetch orthologs
    print("\n3. Fetching orthologs from Ensembl REST API...")
    cache = fetch_all_orthologs(gene_ids, gene_names, cache)
    save_cache(cache, CACHE_FILE)
    print(f"   Done. {len(cache)} genes in cache.")

    # 4. Build gene presence/absence matrix
    print("\n4. Building gene x species matrix...")
    gene_df = build_gene_matrix(gene_ids, gene_names, cache)
    gene_df.to_csv(GENE_OUT, index=False)
    print(f"   Saved to {GENE_OUT}")

    # Count orthologs per species
    sp_cols = [SPECIES_LABELS[sp] for sp in SPECIES]
    print("\n   Orthologs per species:")
    for sp in SPECIES:
        col = SPECIES_LABELS[sp]
        count = gene_df[col].sum()
        mya = DIVERGENCE_MYA[sp]
        print(f"     {col:20s} ({mya:>5d} Mya): {count:3d}/{len(gene_ids)} genes")

    # 5. Build loop both-present matrix
    print("\n5. Building loop x species both-present matrix...")
    loop_df = build_loop_matrix(loops, cache)
    loop_df.to_csv(LOOP_OUT, index=False)
    print(f"   Saved to {LOOP_OUT}")

    # Count conserved loops per species
    print("\n   Conserved loops per species:")
    for sp in SPECIES:
        col = SPECIES_LABELS[sp]
        count = loop_df[col].sum()
        mya = DIVERGENCE_MYA[sp]
        print(f"     {col:20s} ({mya:>5d} Mya): {count:3d}/{len(loops)} loops")

    # 6. Validation checks
    print("\n6. Validation checks:")

    def check_gene(label, gid, species):
        present = species in cache.get(gid, {}).get("orthologs", {})
        ortho_list = cache.get(gid, {}).get("orthologs", {}).get(species, [])
        ids = [o["target_id"] for o in ortho_list] if ortho_list else []
        return present, ids

    # CDK2/CDKN1A in yeast
    cdk2_yeast, cdk2_ids = check_gene("CDK2", "ENSG00000123374", "saccharomyces_cerevisiae")
    cdkn1a_yeast, _ = check_gene("CDKN1A", "ENSG00000124762", "saccharomyces_cerevisiae")
    print(f"   CDK2 in S.cerevisiae:   {cdk2_yeast} {cdk2_ids}")
    print(f"   CDKN1A in S.cerevisiae: {cdkn1a_yeast}")

    # RELA/NFKBIA absent in yeast (NF-κB is metazoan)
    rela_yeast, _ = check_gene("RELA", "ENSG00000173039", "saccharomyces_cerevisiae")
    nfkbia_yeast, _ = check_gene("NFKBIA", "ENSG00000100906", "saccharomyces_cerevisiae")
    print(f"   RELA in S.cerevisiae:   {rela_yeast} (expected: False)")
    print(f"   NFKBIA in S.cerevisiae: {nfkbia_yeast} (expected: False)")

    # TP53 and MDM2 in C. elegans
    tp53_ce, tp53_ids = check_gene("TP53", "ENSG00000141510", "caenorhabditis_elegans")
    mdm2_ce, _ = check_gene("MDM2", "ENSG00000135679", "caenorhabditis_elegans")
    print(f"   TP53 in C.elegans:      {tp53_ce} {tp53_ids}")
    print(f"   MDM2 in C.elegans:      {mdm2_ce} (expected: False -> loop absent)")

    # MTOR in S. pombe (should be present — TOR is ancient)
    mtor_sp, mtor_ids = check_gene("MTOR", "ENSG00000198793", "schizosaccharomyces_pombe")
    print(f"   MTOR in S.pombe:        {mtor_sp} {mtor_ids}")

    # HSPA5 (BiP) in yeast (should be present — ancient chaperone)
    hspa5_yeast, hspa5_ids = check_gene("HSPA5", "ENSG00000044574", "saccharomyces_cerevisiae")
    print(f"   HSPA5 in S.cerevisiae:  {hspa5_yeast} {hspa5_ids}")

    # Summary
    print("\n" + "=" * 70)
    total_loops = len(loops)
    print(f"Summary:")
    for sp in SPECIES:
        col = SPECIES_LABELS[sp]
        gene_count = gene_df[col].sum()
        loop_count = loop_df[col].sum()
        mya = DIVERGENCE_MYA[sp]
        print(f"  {col:20s} ({mya:>5d} Mya): {gene_count:2d}/43 genes, {loop_count:2d}/{total_loops} loops")

    # Conservation gradient
    print(f"\nConservation gradient:")
    vertebrate_conserved = loop_df["D.rerio"].sum()
    metazoan_conserved = loop_df["D.melanogaster"].sum()
    yeast_conserved = loop_df["S.cerevisiae"].sum()
    print(f"  Vertebrate (D. rerio, 435 Mya):       {vertebrate_conserved}/{total_loops} loops")
    print(f"  Metazoan (D. melanogaster, 797 Mya):   {metazoan_conserved}/{total_loops} loops")
    print(f"  Eukaryotic (S. cerevisiae, 1105 Mya):  {yeast_conserved}/{total_loops} loops")
    print("=" * 70)


if __name__ == "__main__":
    main()
