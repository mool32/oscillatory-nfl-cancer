#!/usr/bin/env python3
"""
Phase 2.2a: Deduplicate, cluster, and analyze NFL motifs.

Steps:
1. Resolve KEGG aliases to HGNC standard symbols
2. Deduplicate identical gene sets
3. Cluster by Jaccard similarity > 0.5 → independent modules
4. Merge with 22 manually curated loops (Variant B: KEGG + manual)
5. Hub gene identification
6. Cancer Gene Census enrichment test
"""

import csv
import json
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from scipy import stats
from itertools import combinations

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")

# ---------------------------------------------------------------------------
# KEGG alias → HGNC standard symbol mapping
# Built from KEGG graphics names observed in our data
# ---------------------------------------------------------------------------
ALIAS_MAP = {
    # STAT/JAK family
    "IMD31A": "STAT3", "IMD31B": "STAT3",  # STAT3 aliases
    "AISIMD": "SOCS3",  # SOCS3 alias
    "CANDF7": "STAT2",
    "STAT91": "STAT1",  # STAT1 alias
    "MGF": "STAT5B",
    "EBP-1": "RELA",  # RELA alias
    "AOVD2": "SMAD7",  # SMAD7 alias
    "BSP-1": "SMAD5", "BSP1": "SMAD5",
    "RIGUI": "PER2",
    "MKP-5": "DUSP1",  # DUSP mapping

    # JAK family → collapse to family representative
    "JAK1": "JAK1", "JAK1B": "JAK1", "JTK10": "TYK2", "JAK2": "JAK2",

    # IRS family
    "IRS1": "IRS1", "HIRS-1": "IRS1", "IRS2": "IRS2",

    # Cyclin E
    "CCNE": "CCNE1", "CCNE1": "CCNE1",

    # PI3K subunits → collapse to pathway level
    "P3R3URF-PIK3R3": "PIK3R3", "LOC110117498-PIK3R3": "PIK3CA",
    "CCM4": "PIK3CA", "CLAPO": "PIK3CA",
    "5295": "PIK3CA", "5296": "PIK3R1", "8503": "PIK3R3",

    # CDK inhibitors
    "CDKN1B": "CDKN1B", "CDKN4": "CDKN2A",

    # MAPK family
    "MAPK1": "MAPK1", "ERK": "MAPK1",
    "MAPK14": "MAPK14", "CSBP": "MAPK14", "CSBP1": "MAPK14", "CSBP2": "MAPK14",

    # SOCS family
    "SOCS1": "SOCS1", "SOCS4": "SOCS4", "SOCS7": "SOCS7",
    "CISH": "CISH", "9306": "SOCS6", "9655": "SOCS5",

    # TGF-β receptors
    "TGFBR1": "TGFBR1", "AAT5": "TGFBR1",
    "BMPR1A": "BMPR1A", "BMPR1B": "BMPR1B", "BMPR2": "BMPR2",

    # Others
    "E2F-1": "E2F1", "DILC": "E2F2", "RBAP1": "E2F3",
    "TFDP1": "TFDP1",
    "FOXO6": "FOXO6", "FOXO1": "FOXO1", "FOXO3": "FOXO3", "FOXO4": "FOXO4",
    "PMAIP1": "PMAIP1", "SIVA1": "SIVA1", "BID": "BID", "BBC3": "BBC3",
    "AXIN1": "AXIN1", "AXIN": "AXIN1",
    "GLI": "GLI2", "GLI3": "GLI3",
    "TPTEP2-CSNK1E": "CSNK1E", "LOC400927-CSNK1E": "CSNK1E",
    "M-CSF-R": "CSF1R",
    "LEPR": "LEPR", "PRLR": "PRLR",
    "TNFAIP3": "TNFAIP3", "TRAF6": "TRAF6",
    "PTK2": "PTK2",
    "INSR": "INSR",
    "MPPH": "AKT1", "MPPH2": "AKT2", "AAKG2": "PRKAG2",
}

# Cancer Gene Census genes (core set — from COSMIC CGC)
CGC_GENES = {
    "TP53", "MDM2", "RB1", "E2F1", "MYC", "PTEN", "PIK3CA",
    "CTNNB1", "AXIN1", "AXIN2", "AKT1", "AKT2", "MTOR",
    "STAT3", "STAT5B", "JAK1", "JAK2", "NOTCH1",
    "SMAD3", "SMAD7", "TGFBR1", "BMPR1A",
    "HIF1A", "VHL", "CDK2", "CDKN1A", "CDKN1B", "CDKN2A",
    "CCNE1", "NFE2L2", "RELA", "NFKBIA",
    "GLI1", "GLI2", "GLI3", "SUFU",
    "SOCS1", "CSF1R", "BCL2", "FOXO1", "FOXO3",
    "MAPK1", "TRAF6", "BID",
    "YAP1", "FBXW7", "NFATC1",
}


def load_motifs():
    """Load motifs from CSV."""
    motifs = []
    with open(DATA_DIR / "kegg_short_nfl_motifs.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            motifs.append({
                "id": row["motif_id"],
                "genes_raw": row["genes"].split("|"),
                "length": int(row["length"]),
                "n_inh": int(row["n_inh"]),
                "n_pw": int(row["n_pathways"]),
            })
    return motifs


def resolve_aliases(motifs):
    """Resolve KEGG aliases to HGNC symbols."""
    resolved = []
    for m in motifs:
        genes_hgnc = []
        for g in m["genes_raw"]:
            hgnc = ALIAS_MAP.get(g, g)  # keep original if not in map
            genes_hgnc.append(hgnc)
        m["genes_hgnc"] = genes_hgnc
        m["gene_set"] = frozenset(genes_hgnc)
        resolved.append(m)
    return resolved


def deduplicate(motifs):
    """Collapse identical gene sets."""
    groups = defaultdict(list)
    for m in motifs:
        groups[m["gene_set"]].append(m)

    unique = []
    for gene_set, members in groups.items():
        # Take the one with most pathway evidence
        best = max(members, key=lambda x: x["n_pw"])
        best["n_variants"] = len(members)
        best["all_raw_names"] = [m["genes_raw"] for m in members]
        unique.append(best)

    return unique


def jaccard(s1, s2):
    """Jaccard similarity between two sets."""
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0


def cluster_by_jaccard(motifs, threshold=0.5):
    """
    Single-linkage clustering by Jaccard similarity.
    Two motifs with J > threshold are in the same cluster.
    Each cluster = one independent module.
    """
    n = len(motifs)
    # Union-Find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i in range(n):
        for j in range(i + 1, n):
            if jaccard(motifs[i]["gene_set"], motifs[j]["gene_set"]) > threshold:
                union(i, j)

    # Group by cluster
    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)

    return clusters


def identify_hubs(motifs):
    """Find genes present in the most motifs."""
    gene_counts = Counter()
    for m in motifs:
        for g in m["gene_set"]:
            gene_counts[g] += 1
    return gene_counts


def cgc_enrichment(motifs, gene_counts):
    """Test if hub genes are enriched for Cancer Gene Census."""
    # All unique genes in motifs
    all_genes = set()
    for m in motifs:
        all_genes.update(m["gene_set"])

    n_total = len(all_genes)
    n_cgc = len(all_genes & CGC_GENES)
    n_non_cgc = n_total - n_cgc

    # Hub genes = those in ≥3 motifs
    hub_genes = {g for g, c in gene_counts.items() if c >= 3}
    n_hub = len(hub_genes)
    n_hub_cgc = len(hub_genes & CGC_GENES)
    n_hub_non_cgc = n_hub - n_hub_cgc

    # Fisher exact
    # 2×2: hub_cgc, hub_non_cgc, non_hub_cgc, non_hub_non_cgc
    non_hub_cgc = n_cgc - n_hub_cgc
    non_hub_non_cgc = n_non_cgc - n_hub_non_cgc

    table = [[n_hub_cgc, n_hub_non_cgc], [non_hub_cgc, non_hub_non_cgc]]
    OR, p_fisher = stats.fisher_exact(table)

    return {
        "n_total_genes": n_total,
        "n_cgc": n_cgc,
        "n_hub": n_hub,
        "n_hub_cgc": n_hub_cgc,
        "hub_genes": sorted(hub_genes),
        "hub_cgc_genes": sorted(hub_genes & CGC_GENES),
        "OR": OR,
        "p_fisher": p_fisher,
        "table": table,
    }


def build_module_table(clusters, motifs):
    """Build table of independent modules with representative motif."""
    modules = []
    for cid, member_indices in sorted(clusters.items(), key=lambda x: -len(x[1])):
        members = [motifs[i] for i in member_indices]

        # All genes across cluster
        all_genes = set()
        for m in members:
            all_genes.update(m["gene_set"])

        # Representative = largest pathway evidence
        rep = max(members, key=lambda x: x["n_pw"])

        # Pathway label
        pathway_label = infer_pathway(all_genes)

        modules.append({
            "module_id": len(modules) + 1,
            "n_motifs": len(members),
            "n_unique_genes": len(all_genes),
            "all_genes": sorted(all_genes),
            "representative": rep["genes_hgnc"],
            "length": rep["length"],
            "n_pw": rep["n_pw"],
            "pathway_label": pathway_label,
            "n_cgc": len(all_genes & CGC_GENES),
        })

    return modules


def infer_pathway(genes):
    """Infer pathway name from gene set."""
    pathway_markers = {
        "JAK-STAT/SOCS": {"STAT3", "SOCS3", "JAK1", "JAK2", "STAT1", "STAT5B", "SOCS1", "CISH"},
        "p53/MDM2": {"TP53", "MDM2"},
        "Rb/E2F": {"RB1", "E2F1", "CCNE1", "TFDP1", "E2F2", "E2F3"},
        "NF-κB": {"RELA", "NFKBIA", "TNFAIP3", "TRAF6"},
        "TGF-β/SMAD": {"SMAD3", "SMAD7", "SMAD2", "SMAD4", "SMAD5", "TGFBR1"},
        "mTOR/IRS": {"MTOR", "IRS1", "IRS2", "PIK3CA", "PIK3R3", "PIK3R1"},
        "MAPK/ERK": {"MAPK1", "MAPK14", "DUSP1"},
        "Cell Cycle/CDK": {"CDK2", "CDKN1A", "CDKN1B", "CDKN2A", "FOXO1", "FOXO3", "FOXO4", "FOXO6"},
        "PI3K/PTEN": {"PTEN", "PTK2", "PIK3CA", "PIK3R3"},
        "Hedgehog/GLI": {"GLI1", "GLI2", "GLI3", "SUFU"},
        "p53/Apoptosis": {"TP53", "BCL2", "PMAIP1", "SIVA1", "BID", "BBC3"},
        "Wnt/β-catenin": {"CTNNB1", "CDH1", "AXIN1"},
        "Hippo/YAP": {"YAP1", "AXIN1", "CSNK1E"},
        "BMP/SMAD": {"BMPR1A", "BMPR1B", "BMPR2", "SMAD3", "MAPK1"},
    }

    best_match = "Unknown"
    best_overlap = 0
    for name, markers in pathway_markers.items():
        overlap = len(genes & markers)
        if overlap > best_overlap:
            best_overlap = overlap
            best_match = name

    return best_match


def main():
    print("=" * 70)
    print("Phase 2.2a: NFL Motif Deduplication, Clustering, and Hub Analysis")
    print("=" * 70)

    # 1. Load
    motifs = load_motifs()
    n_2node = sum(1 for m in motifs if m["length"] == 2)
    n_3node = sum(1 for m in motifs if m["length"] == 3)
    print(f"\nLoaded {len(motifs)} motifs ({n_2node} two-node, {n_3node} three-node)")

    # 2. Resolve aliases
    motifs = resolve_aliases(motifs)

    # Show alias resolution stats
    n_resolved = sum(1 for m in motifs for g in m["genes_raw"] if g in ALIAS_MAP)
    print(f"Resolved {n_resolved} alias references")

    # 3. Deduplicate
    unique = deduplicate(motifs)
    n_2u = sum(1 for m in unique if m["length"] == 2)
    n_3u = sum(1 for m in unique if m["length"] == 3)
    print(f"\nAfter deduplication: {len(unique)} unique motifs ({n_2u} two-node, {n_3u} three-node)")
    print(f"  Collapsed {len(motifs) - len(unique)} duplicate gene sets")

    # Show what collapsed
    multi = [m for m in unique if m["n_variants"] > 1]
    print(f"\n  {len(multi)} motifs had multiple KEGG variants:")
    for m in sorted(multi, key=lambda x: -x["n_variants"])[:10]:
        print(f"    {' | '.join(m['genes_hgnc'])} ({m['n_variants']} variants)")

    # 4. Cluster by Jaccard > 0.5
    print(f"\n--- Clustering (Jaccard > 0.5, single-linkage) ---")
    clusters = cluster_by_jaccard(unique, threshold=0.5)
    print(f"  {len(unique)} unique motifs → {len(clusters)} independent modules")

    # 5. Build module table
    modules = build_module_table(clusters, unique)

    print(f"\n{'='*70}")
    print(f"{'#':>3} {'Pathway':<25} {'Motifs':>6} {'Genes':>5} {'CGC':>4} {'Representative'}")
    print("-" * 90)
    for mod in modules:
        rep_str = " → ".join(mod["representative"][:4])
        print(f"{mod['module_id']:3d} {mod['pathway_label']:<25} "
              f"{mod['n_motifs']:6d} {mod['n_unique_genes']:5d} {mod['n_cgc']:4d} "
              f"{rep_str}")

    # 6. Hub genes
    print(f"\n{'='*70}")
    print("Hub Gene Analysis")
    print("=" * 70)
    gene_counts = identify_hubs(unique)
    print(f"\nTop 20 genes by motif participation:")
    for gene, count in gene_counts.most_common(20):
        cgc_mark = " [CGC]" if gene in CGC_GENES else ""
        print(f"  {gene:15s}: {count:3d} motifs{cgc_mark}")

    # 7. CGC enrichment
    print(f"\n{'='*70}")
    print("Cancer Gene Census Enrichment")
    print("=" * 70)
    cgc_result = cgc_enrichment(unique, gene_counts)
    print(f"\nAll genes in NFL motifs: {cgc_result['n_total_genes']}")
    print(f"  CGC genes: {cgc_result['n_cgc']} ({100*cgc_result['n_cgc']/cgc_result['n_total_genes']:.1f}%)")
    print(f"\nHub genes (≥3 motifs): {cgc_result['n_hub']}")
    print(f"  Hub CGC genes: {cgc_result['n_hub_cgc']} ({100*cgc_result['n_hub_cgc']/cgc_result['n_hub']:.1f}%)")
    print(f"\nFisher exact test (hub × CGC):")
    print(f"  2×2 table: {cgc_result['table']}")
    print(f"  OR = {cgc_result['OR']:.2f}")
    print(f"  p = {cgc_result['p_fisher']:.4f}")

    if cgc_result["p_fisher"] < 0.05:
        print(f"  → SIGNIFICANT: hub genes are enriched for cancer genes")
    else:
        print(f"  → Not significant at α=0.05")

    print(f"\nHub CGC genes: {', '.join(cgc_result['hub_cgc_genes'])}")

    # 8. Loop-loop overlap matrix (for modules)
    print(f"\n{'='*70}")
    print("Module-Module Overlap Matrix")
    print("=" * 70)
    n_mod = len(modules)
    overlap_matrix = np.zeros((n_mod, n_mod))
    for i in range(n_mod):
        for j in range(n_mod):
            gi = set(modules[i]["all_genes"])
            gj = set(modules[j]["all_genes"])
            overlap_matrix[i, j] = jaccard(gi, gj)

    # Print significant overlaps (> 0)
    cross_module_links = []
    for i in range(n_mod):
        for j in range(i+1, n_mod):
            if overlap_matrix[i, j] > 0:
                shared = set(modules[i]["all_genes"]) & set(modules[j]["all_genes"])
                cross_module_links.append({
                    "mod_i": modules[i]["module_id"],
                    "mod_j": modules[j]["module_id"],
                    "pw_i": modules[i]["pathway_label"],
                    "pw_j": modules[j]["pathway_label"],
                    "jaccard": overlap_matrix[i, j],
                    "shared_genes": sorted(shared),
                })

    cross_module_links.sort(key=lambda x: -x["jaccard"])
    print(f"\nCross-module gene sharing ({len(cross_module_links)} pairs with overlap):")
    for link in cross_module_links[:20]:
        print(f"  {link['pw_i']:20s} × {link['pw_j']:20s}  "
              f"J={link['jaccard']:.3f}  shared: {', '.join(link['shared_genes'][:5])}")

    # 9. Save everything
    save_results(modules, gene_counts, cgc_result, cross_module_links, unique)
    print("\nDone!")


def save_results(modules, gene_counts, cgc_result, cross_links, unique_motifs):
    """Save all results."""
    # Module table
    with open(DATA_DIR / "nfl_modules.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["module_id", "pathway_label", "n_motifs", "n_genes",
                     "n_cgc", "all_genes", "representative"])
        for mod in modules:
            w.writerow([mod["module_id"], mod["pathway_label"],
                        mod["n_motifs"], mod["n_unique_genes"], mod["n_cgc"],
                        "|".join(mod["all_genes"]),
                        "|".join(mod["representative"])])

    # Hub genes
    with open(DATA_DIR / "hub_genes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "n_motifs", "is_cgc"])
        for gene, count in gene_counts.most_common():
            w.writerow([gene, count, gene in CGC_GENES])

    # Cross-module links
    with open(DATA_DIR / "module_overlaps.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mod_i", "mod_j", "pw_i", "pw_j", "jaccard", "shared_genes"])
        for link in cross_links:
            w.writerow([link["mod_i"], link["mod_j"], link["pw_i"], link["pw_j"],
                        f"{link['jaccard']:.4f}", "|".join(link["shared_genes"])])

    # Unique motifs (filtered)
    with open(DATA_DIR / "unique_nfl_motifs.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genes_hgnc", "length", "n_inh", "n_pw", "n_variants"])
        for m in unique_motifs:
            w.writerow(["|".join(m["genes_hgnc"]), m["length"],
                        m["n_inh"], m["n_pw"], m["n_variants"]])

    # Summary JSON
    summary = {
        "n_raw_motifs": sum(m["n_variants"] for m in unique_motifs),
        "n_unique_motifs": len(unique_motifs),
        "n_modules": len([m for m in modules]),
        "n_hub_genes": cgc_result["n_hub"],
        "cgc_enrichment_OR": cgc_result["OR"],
        "cgc_enrichment_p": cgc_result["p_fisher"],
    }
    with open(DATA_DIR / "nfl_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved: nfl_modules.csv, hub_genes.csv, module_overlaps.csv, unique_nfl_motifs.csv")


if __name__ == "__main__":
    main()
