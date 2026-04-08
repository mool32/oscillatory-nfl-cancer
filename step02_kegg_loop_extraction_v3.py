#!/usr/bin/env python3
"""
Phase 2.1 v3 — KEGG negative feedback loops via bounded DFS from seed genes.

Strategy: For each seed gene, do bounded DFS (depth ≤ 6) looking for paths
that return to the starting node. This avoids nx.simple_cycles entirely.
Memory-efficient and guaranteed to terminate.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
import csv
import random
import json
import time
import sys

# ---------------------------------------------------------------------------
KGML_DIR = Path("/Users/teo/Desktop/research/oscilatory/data/kegg_kgml")
OUTPUT_DIR = Path("/Users/teo/Desktop/research/paper2/data")
MAX_DEPTH = 6  # max loop length

SEED_GENES = {
    "RELA", "NFKBIA", "TP53", "MDM2", "MAPK1", "DUSP1",
    "CTNNB1", "AXIN2", "NOTCH1", "FBXW7", "ARNTL", "PER2",
    "YAP1", "LATS1", "MTOR", "TSC1", "STAT3", "SOCS3",
    "SMAD3", "SMAD7", "NFE2L2", "KEAP1", "GLI1", "SUFU",
    "ERN1", "HSPA5", "PLCG1", "ATP2A2", "CDK2", "CDKN1A",
    "E2F1", "RB1", "HIF1A", "VHL", "NFATC1", "RCAN1",
    "PIK3CA", "PTEN", "MYC", "TCF3", "ID1", "IL6", "NFKBIZ",
}

ACTIVATION_TYPES = {"activation", "expression", "phosphorylation",
                    "ubiquitination", "glycosylation", "methylation", "indirect effect"}
INHIBITION_TYPES = {"inhibition", "repression", "dephosphorylation", "dissociation"}


def build_global_graph():
    """Parse all KGML → adjacency dict (faster than NetworkX for DFS)."""
    kgml_files = sorted(KGML_DIR.glob("*.kgml"))
    print(f"Found {len(kgml_files)} KGML files")

    # adjacency: node -> [(neighbor, sign, pathway_id), ...]
    adj = defaultdict(list)
    global_sym = {}
    pw_names = {}
    edge_set = set()  # for dedup

    for kf in kgml_files:
        try:
            tree = ET.parse(str(kf))
            root = tree.getroot()
            pid = root.get("name", "").replace("path:", "")
            ptitle = root.get("title", "")
            pw_names[pid] = ptitle

            entry_map = {}
            for entry in root.findall("entry"):
                eid = entry.get("id")
                if entry.get("type") == "gene":
                    kegg_ids = [x.replace("hsa:", "") for x in entry.get("name", "").split()]
                    graphics = entry.find("graphics")
                    display = graphics.get("name", "") if graphics is not None else ""
                    symbols = []
                    if display:
                        for p in display.rstrip(".").split(","):
                            s = p.strip().split("/")[0].strip()
                            if s and len(s) < 20:
                                symbols.append(s)
                    entry_map[eid] = kegg_ids
                    for i, kid in enumerate(kegg_ids):
                        if i < len(symbols):
                            global_sym[kid] = symbols[i]
                        elif symbols:
                            global_sym[kid] = symbols[0]

            for entry in root.findall("entry"):
                eid = entry.get("id")
                if entry.get("type") == "group":
                    comps = []
                    for c in entry.findall("component"):
                        cid = c.get("id")
                        if cid in entry_map:
                            comps.extend(entry_map[cid])
                    if comps:
                        entry_map[eid] = comps

            for rel in root.findall("relation"):
                e1, e2 = rel.get("entry1"), rel.get("entry2")
                if e1 not in entry_map or e2 not in entry_map:
                    continue
                sign = 0
                for st in rel.findall("subtype"):
                    sn = st.get("name")
                    if sn in ACTIVATION_TYPES: sign = max(sign, 1)
                    elif sn in INHIBITION_TYPES: sign = -1
                if sign == 0: sign = 1

                for s in entry_map[e1]:
                    for t in entry_map[e2]:
                        edge_key = (s, t)
                        if edge_key not in edge_set:
                            edge_set.add(edge_key)
                            adj[s].append((t, sign, pid))
                        # Update sign to inhibition if new evidence
                        # (handled by storing multiple entries, pick worst)
        except Exception as e:
            pass

    # Consolidate: for each (s,t), determine final sign
    # If any pathway says inhibition, it's inhibition
    consolidated = defaultdict(list)  # node -> [(neighbor, sign)]
    edge_info = {}  # (s,t) -> {sign, pathways}

    for src, neighbors in adj.items():
        for tgt, sign, pid in neighbors:
            key = (src, tgt)
            if key not in edge_info:
                edge_info[key] = {"sign": sign, "pathways": {pid}}
            else:
                edge_info[key]["pathways"].add(pid)
                if sign == -1:
                    edge_info[key]["sign"] = -1

    # Build final adjacency
    final_adj = defaultdict(list)
    for (src, tgt), info in edge_info.items():
        final_adj[src].append((tgt, info["sign"]))

    n_nodes = len(set(list(final_adj.keys()) + [t for ns in final_adj.values() for t, _ in ns]))
    n_edges = len(edge_info)
    n_inh = sum(1 for v in edge_info.values() if v["sign"] == -1)
    print(f"Global graph: {n_nodes} nodes, {n_edges} edges ({n_inh} inhibitory)")

    return final_adj, edge_info, global_sym, pw_names


def find_seed_ids(global_sym):
    """Map SEED_GENES symbols → KEGG node IDs."""
    sym_to_ids = defaultdict(set)
    for kid, sym in global_sym.items():
        sym_to_ids[sym.upper()].add(kid)

    seed_ids = set()
    found = set()
    for gene in SEED_GENES:
        if gene.upper() in sym_to_ids:
            seed_ids.update(sym_to_ids[gene.upper()])
            found.add(gene)

    print(f"Seed genes in graph: {len(found)}/43")
    missing = SEED_GENES - found
    if missing:
        print(f"  Missing: {', '.join(sorted(missing))}")
    return seed_ids, found


def dfs_find_cycles(adj, start, max_depth):
    """
    Bounded DFS from `start` to find all simple cycles passing through it.
    Returns list of cycles (each cycle = list of node IDs).
    """
    cycles = []

    # Stack: (current_node, path, visited_set)
    stack = [(start, [start], {start})]

    while stack:
        node, path, visited = stack.pop()

        if len(path) > max_depth:
            continue

        for neighbor, sign in adj.get(node, []):
            if neighbor == start and len(path) >= 2:
                # Found a cycle!
                cycles.append(list(path))
            elif neighbor not in visited and len(path) < max_depth:
                stack.append((neighbor, path + [neighbor], visited | {neighbor}))

    return cycles


def main():
    random.seed(42)
    print("=" * 60)
    print("Phase 2.1 v3: KEGG Negative Feedback Loops (bounded DFS)")
    print("=" * 60)

    # Build graph
    final_adj, edge_info, global_sym, pw_names = build_global_graph()

    # Seed nodes
    seed_ids, found_genes = find_seed_ids(global_sym)

    # Find cycles from each seed
    all_loops = {}  # frozenset(cycle) -> info
    t0 = time.time()

    for i, seed in enumerate(sorted(seed_ids)):
        sym = global_sym.get(seed, seed)

        cycles = dfs_find_cycles(final_adj, seed, MAX_DEPTH)

        for cycle in cycles:
            # Count inhibitory edges
            n_inh = 0
            pws = set()
            for j in range(len(cycle)):
                src = cycle[j]
                tgt = cycle[(j + 1) % len(cycle)]
                key = (src, tgt)
                if key in edge_info:
                    if edge_info[key]["sign"] == -1:
                        n_inh += 1
                    pws.update(edge_info[key]["pathways"])

            # Negative feedback = odd inhibitory edges
            if n_inh % 2 == 1:
                dedup = frozenset(cycle)
                if dedup not in all_loops:
                    symbols = [global_sym.get(n, n) for n in cycle]
                    seeds_in = [global_sym.get(n, n) for n in cycle if n in seed_ids]
                    pw_list = [f"{p}:{pw_names.get(p, '?')}" for p in pws]

                    all_loops[dedup] = {
                        "cycle_nodes": cycle,
                        "gene_symbols": symbols,
                        "loop_length": len(cycle),
                        "n_inh": n_inh,
                        "n_act": len(cycle) - n_inh,
                        "pathways": pw_list,
                        "seeds_in_loop": seeds_in,
                    }

        if (i + 1) % 20 == 0 or (i + 1) == len(seed_ids):
            el = time.time() - t0
            print(f"  [{i+1}/{len(seed_ids)}] {sym:10s}: {len(cycles):6d} cycles, "
                  f"unique_nfl={len(all_loops):5d} ({el:.1f}s)")

    elapsed = time.time() - t0
    loops = sorted(all_loops.values(), key=lambda x: (x["loop_length"], x["gene_symbols"]))

    print(f"\n{'='*60}")
    print(f"Total unique negative feedback loops: {len(loops)}")
    print(f"Time: {elapsed:.1f}s")

    # Length distribution
    length_dist = defaultdict(int)
    for l in loops:
        length_dist[l["loop_length"]] += 1
    print("\nLength distribution:")
    for le in sorted(length_dist):
        print(f"  Length {le}: {length_dist[le]}")

    # Validation
    print("\n=== Validation ===")
    known = [
        ("NF-kB/NFKBIA",  {"RELA", "NFKBIA"}),
        ("p53/MDM2",       {"TP53", "MDM2"}),
        ("ERK/DUSP1",      {"MAPK1", "DUSP1"}),
        ("Wnt/AXIN2",      {"CTNNB1", "AXIN2"}),
        ("mTOR/TSC1",      {"MTOR", "TSC1"}),
        ("STAT3/SOCS3",    {"STAT3", "SOCS3"}),
        ("E2F1/RB1",       {"E2F1", "RB1"}),
        ("HIF1A/VHL",      {"HIF1A", "VHL"}),
    ]
    for name, targets in known:
        found = False
        example = None
        for loop in loops:
            syms = {s.upper() for s in loop["gene_symbols"]}
            if all(any(s.startswith(t) for s in syms) for t in targets):
                found = True
                example = loop
                break
        mark = "✓" if found else "✗"
        if example:
            path = " → ".join(example["gene_symbols"][:8])
            print(f"  {mark} {name:20s} (len={example['loop_length']}: {path})")
        else:
            print(f"  {mark} {name}")

    # Random sample
    print(f"\n=== 10 Random Loops ===")
    sample = random.sample(loops, min(10, len(loops)))
    for i, l in enumerate(sample, 1):
        path = " → ".join(l["gene_symbols"][:8])
        print(f"  [{i}] len={l['loop_length']}, inh={l['n_inh']}: {path}")

    # Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = OUTPUT_DIR / "kegg_negative_feedback_loops.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["loop_id", "gene_members", "loop_length", "n_inh", "n_act",
                     "pathways", "seed_genes"])
        for i, l in enumerate(loops, 1):
            w.writerow([f"NFL_{i:04d}", "|".join(l["gene_symbols"]),
                        l["loop_length"], l["n_inh"], l["n_act"],
                        "|".join(l["pathways"]), "|".join(l["seeds_in_loop"])])

    summary = {
        "n_loops": len(loops),
        "length_dist": dict(length_dist),
        "max_depth": MAX_DEPTH,
        "n_seed_genes_found": len(found_genes),
        "elapsed_s": round(elapsed, 1),
    }
    with open(OUTPUT_DIR / "kegg_loop_stats.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Text summary
    with open(OUTPUT_DIR / "kegg_loop_summary.txt", "w") as f:
        f.write(f"KEGG NFL Extraction v3 (bounded DFS, max_depth={MAX_DEPTH})\n")
        f.write(f"Unique negative feedback loops: {len(loops)}\n")
        f.write(f"Time: {elapsed:.1f}s\n\n")
        for le in sorted(length_dist):
            f.write(f"Length {le}: {length_dist[le]}\n")

    print(f"\nSaved {len(loops)} loops → {csv_path}")
    print("Done!")


if __name__ == "__main__":
    main()
