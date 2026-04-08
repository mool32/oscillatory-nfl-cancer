#!/usr/bin/env python3
"""
Phase 2.1 v2 — Extract negative feedback loops from KEGG KGML files.

Strategy: instead of enumerating ALL cycles on the global graph (intractable),
we use a targeted approach:
  1. Build global signed DiGraph from all KGML files
  2. For each of our 43 seed genes: extract k-hop ego subgraph (k=5)
  3. Find simple cycles within each ego subgraph (tractable)
  4. Filter for negative feedback (odd inhibitory edges)
  5. Deduplicate and merge

This captures all loops of length ≤10 that pass through at least one seed gene,
which is exactly what we need for Paper 2.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
import csv
import random
import sys
import time

import networkx as nx

# ---------------------------------------------------------------------------
KGML_DIR = Path("/Users/teo/Desktop/research/oscilatory/data/kegg_kgml")
OUTPUT_DIR = Path("/Users/teo/Desktop/research/paper2/data")
MAX_CYCLE_LENGTH = 8  # conservative
EGO_RADIUS = 5        # k-hop neighborhood

# Our 43 seed genes (from Phase 1.1)
SEED_GENES = {
    "RELA", "NFKBIA", "TP53", "MDM2", "MAPK1", "DUSP1",
    "CTNNB1", "AXIN2", "NOTCH1", "FBXW7", "ARNTL", "PER2",
    "YAP1", "LATS1", "MTOR", "TSC1", "STAT3", "SOCS3",
    "SMAD3", "SMAD7", "NFE2L2", "KEAP1", "GLI1", "SUFU",
    "ERN1", "HSPA5", "PLCG1", "ATP2A2", "CDK2", "CDKN1A",
    "E2F1", "RB1", "HIF1A", "VHL", "NFATC1", "RCAN1",
    "PIK3CA", "PTEN", "MYC", "TCF3", "ID1", "IL6", "NFKBIZ",
}

ACTIVATION_TYPES = {
    "activation", "expression", "phosphorylation",
    "ubiquitination", "glycosylation", "methylation",
    "indirect effect",
}
INHIBITION_TYPES = {
    "inhibition", "repression", "dephosphorylation",
    "dissociation",
}


def parse_kgml(kgml_path):
    """Parse a single KGML file. Returns edges and symbol mappings."""
    tree = ET.parse(str(kgml_path))
    root = tree.getroot()
    pathway_id = root.get("name", "").replace("path:", "")
    pathway_title = root.get("title", "")

    entry_map = {}
    id_to_symbols = {}

    for entry in root.findall("entry"):
        eid = entry.get("id")
        etype = entry.get("type")
        ename = entry.get("name", "")

        if etype == "gene":
            kegg_ids = [x.replace("hsa:", "") for x in ename.split()]
            graphics = entry.find("graphics")
            display = graphics.get("name", "") if graphics is not None else ""

            symbols = []
            if display:
                for p in display.rstrip(".").split(","):
                    sym = p.strip().split("/")[0].strip()
                    if sym and len(sym) < 20 and not sym.startswith("MGC:"):
                        symbols.append(sym)

            entry_map[eid] = {"kegg_ids": kegg_ids, "symbols": symbols}

            if len(symbols) >= len(kegg_ids):
                for i, kid in enumerate(kegg_ids):
                    id_to_symbols[kid] = symbols[i] if i < len(symbols) else symbols[0]
            elif symbols:
                for kid in kegg_ids:
                    id_to_symbols[kid] = symbols[0]

    # Groups
    for entry in root.findall("entry"):
        eid = entry.get("id")
        if entry.get("type") == "group":
            comps = []
            for c in entry.findall("component"):
                cid = c.get("id")
                if cid in entry_map:
                    comps.extend(entry_map[cid]["kegg_ids"])
            if comps:
                entry_map[eid] = {"kegg_ids": comps, "symbols": [id_to_symbols.get(k, k) for k in comps]}

    # Relations -> edges
    edges = []
    for rel in root.findall("relation"):
        e1, e2 = rel.get("entry1"), rel.get("entry2")
        if e1 not in entry_map or e2 not in entry_map:
            continue

        sign = 0
        for st in rel.findall("subtype"):
            sn = st.get("name")
            if sn in ACTIVATION_TYPES:
                sign = max(sign, 1)
            elif sn in INHIBITION_TYPES:
                sign = -1

        if sign == 0:
            sign = 1

        for s in entry_map[e1]["kegg_ids"]:
            for t in entry_map[e2]["kegg_ids"]:
                edges.append((s, t, sign))

    return pathway_id, pathway_title, edges, id_to_symbols


def build_global_graph():
    """Parse all KGML and build global DiGraph."""
    kgml_files = sorted(KGML_DIR.glob("*.kgml"))
    print(f"Found {len(kgml_files)} KGML files")

    G = nx.DiGraph()
    global_sym = {}
    pw_names = {}

    for kf in kgml_files:
        try:
            pid, ptitle, edges, id2sym = parse_kgml(kf)
            pw_names[pid] = ptitle
            global_sym.update(id2sym)

            for s, t, sign in edges:
                if G.has_edge(s, t):
                    G[s][t]["pathways"].add(pid)
                    if sign == -1:
                        G[s][t]["sign"] = -1
                else:
                    G.add_edge(s, t, sign=sign, pathways={pid})

                for n in [s, t]:
                    if n in id2sym:
                        G.nodes[n]["symbol"] = id2sym[n]
        except Exception as e:
            print(f"  WARN: {kf.name}: {e}")

    n_act = sum(1 for _, _, d in G.edges(data=True) if d["sign"] == 1)
    n_inh = sum(1 for _, _, d in G.edges(data=True) if d["sign"] == -1)
    print(f"Global graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  Activation: {n_act}, Inhibition: {n_inh}")

    return G, global_sym, pw_names


def find_seed_node_ids(G, global_sym):
    """Map seed gene symbols to KEGG node IDs in the graph."""
    # Reverse map: symbol -> node_id(s)
    sym_to_nodes = defaultdict(set)
    for nid in G.nodes():
        sym = global_sym.get(nid, G.nodes[nid].get("symbol", ""))
        if sym:
            sym_to_nodes[sym.upper()].add(nid)

    seed_ids = set()
    found = set()
    missing = set()

    for gene in SEED_GENES:
        key = gene.upper()
        if key in sym_to_nodes:
            seed_ids.update(sym_to_nodes[key])
            found.add(gene)
        else:
            missing.add(gene)

    print(f"\nSeed genes found in KEGG graph: {len(found)}/43")
    if missing:
        print(f"  Missing: {', '.join(sorted(missing))}")

    return seed_ids, found, missing


def find_loops_from_seeds(G, seed_ids, global_sym, pw_names):
    """
    For each seed gene, extract ego subgraph and find cycles.
    This is tractable because ego subgraphs are small (~100-500 nodes).
    """
    all_loops = {}  # dedup_key -> loop_info
    total_cycles = 0

    print(f"\nSearching for negative feedback loops from {len(seed_ids)} seed nodes...")
    print(f"  Ego radius: {EGO_RADIUS} hops, max cycle length: {MAX_CYCLE_LENGTH}")

    t0 = time.time()

    for i, seed in enumerate(sorted(seed_ids)):
        sym = global_sym.get(seed, seed)

        # Get ego subgraph (undirected distance for neighborhood, but keep directed edges)
        # Use BFS on undirected view to find k-hop neighborhood
        G_undir = G.to_undirected(as_view=True)
        try:
            neighbors = set(nx.single_source_shortest_path_length(G_undir, seed, cutoff=EGO_RADIUS).keys())
        except nx.NetworkXError:
            continue

        if len(neighbors) < 2:
            continue

        # Extract directed subgraph
        subG = G.subgraph(neighbors).copy()

        # Remove self-loops for cycle finding (we'll count them separately)
        self_loops_here = [(u, v) for u, v in subG.edges() if u == v]
        subG.remove_edges_from(self_loops_here)

        if subG.number_of_edges() == 0:
            continue

        # Find cycles in this subgraph
        try:
            cycles = list(nx.simple_cycles(subG, length_bound=MAX_CYCLE_LENGTH))
        except TypeError:
            # Older NetworkX
            cycles = [c for c in nx.simple_cycles(subG) if len(c) <= MAX_CYCLE_LENGTH]

        # Filter for negative feedback
        for cycle in cycles:
            n_inh = 0
            edge_pws = set()

            for j in range(len(cycle)):
                src = cycle[j]
                tgt = cycle[(j + 1) % len(cycle)]
                if G.has_edge(src, tgt):
                    ed = G[src][tgt]
                    if ed["sign"] == -1:
                        n_inh += 1
                    edge_pws.update(ed.get("pathways", set()))

            if n_inh % 2 == 1:  # negative feedback
                key = tuple(sorted(cycle))
                if key not in all_loops:
                    symbols = [global_sym.get(n, n) for n in cycle]
                    pw_list = [f"{p}:{pw_names.get(p, p)}" for p in edge_pws]

                    # Check which seed genes are in this loop
                    seeds_in_loop = [global_sym.get(n, n) for n in cycle if n in seed_ids]

                    all_loops[key] = {
                        "cycle_nodes": cycle,
                        "gene_symbols": symbols,
                        "loop_length": len(cycle),
                        "n_inhibitory_edges": n_inh,
                        "n_activating_edges": len(cycle) - n_inh,
                        "pathway_memberships": pw_list,
                        "seeds_in_loop": seeds_in_loop,
                    }

        total_cycles += len(cycles)

        if (i + 1) % 10 == 0 or (i + 1) == len(seed_ids):
            elapsed = time.time() - t0
            print(f"  [{i+1}/{len(seed_ids)}] {sym:10s} ego={len(neighbors):5d} nodes, "
                  f"cycles={len(cycles):6d}, neg_fb_unique={len(all_loops):5d} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nTotal cycles examined: {total_cycles}")
    print(f"Unique negative feedback loops: {len(all_loops)}")
    print(f"Time: {elapsed:.1f}s")

    # Also count autorepressive self-loops with inhibitory edge
    self_loops = []
    for nid in seed_ids:
        if G.has_edge(nid, nid) and G[nid][nid]["sign"] == -1:
            sym = global_sym.get(nid, nid)
            pws = G[nid][nid].get("pathways", set())
            pw_list = [f"{p}:{pw_names.get(p, p)}" for p in pws]
            self_loops.append({
                "cycle_nodes": [nid],
                "gene_symbols": [sym],
                "loop_length": 1,
                "n_inhibitory_edges": 1,
                "n_activating_edges": 0,
                "pathway_memberships": pw_list,
                "seeds_in_loop": [sym],
            })
    print(f"Self-inhibitory loops (seed genes): {len(self_loops)}")

    loops_list = list(all_loops.values()) + self_loops
    loops_list.sort(key=lambda x: (x["loop_length"], x["gene_symbols"]))
    return loops_list, elapsed


def validate_known_loops(loops, global_sym):
    """Validate against known negative feedback pairs."""
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
        ("MYC/FBXW7",      {"MYC", "FBXW7"}),
    ]

    for name, targets in known:
        found = False
        example = None
        for loop in loops:
            syms_upper = {s.upper() for s in loop["gene_symbols"]}
            if targets.issubset(syms_upper):
                found = True
                example = loop
                break
            # Also check prefix match (gene families)
            if all(any(s.startswith(t) for s in syms_upper) for t in targets):
                found = True
                example = loop
                break

        status = "✓ FOUND" if found else "✗ MISSING"
        if example:
            path = " → ".join(example["gene_symbols"][:8])
            extra = "..." if len(example["gene_symbols"]) > 8 else ""
            print(f"  {status} {name:20s} (len={example['loop_length']}: {path}{extra})")
        else:
            print(f"  {status} {name}")


def save_results(loops, G, global_sym, pw_names, elapsed):
    """Save CSV + summary."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # CSV
    csv_path = OUTPUT_DIR / "kegg_negative_feedback_loops.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["loop_id", "gene_members", "loop_length", "n_inh", "n_act",
                     "pathways", "seed_genes_in_loop"])
        for i, loop in enumerate(loops, 1):
            w.writerow([
                f"NFL_{i:04d}",
                "|".join(loop["gene_symbols"]),
                loop["loop_length"],
                loop["n_inhibitory_edges"],
                loop["n_activating_edges"],
                "|".join(loop["pathway_memberships"]),
                "|".join(loop["seeds_in_loop"]),
            ])
    print(f"\nSaved {len(loops)} loops → {csv_path}")

    # Summary
    length_dist = defaultdict(int)
    for l in loops:
        length_dist[l["loop_length"]] += 1

    summary_path = OUTPUT_DIR / "kegg_loop_summary.txt"
    with open(summary_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("KEGG Negative Feedback Loop Extraction (v2 — ego-centric)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Global graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges\n")
        f.write(f"Seed genes: 43, found in graph: {sum(1 for n in G.nodes() if global_sym.get(n, G.nodes[n].get('symbol', '')).upper() in {g.upper() for g in SEED_GENES})}\n")
        f.write(f"Ego radius: {EGO_RADIUS}, max cycle length: {MAX_CYCLE_LENGTH}\n")
        f.write(f"Time: {elapsed:.1f}s\n\n")
        f.write(f"Total unique negative feedback loops: {len(loops)}\n\n")
        f.write("Length distribution:\n")
        for le in sorted(length_dist):
            f.write(f"  Length {le}: {length_dist[le]}\n")

        # Top pathways
        pw_counts = defaultdict(int)
        for loop in loops:
            for pw in loop["pathway_memberships"]:
                pid = pw.split(":")[0]
                pw_counts[pid] += 1

        f.write("\nTop 20 pathways by loop count:\n")
        for pid, cnt in sorted(pw_counts.items(), key=lambda x: -x[1])[:20]:
            f.write(f"  {pid} ({pw_names.get(pid, '?')}): {cnt}\n")

    print(f"Saved summary → {summary_path}")

    # Also save loop statistics for Paper 2
    stats = {
        "n_total_loops": len(loops),
        "n_self_loops": length_dist.get(1, 0),
        "length_distribution": dict(length_dist),
        "graph_nodes": G.number_of_nodes(),
        "graph_edges": G.number_of_edges(),
    }
    import json
    with open(OUTPUT_DIR / "kegg_loop_stats.json", "w") as f:
        json.dump(stats, f, indent=2)


def main():
    random.seed(42)
    print("=" * 60)
    print("Phase 2.1 v2: KEGG Negative Feedback Loop Extraction")
    print("(Ego-centric approach — tractable on large graphs)")
    print("=" * 60)

    # Build graph
    G, global_sym, pw_names = build_global_graph()

    # Find seed nodes
    seed_ids, found, missing = find_seed_node_ids(G, global_sym)

    # Find loops
    loops, elapsed = find_loops_from_seeds(G, seed_ids, global_sym, pw_names)

    if not loops:
        print("ERROR: No loops found!")
        sys.exit(1)

    # Validate
    validate_known_loops(loops, global_sym)

    # Random sample
    print(f"\n=== 10 Random Loops ===")
    sample = random.sample(loops, min(10, len(loops)))
    for i, loop in enumerate(sample, 1):
        syms = " → ".join(loop["gene_symbols"][:8])
        extra = "..." if len(loop["gene_symbols"]) > 8 else ""
        print(f"  [{i}] len={loop['loop_length']}, inh={loop['n_inhibitory_edges']}: {syms}{extra}")

    # Save
    save_results(loops, G, global_sym, pw_names, elapsed)
    print("\nDone!")


if __name__ == "__main__":
    main()
