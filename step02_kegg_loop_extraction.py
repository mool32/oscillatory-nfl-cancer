#!/usr/bin/env python3
"""
Phase 2.1 — Extract ALL negative feedback loops from KEGG KGML files.

Parses all KGML files into a global signaling DiGraph, finds simple cycles,
and filters for negative feedback loops (odd number of inhibitory edges).

Output:
  - data/kegg_negative_feedback_loops.csv
  - data/kegg_loop_summary.txt
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
# Configuration
# ---------------------------------------------------------------------------
KGML_DIR = Path("/Users/teo/Desktop/research/oscilatory/data/kegg_kgml")
OUTPUT_DIR = Path("/Users/teo/Desktop/research/paper2/data")
MAX_CYCLE_LENGTH = 10  # Will reduce if too slow

# Edge sign classification (from existing parser)
ACTIVATION_TYPES = {
    "activation", "expression", "phosphorylation",
    "ubiquitination", "glycosylation", "methylation",
    "indirect effect",
}
INHIBITION_TYPES = {
    "inhibition", "repression", "dephosphorylation",
    "dissociation",
}


# ---------------------------------------------------------------------------
# KGML Parsing (adapted from kegg_extract.py, enhanced for symbol mapping)
# ---------------------------------------------------------------------------

def parse_kgml(kgml_path: Path) -> tuple:
    """
    Parse a single KGML file into edges and node metadata.

    Returns:
        pathway_id: str (e.g. 'hsa04064')
        pathway_title: str
        edges: list of (source_id, target_id, sign, subtypes)
        id_to_symbols: dict mapping KEGG gene ID -> list of HGNC symbols
        id_to_display: dict mapping KEGG gene ID -> display name
    """
    tree = ET.parse(str(kgml_path))
    root = tree.getroot()

    pathway_id = root.get("name", "").replace("path:", "")
    pathway_title = root.get("title", "")

    # Build entry_id -> gene info mapping
    entry_map = {}  # entry XML id -> {type, kegg_ids, display_name}
    id_to_symbols = {}  # kegg gene id (str) -> first symbol from display name
    id_to_display = {}

    for entry in root.findall("entry"):
        entry_id = entry.get("id")
        entry_type = entry.get("type")
        entry_name = entry.get("name", "")

        if entry_type == "gene":
            # KEGG IDs like "hsa:598" -> "598"
            kegg_ids = [x.replace("hsa:", "") for x in entry_name.split()]

            graphics = entry.find("graphics")
            display_name = graphics.get("name", "") if graphics is not None else ""

            # Extract HGNC symbols from display name
            # Format: "BCL2L1, BCL-XL/S, BCL2L, ..." or "IRAK1, pelle..."
            # For multi-gene entries: "TRAF2, MGC:45012..." but the entry has multiple hsa: ids
            # The first token before comma is typically the primary symbol
            symbols = []
            if display_name:
                # Split by "..." first (multi-gene entries truncate)
                clean = display_name.rstrip(".")
                parts = [p.strip() for p in clean.split(",")]
                # For multi-gene entries, map each kegg_id to its primary symbol
                # KEGG often lists them in order matching the hsa: ids
                for p in parts:
                    # Take first word (the gene symbol)
                    sym = p.strip().split("/")[0].strip()
                    if sym and not sym.startswith("MGC:") and len(sym) < 20:
                        symbols.append(sym)

            entry_map[entry_id] = {
                "type": "gene",
                "kegg_ids": kegg_ids,
                "symbols": symbols,
                "display_name": display_name,
            }

            # Map kegg IDs to symbols
            # If we have matching counts, map 1:1; otherwise use first symbol for all
            if len(symbols) >= len(kegg_ids):
                for i, kid in enumerate(kegg_ids):
                    id_to_symbols[kid] = symbols[i] if i < len(symbols) else symbols[0]
                    id_to_display[kid] = display_name
            else:
                primary = symbols[0] if symbols else kegg_ids[0]
                for kid in kegg_ids:
                    id_to_symbols[kid] = primary
                    id_to_display[kid] = display_name

    # Handle group entries
    for entry in root.findall("entry"):
        entry_id = entry.get("id")
        entry_type = entry.get("type")
        if entry_type == "group":
            components = []
            for comp in entry.findall("component"):
                comp_id = comp.get("id")
                if comp_id in entry_map:
                    components.extend(entry_map[comp_id]["kegg_ids"])
            if components:
                entry_map[entry_id] = {
                    "type": "group",
                    "kegg_ids": components,
                    "symbols": [id_to_symbols.get(k, k) for k in components],
                    "display_name": "",
                }

    # Parse relations -> edges
    edges = []
    for relation in root.findall("relation"):
        entry1_id = relation.get("entry1")
        entry2_id = relation.get("entry2")

        if entry1_id not in entry_map or entry2_id not in entry_map:
            continue

        # Determine sign from subtypes
        subtypes = []
        sign = 0
        for subtype in relation.findall("subtype"):
            st_name = subtype.get("name")
            subtypes.append(st_name)
            if st_name in ACTIVATION_TYPES:
                sign = max(sign, 1)
            elif st_name in INHIBITION_TYPES:
                sign = -1  # inhibition overrides

        if sign == 0:
            sign = 1  # default to activation

        sources = entry_map[entry1_id]["kegg_ids"]
        targets = entry_map[entry2_id]["kegg_ids"]

        for s in sources:
            for t in targets:
                edges.append((s, t, sign, subtypes))

    return pathway_id, pathway_title, edges, id_to_symbols, id_to_display


def build_global_graph(kgml_dir: Path):
    """
    Parse all KGML files and merge into a single global DiGraph.

    Nodes: KEGG gene IDs (will be labeled with HGNC symbols)
    Edges carry: sign (+1/-1), source_pathways (set), subtypes
    """
    kgml_files = sorted(kgml_dir.glob("*.kgml"))
    print(f"Found {len(kgml_files)} KGML files")

    G = nx.DiGraph()
    global_id_to_symbol = {}
    pathway_names = {}

    n_parsed = 0
    n_failed = 0

    for kf in kgml_files:
        try:
            pid, ptitle, edges, id2sym, id2disp = parse_kgml(kf)
            pathway_names[pid] = ptitle
            global_id_to_symbol.update(id2sym)

            for src, tgt, sign, subtypes in edges:
                if src == tgt:
                    # Self-loops: still add them, mark specially
                    pass

                if G.has_edge(src, tgt):
                    # Edge exists: update pathway membership
                    G[src][tgt]["source_pathways"].add(pid)
                    # If conflicting signs from different pathways,
                    # keep the one with inhibition (conservative for finding neg loops)
                    if sign == -1:
                        G[src][tgt]["sign"] = -1
                else:
                    G.add_edge(src, tgt,
                               sign=sign,
                               source_pathways={pid},
                               subtypes=subtypes)

                # Ensure nodes exist with symbol info
                for nid in [src, tgt]:
                    if nid in id2sym:
                        G.nodes[nid]["symbol"] = id2sym[nid]
                    if nid in id2disp:
                        G.nodes[nid]["display_name"] = id2disp[nid]

            n_parsed += 1
        except Exception as e:
            print(f"  WARNING: Failed to parse {kf.name}: {e}")
            n_failed += 1

    print(f"Parsed {n_parsed} pathways ({n_failed} failed)")
    print(f"Global graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Count edge types
    n_act = sum(1 for _, _, d in G.edges(data=True) if d.get("sign", 1) == 1)
    n_inh = sum(1 for _, _, d in G.edges(data=True) if d.get("sign", 1) == -1)
    n_self = sum(1 for u, v in G.edges() if u == v)
    print(f"  Activation edges: {n_act}")
    print(f"  Inhibition edges: {n_inh}")
    print(f"  Self-loops: {n_self}")

    # Connected components
    n_components = nx.number_weakly_connected_components(G)
    largest_cc = max(nx.weakly_connected_components(G), key=len)
    print(f"  Weakly connected components: {n_components}")
    print(f"  Largest component: {len(largest_cc)} nodes")

    return G, global_id_to_symbol, pathway_names


# ---------------------------------------------------------------------------
# Cycle Finding
# ---------------------------------------------------------------------------

def find_negative_feedback_loops(G: nx.DiGraph, id_to_symbol: dict,
                                  pathway_names: dict, max_length: int = 10):
    """
    Find all simple cycles up to max_length, filter for negative feedback.

    A negative feedback loop has an ODD number of inhibitory edges.
    """
    print(f"\nFinding simple cycles with length <= {max_length}...")
    print(f"Graph size: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    t0 = time.time()

    # Use length_bound parameter (NetworkX >= 3.1)
    try:
        all_cycles = list(nx.simple_cycles(G, length_bound=max_length))
    except TypeError:
        # Older NetworkX without length_bound
        print("  WARNING: NetworkX version doesn't support length_bound.")
        print("  Using Johnson's algorithm and filtering by length post-hoc.")
        print("  This may be slow...")
        all_cycles = [c for c in nx.simple_cycles(G) if len(c) <= max_length]

    elapsed = time.time() - t0
    print(f"  Found {len(all_cycles)} total simple cycles in {elapsed:.1f}s")

    # Filter for negative feedback loops (odd number of inhibitory edges)
    neg_loops = []

    for cycle in all_cycles:
        n_inh = 0
        n_act = 0
        edge_pathways = set()
        is_self = len(cycle) == 1

        for i in range(len(cycle)):
            src = cycle[i]
            tgt = cycle[(i + 1) % len(cycle)]

            if G.has_edge(src, tgt):
                edata = G[src][tgt]
                sign = edata.get("sign", 1)
                if sign == -1:
                    n_inh += 1
                else:
                    n_act += 1
                edge_pathways.update(edata.get("source_pathways", set()))
            else:
                # Edge missing (shouldn't happen in a cycle)
                n_act += 1

        # Negative feedback = odd number of inhibitory edges
        if n_inh % 2 == 1:
            # Get gene symbols
            symbols = [id_to_symbol.get(nid, nid) for nid in cycle]

            # Deduplication key: sorted tuple of members
            dedup_key = tuple(sorted(cycle))

            # Pathway memberships with names
            pw_list = []
            for pw in edge_pathways:
                name = pathway_names.get(pw, pw)
                pw_list.append(f"{pw}:{name}")

            neg_loops.append({
                "cycle_nodes": cycle,
                "gene_symbols": symbols,
                "loop_length": len(cycle),
                "n_inhibitory_edges": n_inh,
                "n_activating_edges": n_act,
                "pathway_memberships": pw_list,
                "is_self_loop": is_self,
                "dedup_key": dedup_key,
            })

    print(f"  Negative feedback loops (before dedup): {len(neg_loops)}")

    # Deduplicate: same cycle found from different starting nodes
    seen_keys = set()
    unique_loops = []
    for loop in neg_loops:
        key = loop["dedup_key"]
        if key not in seen_keys:
            seen_keys.add(key)
            unique_loops.append(loop)

    print(f"  Negative feedback loops (after dedup): {len(unique_loops)}")

    return unique_loops


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_known_loops(loops: list, id_to_symbol: dict):
    """Check that known negative feedback loops are found."""
    print("\n=== Validation: Known Loops ===")

    # Build symbol -> kegg_id reverse map
    sym_to_id = {}
    for kid, sym in id_to_symbol.items():
        sym_to_id[sym.upper()] = kid

    # Known negative feedback loops to check
    known_pairs = [
        ("NF-kB / NFKBIA", ["RELA", "NFKBIA"]),
        ("p53 / MDM2", ["TP53", "MDM2"]),
        ("ERK / DUSP", ["MAPK1", "DUSP"]),  # DUSP family - partial match
    ]

    all_loop_symbols_upper = []
    for loop in loops:
        syms_upper = set(s.upper() for s in loop["gene_symbols"])
        all_loop_symbols_upper.append(syms_upper)

    for name, targets in known_pairs:
        found = False
        matching_loops = []
        for i, syms in enumerate(all_loop_symbols_upper):
            # Check if all target genes (or prefix matches for families) are in the loop
            all_match = True
            for t in targets:
                if not any(s.startswith(t) for s in syms):
                    all_match = False
                    break
            if all_match:
                found = True
                matching_loops.append(loops[i])

        status = "FOUND" if found else "NOT FOUND"
        print(f"  {name}: {status}", end="")
        if matching_loops:
            # Show first match
            m = matching_loops[0]
            print(f" (example: {' -> '.join(m['gene_symbols'][:6])}{'...' if len(m['gene_symbols']) > 6 else ''})")
        else:
            print()
            # Debug: check if the genes exist at all
            for t in targets:
                matches = [s for s in sym_to_id if s.startswith(t)]
                if matches:
                    print(f"    Gene {t} exists in graph as: {matches[:5]}")
                else:
                    print(f"    Gene {t} NOT found in graph symbol map")


def print_random_loops(loops: list, n: int = 10):
    """Print n random loops for manual inspection."""
    print(f"\n=== {n} Random Loops for Inspection ===")
    sample = random.sample(loops, min(n, len(loops)))
    for i, loop in enumerate(sample, 1):
        syms = " -> ".join(loop["gene_symbols"])
        pw = "; ".join(loop["pathway_memberships"][:3])
        if len(loop["pathway_memberships"]) > 3:
            pw += f" (+{len(loop['pathway_memberships']) - 3} more)"
        print(f"  [{i}] Length={loop['loop_length']}, "
              f"inh={loop['n_inhibitory_edges']}, act={loop['n_activating_edges']}")
        print(f"      {syms}")
        print(f"      Pathways: {pw}")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_results(loops: list, output_dir: Path, G: nx.DiGraph,
                 id_to_symbol: dict, pathway_names: dict, elapsed: float):
    """Save CSV and summary."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # CSV
    csv_path = output_dir / "kegg_negative_feedback_loops.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "loop_id", "gene_members", "loop_length",
            "n_inhibitory_edges", "n_activating_edges",
            "pathway_memberships", "is_self_loop"
        ])
        for i, loop in enumerate(loops, 1):
            writer.writerow([
                f"NFL_{i:04d}",
                "|".join(loop["gene_symbols"]),
                loop["loop_length"],
                loop["n_inhibitory_edges"],
                loop["n_activating_edges"],
                "|".join(loop["pathway_memberships"]),
                loop["is_self_loop"],
            ])
    print(f"\nSaved {len(loops)} loops to {csv_path}")

    # Summary stats
    summary_path = output_dir / "kegg_loop_summary.txt"

    lengths = [l["loop_length"] for l in loops]
    n_self = sum(1 for l in loops if l["is_self_loop"])
    length_dist = defaultdict(int)
    for le in lengths:
        length_dist[le] += 1

    # Pathway participation
    pw_counts = defaultdict(int)
    for loop in loops:
        for pw in loop["pathway_memberships"]:
            pw_id = pw.split(":")[0]
            pw_counts[pw_id] += 1

    with open(summary_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("KEGG Negative Feedback Loop Extraction Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Global graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges\n")
        n_act = sum(1 for _, _, d in G.edges(data=True) if d.get("sign", 1) == 1)
        n_inh = sum(1 for _, _, d in G.edges(data=True) if d.get("sign", 1) == -1)
        f.write(f"  Activation edges: {n_act}\n")
        f.write(f"  Inhibition edges: {n_inh}\n")
        n_comp = nx.number_weakly_connected_components(G)
        f.write(f"  Weakly connected components: {n_comp}\n\n")

        f.write(f"Total negative feedback loops: {len(loops)}\n")
        f.write(f"  Self-loops (gene inhibits itself): {n_self}\n")
        f.write(f"  Multi-gene loops: {len(loops) - n_self}\n\n")

        f.write("Loop length distribution:\n")
        for le in sorted(length_dist.keys()):
            f.write(f"  Length {le}: {length_dist[le]} loops\n")

        f.write(f"\nMax cycle length searched: {MAX_CYCLE_LENGTH}\n")
        f.write(f"Cycle finding time: {elapsed:.1f}s\n\n")

        f.write("Top 20 pathways by loop participation:\n")
        for pw_id, count in sorted(pw_counts.items(), key=lambda x: -x[1])[:20]:
            pw_name = pathway_names.get(pw_id, pw_id)
            f.write(f"  {pw_id} ({pw_name}): {count} loops\n")

        # Count assessment
        f.write(f"\n--- Quality Assessment ---\n")
        if len(loops) < 20:
            f.write("WARNING: <20 loops found. Criteria may be too strict.\n")
        elif len(loops) > 2000:
            f.write("WARNING: >2000 loops found. Criteria may be too loose.\n")
        else:
            f.write(f"Loop count ({len(loops)}) is in expected range [20-2000].\n")

    print(f"Saved summary to {summary_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    random.seed(42)
    print("=" * 60)
    print("Phase 2.1: KEGG Negative Feedback Loop Extraction")
    print("=" * 60)

    # Step 1: Build global graph
    print("\n--- Step 1: Parse KGML files and build global graph ---")
    G, id_to_symbol, pathway_names = build_global_graph(KGML_DIR)

    # Step 2: Find cycles
    print("\n--- Step 2: Find negative feedback loops ---")
    max_len = MAX_CYCLE_LENGTH

    # If graph is very large, start with shorter cycles
    if G.number_of_nodes() > 5000:
        print(f"WARNING: Large graph ({G.number_of_nodes()} nodes). Starting with max_length=6")
        max_len = 6

    t0 = time.time()
    loops = find_negative_feedback_loops(G, id_to_symbol, pathway_names, max_length=max_len)
    elapsed = time.time() - t0

    # If too slow and got zero results, try smaller
    if elapsed > 300 and len(loops) == 0:
        print(f"\nWARNING: Cycle finding took {elapsed:.0f}s. Reducing max_length to 5.")
        max_len = 5
        t0 = time.time()
        loops = find_negative_feedback_loops(G, id_to_symbol, pathway_names, max_length=max_len)
        elapsed = time.time() - t0

    if len(loops) == 0:
        print("\nERROR: No negative feedback loops found! Check graph connectivity.")
        # Debug info
        print(f"Nodes with self-loops: {sum(1 for u, v in G.edges() if u == v)}")
        print(f"Nodes with in+out degree: {sum(1 for n in G.nodes() if G.in_degree(n) > 0 and G.out_degree(n) > 0)}")
        sys.exit(1)

    # Step 3: Sort by loop length
    loops.sort(key=lambda x: (x["loop_length"], x["gene_symbols"]))

    # Step 4: Validate
    validate_known_loops(loops, id_to_symbol)

    # Step 5: Print random sample
    print_random_loops(loops, n=10)

    # Step 6: Save
    print("\n--- Step 3: Save results ---")
    save_results(loops, OUTPUT_DIR, G, id_to_symbol, pathway_names, elapsed)

    print("\nDone!")


if __name__ == "__main__":
    main()
