#!/usr/bin/env python3
"""
Phase 2.1 final — Validate our 22 feedback loops in KEGG + discover nearby loops.

Two-pronged approach:
1. VALIDATION: For each of 22 loops, find shortest path A→B (activation)
   and B→A (or B⊣A, inhibition) in the KEGG graph.
2. DISCOVERY: For each seed gene, find all 2-node and 3-node negative feedback
   motifs (tractable, specific, interpretable).
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict, deque
import csv
import json
import time

KGML_DIR = Path("/Users/teo/Desktop/research/oscilatory/data/kegg_kgml")
OUTPUT_DIR = Path("/Users/teo/Desktop/research/paper2/data")

ACTIVATION_TYPES = {"activation", "expression", "phosphorylation",
                    "ubiquitination", "glycosylation", "methylation", "indirect effect"}
INHIBITION_TYPES = {"inhibition", "repression", "dephosphorylation", "dissociation"}

# Our 22 loops from Phase 1.1
LOOPS_22 = [
    {"id": 1,  "pathway": "NF-κB / IκBα",    "act": "RELA",    "inh": "NFKBIA"},
    {"id": 2,  "pathway": "p53 / Mdm2",       "act": "TP53",    "inh": "MDM2"},
    {"id": 3,  "pathway": "ERK / DUSP",        "act": "MAPK1",   "inh": "DUSP1"},
    {"id": 4,  "pathway": "Wnt / APC-Axin",    "act": "CTNNB1",  "inh": "AXIN2"},
    {"id": 5,  "pathway": "Notch / FBXW7",     "act": "NOTCH1",  "inh": "FBXW7"},
    {"id": 6,  "pathway": "Circadian",         "act": "ARNTL",   "inh": "PER2"},
    {"id": 7,  "pathway": "Hippo / LATS",      "act": "YAP1",    "inh": "LATS1"},
    {"id": 8,  "pathway": "mTOR / TSC",        "act": "MTOR",    "inh": "TSC1"},
    {"id": 9,  "pathway": "JAK-STAT / SOCS",   "act": "STAT3",   "inh": "SOCS3"},
    {"id": 10, "pathway": "TGF-β / SMAD6-7",   "act": "SMAD3",   "inh": "SMAD7"},
    {"id": 11, "pathway": "NRF2 / KEAP1",      "act": "NFE2L2",  "inh": "KEAP1"},
    {"id": 12, "pathway": "Hedgehog / SUFU",    "act": "GLI1",    "inh": "SUFU"},
    {"id": 13, "pathway": "UPR / IRE1-BiP",    "act": "ERN1",    "inh": "HSPA5"},
    {"id": 14, "pathway": "Calcium / SERCA",   "act": "PLCG1",   "inh": "ATP2A2"},
    {"id": 15, "pathway": "Cell Cycle CKI/CDK","act": "CDK2",    "inh": "CDKN1A"},
    {"id": 16, "pathway": "Rb / E2F",          "act": "E2F1",    "inh": "RB1"},
    {"id": 17, "pathway": "HIF / VHL",         "act": "HIF1A",   "inh": "VHL"},
    {"id": 18, "pathway": "NFAT / RCAN",       "act": "NFATC1",  "inh": "RCAN1"},
    {"id": 19, "pathway": "PI3K-AKT / PTEN",   "act": "PIK3CA",  "inh": "PTEN"},
    {"id": 20, "pathway": "Myc / FBXW7",       "act": "MYC",     "inh": "FBXW7"},
    {"id": 21, "pathway": "Id / bHLH",         "act": "TCF3",    "inh": "ID1"},
    {"id": 22, "pathway": "IκBζ / IL-6",       "act": "IL6",     "inh": "NFKBIZ"},
]


def build_graph():
    """Build global graph with better symbol resolution."""
    kgml_files = sorted(KGML_DIR.glob("*.kgml"))
    print(f"Parsing {len(kgml_files)} KGML files...")

    # edge_info: (src_kegg_id, tgt_kegg_id) -> {sign, pathways}
    edge_info = {}
    global_sym = {}  # kegg_id -> primary symbol
    # Also build alias map: multiple symbols per KEGG id
    all_aliases = defaultdict(set)  # kegg_id -> set of all symbols seen
    pw_names = {}

    for kf in kgml_files:
        try:
            tree = ET.parse(str(kf))
            root = tree.getroot()
            pid = root.get("name", "").replace("path:", "")
            pw_names[pid] = root.get("title", "")

            entry_map = {}
            for entry in root.findall("entry"):
                eid = entry.get("id")
                if entry.get("type") != "gene":
                    if entry.get("type") == "group":
                        comps = []
                        for c in entry.findall("component"):
                            cid = c.get("id")
                            if cid in entry_map:
                                comps.extend(entry_map[cid])
                        if comps:
                            entry_map[eid] = comps
                    continue

                kegg_ids = [x.replace("hsa:", "") for x in entry.get("name", "").split()]
                entry_map[eid] = kegg_ids

                graphics = entry.find("graphics")
                if graphics is not None:
                    display = graphics.get("name", "")
                    if display:
                        parts = [p.strip() for p in display.rstrip(".").split(",")]
                        for i, kid in enumerate(kegg_ids):
                            if i < len(parts):
                                sym = parts[i].strip().split("/")[0].strip()
                                if sym:
                                    if kid not in global_sym:
                                        global_sym[kid] = sym
                                    all_aliases[kid].add(sym.upper())
                                    # Also add all slash variants
                                    for variant in parts[i].strip().split("/"):
                                        v = variant.strip()
                                        if v:
                                            all_aliases[kid].add(v.upper())

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
                        key = (s, t)
                        if key not in edge_info:
                            edge_info[key] = {"sign": sign, "pathways": {pid}}
                        else:
                            edge_info[key]["pathways"].add(pid)
                            if sign == -1:
                                edge_info[key]["sign"] = -1
        except:
            pass

    # Build adjacency
    adj = defaultdict(list)
    for (s, t), info in edge_info.items():
        adj[s].append((t, info["sign"]))

    # Reverse symbol map
    sym_to_ids = defaultdict(set)
    for kid, aliases in all_aliases.items():
        for alias in aliases:
            sym_to_ids[alias].add(kid)
    for kid, sym in global_sym.items():
        sym_to_ids[sym.upper()].add(kid)

    n_nodes = len(set(list(adj.keys()) + [t for ns in adj.values() for t, _ in ns]))
    print(f"Graph: {n_nodes} nodes, {len(edge_info)} edges")
    print(f"Symbol map: {len(global_sym)} gene IDs, {len(sym_to_ids)} unique symbols")

    return adj, edge_info, global_sym, sym_to_ids, pw_names


def bfs_shortest_path(adj, src_ids, tgt_ids, max_depth=6):
    """BFS from any src_id to any tgt_id, return shortest path + signs."""
    queue = deque()
    visited = set()

    for sid in src_ids:
        queue.append((sid, [sid]))
        visited.add(sid)

    while queue:
        node, path = queue.popleft()
        if len(path) > max_depth + 1:
            break

        for neighbor, sign in adj.get(node, []):
            if neighbor in tgt_ids:
                return path + [neighbor]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None


def get_path_info(path, edge_info, global_sym):
    """Get sign and pathway info for a path."""
    signs = []
    pws = set()
    for i in range(len(path) - 1):
        key = (path[i], path[i+1])
        if key in edge_info:
            signs.append(edge_info[key]["sign"])
            pws.update(edge_info[key]["pathways"])
        else:
            signs.append(0)

    symbols = [global_sym.get(n, n) for n in path]
    return symbols, signs, pws


def validate_22_loops(adj, edge_info, global_sym, sym_to_ids, pw_names):
    """For each of 22 loops, check KEGG evidence for the feedback circuit."""
    print("\n" + "=" * 70)
    print("VALIDATION: 22 feedback loops in KEGG")
    print("=" * 70)

    results = []

    for loop in LOOPS_22:
        act_gene = loop["act"]
        inh_gene = loop["inh"]

        act_ids = sym_to_ids.get(act_gene.upper(), set())
        inh_ids = sym_to_ids.get(inh_gene.upper(), set())

        result = {
            "loop_id": loop["id"],
            "pathway": loop["pathway"],
            "activator": act_gene,
            "inhibitor": inh_gene,
            "act_in_kegg": len(act_ids) > 0,
            "inh_in_kegg": len(inh_ids) > 0,
            "act_to_inh_path": None,
            "inh_to_act_path": None,
            "act_to_inh_len": None,
            "inh_to_act_len": None,
            "has_inhibitory_edge": False,
            "full_loop_in_kegg": False,
        }

        if not act_ids or not inh_ids:
            status = "MISSING GENE"
            if not act_ids:
                status += f" ({act_gene})"
            if not inh_ids:
                status += f" ({inh_gene})"
            print(f"  [{loop['id']:2d}] {loop['pathway']:30s}  {status}")
            results.append(result)
            continue

        # Forward: activator → inhibitor (should be activation path)
        fwd_path = bfs_shortest_path(adj, act_ids, inh_ids, max_depth=6)
        # Reverse: inhibitor → activator (should contain inhibition)
        rev_path = bfs_shortest_path(adj, inh_ids, act_ids, max_depth=6)

        # Also check direct edges
        direct_act_to_inh = any((a, b) in edge_info for a in act_ids for b in inh_ids)
        direct_inh_to_act = any((b, a) in edge_info for a in act_ids for b in inh_ids)

        # Check if inhibitory edge exists anywhere in the loop
        has_inh = False
        if rev_path:
            _, signs, _ = get_path_info(rev_path, edge_info, global_sym)
            has_inh = -1 in signs

        # Also check direct inhibition
        if not has_inh:
            for b in inh_ids:
                for a in act_ids:
                    if (b, a) in edge_info and edge_info[(b, a)]["sign"] == -1:
                        has_inh = True
                        break

        result["has_inhibitory_edge"] = has_inh

        if fwd_path:
            syms, signs, pws = get_path_info(fwd_path, edge_info, global_sym)
            result["act_to_inh_path"] = " → ".join(syms)
            result["act_to_inh_len"] = len(fwd_path) - 1
        if rev_path:
            syms, signs, pws = get_path_info(rev_path, edge_info, global_sym)
            result["inh_to_act_path"] = " → ".join(syms)
            result["inh_to_act_len"] = len(rev_path) - 1

        result["full_loop_in_kegg"] = (fwd_path is not None) and (rev_path is not None) and has_inh

        # Status
        if result["full_loop_in_kegg"]:
            status = f"✓ FULL LOOP (fwd={result['act_to_inh_len']}, rev={result['inh_to_act_len']})"
        elif fwd_path and rev_path:
            status = f"~ PARTIAL (paths exist, no inhibition annotated)"
        elif fwd_path or rev_path:
            direction = "fwd" if fwd_path else "rev"
            status = f"~ ONE-WAY ({direction} only)"
        else:
            status = "✗ NO PATHS"

        print(f"  [{loop['id']:2d}] {loop['pathway']:30s}  {status}")
        if fwd_path:
            syms, _, _ = get_path_info(fwd_path, edge_info, global_sym)
            print(f"       fwd: {' → '.join(syms)}")
        if rev_path:
            syms, signs, _ = get_path_info(rev_path, edge_info, global_sym)
            sign_str = " → ".join(f"{'⊣' if s==-1 else '→'} {syms[i+1]}" for i, s in enumerate(signs))
            print(f"       rev: {syms[0]} {sign_str}")

        results.append(result)

    # Summary
    n_full = sum(1 for r in results if r["full_loop_in_kegg"])
    n_partial = sum(1 for r in results if not r["full_loop_in_kegg"] and
                    (r["act_to_inh_path"] or r["inh_to_act_path"]))
    n_missing = sum(1 for r in results if not r["act_in_kegg"] or not r["inh_in_kegg"])
    n_none = 22 - n_full - n_partial - n_missing

    print(f"\n--- Summary ---")
    print(f"  Full loop in KEGG:    {n_full}/22")
    print(f"  Partial evidence:     {n_partial}/22")
    print(f"  Gene missing in KEGG: {n_missing}/22")
    print(f"  No KEGG evidence:     {n_none}/22")

    return results


def discover_short_nfl(adj, edge_info, global_sym, sym_to_ids, pw_names):
    """Find ALL 2-node and 3-node negative feedback loops involving seed genes."""
    print("\n" + "=" * 70)
    print("DISCOVERY: Short (2-3 node) negative feedback motifs")
    print("=" * 70)

    all_seeds = set()
    for loop in LOOPS_22:
        all_seeds.update(sym_to_ids.get(loop["act"].upper(), set()))
        all_seeds.update(sym_to_ids.get(loop["inh"].upper(), set()))

    # 2-node loops: A→B and B⊣A (or A⊣B and B→A)
    two_node = []
    checked = set()
    for a in all_seeds:
        for b, sign_ab in adj.get(a, []):
            if (a, b) in checked or a == b:
                continue
            checked.add((a, b))
            # Check reverse
            if (b, a) in edge_info:
                sign_ba = edge_info[(b, a)]["sign"]
                # Negative feedback = exactly one inhibition
                if (sign_ab == -1) != (sign_ba == -1):  # XOR
                    sym_a = global_sym.get(a, a)
                    sym_b = global_sym.get(b, b)
                    pws = edge_info[(a, b)]["pathways"] | edge_info[(b, a)]["pathways"]
                    two_node.append({
                        "genes": [sym_a, sym_b],
                        "length": 2,
                        "n_inh": 1,
                        "pathways": [f"{p}:{pw_names.get(p,'?')}" for p in pws],
                    })

    print(f"\n2-node NFL motifs: {len(two_node)}")
    for m in two_node[:20]:
        print(f"  {m['genes'][0]} ↔ {m['genes'][1]}  ({len(m['pathways'])} pathways)")

    # 3-node loops: A→B→C→A with odd inhibition
    three_node = []
    seen = set()
    for a in all_seeds:
        for b, s_ab in adj.get(a, []):
            if b == a:
                continue
            for c, s_bc in adj.get(b, []):
                if c == a or c == b:
                    continue
                if (c, a) in edge_info:
                    s_ca = edge_info[(c, a)]["sign"]
                    n_inh = sum(1 for s in [s_ab, s_bc, s_ca] if s == -1)
                    if n_inh % 2 == 1:
                        key = frozenset([a, b, c])
                        if key not in seen:
                            seen.add(key)
                            sym_a = global_sym.get(a, a)
                            sym_b = global_sym.get(b, b)
                            sym_c = global_sym.get(c, c)
                            pws = set()
                            for e in [(a,b), (b,c), (c,a)]:
                                if e in edge_info:
                                    pws.update(edge_info[e]["pathways"])
                            three_node.append({
                                "genes": [sym_a, sym_b, sym_c],
                                "length": 3,
                                "n_inh": n_inh,
                                "pathways": [f"{p}:{pw_names.get(p,'?')}" for p in pws],
                            })

    print(f"\n3-node NFL motifs: {len(three_node)}")
    for m in three_node[:30]:
        print(f"  {' → '.join(m['genes'])} (inh={m['n_inh']}, {len(m['pathways'])} pw)")

    return two_node, three_node


def save_all(validation_results, two_node, three_node):
    """Save everything."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Validation CSV
    with open(OUTPUT_DIR / "kegg_loop_validation.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["loop_id", "pathway", "activator", "inhibitor",
                     "act_in_kegg", "inh_in_kegg",
                     "fwd_path", "fwd_len", "rev_path", "rev_len",
                     "has_inhibition", "full_loop"])
        for r in validation_results:
            w.writerow([r["loop_id"], r["pathway"], r["activator"], r["inhibitor"],
                        r["act_in_kegg"], r["inh_in_kegg"],
                        r["act_to_inh_path"] or "", r["act_to_inh_len"] or "",
                        r["inh_to_act_path"] or "", r["inh_to_act_len"] or "",
                        r["has_inhibitory_edge"], r["full_loop_in_kegg"]])

    # Short motifs CSV
    all_motifs = two_node + three_node
    with open(OUTPUT_DIR / "kegg_short_nfl_motifs.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["motif_id", "genes", "length", "n_inh", "n_pathways"])
        for i, m in enumerate(all_motifs, 1):
            w.writerow([f"M_{i:04d}", "|".join(m["genes"]), m["length"],
                        m["n_inh"], len(m["pathways"])])

    print(f"\nSaved validation → {OUTPUT_DIR / 'kegg_loop_validation.csv'}")
    print(f"Saved {len(all_motifs)} short motifs → {OUTPUT_DIR / 'kegg_short_nfl_motifs.csv'}")


def main():
    print("=" * 70)
    print("Phase 2.1: KEGG Loop Validation + Short Motif Discovery")
    print("=" * 70)

    adj, edge_info, global_sym, sym_to_ids, pw_names = build_graph()

    # Check which seed genes are found
    print("\nSeed gene resolution:")
    for loop in LOOPS_22:
        for gene in [loop["act"], loop["inh"]]:
            ids = sym_to_ids.get(gene.upper(), set())
            if ids:
                syms = [global_sym.get(i, i) for i in ids]
                print(f"  {gene:10s} → {len(ids)} KEGG IDs ({', '.join(syms[:3])})")
            else:
                print(f"  {gene:10s} → NOT FOUND")

    validation = validate_22_loops(adj, edge_info, global_sym, sym_to_ids, pw_names)
    two, three = discover_short_nfl(adj, edge_info, global_sym, sym_to_ids, pw_names)
    save_all(validation, two, three)
    print("\nDone!")


if __name__ == "__main__":
    main()
