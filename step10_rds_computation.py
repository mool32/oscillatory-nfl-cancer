#!/usr/bin/env python3
"""
Reversibility Dependency Score (RDS) computation.

For each gene G in NFL motifs:
  1. Find all negative feedback loops containing G
  2. For each loop: check if G is irreplaceable (no alternative feedback path)
  3. Determine position sign: rise (+1) or recovery (−1)
  4. RDS(G) = Σ_l Sign(G,l) × R(G,l)

Then test predictions against CGC, oncogene/TSG classification, drug targets.
"""

import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from scipy import stats

DATA_DIR = Path("/Users/teo/Desktop/research/paper2/data")
KGML_DIR = Path("/Users/teo/Desktop/research/oscilatory/data/kegg_kgml")

# Alias resolution (from step05)
ALIAS_MAP = {
    "IMD31A": "STAT3", "IMD31B": "STAT3", "AISIMD": "SOCS3",
    "CANDF7": "STAT2", "STAT91": "STAT1", "MGF": "STAT5B",
    "EBP-1": "RELA", "AOVD2": "SMAD7", "BSP-1": "SMAD5", "BSP1": "SMAD5",
    "JAK1B": "JAK1", "JTK10": "TYK2",
    "HIRS-1": "IRS1", "CCNE": "CCNE1",
    "P3R3URF-PIK3R3": "PIK3R3", "LOC110117498-PIK3R3": "PIK3CA",
    "CCM4": "PIK3CA", "CLAPO": "PIK3CA",
    "5295": "PIK3CA", "5296": "PIK3R1", "8503": "PIK3R3",
    "CDKN4": "CDKN2A", "ERK": "MAPK1",
    "MAPK14": "MAPK14", "CSBP": "MAPK14", "CSBP1": "MAPK14", "CSBP2": "MAPK14",
    "CISH": "CISH", "9306": "SOCS6", "9655": "SOCS5",
    "TGFBR1": "TGFBR1", "AAT5": "TGFBR1",
    "E2F-1": "E2F1", "DILC": "E2F2", "RBAP1": "E2F3",
    "AXIN": "AXIN1", "GLI": "GLI2",
    "TPTEP2-CSNK1E": "CSNK1E", "LOC400927-CSNK1E": "CSNK1E",
    "M-CSF-R": "CSF1R", "TNFAIP3": "TNFAIP3", "TRAF6": "TRAF6",
    "MPPH": "AKT1", "MPPH2": "AKT2",
    "CVID12": "NFKBIA",
}

ACTIVATION_TYPES = {"activation", "expression", "phosphorylation",
                    "ubiquitination", "glycosylation", "methylation", "indirect effect"}
INHIBITION_TYPES = {"inhibition", "repression", "dephosphorylation", "dissociation"}

# Rise/recovery classification from Paper 1
# +1 = rise arm (GoF → stuck ON → oncogene prediction)
# -1 = recovery arm (LoF → brake removed → TSG prediction)
POSITION_SIGN = {
    # NF-κB
    "RELA": +1, "NFKB1": +1, "RELB": +1,
    "NFKBIA": -1, "TNFAIP3": -1, "TRAF6": +1,
    # p53
    "TP53": -1, "MDM2": +1, "MDM4": +1, "ATM": -1, "CHEK2": -1,
    # ERK/MAPK
    "MAPK1": +1, "MAPK3": +1, "MAPK14": +1,
    "DUSP1": -1,
    "SOCS1": -1, "SOCS3": -1, "SOCS4": -1, "SOCS5": -1, "SOCS6": -1, "SOCS7": -1,
    "CISH": -1,
    # Wnt
    "CTNNB1": +1, "CDH1": -1, "AXIN1": -1, "AXIN2": -1,
    # Notch
    "NOTCH1": +1, "FBXW7": -1,
    # Circadian
    "CLOCK": +1, "ARNTL": +1, "PER2": -1, "CRY1": -1,
    # Hippo
    "YAP1": +1, "LATS1": -1, "LATS2": -1, "STK3": -1,
    # mTOR
    "MTOR": +1, "TSC1": -1, "TSC2": -1, "IRS1": +1, "IRS2": +1,
    "PIK3CA": +1, "PIK3R1": -1, "PIK3R3": +1, "PTEN": -1,
    "AKT1": +1, "AKT2": +1, "AKT3": +1, "PTK2": +1,
    # JAK-STAT
    "JAK1": +1, "JAK2": +1, "TYK2": +1,
    "STAT1": +1, "STAT2": +1, "STAT3": +1, "STAT5A": +1, "STAT5B": +1,
    "LEPR": +1, "CSF1R": +1, "PRLR": +1,
    # TGF-β
    "TGFBR1": +1, "SMAD2": +1, "SMAD3": +1, "SMAD4": +1,
    "SMAD1": +1, "SMAD5": +1,
    "SMAD7": -1,
    # BMP receptors
    "BMPR1A": +1, "BMPR1B": +1, "BMPR2": +1,
    # NRF2
    "NFE2L2": +1, "KEAP1": -1,
    # Hedgehog
    "GLI1": +1, "GLI2": +1, "GLI3": -1, "SUFU": -1,
    # UPR
    "ERN1": +1, "HSPA5": -1,
    # Calcium
    "PLCG1": +1, "ATP2A2": -1,
    # Cell cycle
    "CDK2": +1, "CDKN1A": -1, "CDKN1B": -1, "CDKN2A": -1,
    "CCNE1": +1, "E2F1": +1, "E2F2": +1, "E2F3": +1, "TFDP1": +1,
    "RB1": -1, "RBL1": -1, "RBL2": -1,
    # FOXO
    "FOXO1": -1, "FOXO3": -1, "FOXO4": -1, "FOXO6": -1,
    # p53-apoptosis
    "BCL2": +1, "PMAIP1": -1, "SIVA1": -1, "BID": -1, "BBC3": -1, "TP73": -1,
    # NFAT
    "NFATC1": +1, "RCAN1": -1,
    # IL-6/NF-κB
    "IL6": +1, "NFKBIZ": -1,
    # MYC
    "MYC": +1,
    # HIF
    "HIF1A": +1, "VHL": -1,
    # Hippo cross
    "CSNK1E": -1,
}

# Full CGC
def load_cgc():
    genes = set()
    with open(DATA_DIR / "cosmic_cgc_full.txt") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    return genes

# CGC classification (oncogene/TSG/both) — manual from COSMIC
CGC_ROLE = {
    "AKT1": "oncogene", "AKT2": "oncogene", "AKT3": "oncogene",
    "AXIN1": "TSG", "BCL2": "oncogene", "BMPR1A": "TSG",
    "CCNE1": "oncogene", "CDH1": "TSG", "CDK2": "oncogene",
    "CDKN1A": "TSG", "CDKN1B": "TSG", "CDKN2A": "TSG",
    "CSF1R": "oncogene", "CTNNB1": "oncogene",
    "E2F1": "oncogene",
    "FOXO1": "both", "FOXO3": "TSG", "FOXO4": "oncogene",
    "GLI1": "oncogene", "GLI2": "oncogene", "GLI3": "TSG",
    "JAK1": "oncogene", "JAK2": "oncogene",
    "MAPK1": "oncogene", "MDM2": "oncogene", "MTOR": "oncogene",
    "PIK3CA": "oncogene", "PIK3R1": "TSG",
    "PTEN": "TSG", "RB1": "TSG", "RELA": "oncogene",
    "SMAD2": "TSG", "SMAD3": "TSG", "SMAD4": "TSG",
    "SOCS1": "TSG", "STAT3": "oncogene", "STAT5B": "oncogene",
    "TGFBR1": "TSG", "TNFAIP3": "TSG", "TP53": "TSG",
}


def build_kegg_graph():
    """Parse KGML → adjacency with HGNC symbols and edge signs."""
    kgml_files = sorted(KGML_DIR.glob("*.kgml"))
    print(f"Parsing {len(kgml_files)} KGML files...")

    # kegg_id → HGNC symbol
    sym = {}
    # Edges: (src_hgnc, tgt_hgnc) → sign
    edges = {}

    for kf in kgml_files:
        try:
            tree = ET.parse(str(kf))
            root = tree.getroot()
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
                            sym[kid] = symbols[i]

            for entry in root.findall("entry"):
                if entry.get("type") == "group":
                    comps = []
                    for c in entry.findall("component"):
                        cid = c.get("id")
                        if cid in entry_map:
                            comps.extend(entry_map[cid])
                    if comps:
                        entry_map[entry.get("id")] = comps

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
                        s_hgnc = ALIAS_MAP.get(sym.get(s, s), sym.get(s, s))
                        t_hgnc = ALIAS_MAP.get(sym.get(t, t), sym.get(t, t))
                        key = (s_hgnc, t_hgnc)
                        if key not in edges:
                            edges[key] = sign
                        elif sign == -1:
                            edges[key] = -1  # inhibition overrides
        except:
            pass

    # Build adjacency
    adj = defaultdict(list)
    for (s, t), sign in edges.items():
        adj[s].append((t, sign))

    all_nodes = set(s for s, _ in edges) | set(t for _, t in edges)
    print(f"  {len(all_nodes)} nodes, {len(edges)} edges")
    return adj, edges, all_nodes


def load_nfl_motifs():
    """Load unique NFL motifs."""
    motifs = []
    with open(DATA_DIR / "unique_nfl_motifs.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes = row["genes_hgnc"].split("|")
            motifs.append({
                "genes": genes,
                "gene_set": frozenset(genes),
                "length": int(row["length"]),
                "n_inh": int(row["n_inh"]),
                "n_pw": int(row["n_pw"]),
            })
    return motifs


def check_irreplaceability(gene, loop_genes, adj, edges):
    """
    Check if gene G is irreplaceable in a feedback loop.

    G is irreplaceable if removing it from the network destroys
    the negative feedback path among remaining loop members.

    For a loop [A, B, C] where we remove B:
    Check if there's still a path A → ... → C with the same
    net inhibitory effect, using only edges NOT through B.
    """
    loop_list = list(loop_genes)
    if gene not in loop_list:
        return True  # shouldn't happen

    # For simplicity with 2-3 node loops:
    if len(loop_list) == 2:
        # Two-node loop: A ↔ B. Removing either breaks the loop entirely.
        return True

    if len(loop_list) == 3:
        # Three-node loop: A → B → C ⊣ A (or variants).
        # Removing B: need A → ? → C without B. Check 1-hop alternatives.
        remaining = [g for g in loop_list if g != gene]
        g1, g2 = remaining

        # Check if there's any path from g1 to g2 (or g2 to g1)
        # that doesn't go through gene, with at most 2 hops
        for mid_gene, sign in adj.get(g1, []):
            if mid_gene == gene:
                continue
            for tgt, sign2 in adj.get(mid_gene, []):
                if tgt == g2:
                    # Alternative path exists: g1 → mid → g2
                    # Check if net sign preserves negative feedback
                    return False  # There's an alternative — gene is replaceable

        for mid_gene, sign in adj.get(g2, []):
            if mid_gene == gene:
                continue
            for tgt, sign2 in adj.get(mid_gene, []):
                if tgt == g1:
                    return False

        return True  # No alternative found — gene is irreplaceable

    # For longer loops, conservative: assume irreplaceable
    return True


def compute_rds(motifs, adj, edges):
    """Compute RDS for all genes in NFL motifs."""
    gene_loops = defaultdict(list)  # gene → list of loops

    for m in motifs:
        for g in m["genes"]:
            gene_loops[g].append(m)

    rds = {}
    details = {}

    for gene in sorted(gene_loops.keys()):
        loops = gene_loops[gene]
        total_rds = 0
        n_irreplaceable = 0
        n_total = len(loops)
        sign = POSITION_SIGN.get(gene, 0)

        loop_details = []
        for loop in loops:
            irreplaceable = check_irreplaceability(gene, loop["gene_set"], adj, edges)
            r = 1 if irreplaceable else 0
            n_irreplaceable += r
            contribution = sign * r
            total_rds += contribution
            loop_details.append({
                "genes": loop["genes"],
                "irreplaceable": irreplaceable,
                "contribution": contribution,
            })

        rds[gene] = {
            "RDS": total_rds,
            "abs_RDS": abs(total_rds),
            "n_loops": n_total,
            "n_irreplaceable": n_irreplaceable,
            "sign": sign,
            "sign_label": "rise" if sign > 0 else "recovery" if sign < 0 else "unknown",
            "loops": loop_details,
        }

    return rds


def test_predictions(rds, cgc):
    """Test all 5 predictions."""
    print(f"\n{'='*70}")
    print("PREDICTION TESTS")
    print("=" * 70)

    genes = sorted(rds.keys())
    abs_rds = np.array([rds[g]["abs_RDS"] for g in genes])
    signed_rds = np.array([rds[g]["RDS"] for g in genes])
    is_cgc = np.array([g in cgc for g in genes])
    n_loops = np.array([rds[g]["n_loops"] for g in genes])

    # ===================================================================
    # PREDICTION 1: |RDS| correlates with CGC membership
    # ===================================================================
    print(f"\n--- Prediction 1: |RDS| > 0 → CGC membership ---")
    rds_pos = abs_rds > 0
    rds_zero = abs_rds == 0

    cgc_in_pos = sum(is_cgc[rds_pos])
    total_pos = sum(rds_pos)
    cgc_in_zero = sum(is_cgc[rds_zero])
    total_zero = sum(rds_zero)

    print(f"  |RDS| > 0: {cgc_in_pos}/{total_pos} CGC ({100*cgc_in_pos/total_pos:.1f}%)" if total_pos > 0 else "")
    print(f"  |RDS| = 0: {cgc_in_zero}/{total_zero} CGC ({100*cgc_in_zero/total_zero:.1f}%)" if total_zero > 0 else "")

    if total_pos > 0 and total_zero > 0:
        table1 = [[cgc_in_pos, total_pos - cgc_in_pos],
                   [cgc_in_zero, total_zero - cgc_in_zero]]
        OR1, p1 = stats.fisher_exact(table1)
        print(f"  Fisher exact: OR = {OR1:.2f}, p = {p1:.4f}")

    # Point-biserial correlation |RDS| vs CGC
    rpb, p_pb = stats.pointbiserialr(is_cgc.astype(int), abs_rds)
    print(f"  Point-biserial r(|RDS|, CGC) = {rpb:.3f}, p = {p_pb:.4f}")

    # Spearman
    rho1, p_rho1 = stats.spearmanr(abs_rds, is_cgc.astype(int))
    print(f"  Spearman ρ(|RDS|, CGC) = {rho1:.3f}, p = {p_rho1:.4f}")

    # ===================================================================
    # PREDICTION 2: sign(RDS) predicts oncogene vs TSG
    # ===================================================================
    print(f"\n--- Prediction 2: RDS sign → oncogene/TSG ---")

    correct = 0
    total_classified = 0
    og_rds = []
    tsg_rds = []
    both_rds = []

    for g in genes:
        if g in CGC_ROLE:
            role = CGC_ROLE[g]
            s_rds = rds[g]["RDS"]
            if role == "oncogene":
                og_rds.append(s_rds)
                if s_rds > 0: correct += 1
                total_classified += 1
            elif role == "TSG":
                tsg_rds.append(s_rds)
                if s_rds < 0: correct += 1
                total_classified += 1
            elif role == "both":
                both_rds.append(s_rds)
                total_classified += 1

    print(f"  Classified genes: {total_classified}")
    print(f"  Correct sign predictions: {correct}/{total_classified - len(both_rds)} ({100*correct/(total_classified - len(both_rds)):.1f}%)" if total_classified > len(both_rds) else "")

    if og_rds and tsg_rds:
        print(f"\n  Oncogenes (n={len(og_rds)}): mean RDS = {np.mean(og_rds):+.2f}, median = {np.median(og_rds):+.1f}")
        print(f"  TSGs (n={len(tsg_rds)}): mean RDS = {np.mean(tsg_rds):+.2f}, median = {np.median(tsg_rds):+.1f}")
        if both_rds:
            print(f"  Both (n={len(both_rds)}): mean RDS = {np.mean(both_rds):+.2f}, median = {np.median(both_rds):+.1f}")

        U, p_u = stats.mannwhitneyu(og_rds, tsg_rds, alternative='greater')
        print(f"  Mann-Whitney (OG > TSG): U = {U:.0f}, p = {p_u:.4f}")

    # ===================================================================
    # PREDICTION 3: |RDS| vs loop count (RDS adds value beyond counting)
    # ===================================================================
    print(f"\n--- Prediction 3: |RDS| outperforms loop count ---")

    rpb_rds, p_rds = stats.pointbiserialr(is_cgc.astype(int), abs_rds)
    rpb_loops, p_loops = stats.pointbiserialr(is_cgc.astype(int), n_loops)
    print(f"  r(|RDS|, CGC) = {rpb_rds:.3f}, p = {p_rds:.4f}")
    print(f"  r(n_loops, CGC) = {rpb_loops:.3f}, p = {p_loops:.4f}")
    if abs(rpb_rds) > abs(rpb_loops):
        print(f"  → |RDS| is a better predictor than loop count alone")
    else:
        print(f"  → Loop count is sufficient (irreplaceability doesn't add value)")

    # ===================================================================
    # PREDICTION 5: |RDS| ≈ 0 for "both" genes
    # ===================================================================
    print(f"\n--- Prediction 5: 'Both' genes have |RDS| closer to 0 ---")
    og_abs = [abs(r) for r in og_rds]
    tsg_abs = [abs(r) for r in tsg_rds]
    both_abs = [abs(r) for r in both_rds]

    if og_abs:
        print(f"  Oncogenes: mean |RDS| = {np.mean(og_abs):.2f}")
    if tsg_abs:
        print(f"  TSGs: mean |RDS| = {np.mean(tsg_abs):.2f}")
    if both_abs:
        print(f"  Both: mean |RDS| = {np.mean(both_abs):.2f}")
        if len(both_abs) >= 2:
            # Compare both vs (og+tsg)
            single_role = og_abs + tsg_abs
            U5, p5 = stats.mannwhitneyu(single_role, both_abs, alternative='greater')
            print(f"  Mann-Whitney (single > both): U = {U5:.0f}, p = {p5:.4f}")

    return {
        "rpb_rds": rpb_rds,
        "p_rds": p_rds,
        "rpb_loops": rpb_loops,
        "p_loops": p_loops,
        "rho": rho1,
        "p_rho": p_rho1,
    }


def main():
    print("=" * 70)
    print("Reversibility Dependency Score (RDS) Computation")
    print("=" * 70)

    # Load data
    adj, edges, all_nodes = build_kegg_graph()
    motifs = load_nfl_motifs()
    cgc = load_cgc()

    print(f"\nNFL motifs: {len(motifs)}")
    print(f"Unique genes: {len(set(g for m in motifs for g in m['genes']))}")

    # Compute RDS
    rds = compute_rds(motifs, adj, edges)

    # Print RDS table
    print(f"\n{'='*70}")
    print(f"{'Gene':15s} {'RDS':>5s} {'|RDS|':>6s} {'Loops':>6s} {'Irrepl':>7s} {'Sign':>8s} {'CGC':>4s} {'Role':>10s}")
    print("-" * 70)

    for gene in sorted(rds.keys(), key=lambda g: -abs(rds[g]["RDS"])):
        r = rds[gene]
        cgc_mark = "✓" if gene in cgc else ""
        role = CGC_ROLE.get(gene, "")
        print(f"{gene:15s} {r['RDS']:+5d} {r['abs_RDS']:6d} {r['n_loops']:6d} "
              f"{r['n_irreplaceable']:7d} {r['sign_label']:>8s} {cgc_mark:>4s} {role:>10s}")

    # Test predictions
    test_results = test_predictions(rds, cgc)

    # Save
    save_rds(rds, test_results)
    print("\nDone!")


def save_rds(rds, test_results):
    """Save RDS table and test results."""
    # CSV
    with open(DATA_DIR / "rds_scores.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "RDS", "abs_RDS", "n_loops", "n_irreplaceable",
                     "sign", "sign_label"])
        for gene in sorted(rds.keys(), key=lambda g: -abs(rds[g]["RDS"])):
            r = rds[gene]
            w.writerow([gene, r["RDS"], r["abs_RDS"], r["n_loops"],
                        r["n_irreplaceable"], r["sign"], r["sign_label"]])

    # JSON with full details
    json_rds = {}
    for gene, r in rds.items():
        json_rds[gene] = {
            "RDS": r["RDS"],
            "abs_RDS": r["abs_RDS"],
            "n_loops": r["n_loops"],
            "n_irreplaceable": r["n_irreplaceable"],
            "sign": r["sign"],
            "sign_label": r["sign_label"],
        }

    with open(DATA_DIR / "rds_full.json", "w") as f:
        json.dump({"rds": json_rds, "test_results": {
            k: float(v) if isinstance(v, (float, np.floating)) else v
            for k, v in test_results.items()
        }}, f, indent=2)

    print(f"Saved: rds_scores.csv, rds_full.json")


if __name__ == "__main__":
    main()
