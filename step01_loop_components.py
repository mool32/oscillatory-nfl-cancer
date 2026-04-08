#!/usr/bin/env python3
"""
Phase 1.1 — Core (activator, inhibitor) pairs for oscillatory feedback loops.

Combines 14 oscillatory pathways from cancer analysis + 19 negative feedback
circuits from entropy analysis into a unified table of ~22 unique loops.
Fetches Ensembl gene IDs via REST API.

Output: data/loop_components.csv
"""

import csv
import json
import os
import time
import urllib.request
import urllib.error

# ── Define the 22 unique feedback loops ──────────────────────────────────────

LOOPS = [
    # From cancer report (14) — all also have entropy-paper counterparts except
    # Circadian, UPR, Calcium which are cancer-only
    {
        "pathway_name": "NF-κB / IκBα",
        "activator_gene": "RELA",
        "inhibitor_gene": "NFKBIA",
        "feedback_type": "SEQ",
        "mechanism": "sequestration of NF-κB in cytoplasm",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "15516999",  # Hoffmann et al. 2002 / Hayden & Ghosh 2008
    },
    {
        "pathway_name": "p53 / Mdm2",
        "activator_gene": "TP53",
        "inhibitor_gene": "MDM2",
        "feedback_type": "ENZ",
        "mechanism": "E3 ubiquitin ligase targeting p53 for degradation",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "15838523",  # Harris & Levine 2005
    },
    {
        "pathway_name": "ERK / DUSP",
        "activator_gene": "MAPK1",
        "inhibitor_gene": "DUSP1",
        "feedback_type": "ENZ",
        "mechanism": "dual-specificity phosphatase dephosphorylating ERK",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "17496910",  # Caunt & Keyse 2013 / Santos et al. 2007
    },
    {
        "pathway_name": "Wnt / APC-Axin",
        "activator_gene": "CTNNB1",
        "inhibitor_gene": "AXIN2",
        "feedback_type": "SEQ",
        "mechanism": "scaffold for β-catenin destruction complex",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "12717450",  # Lustig et al. 2002
    },
    {
        "pathway_name": "Notch / FBXW7",
        "activator_gene": "NOTCH1",
        "inhibitor_gene": "FBXW7",
        "feedback_type": "ENZ",
        "mechanism": "F-box E3 ligase targeting NICD for proteasomal degradation",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "18023131",  # O'Neil et al. 2007
    },
    {
        "pathway_name": "Circadian (CLOCK-BMAL1 / PER-CRY)",
        "activator_gene": "ARNTL",
        "inhibitor_gene": "PER2",
        "feedback_type": "SEQ",
        "mechanism": "transcriptional repression of CLOCK-BMAL1 by PER/CRY complex",
        "in_cancer_report": "yes",
        "in_entropy_paper": "no",
        "literature_PMID": "12397363",  # Reppert & Weaver 2002
    },
    {
        "pathway_name": "Hippo / LATS",
        "activator_gene": "YAP1",
        "inhibitor_gene": "LATS1",
        "feedback_type": "ENZ",
        "mechanism": "kinase phosphorylating YAP for cytoplasmic retention/degradation",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "25049502",  # Yu et al. 2015 / Moroishi et al. 2015
    },
    {
        "pathway_name": "mTOR / TSC",
        "activator_gene": "MTOR",
        "inhibitor_gene": "TSC1",
        "feedback_type": "ENZ",
        "mechanism": "GAP complex inhibiting Rheb GTPase upstream of mTORC1",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "15314020",  # Inoki et al. 2003
    },
    {
        "pathway_name": "JAK-STAT / SOCS",
        "activator_gene": "STAT3",
        "inhibitor_gene": "SOCS3",
        "feedback_type": "ENZ",
        "mechanism": "pseudosubstrate inhibition of JAK + E3 ligase adapter",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "11909529",  # Yoshimura et al. 2007
    },
    {
        "pathway_name": "TGF-β / SMAD6-7",
        "activator_gene": "SMAD3",
        "inhibitor_gene": "SMAD7",
        "feedback_type": "SEQ",
        "mechanism": "inhibitory SMAD blocks receptor phosphorylation of R-SMADs",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "9335507",  # Nakao et al. 1997
    },
    {
        "pathway_name": "NRF2 / KEAP1",
        "activator_gene": "NFE2L2",
        "inhibitor_gene": "KEAP1",
        "feedback_type": "SEQ",
        "mechanism": "cytoplasmic sequestration and CUL3-mediated ubiquitination",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "15001544",  # Itoh et al. 2004
    },
    {
        "pathway_name": "Hedgehog / SUFU",
        "activator_gene": "GLI1",
        "inhibitor_gene": "SUFU",
        "feedback_type": "ENZ",
        "mechanism": "cytoplasmic sequestration/degradation of GLI transcription factors",
        "in_cancer_report": "yes",
        "in_entropy_paper": "yes",
        "literature_PMID": "17360644",  # Svard et al. 2006
    },
    {
        "pathway_name": "UPR / IRE1-BiP",
        "activator_gene": "ERN1",
        "inhibitor_gene": "HSPA5",
        "feedback_type": "SEQ",
        "mechanism": "BiP chaperone binds IRE1 lumenal domain, repressing activation",
        "in_cancer_report": "yes",
        "in_entropy_paper": "no",
        "literature_PMID": "15722559",  # Bertolotti et al. 2000
    },
    {
        "pathway_name": "Calcium / SERCA",
        "activator_gene": "PLCG1",
        "inhibitor_gene": "ATP2A2",
        "feedback_type": "ENZ",
        "mechanism": "SERCA pump restores ER Ca²⁺, terminating cytoplasmic signal",
        "in_cancer_report": "yes",
        "in_entropy_paper": "no",
        "literature_PMID": "12559096",  # Bhatt et al. 2000 / Berridge et al. 2003
    },
    # ── Entropy-paper only loops (8 additional) ─────────────────────────────
    {
        "pathway_name": "Cell Cycle CKI / CDK",
        "activator_gene": "CDK2",
        "inhibitor_gene": "CDKN1A",
        "feedback_type": "SEQ",
        "mechanism": "stoichiometric CDK inhibitor (p21) binding and inactivation",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "8259215",  # Harper et al. 1993
    },
    {
        "pathway_name": "Rb / E2F",
        "activator_gene": "E2F1",
        "inhibitor_gene": "RB1",
        "feedback_type": "SEQ",
        "mechanism": "Rb sequesters E2F, blocking transcriptional activation",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "8259210",  # Weinberg 1995 / Chellappan et al. 1991
    },
    {
        "pathway_name": "HIF / VHL",
        "activator_gene": "HIF1A",
        "inhibitor_gene": "VHL",
        "feedback_type": "SEQ",
        "mechanism": "VHL E3 ligase targets hydroxylated HIF-1α for degradation",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "10205047",  # Maxwell et al. 1999
    },
    {
        "pathway_name": "NFAT / RCAN",
        "activator_gene": "NFATC1",
        "inhibitor_gene": "RCAN1",
        "feedback_type": "SEQ",
        "mechanism": "RCAN1 inhibits calcineurin phosphatase, blocking NFAT dephosphorylation",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "11331884",  # Rothermel et al. 2001
    },
    {
        "pathway_name": "PI3K-AKT / PTEN",
        "activator_gene": "PIK3CA",
        "inhibitor_gene": "PTEN",
        "feedback_type": "ENZ",
        "mechanism": "lipid phosphatase dephosphorylating PIP3 to PIP2",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "9468138",  # Maehama & Dixon 1998
    },
    {
        "pathway_name": "Myc / FBXW7",
        "activator_gene": "MYC",
        "inhibitor_gene": "FBXW7",
        "feedback_type": "ENZ",
        "mechanism": "SCF-FBXW7 E3 ligase ubiquitinates phospho-Myc for proteasomal degradation",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "15580293",  # Welcker et al. 2004
    },
    {
        "pathway_name": "Id / bHLH",
        "activator_gene": "TCF3",
        "inhibitor_gene": "ID1",
        "feedback_type": "SEQ",
        "mechanism": "Id1 sequesters E2A/TCF3 bHLH factors via dominant-negative heterodimerization",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "2289870",  # Benezra et al. 1990
    },
    {
        "pathway_name": "IκBζ / IL-6",
        "activator_gene": "IL6",
        "inhibitor_gene": "NFKBIZ",
        "feedback_type": "SEQ",
        "mechanism": "IκBζ modulates NF-κB target gene selectivity, dampening IL-6 loop",
        "in_cancer_report": "no",
        "in_entropy_paper": "yes",
        "literature_PMID": "15226429",  # Yamamoto et al. 2004
    },
]


# ── Ensembl REST API lookup ──────────────────────────────────────────────────

ENSEMBL_API = "https://rest.ensembl.org"

# Hardcoded fallback IDs for the most critical genes (in case API is slow/down)
ENSEMBL_FALLBACK = {
    "RELA":    "ENSG00000173039",
    "NFKBIA":  "ENSG00000100906",
    "TP53":    "ENSG00000141510",
    "MDM2":    "ENSG00000135679",
    "MAPK1":   "ENSG00000100030",
    "DUSP1":   "ENSG00000120129",
    "CTNNB1":  "ENSG00000168036",
    "AXIN2":   "ENSG00000168646",
    "NOTCH1":  "ENSG00000148400",
    "FBXW7":   "ENSG00000112270",
    "ARNTL":   "ENSG00000133794",
    "PER2":    "ENSG00000132326",
    "YAP1":    "ENSG00000137693",
    "LATS1":   "ENSG00000131023",
    "MTOR":    "ENSG00000198793",
    "TSC1":    "ENSG00000165699",
    "STAT3":   "ENSG00000168610",
    "SOCS3":   "ENSG00000120833",
    "SMAD3":   "ENSG00000166949",
    "SMAD7":   "ENSG00000101665",
    "NFE2L2":  "ENSG00000116044",
    "KEAP1":   "ENSG00000079999",
    "GLI1":    "ENSG00000111087",
    "SUFU":    "ENSG00000107882",
    "ERN1":    "ENSG00000178607",
    "HSPA5":   "ENSG00000044574",
    "PLCG1":   "ENSG00000124181",
    "ATP2A2":  "ENSG00000174437",
    "CDK2":    "ENSG00000123572",
    "CDKN1A":  "ENSG00000124762",
    "E2F1":    "ENSG00000101412",
    "RB1":     "ENSG00000139687",
    "HIF1A":   "ENSG00000100644",
    "VHL":     "ENSG00000134086",
    "NFATC1":  "ENSG00000131196",
    "RCAN1":   "ENSG00000159200",
    "PIK3CA":  "ENSG00000121879",
    "PTEN":    "ENSG00000171862",
    "MYC":     "ENSG00000136997",
    "TCF3":    "ENSG00000071564",
    "ID1":     "ENSG00000125968",
    "IL6":     "ENSG00000136244",
    "NFKBIZ":  "ENSG00000144802",
}


def fetch_ensembl_id(symbol: str, retries: int = 2) -> str:
    """Fetch Ensembl gene ID for a HGNC symbol via REST API."""
    url = f"{ENSEMBL_API}/xrefs/symbol/homo_sapiens/{symbol}?content-type=application/json"
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())
            for entry in data:
                if entry.get("type") == "gene" and entry.get("id", "").startswith("ENSG"):
                    return entry["id"]
            # If no gene type found, return first ENSG hit
            for entry in data:
                if entry.get("id", "").startswith("ENSG"):
                    return entry["id"]
            return ENSEMBL_FALLBACK.get(symbol, "NOT_FOUND")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            if attempt < retries:
                time.sleep(1)
                continue
            print(f"  WARNING: API failed for {symbol}: {e}")
            return ENSEMBL_FALLBACK.get(symbol, "API_FAILED")


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "data")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "loop_components.csv")

    # Collect all unique gene symbols
    all_symbols = set()
    for loop in LOOPS:
        all_symbols.add(loop["activator_gene"])
        all_symbols.add(loop["inhibitor_gene"])

    print(f"Fetching Ensembl IDs for {len(all_symbols)} unique genes...")
    ensembl_map = {}
    for i, sym in enumerate(sorted(all_symbols), 1):
        print(f"  [{i}/{len(all_symbols)}] {sym}...", end=" ", flush=True)
        eid = fetch_ensembl_id(sym)
        ensembl_map[sym] = eid
        print(eid)
        time.sleep(0.15)  # polite rate limiting

    # Write CSV
    fieldnames = [
        "loop_id", "pathway_name",
        "activator_gene", "activator_ensembl",
        "inhibitor_gene", "inhibitor_ensembl",
        "feedback_type", "mechanism",
        "in_cancer_report", "in_entropy_paper",
        "literature_PMID",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for idx, loop in enumerate(LOOPS, 1):
            row = {
                "loop_id": idx,
                "pathway_name": loop["pathway_name"],
                "activator_gene": loop["activator_gene"],
                "activator_ensembl": ensembl_map.get(loop["activator_gene"], ""),
                "inhibitor_gene": loop["inhibitor_gene"],
                "inhibitor_ensembl": ensembl_map.get(loop["inhibitor_gene"], ""),
                "feedback_type": loop["feedback_type"],
                "mechanism": loop["mechanism"],
                "in_cancer_report": loop["in_cancer_report"],
                "in_entropy_paper": loop["in_entropy_paper"],
                "literature_PMID": loop["literature_PMID"],
            }
            writer.writerow(row)

    print(f"\nWrote {len(LOOPS)} loops to {out_path}")

    # Verification: print summary
    api_ok = sum(1 for v in ensembl_map.values() if v.startswith("ENSG"))
    api_fail = sum(1 for v in ensembl_map.values() if not v.startswith("ENSG"))
    print(f"Ensembl IDs resolved: {api_ok}/{len(ensembl_map)}")
    if api_fail:
        failed = [k for k, v in ensembl_map.items() if not v.startswith("ENSG")]
        print(f"  Failed: {failed}")


if __name__ == "__main__":
    main()
