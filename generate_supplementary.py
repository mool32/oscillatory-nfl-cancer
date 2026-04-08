#!/usr/bin/env python3
"""
Generate Supplementary Tables S1-S11 for Paper 2.
"""

import json
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict

DATA = Path("/Users/teo/Desktop/research/paper2/data")
SUPP = Path("/Users/teo/Desktop/research/paper2/supplementary")
SUPP.mkdir(exist_ok=True)


def load_json(name):
    with open(DATA / name) as f:
        return json.load(f)


def load_csv(name):
    rows = []
    with open(DATA / name) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


# ── Table S1: 22 core feedback loops ──────────────────────────────────
def table_s1():
    loops = load_csv('loop_components.csv')
    with open(SUPP / 'table_s1_core_loops.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Loop ID', 'Pathway', 'Activator Gene', 'Activator Ensembl',
                        'Inhibitor Gene', 'Inhibitor Ensembl', 'Feedback Type',
                        'Mechanism', 'In Cancer Report', 'In Entropy Paper', 'Literature PMID'])
        for row in loops:
            writer.writerow([
                row.get('loop_id', ''),
                row.get('pathway_name', ''),
                row.get('activator_gene', ''),
                row.get('activator_ensembl', ''),
                row.get('inhibitor_gene', ''),
                row.get('inhibitor_ensembl', ''),
                row.get('feedback_type', ''),
                row.get('mechanism', ''),
                row.get('in_cancer_report', ''),
                row.get('in_entropy_paper', ''),
                row.get('literature_PMID', ''),
            ])
    print("  S1 done (core loops)")


# ── Table S2: Ortholog conservation matrix ────────────────────────────
def table_s2():
    orthologs = load_csv('loop_ortholog_matrix.csv')
    with open(SUPP / 'table_s2_ortholog_matrix.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        species = ['M.musculus','R.norvegicus','C.l.familiaris','G.gallus',
                   'X.tropicalis','D.rerio','C.intestinalis','S.purpuratus',
                   'D.melanogaster','C.elegans','N.vectensis',
                   'S.cerevisiae','S.pombe','A.thaliana','D.discoideum']
        writer.writerow(['Loop ID', 'Pathway', 'Activator', 'Inhibitor'] + species + ['N Species'])
        for row in orthologs:
            writer.writerow([
                row.get('loop_id', ''),
                row.get('pathway_name', ''),
                row.get('activator', ''),
                row.get('inhibitor', ''),
            ] + [row.get(sp, '0') for sp in species] + [row.get('n_species_both', '')])
    print("  S2 done (ortholog matrix)")


# ── Table S3: 128 unique NFL motifs ──────────────────────────────────
def table_s3():
    motifs = load_csv('unique_nfl_motifs.csv')
    with open(SUPP / 'table_s3_unique_motifs.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Motif Genes (HGNC)', 'Length', 'N Inhibitory Edges',
                        'N Pathways', 'N Variants'])
        for row in motifs:
            writer.writerow([
                row.get('genes_hgnc', ''),
                row.get('length', ''),
                row.get('n_inh', ''),
                row.get('n_pw', ''),
                row.get('n_variants', ''),
            ])
    print(f"  S3 done ({len(motifs)} unique motifs)")


# ── Table S4: 59 non-redundant modules ───────────────────────────────
def table_s4():
    modules = load_csv('nfl_modules.csv')
    with open(SUPP / 'table_s4_modules.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Module ID', 'Pathway Label', 'N Motifs', 'N Genes',
                        'N CGC', 'All Genes', 'Representative Genes'])
        for row in modules:
            writer.writerow([
                row.get('module_id', ''),
                row.get('pathway_label', ''),
                row.get('n_motifs', ''),
                row.get('n_genes', ''),
                row.get('n_cgc', ''),
                row.get('all_genes', ''),
                row.get('representative', ''),
            ])
    print(f"  S4 done ({len(modules)} modules)")


# ── Table S5: Community detection at 3 resolutions ───────────────────
def table_s5():
    comm = load_json('community_detection_results.json')
    with open(SUPP / 'table_s5_community_detection.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Resolution', 'N Communities', 'Modularity Q', 'ARI vs KEGG'])
        for res in ['0.5', '1.0', '1.5']:
            d = comm['louvain_by_resolution'][res]
            writer.writerow([res, d['n_communities'], f"{d['Q']:.3f}", f"{d['ARI']:.3f}"])
        writer.writerow([])
        writer.writerow(['--- Resolution 1.0 Community Details ---'])
        writer.writerow(['Community', 'N Motifs', 'N Genes', 'Dominant KEGG Pathway', 'Top Genes'])
        for c in comm['resolution_1_0']['communities']:
            top = c.get('top_genes_by_frequency', [])[:5]
            writer.writerow([
                c['community_idx'],
                c['n_motifs'],
                c['n_genes'],
                c.get('dominant_kegg_pathway', c.get('auto_label', '')),
                '|'.join(top) if top else '',
            ])
    print("  S5 done (community detection)")


# ── Table S6: Biological validation ──────────────────────────────────
def table_s6():
    val = load_csv('bio_validation_20motifs.csv')
    with open(SUPP / 'table_s6_bio_validation.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Motif Genes', 'Verdict', 'Notes',
                        'Edge1', 'Edge1 Valid', 'Edge2', 'Edge2 Valid',
                        'Edge3', 'Edge3 Valid',
                        'Feedback Described', 'Could Oscillate'])
        for row in val:
            writer.writerow([
                row.get('motif_genes', ''),
                row.get('verdict', ''),
                row.get('notes', ''),
                row.get('edge1', ''),
                row.get('edge1_valid', ''),
                row.get('edge2', ''),
                row.get('edge2_valid', ''),
                row.get('edge3', ''),
                row.get('edge3_valid', ''),
                row.get('feedback_described', ''),
                row.get('could_oscillate', ''),
            ])
    print(f"  S6 done ({len(val)} validated motifs)")


# ── Table S7: KEGG signaling gene list with NFL/CGC status ───────────
def table_s7():
    # We don't have the full 16760 gene list as a single file.
    # Create a summary from available data.
    cgc = load_json('cgc_enrichment_CORRECTED.json')
    with open(SUPP / 'table_s7_kegg_gene_summary.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'N Genes', 'N CGC', 'CGC %'])
        writer.writerow(['NFL genes', cgc['nfl_genes'], cgc['nfl_cgc'], f"{cgc['nfl_cgc_pct']:.1f}"])
        writer.writerow(['Non-NFL KEGG genes', cgc['kegg_genes'] - cgc['nfl_genes'],
                        cgc['non_nfl_cgc'], f"{cgc['non_nfl_cgc_pct']:.1f}"])
        writer.writerow(['Total KEGG signaling', cgc['kegg_genes'],
                        cgc['nfl_cgc'] + cgc['non_nfl_cgc'], ''])
        writer.writerow([])
        writer.writerow(['Enrichment Test', 'Fold/OR', 'p-value', ''])
        writer.writerow(['Test 1: NFL vs Genome', f"{cgc['test1_fold']:.1f}x", f"{cgc['test1_p']:.2e}", ''])
        writer.writerow(['Test 2: NFL vs KEGG', f"{cgc['test2_fold']:.1f}x", f"{cgc['test2_p']:.2e}", ''])
        writer.writerow(['Test 3: NFL vs non-NFL', f"OR={cgc['test3_OR']:.1f}", f"{cgc['test3_p']:.2e}", ''])
    print("  S7 done (KEGG gene summary)")


# ── Table S8: Irreversible Authority scores ──────────────────────────
def table_s8():
    vuln = load_json('vulnerability_metric.json')
    ia = load_json('irreversible_authority.json')

    IA_SCORES = {
        "NF-κB": 2, "ERK/MAPK": 2, "JAK-STAT": 2, "p53": 3, "Wnt": 3,
        "Notch": 3, "Hippo": 2, "TGF-β": 2, "mTOR": 1, "Calcium": 1,
        "Cell Cycle": 2, "Circadian": 0, "NRF2": 0, "PI3K/PTEN": 2,
        "AMPK": 1, "SREBP": 0, "ATR/CHK1": 2, "Rho/ROCK": 1,
        "PPAR/LXR": 0, "Autophagy": 1,
    }

    IA_DECISIONS = {
        "p53": "apoptosis, senescence, DNA repair commitment",
        "Wnt": "terminal differentiation, stem cell fate, EMT",
        "Notch": "lateral inhibition, lineage commitment, angiogenesis",
        "NF-κB": "inflammatory commitment, survival decision",
        "ERK/MAPK": "proliferation vs differentiation",
        "JAK-STAT": "immune activation, differentiation",
        "Hippo": "organ size, contact inhibition",
        "TGF-β": "EMT, growth arrest",
        "Cell Cycle": "mitotic commitment, cytokinesis",
        "PI3K/PTEN": "survival, metabolic commitment",
        "ATR/CHK1": "replication checkpoint, fork restart",
        "mTOR": "growth commitment",
        "AMPK": "metabolic switch",
        "Calcium": "NFAT activation",
        "Rho/ROCK": "cytoskeletal commitment",
        "Autophagy": "autophagic flux",
        "Circadian": "(reversible oscillation)",
        "NRF2": "(reversible stress response)",
        "SREBP": "(reversible lipid regulation)",
        "PPAR/LXR": "(reversible transcription)",
    }

    with open(SUPP / 'table_s8_irreversible_authority.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Module', 'IA Score', 'Authority', 'N Genes', 'N CGC',
                        'CGC Fraction', 'Irreversible Decisions'])
        for m in vuln['module_data']:
            name = m['module']
            writer.writerow([
                name,
                IA_SCORES.get(name, ''),
                m['authority'],
                m['n_genes'],
                m['n_cgc'],
                f"{m['cgc_frac']:.3f}",
                IA_DECISIONS.get(name, ''),
            ])
        writer.writerow([])
        writer.writerow(['Metric', 'Spearman rho', 'p-value', '', '', '', ''])
        writer.writerow(['IA', f"{ia['ia_spearman_rho']:.3f}", f"{ia['ia_spearman_p']}", '', '', '', ''])
        writer.writerow(['Authority', f"{ia['auth_spearman_rho']:.3f}", f"{ia['auth_spearman_p']}", '', '', '', ''])
    print("  S8 done (IA scores)")


# ── Table S9: Eigenspace cancer candidates ───────────────────────────
def table_s9():
    eigen = load_json('eigenspace_cancer.json')
    notch = load_json('notch_validation.json')

    candidates = eigen['step4']['top_candidates']
    validated = set(notch.get('validated', []))

    with open(SUPP / 'table_s9_cancer_candidates.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Rank', 'Gene', 'IntOGen Validated', 'Source'])
        for i, gene in enumerate(candidates, 1):
            writer.writerow([i, gene, 'Yes' if gene in validated else 'No',
                           'Eigenspace proximity to CGC centroid'])
        writer.writerow([])
        writer.writerow(['Summary', '', '', ''])
        writer.writerow(['Total candidates', len(candidates), '', ''])
        writer.writerow(['IntOGen validated', notch['n_intogen'], '', ''])
        writer.writerow(['Precision', f"{notch['precision_intogen']:.0%}", '', ''])
        writer.writerow(['Enrichment p', f"{notch['intogen_enrichment_p']}", '', ''])
    print(f"  S9 done ({len(candidates)} candidates)")


# ── Table S10: Tissue-specific eigendecomposition ────────────────────
def table_s10():
    tissue = load_json('tissue_specific_eigen.json')

    with open(SUPP / 'table_s10_tissue_eigen.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Tissue Class', 'N Cell Types', 'N Kaiser', 'PC1 Var%',
                        'PC1 Sim to Global', 'PC1 Top Modules'])
        # Global first
        g = tissue['global']
        writer.writerow(['GLOBAL (all)', g['n_cell_types'], g['n_kaiser'],
                        f"{g['variance_pct'][0]:.1f}", '1.000', '(reference)'])
        for t_name, t_data in sorted(tissue['tissues'].items()):
            pc1_sim = t_data['pcs'][0]['sim_global'] if t_data['pcs'] else ''
            top_mods = ', '.join([f"{m}({v:.2f})" for m, v in t_data['pcs'][0]['top_modules']])  if t_data['pcs'] else ''
            writer.writerow([
                t_name,
                t_data['n_cell_types'],
                t_data['n_kaiser'],
                f"{t_data['variance_pct'][0]:.1f}",
                f"{pc1_sim:.3f}" if isinstance(pc1_sim, float) else pc1_sim,
                top_mods,
            ])
        writer.writerow([])
        writer.writerow(['Verdict', tissue['verdict'], '', '', '', ''])
        writer.writerow(['Mean PC1 similarity', f"{tissue['pc1_mean_similarity']:.3f}", '', '', '', ''])
    print("  S10 done (tissue eigen)")


# ── Table S11: MOCA developmental trajectories ──────────────────────
def table_s11():
    moca = load_json('moca_differentiation.json')

    with open(SUPP / 'table_s11_moca_trajectories.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Trajectory', 'Cancer Relevance', 'rho(day, E_intra)', 'p-value',
                        'Direction', 'E_intra E9.5', 'E_intra E13.5', 'IA Score (adult)'])

        for traj_name, traj_data in sorted(moca['trajectory_stats'].items()):
            eigen_pos = moca['eigenspace_positions'].get(traj_name, {})
            e_means = traj_data.get('E_intra_means', [])
            writer.writerow([
                traj_name,
                traj_data.get('cancer_relevance', ''),
                f"{traj_data['rho_day_eintra']:.3f}" if traj_data['rho_day_eintra'] else '',
                f"{traj_data['p_day_eintra']:.4f}" if traj_data['p_day_eintra'] else '',
                traj_data.get('direction', ''),
                f"{e_means[0]:.3f}" if e_means else '',
                f"{e_means[-1]:.3f}" if e_means else '',
                f"{eigen_pos.get('ia_score', '')}" if eigen_pos.get('ia_score') else '',
            ])
        writer.writerow([])
        writer.writerow(['Summary', '', '', '', '', '', '', ''])
        writer.writerow(['N trajectories', moca['n_trajectories'], '', '', '', '', '', ''])
        writer.writerow(['All decreasing', len(moca['decreasing_trajectories']), '', '', '', '', '', ''])
        writer.writerow(['Dev-aging replication rho',
                        f"{moca['dev_aging_replication']['rho']:.3f}" if moca['dev_aging_replication']['rho'] else '',
                        '', '', '', '', '', ''])
    print("  S11 done (MOCA trajectories)")


# ── Table: Loop ages (bonus, used in Fig 4) ──────────────────────────
def table_ages():
    ages = load_csv('loop_ages.csv')
    with open(SUPP / 'table_loop_ages.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Loop ID', 'Pathway', 'Oldest Species', 'Age (Mya)',
                        'N Species', 'Confidence', 'Asymmetry (A)', 'Feedback Type', 'Cluster'])
        for row in ages:
            writer.writerow([
                row.get('loop_id', ''),
                row.get('pathway', ''),
                row.get('oldest_species', ''),
                row.get('age_mya', ''),
                row.get('n_species_both', ''),
                row.get('confidence', ''),
                row.get('A', ''),
                row.get('feedback_type', ''),
                row.get('cluster', ''),
            ])
    print(f"  Loop ages table done ({len(ages)} loops)")


if __name__ == "__main__":
    print("Generating Supplementary Tables...")
    table_s1()
    table_s2()
    table_s3()
    table_s4()
    table_s5()
    table_s6()
    table_s7()
    table_s8()
    table_s9()
    table_s10()
    table_s11()
    table_ages()
    print(f"\nAll tables saved to {SUPP}/")
    import os
    for f in sorted(os.listdir(SUPP)):
        size = os.path.getsize(SUPP / f)
        print(f"  {f}: {size:,} bytes")
