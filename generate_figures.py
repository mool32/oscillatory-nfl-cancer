#!/usr/bin/env python3
"""
Generate all 7 main figures + supplementary for Paper 2 manuscript.
"""

import json
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker

# ── Style ─────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.linewidth': 0.8,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

DATA = Path("/Users/teo/Desktop/research/paper2/data")
FIGS = Path("/Users/teo/Desktop/research/paper2/figures")
FIGS.mkdir(exist_ok=True)

# Color palette
C_NFL = '#2166AC'       # blue - NFL genes
C_NON = '#B2182B'       # red - non-NFL
C_CGC = '#E69F00'       # gold - CGC
C_NONCGC = '#999999'    # gray
C_HIGH = '#D73027'      # red for high cancer
C_MED = '#FC8D59'       # orange for medium
C_LOW = '#91BFDB'       # light blue for low
C_ACCENT = '#4DAF4A'    # green accent
C_BG = '#F7F7F7'


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


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1: NFL Extraction Pipeline
# ══════════════════════════════════════════════════════════════════════════
def figure1():
    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # Panel A: Pipeline schematic (text-based)
    ax = fig.add_subplot(gs[0, 0])
    steps = [
        ("159 KEGG\nKGML files", 0.9),
        ("Signed\nDiGraphs", 0.72),
        ("Depth-limited\nDFS (d\u22643)", 0.54),
        ("228 raw\nNFL motifs", 0.36),
        ("Alias resolution\n+ dedup", 0.18),
        ("128 unique\nNFLs", 0.0),
    ]
    for i, (label, y) in enumerate(steps):
        color = C_NFL if i in [3, 5] else '#DDDDDD'
        ax.add_patch(plt.Rectangle((0.15, y), 0.7, 0.14,
                     facecolor=color, edgecolor='black', linewidth=0.8, alpha=0.7))
        ax.text(0.5, y + 0.07, label, ha='center', va='center', fontsize=7.5,
                fontweight='bold' if i in [3, 5] else 'normal',
                color='white' if i in [3, 5] else 'black')
        if i < len(steps) - 1:
            ax.annotate('', xy=(0.5, steps[i+1][1] + 0.14), xytext=(0.5, y),
                       arrowprops=dict(arrowstyle='->', color='black', lw=1))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.05)
    ax.axis('off')
    ax.set_title('A  Extraction pipeline', loc='left', fontweight='bold')

    # Panel B: Example motifs
    ax = fig.add_subplot(gs[0, 1])
    # Two-node: p53 <-> MDM2
    circle_kw = dict(radius=0.12, facecolor='white', edgecolor='black', linewidth=1.2)
    ax.add_patch(plt.Circle((0.25, 0.75), **circle_kw))
    ax.text(0.25, 0.75, 'TP53', ha='center', va='center', fontsize=8, fontweight='bold')
    ax.add_patch(plt.Circle((0.75, 0.75), **circle_kw))
    ax.text(0.75, 0.75, 'MDM2', ha='center', va='center', fontsize=8, fontweight='bold')
    ax.annotate('', xy=(0.62, 0.78), xytext=(0.38, 0.78),
               arrowprops=dict(arrowstyle='->', color=C_ACCENT, lw=2))
    ax.text(0.5, 0.83, '+', ha='center', color=C_ACCENT, fontsize=10, fontweight='bold')
    ax.annotate('', xy=(0.38, 0.72), xytext=(0.62, 0.72),
               arrowprops=dict(arrowstyle='-|>', color=C_NON, lw=2))
    ax.text(0.5, 0.66, '\u2013', ha='center', color=C_NON, fontsize=14, fontweight='bold')
    ax.text(0.5, 0.58, '2-node NFL', ha='center', fontsize=8, style='italic')

    # Three-node: E2F1 -> CCNE1 -> RB1 -| E2F1
    for (x, y, name) in [(0.25, 0.25, 'E2F1'), (0.75, 0.25, 'CCNE1'), (0.5, 0.0, 'RB1')]:
        ax.add_patch(plt.Circle((x, y), **circle_kw))
        ax.text(x, y, name, ha='center', va='center', fontsize=7, fontweight='bold')
    ax.annotate('', xy=(0.62, 0.28), xytext=(0.38, 0.28),
               arrowprops=dict(arrowstyle='->', color=C_ACCENT, lw=1.5))
    ax.annotate('', xy=(0.55, 0.12), xytext=(0.70, 0.16),
               arrowprops=dict(arrowstyle='->', color=C_ACCENT, lw=1.5))
    ax.annotate('', xy=(0.30, 0.16), xytext=(0.45, 0.12),
               arrowprops=dict(arrowstyle='-|>', color=C_NON, lw=1.5))
    ax.text(0.5, -0.12, '3-node NFL', ha='center', fontsize=8, style='italic')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.2, 1.0)
    ax.axis('off')
    ax.set_title('B  Example NFL motifs', loc='left', fontweight='bold')

    # Panel C: Deduplication cascade (bar chart)
    ax = fig.add_subplot(gs[1, 0])
    stages = ['Raw\nmotifs', 'Unique\n(dedup)', 'Modules\n(Jaccard)', 'Pathway\ngroups']
    counts = [228, 128, 59, 14]
    colors = ['#AAAAAA', '#7FBADC', C_NFL, '#1A4472']
    bars = ax.bar(range(4), counts, color=colors, edgecolor='black', linewidth=0.5, width=0.65)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(count), ha='center', fontweight='bold', fontsize=10)
    ax.set_xticks(range(4))
    ax.set_xticklabels(stages, fontsize=8)
    ax.set_ylabel('Count')
    ax.set_ylim(0, 270)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('C  Deduplication cascade', loc='left', fontweight='bold')

    # Panel D: Biological validation pie
    ax = fig.add_subplot(gs[1, 1])
    sizes = [7, 9, 4]
    labels = ['Full validation\n(35%)', 'Partial support\n(45%)', 'Artifact\n(20%)']
    colors_pie = [C_ACCENT, C_CGC, '#CCCCCC']
    wedges, texts = ax.pie(sizes, labels=labels, colors=colors_pie,
                            startangle=90, counterclock=False,
                            wedgeprops=dict(edgecolor='black', linewidth=0.5))
    for t in texts:
        t.set_fontsize(8)
    ax.text(0, -1.4, 'N = 20 randomly sampled motifs', ha='center', fontsize=8, style='italic')
    ax.set_title('D  Biological validation', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig1_pipeline.pdf')
    fig.savefig(FIGS / 'fig1_pipeline.png')
    plt.close(fig)
    print("  Fig 1 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2: CGC Enrichment (MAIN RESULT)
# ══════════════════════════════════════════════════════════════════════════
def figure2():
    cgc = load_json('cgc_enrichment_CORRECTED.json')

    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # Panel A: Three-test waterfall
    ax = fig.add_subplot(gs[0, 0])
    tests = ['vs Genome\n(Test 1)', 'vs KEGG\n(Test 2)', 'vs non-NFL\n(Test 3)']
    values = [cgc['test1_fold'], cgc['test2_fold'], cgc['test3_OR']]
    p_vals = [cgc['test1_p'], cgc['test2_p'], cgc['test3_p']]
    colors_bar = ['#91BFDB', C_NFL, '#1A4472']
    bars = ax.bar(range(3), values, color=colors_bar, edgecolor='black', linewidth=0.5, width=0.6)
    for i, (bar, val, p) in enumerate(zip(bars, values, p_vals)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                f'{val:.1f}x', ha='center', fontweight='bold', fontsize=10)
        # p-value annotation
        p_str = f'p = {p:.0e}' if p > 0 else 'p < 1e-44'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 6,
                p_str, ha='center', fontsize=6.5, color='#555555')
    ax.set_xticks(range(3))
    ax.set_xticklabels(tests, fontsize=8)
    ax.set_ylabel('Enrichment (fold or OR)')
    ax.set_ylim(0, 75)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('A  Three nested enrichment tests', loc='left', fontweight='bold')

    # Panel B: 2x2 table
    ax = fig.add_subplot(gs[0, 1])
    table = cgc['test3_table']
    ax.axis('off')
    cell_text = [
        [str(table[0][0]), str(table[0][1]), str(sum(table[0]))],
        [str(table[1][0]), str(table[1][1].replace(',','') if isinstance(table[1][1], str) else table[1][1]),
         str(sum(table[1]))],
    ]
    # Fix: table values are ints
    cell_text = [
        [str(table[0][0]), str(table[0][1]), str(table[0][0]+table[0][1])],
        [str(table[1][0]), str(table[1][1]), str(table[1][0]+table[1][1])],
    ]
    col_labels = ['CGC', 'Not CGC', 'Total']
    row_labels = ['NFL gene', 'Non-NFL']

    tbl = ax.table(cellText=cell_text, colLabels=col_labels, rowLabels=row_labels,
                   loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.8)
    # Color the key cell
    tbl[1, 0].set_facecolor('#FFE0B2')
    tbl[1, 0].set_text_props(fontweight='bold')
    ax.text(0.5, 0.15, f'OR = {cgc["test3_OR"]:.1f},  p = 9\u00d710\u207b\u2074\u2074',
            ha='center', transform=ax.transAxes, fontsize=10, fontweight='bold', color=C_NON)
    ax.set_title('B  Critical test: NFL vs non-NFL', loc='left', fontweight='bold')

    # Panel C: Convergent validation
    ax = fig.add_subplot(gs[1, 0])
    methods = ['Paper 1\n(temporal)', 'Paper 2\n(topological)']
    vals = [cgc['paper1_fold'], cgc['paper2_fold']]
    bars = ax.bar(range(2), vals, color=[C_CGC, C_NFL], edgecolor='black', linewidth=0.5, width=0.5)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{val:.1f}x', ha='center', fontweight='bold', fontsize=11)
    ax.set_xticks(range(2))
    ax.set_xticklabels(methods, fontsize=9)
    ax.set_ylabel('CGC enrichment vs genome')
    ax.set_ylim(0, 16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=1, color='black', ls='--', lw=0.5, alpha=0.3)
    ax.set_title('C  Convergent validation', loc='left', fontweight='bold')

    # Panel D: Venn diagram
    ax = fig.add_subplot(gs[1, 1])
    # Simple 2-circle Venn
    from matplotlib.patches import Circle as MplCircle
    c1 = MplCircle((0.38, 0.5), 0.3, facecolor=C_CGC, alpha=0.4, edgecolor='black', lw=1.2)
    c2 = MplCircle((0.62, 0.5), 0.3, facecolor=C_NFL, alpha=0.4, edgecolor='black', lw=1.2)
    ax.add_patch(c1)
    ax.add_patch(c2)
    ax.text(0.25, 0.5, f'{cgc["venn_paper1_only"]}', ha='center', va='center',
            fontsize=16, fontweight='bold')
    ax.text(0.5, 0.5, f'{cgc["venn_both"]}', ha='center', va='center',
            fontsize=16, fontweight='bold')
    ax.text(0.75, 0.5, f'{cgc["venn_nfl_only"]}', ha='center', va='center',
            fontsize=16, fontweight='bold')
    ax.text(0.25, 0.12, 'Paper 1\nonly', ha='center', fontsize=8)
    ax.text(0.75, 0.12, 'NFL\nonly', ha='center', fontsize=8)
    ax.text(0.5, 0.88, 'Shared', ha='center', fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('D  CGC genes by approach', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig2_enrichment.pdf')
    fig.savefig(FIGS / 'fig2_enrichment.png')
    plt.close(fig)
    print("  Fig 2 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3: Modularity and Cross-Talk
# ══════════════════════════════════════════════════════════════════════════
def figure3():
    comm = load_json('community_detection_results.json')

    fig = plt.figure(figsize=(7.5, 4.5))
    gs = GridSpec(1, 3, wspace=0.35, width_ratios=[1.2, 1, 0.8])

    # Panel A: Community sizes
    ax = fig.add_subplot(gs[0, 0])
    communities = comm['resolution_1_0']['communities']
    comm_sorted = sorted(communities, key=lambda c: c['n_motifs'], reverse=True)
    labels_comm = []
    sizes_comm = []
    for c in comm_sorted:
        lbl = c.get('auto_label', c.get('dominant_kegg_pathway', f"C{c['community_idx']}"))
        if len(lbl) > 20:
            lbl = lbl[:18] + '..'
        labels_comm.append(lbl)
        sizes_comm.append(c['n_motifs'])

    cmap = plt.cm.Set3(np.linspace(0, 1, len(labels_comm)))
    bars = ax.barh(range(len(labels_comm)), sizes_comm, color=cmap, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(labels_comm)))
    ax.set_yticklabels(labels_comm, fontsize=7)
    ax.set_xlabel('Number of motifs')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Add Q annotation
    Q = comm['resolution_1_0']['modularity_Q']
    ax.text(0.95, 0.95, f'Q = {Q:.3f}', transform=ax.transAxes, ha='right', va='top',
            fontsize=9, fontweight='bold', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    ax.set_title('A  Louvain communities', loc='left', fontweight='bold')

    # Panel B: ARI across resolutions
    ax = fig.add_subplot(gs[0, 1])
    resolutions = sorted(comm['louvain_by_resolution'].keys(), key=float)
    aris = [comm['louvain_by_resolution'][r]['ARI'] for r in resolutions]
    qs = [comm['louvain_by_resolution'][r]['Q'] for r in resolutions]
    x = [float(r) for r in resolutions]

    ax.plot(x, aris, 'o-', color=C_NFL, lw=2, markersize=8, label='ARI', zorder=3)
    ax2 = ax.twinx()
    ax2.plot(x, qs, 's--', color=C_CGC, lw=2, markersize=8, label='Q', zorder=3)
    ax.set_xlabel('Louvain resolution')
    ax.set_ylabel('ARI', color=C_NFL)
    ax2.set_ylabel('Modularity Q', color=C_CGC)
    ax.set_ylim(0, 0.7)
    ax2.set_ylim(0, 0.7)
    ax.spines['top'].set_visible(False)
    ax.set_title('B  Resolution sensitivity', loc='left', fontweight='bold')

    # Panel C: Cross-talk hub genes
    ax = fig.add_subplot(gs[0, 2])
    cross = comm['pathway_overlap_network']['cross_talk_genes']
    # Sort by number of pathways
    hubs = sorted(cross.items(), key=lambda x: len(x[1]), reverse=True)[:8]
    genes_h = [h[0] for h in hubs]
    n_pw = [len(h[1]) for h in hubs]
    colors_h = [C_CGC if g in ['MAPK1','TP53','SMAD2','SMAD3','PIK3CA','STAT3','KRAS','AKT1']
                else C_NONCGC for g in genes_h]
    ax.barh(range(len(genes_h)), n_pw, color=colors_h, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(genes_h)))
    ax.set_yticklabels(genes_h, fontsize=8, fontweight='bold')
    ax.set_xlabel('Pathways bridged')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Legend
    ax.legend(handles=[
        mpatches.Patch(color=C_CGC, label='CGC gene'),
        mpatches.Patch(color=C_NONCGC, label='Non-CGC'),
    ], fontsize=7, loc='lower right')
    ax.set_title('C  Hub genes', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig3_modularity.pdf')
    fig.savefig(FIGS / 'fig3_modularity.png')
    plt.close(fig)
    print("  Fig 3 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4: Evolutionary Conservation
# ══════════════════════════════════════════════════════════════════════════
def figure4():
    ages = load_csv('loop_ages.csv')
    orthologs = load_csv('loop_ortholog_matrix.csv')

    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # Panel A: Ortholog heatmap
    ax = fig.add_subplot(gs[0, :])
    species_order = ['M.musculus','R.norvegicus','C.l.familiaris','G.gallus',
                     'X.tropicalis','D.rerio','C.intestinalis','S.purpuratus',
                     'D.melanogaster','C.elegans','N.vectensis',
                     'S.cerevisiae','S.pombe','A.thaliana','D.discoideum']
    species_short = ['Mouse','Rat','Dog','Chicken','Frog','Zebrafish',
                     'Ciona','Urchin','Fly','Worm','Anemone',
                     'S.cer','S.pom','Arab.','Dicty.']

    # Build matrix
    loop_names = [row['pathway_name'] for row in orthologs]
    mat = np.zeros((len(orthologs), len(species_order)))
    for i, row in enumerate(orthologs):
        for j, sp in enumerate(species_order):
            mat[i, j] = int(row.get(sp, 0))

    im = ax.imshow(mat, cmap='Blues', aspect='auto', interpolation='nearest')
    ax.set_xticks(range(len(species_short)))
    ax.set_xticklabels(species_short, rotation=45, ha='right', fontsize=6.5)
    ax.set_yticks(range(len(loop_names)))
    ax.set_yticklabels([n[:25] for n in loop_names], fontsize=6)
    # Vertical line separating vertebrates/invertebrates
    ax.axvline(x=5.5, color='red', ls='--', lw=1, alpha=0.7)
    ax.axvline(x=10.5, color='red', ls='--', lw=1, alpha=0.7)
    ax.text(2.5, -1.5, 'Vertebrates', ha='center', fontsize=7, color='red')
    ax.text(8, -1.5, 'Invertebrates', ha='center', fontsize=7, color='red')
    ax.text(13, -1.5, 'Unicellular+', ha='center', fontsize=7, color='red')
    ax.set_title('A  Ortholog conservation matrix (22 loops \u00d7 15 species)', loc='left', fontweight='bold')

    # Panel C: Age vs A scatter
    ax = fig.add_subplot(gs[1, 0])
    age_vals = []
    a_vals = []
    names_a = []
    nrf2_idx = None
    for row in ages:
        try:
            age = float(row['age_mya'])
            A = float(row['A'])
            age_vals.append(age)
            a_vals.append(A)
            names_a.append(row['pathway'])
            if 'NRF2' in row['pathway'] or 'Nrf2' in row['pathway']:
                nrf2_idx = len(age_vals) - 1
        except (ValueError, KeyError):
            pass

    ax.scatter(age_vals, a_vals, c=C_NFL, s=50, edgecolor='black', linewidth=0.5, zorder=3)
    if nrf2_idx is not None:
        ax.scatter([age_vals[nrf2_idx]], [a_vals[nrf2_idx]], c='red', s=80,
                   edgecolor='black', linewidth=1, zorder=4, marker='D')
        ax.annotate('NRF2/KEAP1', (age_vals[nrf2_idx], a_vals[nrf2_idx]),
                   xytext=(10, -15), textcoords='offset points', fontsize=7,
                   arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

    # Fit line without NRF2
    from scipy import stats as sp_stats
    ages_no_nrf2 = [a for i, a in enumerate(age_vals) if i != nrf2_idx]
    a_no_nrf2 = [a for i, a in enumerate(a_vals) if i != nrf2_idx]
    if len(ages_no_nrf2) >= 3:
        slope, intercept, r, p, se = sp_stats.linregress(ages_no_nrf2, a_no_nrf2)
        x_fit = np.linspace(min(age_vals), max(age_vals), 100)
        ax.plot(x_fit, slope * x_fit + intercept, '--', color='gray', lw=1, alpha=0.7)
        rho, p_rho = sp_stats.spearmanr(ages_no_nrf2, a_no_nrf2)
        ax.text(0.05, 0.95, f'\u03c1 = {rho:.2f}, p = {p_rho:.3f}\n(without NRF2)',
                transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    ax.set_xlabel('Loop age (Mya)')
    ax.set_ylabel('Waveform asymmetry (A)')
    ax.axhline(y=0.35, color='gray', ls=':', lw=0.8, alpha=0.5)
    ax.text(1550, 0.36, 'Class I/II boundary', fontsize=6.5, color='gray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('C  Age vs waveform asymmetry', loc='left', fontweight='bold')

    # Panel D: Conservation summary bar
    ax = fig.add_subplot(gs[1, 1])
    categories = ['Vertebrate\n(6 spp)', 'Invertebrate\n(5 spp)', 'Unicellular\n(4 spp)']
    # Count loops present in at least one species per category
    vert_spp = ['M.musculus','R.norvegicus','C.l.familiaris','G.gallus','X.tropicalis','D.rerio']
    invert_spp = ['C.intestinalis','S.purpuratus','D.melanogaster','C.elegans','N.vectensis']
    uni_spp = ['S.cerevisiae','S.pombe','A.thaliana','D.discoideum']

    def count_loops_in(species_list):
        count = 0
        for row in orthologs:
            if any(int(row.get(sp, 0)) == 1 for sp in species_list):
                count += 1
        return count

    counts = [count_loops_in(vert_spp), count_loops_in(invert_spp), count_loops_in(uni_spp)]
    pcts = [c/len(orthologs)*100 for c in counts]
    bars = ax.bar(range(3), pcts, color=[C_ACCENT, C_CGC, C_NON],
                  edgecolor='black', linewidth=0.5, width=0.55)
    for bar, pct, cnt in zip(bars, pcts, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{cnt}/{len(orthologs)}\n({pct:.0f}%)', ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(range(3))
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylabel('Loops conserved (%)')
    ax.set_ylim(0, 115)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('D  Conservation gradient', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig4_evolution.pdf')
    fig.savefig(FIGS / 'fig4_evolution.png')
    plt.close(fig)
    print("  Fig 4 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 5: Irreversible Authority
# ══════════════════════════════════════════════════════════════════════════
def figure5():
    ia = load_json('irreversible_authority.json')
    vuln = load_json('vulnerability_metric.json')

    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # Module data from vulnerability_metric
    modules = vuln['module_data']

    # IA scores (from manuscript table)
    IA_SCORES = {
        "NF-\u03baB": 2, "ERK/MAPK": 2, "JAK-STAT": 2, "p53": 3, "Wnt": 3,
        "Notch": 3, "Hippo": 2, "TGF-\u03b2": 2, "mTOR": 1, "Calcium": 1,
        "Cell Cycle": 2, "Circadian": 0, "NRF2": 0, "PI3K/PTEN": 2,
        "AMPK": 1, "SREBP": 0, "ATR/CHK1": 2, "Rho/ROCK": 1,
        "PPAR/LXR": 0, "Autophagy": 1,
    }

    # Panel A: IA vs CGC fraction scatter
    ax = fig.add_subplot(gs[0, 0])
    for m in modules:
        name = m['module']
        ia_val = IA_SCORES.get(name, 0)
        cgc_frac = m['cgc_frac']
        # Jitter IA slightly for visibility
        jitter = np.random.uniform(-0.08, 0.08)
        color = C_HIGH if cgc_frac > 0.4 else (C_MED if cgc_frac > 0.15 else C_LOW)
        ax.scatter(ia_val + jitter, cgc_frac, c=color, s=60, edgecolor='black',
                   linewidth=0.5, zorder=3)
        if cgc_frac > 0.5 or (ia_val >= 3) or name in ['Circadian', 'SREBP']:
            ax.annotate(name, (ia_val + jitter, cgc_frac),
                       xytext=(5, 5), textcoords='offset points', fontsize=6)

    ax.set_xlabel('Irreversible Authority (IA)')
    ax.set_ylabel('CGC fraction')
    ax.set_xticks([0, 1, 2, 3])
    ax.text(0.05, 0.95, f'rho = {ia["ia_spearman_rho"]:.3f}\np < 1e-5',
            transform=ax.transAxes, fontsize=9, va='top', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('A  IA predicts CGC fraction', loc='left', fontweight='bold')

    # Panel B: Dose-response
    ax = fig.add_subplot(gs[0, 1])
    ia_groups = defaultdict(list)
    for m in modules:
        name = m['module']
        ia_val = IA_SCORES.get(name, 0)
        ia_groups[ia_val].append(m['cgc_frac'])

    ia_levels = sorted(ia_groups.keys())
    means = [np.mean(ia_groups[k]) for k in ia_levels]
    sems = [np.std(ia_groups[k]) / np.sqrt(len(ia_groups[k])) if len(ia_groups[k]) > 1 else 0
            for k in ia_levels]
    ns = [len(ia_groups[k]) for k in ia_levels]

    bars = ax.bar(ia_levels, means, yerr=sems, color=[C_LOW, C_MED, C_CGC, C_HIGH],
                  edgecolor='black', linewidth=0.5, width=0.6, capsize=4)
    for x, n in zip(ia_levels, ns):
        ax.text(x, -0.04, f'N={n}', ha='center', fontsize=7, color='gray')
    ax.set_xlabel('IA score')
    ax.set_ylabel('Mean CGC fraction')
    ax.set_xticks(ia_levels)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('B  Dose-response', loc='left', fontweight='bold')

    # Panel C: Three hypotheses comparison
    ax = fig.add_subplot(gs[1, 0])
    hypotheses = ['Irreversible\nAuthority', 'Authority\n(all)', 'IA + Bottleneck', 'Maintenance\ncost']
    rhos = [0.829, 0.752, 0.869, 0.03]
    colors_h = [C_HIGH, C_CGC, '#4DAF4A', C_NONCGC]
    bars = ax.bar(range(4), rhos, color=colors_h, edgecolor='black', linewidth=0.5, width=0.6)
    for bar, rho in zip(bars, rhos):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{rho:.3f}', ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(range(4))
    ax.set_xticklabels(hypotheses, fontsize=7.5)
    ax.set_ylabel('Spearman \u03c1 with CGC fraction')
    ax.set_ylim(0, 1.0)
    ax.axhline(y=0, color='black', lw=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('C  Hypothesis comparison', loc='left', fontweight='bold')

    # Panel D: Binary split
    ax = fig.add_subplot(gs[1, 1])
    irrev_fracs = [m['cgc_frac'] for m in modules if IA_SCORES.get(m['module'], 0) >= 2]
    rev_fracs = [m['cgc_frac'] for m in modules if IA_SCORES.get(m['module'], 0) < 2]

    parts = ax.violinplot([rev_fracs, irrev_fracs], positions=[0, 1], showmeans=True, showmedians=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor([C_LOW, C_HIGH][i])
        pc.set_alpha(0.7)
    parts['cmeans'].set_color('black')
    parts['cbars'].set_color('black')
    parts['cmins'].set_color('black')
    parts['cmaxes'].set_color('black')

    # Add individual points
    for i, data in enumerate([rev_fracs, irrev_fracs]):
        jitter_x = np.random.uniform(-0.05, 0.05, len(data))
        ax.scatter(i + jitter_x, data, c='black', s=20, zorder=3, alpha=0.6)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Reversible\n(IA < 2)', 'Irreversible\n(IA \u2265 2)'], fontsize=9)
    ax.set_ylabel('CGC fraction')

    mean_rev = np.mean(rev_fracs)
    mean_irrev = np.mean(irrev_fracs)
    ax.text(0.5, 0.95, f'p = 0.0003\n({mean_irrev/mean_rev:.1f}\u00d7 enrichment)',
            transform=ax.transAxes, ha='center', va='top', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('D  Binary validation', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig5_authority.pdf')
    fig.savefig(FIGS / 'fig5_authority.png')
    plt.close(fig)
    print("  Fig 5 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 6: Eigenspace Analysis
# ══════════════════════════════════════════════════════════════════════════
def figure6():
    eigen = load_json('eigen20.json')
    tissue = load_json('tissue_specific_eigen.json')
    notch = load_json('notch_validation.json')

    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # Panel A: Variance explained (scree plot)
    ax = fig.add_subplot(gs[0, 0])
    var_exp = eigen['var_explained'][:10]
    cumvar = np.cumsum(var_exp)
    x = range(1, len(var_exp) + 1)
    ax.bar(x, var_exp, color=C_NFL, edgecolor='black', linewidth=0.5, alpha=0.7, label='Individual')
    ax.plot(x, cumvar, 'ro-', markersize=5, lw=1.5, label='Cumulative')
    ax.axhline(y=100/20, color='gray', ls='--', lw=0.8, label='Kaiser threshold')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance explained (%)')
    ax.set_xticks(range(1, 11))
    ax.legend(fontsize=7, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(0.5, 0.5, f'Kaiser = {eigen["kaiser"]} PCs', transform=ax.transAxes,
            fontsize=9, fontweight='bold')
    ax.set_title('A  Eigenspace: 20 modules \u00d7 154 cell types', loc='left', fontweight='bold')

    # Panel B: Tissue-specific PC1 similarity
    ax = fig.add_subplot(gs[0, 1])
    tissues_data = tissue['tissues']
    tissue_names = []
    pc1_sims = []
    for t_name, t_data in sorted(tissues_data.items()):
        tissue_names.append(t_name.replace('_', '\n'))
        # Get PC1 sim_global
        if t_data['pcs']:
            pc1_sims.append(t_data['pcs'][0]['sim_global'])
        else:
            pc1_sims.append(0)

    # Sort by similarity
    order = np.argsort(pc1_sims)[::-1]
    tissue_names = [tissue_names[i] for i in order]
    pc1_sims = [pc1_sims[i] for i in order]

    colors_t = []
    for sim in pc1_sims:
        if sim > 0.85:
            colors_t.append(C_HIGH)
        elif sim > 0.7:
            colors_t.append(C_CGC)
        elif sim > 0.5:
            colors_t.append(C_MED)
        else:
            colors_t.append(C_LOW)

    bars = ax.barh(range(len(tissue_names)), pc1_sims, color=colors_t,
                   edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(tissue_names)))
    ax.set_yticklabels(tissue_names, fontsize=7.5)
    ax.set_xlabel('Cosine similarity to global PC1')
    ax.set_xlim(0, 1.05)
    ax.axvline(x=tissue['pc1_mean_similarity'], color='black', ls='--', lw=1, alpha=0.5)
    ax.text(tissue['pc1_mean_similarity'] + 0.02, len(tissue_names) - 0.5,
            f'mean = {tissue["pc1_mean_similarity"]:.2f}', fontsize=7)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('B  PC1 conservation by tissue', loc='left', fontweight='bold')

    # Panel C: Module pair stability heatmap
    ax = fig.add_subplot(gs[1, 0])
    pairs = tissue['pair_correlations']
    pair_names = list(pairs.keys())
    tissue_keys = sorted(list(pairs[pair_names[0]]['by_tissue'].keys()))

    mat = np.zeros((len(pair_names), len(tissue_keys) + 1))
    for i, pn in enumerate(pair_names):
        mat[i, 0] = pairs[pn]['global']
        for j, tk in enumerate(tissue_keys):
            mat[i, j + 1] = pairs[pn]['by_tissue'].get(tk, 0)

    im = ax.imshow(mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(len(tissue_keys) + 1))
    ax.set_xticklabels(['Global'] + [t[:6] for t in tissue_keys], rotation=45, ha='right', fontsize=6.5)
    ax.set_yticks(range(len(pair_names)))
    pair_labels = [pn.replace('_', ' ~ ') for pn in pair_names]
    ax.set_yticklabels(pair_labels, fontsize=6.5)
    plt.colorbar(im, ax=ax, label='Pearson r', shrink=0.8)
    ax.set_title('C  Module pair stability', loc='left', fontweight='bold')

    # Panel D: Novel cancer gene validation
    ax = fig.add_subplot(gs[1, 1])
    n_val = notch['n_intogen']
    n_total = notch['n_candidates']
    n_not = n_total - n_val

    # Stacked bar
    ax.bar([0], [n_val], color=C_ACCENT, edgecolor='black', linewidth=0.5, label=f'IntOGen validated ({n_val})')
    ax.bar([0], [n_not], bottom=[n_val], color=C_NONCGC, edgecolor='black', linewidth=0.5,
           label=f'Not yet validated ({n_not})')
    ax.set_xlim(-0.8, 0.8)
    ax.set_xticks([0])
    ax.set_xticklabels(['Top 20\ncandidates'], fontsize=9)
    ax.set_ylabel('Number of genes')
    ax.legend(fontsize=8)

    # Add precision
    prec = notch['precision_intogen']
    ax.text(0, n_total + 0.5, f'Precision: {prec:.0%}\n43\u00d7 enrichment',
            ha='center', fontsize=10, fontweight='bold')
    ax.set_ylim(0, n_total + 3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('D  Novel cancer gene candidates', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig6_eigenspace.pdf')
    fig.savefig(FIGS / 'fig6_eigenspace.png')
    plt.close(fig)
    print("  Fig 6 done")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 7: Developmental Trajectories
# ══════════════════════════════════════════════════════════════════════════
def figure7():
    moca = load_json('moca_differentiation.json')

    fig = plt.figure(figsize=(7.5, 7))
    gs = GridSpec(2, 2, hspace=0.4, wspace=0.35)

    # Panel A: E_intra trajectories per lineage
    ax = fig.add_subplot(gs[0, 0])
    trajs = moca['trajectory_stats']
    stages = [9.5, 10.5, 11.5, 12.5, 13.5]

    cancer_colors = {'HIGH': C_HIGH, 'MEDIUM': C_MED, 'LOW': C_LOW, 'UNKNOWN': C_NONCGC}
    for traj_name, traj_data in sorted(trajs.items()):
        means = traj_data['E_intra_means']
        cancer = traj_data['cancer_relevance']
        color = cancer_colors.get(cancer, C_NONCGC)
        alpha = 1.0 if cancer == 'HIGH' else 0.5
        lw = 2.0 if cancer == 'HIGH' else 0.8
        days = stages[:len(means)]
        ax.plot(days, means, '-o', color=color, alpha=alpha, lw=lw, markersize=3,
                label=traj_name if cancer == 'HIGH' else None)

    ax.set_xlabel('Embryonic day')
    ax.set_ylabel('Mean E_intra (Shannon entropy)')
    ax.legend(fontsize=7, loc='upper right', title='HIGH cancer risk', title_fontsize=7)
    # Add legend patches for cancer relevance
    leg2 = ax.legend(handles=[
        mpatches.Patch(color=C_HIGH, label='HIGH'),
        mpatches.Patch(color=C_MED, label='MEDIUM'),
        mpatches.Patch(color=C_LOW, label='LOW'),
    ], fontsize=6.5, loc='lower left', title='Cancer risk', title_fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('A  E_intra decreases during differentiation', loc='left', fontweight='bold')

    # Panel B: E_inter trajectories (from moca_group_entropy)
    ax = fig.add_subplot(gs[0, 1])
    # Load the source data
    entropy_data = []
    moca_file = Path("/Users/teo/Desktop/research/oscilatory/results/embryogenesis/analysis/moca_group_entropy.csv")
    if moca_file.exists():
        with open(moca_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                entropy_data.append(row)

    # Group E_inter by trajectory x day
    traj_inter = defaultdict(lambda: defaultdict(list))
    for row in entropy_data:
        traj_inter[row['trajectory']][float(row['day'])].append(float(row['E_inter']))

    cancer_rel = {
        'Hematopoietic': 'HIGH', 'Epithelial': 'HIGH', 'Hepatic': 'HIGH',
        'Neural': 'MEDIUM', 'Neural_crest': 'MEDIUM', 'Mesenchyme': 'MEDIUM', 'Kidney': 'MEDIUM',
        'Endothelial': 'LOW', 'Eye': 'LOW', 'Lens': 'LOW', 'Notochord': 'LOW',
        'Intermediate_mesoderm': 'LOW', 'Other': 'UNKNOWN',
    }

    for traj_name in sorted(traj_inter.keys()):
        days = sorted(traj_inter[traj_name].keys())
        means = [np.mean(traj_inter[traj_name][d]) for d in days]
        cancer = cancer_rel.get(traj_name, 'UNKNOWN')
        color = cancer_colors.get(cancer, C_NONCGC)
        alpha = 1.0 if cancer == 'HIGH' else 0.5
        lw = 2.0 if cancer == 'HIGH' else 0.8
        ax.plot(days, means, '-s', color=color, alpha=alpha, lw=lw, markersize=3)

    ax.set_xlabel('Embryonic day')
    ax.set_ylabel('Mean E_inter (JSD)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('B  E_inter increases (population diversity)', loc='left', fontweight='bold')

    # Panel C: IA score per trajectory
    ax = fig.add_subplot(gs[1, 0])
    eigen_pos = moca['eigenspace_positions']
    traj_names_sorted = sorted(eigen_pos.keys(), key=lambda t: eigen_pos[t].get('ia_score', 0), reverse=True)
    ia_scores_plot = [eigen_pos[t]['ia_score'] for t in traj_names_sorted]
    cancer_labels = [eigen_pos[t]['cancer_relevance'] for t in traj_names_sorted]
    colors_ia = [cancer_colors.get(c, C_NONCGC) for c in cancer_labels]

    # Filter out trajectories with no data
    valid = [(t, s, c) for t, s, c in zip(traj_names_sorted, ia_scores_plot, colors_ia) if s > 0]
    if valid:
        t_v, s_v, c_v = zip(*valid)
        bars = ax.barh(range(len(t_v)), s_v, color=c_v, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(t_v)))
        ax.set_yticklabels([t.replace('_', ' ') for t in t_v], fontsize=7)
        ax.set_xlabel('IA-weighted module activity')
        ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('C  Adult eigenspace: IA score per lineage', loc='left', fontweight='bold')

    # Panel D: Conceptual model
    ax = fig.add_subplot(gs[1, 1])
    ax.axis('off')
    # Draw arrows showing development vs cancer direction
    # Development: high E_intra -> low E_intra (left to right)
    # Cancer: low IA -> high IA (right to left)

    ax.annotate('', xy=(0.85, 0.75), xytext=(0.15, 0.75),
               arrowprops=dict(arrowstyle='->', color=C_NFL, lw=3))
    ax.text(0.5, 0.82, 'DIFFERENTIATION', ha='center', fontsize=10,
            fontweight='bold', color=C_NFL)
    ax.text(0.12, 0.68, 'Pluripotent\n(high E_intra)', ha='center', fontsize=7)
    ax.text(0.88, 0.68, 'Specialized\n(low E_intra)', ha='center', fontsize=7)

    ax.annotate('', xy=(0.15, 0.35), xytext=(0.85, 0.35),
               arrowprops=dict(arrowstyle='->', color=C_NON, lw=3))
    ax.text(0.5, 0.42, 'CANCER (de-differentiation)', ha='center', fontsize=10,
            fontweight='bold', color=C_NON)
    ax.text(0.88, 0.28, 'Normal cell\n(low IA active)', ha='center', fontsize=7)
    ax.text(0.12, 0.28, 'Cancer cell\n(high IA reactivated)', ha='center', fontsize=7)

    ax.text(0.5, 0.08, 'Development compresses state space\nCancer re-expands it via high-IA modules',
            ha='center', fontsize=8, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.95)
    ax.set_title('D  Model: opposing trajectories', loc='left', fontweight='bold')

    fig.savefig(FIGS / 'fig7_development.pdf')
    fig.savefig(FIGS / 'fig7_development.png')
    plt.close(fig)
    print("  Fig 7 done")


# ══════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating Paper 2 figures...")
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()
    print(f"\nAll figures saved to {FIGS}/")
    print("Files:", list(FIGS.glob("*.pdf")))
