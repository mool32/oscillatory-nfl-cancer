[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![License: CC-BY 4.0](https://img.shields.io/badge/Data%20%26%20Manuscript-CC--BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Preprint](https://img.shields.io/badge/Preprint-Manuscript-orange)](paper/manuscript.md)

# Negative feedback loop architecture as a modular predictor of cancer vulnerability across signaling pathways

**128 negative feedback loops algorithmically extracted from 159 KEGG signaling networks; NFL genes show 59-fold CGC enrichment over non-NFL genes in the same pathways**

Theodor Spiro | [ORCID 0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) | tspiro@vaika.org

📄 **Manuscript:** [`paper/manuscript.md`](paper/manuscript.md) (Markdown source) · [`paper/Spiro_2026_Paper2_preprint.docx`](paper/Spiro_2026_Paper2_preprint.docx) (Word preprint)
🧮 **Pipeline entry-points:** [`step02_kegg_loop_extraction_v3.py`](step02_kegg_loop_extraction_v3.py) (NFL extraction) → [`step06_cgc_enrichment.py`](step06_cgc_enrichment.py) (cancer enrichment) → [`step22_irreversible_authority.py`](step22_irreversible_authority.py) (IA metric)
🧬 **Companion paper:** [Temporal architecture of signaling oscillations predicts cancer gene function](https://github.com/mool32/oscillatory-cancer-framework) (Spiro 2026) — independent temporal-classification approach that this work converges with algorithmically

---

## Brief Summary

Whether membership in negative feedback loops (NFLs) — the topological basis of signaling oscillations — predicts cancer gene status independently of pathway identity has not been tested. We algorithmically extract all short NFLs from 159 KEGG signaling networks and demonstrate:

1. **128 unique NFL motifs organize into 14 pathway-level oscillatory modules.** Algorithmic extraction yields 228 raw motifs (21 two-node + 207 three-node), reducing via alias resolution + Jaccard clustering to 128 unique → 59 modules → 14 pathway groups (NF-κB, p53, ERK/MAPK, Wnt, Notch, Circadian, Hippo, mTOR, JAK-STAT, TGF-β, NRF2, Hedgehog, Calcium, Cell Cycle). Biological-validation rate on a random 20-motif sample: 80% (35% full + 45% partial support).
2. **NFL genes are 59-fold enriched for Cancer Gene Census membership** vs non-NFL genes in the same pathways (Fisher exact *p* = 9 × 10⁻⁴⁴). The enrichment converges with the independent temporal-classification approach of the [companion paper](https://github.com/mool32/oscillatory-cancer-framework) (12.8× vs 12.1×); 12 cancer genes are discovered exclusively through NFL extraction.
3. **Data-driven community detection confirms real modularity** (*Q* = 0.39, ARI = 0.48 with KEGG annotations) — the 14 pathway-level groups are not an artifact of the curation but emerge from the topology itself.
4. **Conservation gradient across 15 species** maps universal eukaryotic NFLs (ERK/DUSP, Ca²⁺/SERCA; ~1500 Mya) to vertebrate-specific circuits (Notch, JAK-STAT; ~435 Mya). All ancient loops exhibit asymmetric waveforms.
5. **Irreversible Authority (IA) predicts cancer gene fraction across modules with Spearman ρ = 0.83** (robust under leave-one-out and reclassification perturbation). IA — the count of irreversible cell-fate decisions a module controls — explains the non-uniform distribution of cancer genes across oscillatory modules: high-IA modules (Cell Cycle, p53) are densely cancer-enriched; low-IA modules (Circadian) are not.
6. **Eigenspace projection of module expression across 154 human cell types identifies 20 novel cancer gene candidates, 13 (65%) validated by IntOGen.** Tissue-specific eigendecomposition shows the primary signaling axis is universal (mean cosine similarity 0.70 across tissues) while higher axes are tissue-dependent — explaining why the same module produces different cancer profiles in different tissues.

These results establish feedback-loop topology as a specific, pathway-independent predictor of cancer vulnerability and suggest oscillatory modules as fundamental evolutionary units of cellular regulation.

## Inputs

| Source | Use | Access |
|---|---|---|
| **KEGG signaling pathways** (KGML) | 159 pathway networks → signed DiGraphs → cycle enumeration | [KEGG REST API](https://www.kegg.jp/kegg/rest/) (academic non-commercial) |
| **COSMIC Cancer Gene Census** v103 | Cancer gene annotations (740 genes) | [cancer.sanger.ac.uk/census](https://cancer.sanger.ac.uk/census) (registration required) |
| **Tabula Sapiens** + literature cell-type expression | Module eigenspace across 154 cell types | Public single-cell atlases |
| **Ensembl Compara** orthologs across 15 species | Evolutionary dating (Mya) per loop | [ensembl.org/info/genome/compara](https://www.ensembl.org/info/genome/compara/) |
| **IntOGen** | Independent validation of 20 novel cancer-gene candidates | [intogen.org](https://www.intogen.org/) |

## Repository structure

```
├── paper/
│   ├── manuscript.md                  # Manuscript source (Markdown)
│   └── Spiro_2026_Paper2_preprint.docx  # Word preprint version
├── step01_loop_components.py          # Pipeline scripts in numbered execution order
├── step02_kegg_loop_extraction*.py    # (3 versions kept for provenance; v3 is canonical)
├── step03_ortholog_mapping.py
├── step04_dating_and_age_vs_A.py
├── step05_nfl_clustering.py
├── step06_cgc_enrichment.py           # MAIN RESULT: 59-fold enrichment, p=9e-44
├── step07_cgc_corrected.py            # FDR/multiple-comparison-corrected version
├── step08_fisher_age_class.py
├── step09_celltype_expression.py
├── step10_rds_computation.py
├── step11_rds_polarity_corrected.py
├── step12_cvs_two_axis.py
├── step13_cognitive_load.py
├── step14_receptor_module_coupling.py
├── step15_confound_checks.py
├── step16_eigendecomposition.py
├── step17_new_modules.py
├── step18_convergence_metabolic.py
├── step19_eigen20.py                  # 20-module eigenspace
├── step20_nfl_cgc_new_modules.py
├── step21_vulnerability_metric.py
├── step22_irreversible_authority.py   # IA metric (Spearman rho=0.83 with cancer fraction)
├── step23_eigenspace_cancer.py        # Cancer-region in eigenspace; 20 novel candidates
├── step24_maintenance_cost.py
├── step25_notch_validation.py
├── step26_redundancy_eigenmode.py
├── step27_mutation_trajectories.py
├── step28_tissue_specific_eigen.py    # Tissue-universal vs tissue-specific axes
├── step29_moca_differentiation.py
├── step30_ia_sensitivity.py
├── step31_age_A_sensitivity.py
├── generate_figures.py                # Renders figures/ from data/
├── generate_supplementary.py          # Renders supplementary/ tables
├── data/                              # 45 intermediate JSON/CSV files (every paper number)
├── figures/                           # 7 publication figures (PDF + PNG, 300 DPI)
├── supplementary/                     # 11 supplementary tables (S1-S11)
├── README.md
└── LICENSE
```

### Figures

| File | Topic |
|---|---|
| `figures/fig1_pipeline.pdf/.png` | NFL extraction pipeline + biological validation |
| `figures/fig2_enrichment.pdf/.png` | **Main result:** cancer gene enrichment (3 nested tests) |
| `figures/fig3_modularity.pdf/.png` | Community detection vs KEGG annotations |
| `figures/fig4_evolution.pdf/.png` | Conservation gradient across 15 species |
| `figures/fig5_authority.pdf/.png` | IA metric vs cancer gene fraction (rho=0.83) |
| `figures/fig6_eigenspace.pdf/.png` | Module eigenspace + 20 novel candidates |
| `figures/fig7_development.pdf/.png` | MOCA developmental differentiation analysis |

### Supplementary tables

11 tables (`supplementary/table_s1_*.csv` through `table_s11_*.csv`): core loops, ortholog matrix, unique motifs, modules, community detection, biological validation, KEGG gene summary, IA values, cancer candidates, tissue-specific eigen, MOCA trajectories.

## Reproducing the analysis

### Environment

```bash
git clone https://github.com/mool32/oscillatory-nfl-cancer.git
cd oscillatory-nfl-cancer
pip install numpy pandas scipy networkx matplotlib python-louvain
```

### Run

Scripts execute in numbered order; each is self-contained. The repository commits **all intermediate `data/` and `supplementary/` files**, so paper numbers reproduce without re-querying KEGG/CGC/Ensembl.

```bash
# Core pipeline (NFL extraction through main result)
python step01_loop_components.py
python step02_kegg_loop_extraction_v3.py
python step03_ortholog_mapping.py
python step04_dating_and_age_vs_A.py
python step05_nfl_clustering.py
python step06_cgc_enrichment.py            # main enrichment result
python step07_cgc_corrected.py

# Architectural analyses (IA, eigenspace, tissue-specific axes)
python step22_irreversible_authority.py
python step23_eigenspace_cancer.py
python step28_tissue_specific_eigen.py

# Render figures and supplementary tables
python generate_figures.py
python generate_supplementary.py
```

To re-run from KEGG/CGC/Ensembl directly, you will need API keys / academic licences (see `data/` schema and `step02_kegg_loop_extraction_v3.py` for source URLs).

## Citation

```bibtex
@article{spiro2026nflcancer,
  author  = {Spiro, Theodor},
  title   = {Negative feedback loop architecture as a modular predictor of cancer vulnerability across signaling pathways},
  journal = {Preprint},
  year    = {2026},
  note    = {Manuscript in preparation. https://github.com/mool32/oscillatory-nfl-cancer}
}
```

## Contact

Theodor Spiro — tspiro@vaika.org

## License

- **Code** (`step*.py`, `generate_*.py`): MIT (see [LICENSE](LICENSE))
- **Data** (`data/*`, `supplementary/*`): CC-BY 4.0, with KEGG and COSMIC CGC attribution requirements honored per their respective licences
- **Figures** (`figures/*`): CC-BY 4.0
- **Manuscript** (`paper/manuscript.md`, `paper/*.docx`): CC-BY 4.0
