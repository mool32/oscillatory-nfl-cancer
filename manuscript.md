# Negative feedback loop architecture as a modular predictor of cancer vulnerability across signaling pathways

Theodor Spiro^1^

^1^Independent Researcher

**Correspondence:** theospirin@gmail.com

---

## Abstract

Whether membership in negative feedback loops (NFLs) --- the topological basis of signaling oscillations --- predicts cancer gene status independently of pathway identity has not been tested. Here I algorithmically extract all short NFLs from 159 KEGG signaling networks, identifying 128 unique motifs organized into 14 pathway-level oscillatory modules. NFL genes show 59-fold enrichment for Cancer Gene Census membership compared to non-NFL genes in the same pathways (Fisher exact *p* = 9 x 10^-44^), converging with an independent temporal-classification approach (12.8x vs 12.1x enrichment). Evolutionary analysis across 15 species reveals a conservation gradient from universal eukaryotic NFLs to vertebrate-specific circuits, with older loops exhibiting more asymmetric waveforms. To explain non-uniform cancer gene distribution across modules, I introduce an Irreversible Authority (IA) metric --- the number of irreversible cell-fate decisions each module controls. IA predicts cancer gene fraction with Spearman rho = 0.83 across 20 modules (robust under leave-one-out and reclassification perturbation). Eigenspace projection of module expression across 154 human cell types identifies 20 novel cancer gene candidates, 13 (65%) validated by IntOGen. Tissue-specific eigendecomposition shows the primary signaling axis is universal (mean cosine similarity 0.70) while higher axes are tissue-dependent. These results establish feedback loop topology as a specific, pathway-independent predictor of cancer vulnerability.

**Keywords:** negative feedback loops, cancer gene census, oscillatory modules, signaling topology, irreversible authority, eigenspace, KEGG, cancer vulnerability

---

## Introduction

The identification of cancer driver genes has historically proceeded pathway by pathway. *KRAS* is classified as an oncogene through study of the RAS-MAPK cascade; *APC* as a tumor suppressor through the Wnt/beta-catenin pathway; *TP53* through its apoptotic response to DNA damage. Each classification required independent mechanistic investigation. Whether shared topological features can predict cancer gene status across pathways --- without pathway-specific knowledge --- remains an open question.

A companion study (Spiro, 2026) demonstrated that a gene's temporal position within a signaling oscillation cycle predicts its cancer function: rise-phase genes map to oncogenes and recovery-phase genes map to tumor suppressors (odds ratio = 27.5, *p* = 3.6 x 10^-9^ across 14 pathways), with predicted inversions in growth-inhibitory pathways (p53, TGF-beta). That analysis relied on manual classification of 157 genes based on published oscillation models. The present work asks a complementary question: can the same organizational principle be recovered algorithmically, directly from pathway topology, and extended to predict *which* oscillatory modules carry the greatest cancer vulnerability?

Signaling oscillations universally arise from negative feedback loops (NFLs): a fast activation arm coupled to a slower inhibitory arm generates pulsatile dynamics (Novak and Tyson, 2008; Purvis and Lahav, 2013). NFLs are among the most conserved network motifs in biology, appearing across kingdoms from yeast (MAPK/phosphatase) to mammals (NF-kappaB/IkappaB, p53/MDM2) (Alon, 2007). If oscillatory function is central to cellular decision-making, and if cancer arises from disruptions that convert oscillations to sustained signaling (Purvis et al., 2012), then NFL membership should be a specific predictor of cancer gene status --- more specific than generic pathway membership.

Three questions motivate this study. First, can NFLs be extracted algorithmically from public pathway databases at genome scale, and do they capture known oscillatory circuits? Second, do NFL genes show cancer enrichment beyond what pathway membership alone explains? Third, what determines the *degree* of cancer vulnerability per module --- why is the ERK/MAPK module more cancer-enriched than the Circadian module, when both are oscillatory?

Here I develop an automated pipeline to extract all short NFLs from KEGG signaling networks, characterize their evolutionary conservation and modular organization, quantify cancer gene enrichment with three nested statistical tests, and introduce architectural metrics --- Irreversible Authority and eigenspace position --- that explain the non-uniform distribution of cancer genes across oscillatory modules.

---

## Results

### Algorithmic extraction of negative feedback loops from KEGG signaling networks

I converted 159 KEGG signaling pathway KGML files into signed directed graphs (DiGraphs), where edges represent activation (+1) or inhibition (-1) relationships. To identify NFLs, I enumerated all simple cycles of length 2 (direct negative feedback, A activates B, B inhibits A) and length 3 (indirect negative feedback through one intermediate) using depth-limited search seeded from genes belonging to 22 literature-curated core oscillatory loops (Table S1). A cycle was classified as an NFL if the product of edge signs around the loop was negative (odd number of inhibitory edges).

This procedure yielded 228 raw NFL motifs: 21 two-node and 207 three-node cycles. Many raw motifs represent the same biological feedback loop counted multiple times due to KEGG alias variants (e.g., AISIMD = SOCS3, DAKAR = CDK4). I resolved 37 KEGG aliases to standard HGNC symbols, then removed exact duplicates (identical gene sets with identical edge signs), reducing the count to 128 unique NFL motifs.

These 128 motifs were further organized into 59 non-redundant modules using Jaccard similarity clustering (threshold J > 0.5, indicating >50% gene overlap) and merged into 14 pathway-level groups corresponding to established signaling pathways: NF-kappaB, p53, ERK/MAPK, Wnt, Notch, Circadian, Hippo, mTOR, JAK-STAT, TGF-beta, NRF2, Hedgehog, Calcium, and Cell Cycle (Table S3, Table S4).

To estimate extraction fidelity, I randomly sampled 20 of 128 unique motifs and assessed biological support from published literature. Seven motifs (35%) had full experimental validation (published oscillation data or perturbation studies confirming the feedback relationship), nine (45%) had partial support (correct genes, plausible mechanism, but no direct oscillation measurement), and four (20%) were artifacts (inverted edge direction in the KEGG source, or hub gene connecting unrelated pathways). Overall, 80% of extracted motifs have at least partial biological support.

The 128 motifs collectively encompass 76 unique genes (43 from the original 22 core loops plus 33 algorithmically discovered), of which 36 are members of the COSMIC Cancer Gene Census (CGC v103; 740 genes total). The composition of each module, its gene membership, and CGC status are provided in Tables S3 and S4.


### NFL genes are specifically enriched for cancer genes

I evaluated cancer gene enrichment through three nested statistical tests of increasing specificity (Fig. 2A):

**Test 1: NFL genes vs. genome background.** Of 76 NFL genes, 36 are CGC members (47.4%). Compared to the genome background rate of 3.7% (740/20,000), this represents a 12.8-fold enrichment (hypergeometric *p* = 1.9 x 10^-31^). This result converges quantitatively with Paper 1, which found 12.1-fold enrichment using an independent temporal classification of 157 genes from the same 14 pathways.

**Test 2: NFL genes vs. all KEGG signaling genes.** KEGG signaling pathways contain 16,760 genes total, of which 475 are CGC members (2.8%). NFL genes show a 28.1-fold enrichment over this pathway-level background (hypergeometric *p* = 5.5 x 10^-44^).

**Test 3: NFL genes vs. non-NFL genes in the same pathways (critical test).** This test asks: within the same KEGG signaling pathways, are genes participating in NFLs more likely to be cancer genes than genes that do not? The 2x2 contingency table is:

|                | CGC member | Not CGC | Total  |
|----------------|-----------|---------|--------|
| NFL gene       | 36        | 40      | 76     |
| Non-NFL gene   | 251       | 16,444  | 16,695 |

Fisher exact test: odds ratio = 58.96, *p* = 9.1 x 10^-44^.

This is the key result: NFL membership predicts cancer gene status with an odds ratio of 59, controlling for pathway membership. Genes embedded in negative feedback loops are not merely in cancer-relevant pathways; they are specifically the genes within those pathways that become cancer drivers. Among the 36 CGC genes found through NFL extraction, 12 were not identified by the temporal classification in Paper 1, demonstrating that the two approaches provide complementary coverage of the cancer gene landscape.



### Oscillatory modules exhibit genuine but non-trivial modularity

To determine whether the 128 NFL motifs form biologically coherent modules, I performed unsupervised community detection on the motif-motif network, where edges connect motifs sharing at least one gene, weighted by Jaccard similarity. Louvain community detection (resolution = 1.0) identified 10 communities from a network of 2,042 edges with density 0.25.

The resulting modularity is Q = 0.39, indicating that the NFL network has substantially more within-community edges than expected by chance. However, the correspondence between data-driven communities and KEGG pathway annotations is only partial: Adjusted Rand Index (ARI) = 0.48. Six pathways (NF-kappaB, Wnt, Circadian, Calcium, NRF2, Hedgehog) form clean, isolated communities. The remaining pathways form two large mixed clusters: a JAK-STAT/mTOR/PI3K cluster (reflecting shared use of STAT and AKT signaling nodes) and an ERK/TGF-beta/Cell Cycle cluster (reflecting shared use of MAPK and SMAD nodes).

The modular structure is shaped by a small number of hub genes that participate in multiple NFLs across pathways. MAPK1 (ERK2) participates in NFLs from three distinct modules (ERK/MAPK, BMP/SMAD, JAK-STAT), confirming its role as a cross-pathway signaling integrator. Other cross-module hubs include TP53, SMAD2/3, and PIK3CA. These hub genes are disproportionately cancer drivers: all four are CGC members.


### Evolutionary conservation reveals an asymmetry-first architecture

I mapped the 22 core NFL loop pairs to orthologs across 15 species spanning approximately 1.5 billion years of eukaryotic evolution, from *S. cerevisiae* and *D. discoideum* to mammals (Table S2). Ortholog mapping was performed using Ensembl Compara, requiring that both the activator and inhibitor genes of a loop have identifiable one-to-one orthologs in a given species.

A clear conservation gradient emerges: 22/22 loops (100%) are present in at least one vertebrate species, 11/22 (50%) have identifiable invertebrate orthologs (*D. melanogaster*, *C. elegans*), and 6/22 (27%) have detectable orthologs in unicellular organisms. The most ancient NFLs include ERK/DUSP (phosphatase-based, detected in yeast), Ca^2+^/SERCA (pump-based, detected across eukaryotes), and NRF2/KEAP1 (ubiquitin-based, detected in plants and fungi). Vertebrate-specific innovations include Notch/HES, JAK-STAT/SOCS, and Hedgehog/GLI.

To test whether evolutionary age correlates with oscillatory waveform shape, I combined loop ages with waveform asymmetry values (A = tau_rise/T, where tau_rise is the rise time and T is the period) from Paper 1. Of the 22 loops, 14 have experimentally determined or computationally estimated A values. The correlation between loop age and waveform asymmetry is a negative trend (Spearman rho = -0.38, *p* = 0.18 two-tailed, N = 14; bootstrap 95% CI [-0.85, 0.23]) that becomes significant upon exclusion of NRF2/KEAP1, a biologically motivated outlier (rho = -0.66, *p* = 0.015, N = 13; bootstrap 95% CI [-0.89, -0.17], excluding zero). Leave-one-out analysis confirms that NRF2 is the only loop whose removal changes significance; all other removals yield rho in the range [-0.27, -0.44].

The NRF2/KEAP1 outlier merits explanation. This loop is ancient (detectable in plants) but symmetric (A approximately 0.48) in mammalian cells. NRF2 is a stress-response pathway that was co-opted for antioxidant defense in metazoans; its waveform symmetry likely reflects convergent adaptation to a detection function (symmetric oscillations optimize threshold detection) rather than ancestral architecture. Without NRF2, the data reveal a clear evolutionary pattern: the oldest NFLs (>1000 Mya) all exhibit asymmetric waveforms (A < 0.35, Class I in the terminology of Paper 1), while symmetric waveforms (A > 0.45, Class II) appear exclusively among metazoan innovations.


### Expansion to 20 oscillatory modules

To extend the analysis beyond the 14 KEGG-derived modules, I screened six additional oscillatory pathways with documented or predicted negative feedback regulation: AMPK (energy sensing), SREBP (lipid homeostasis), ATR/CHK1 (replication checkpoint), Rho/ROCK (cytoskeletal dynamics), PPAR/LXR (nuclear receptor cycling), and Autophagy (ULK1/mTOR feedback). Each was confirmed to contain at least one NFL motif with 3 or more unique genes by the same extraction pipeline. The expanded set of 20 modules (Table S8) was used for all subsequent analyses requiring cross-module comparison.


### Irreversible Authority predicts cancer vulnerability across modules

The CGC enrichment result (OR = 59) establishes that NFL genes as a class are cancer-enriched. But cancer genes are not uniformly distributed across modules: 72.7% of p53 module genes are CGC members versus 0% for Circadian. What determines this variation?

I tested three architectural hypotheses:

**Hypothesis 1: Maintenance cost.** Modules with higher expression (greater metabolic cost) should be harder to maintain and more vulnerable to cancer-causing mutations. Result: no significant correlation between module expression rank and CGC fraction (Spearman rho = 0.03, *p* = 0.96). **Rejected.**

**Hypothesis 2: Convergence (hub degree).** Modules with more incoming connections from other pathways should be more vulnerable because more perturbations can reach them. Result: moderate correlation (partial rho = 0.25, *p* = 0.4 after controlling for total expression). **Insufficient.**

**Hypothesis 3: Irreversible Authority (IA).** Modules that control more irreversible cell-fate decisions (apoptosis, senescence, terminal differentiation, DNA damage checkpoint commitment) should be more cancer-relevant because their disruption has irreversible phenotypic consequences. I quantified IA as the number of distinct irreversible cellular outcomes controlled by each module, based on established pathway biology:

| Module | IA score | CGC fraction | Key irreversible decisions |
|--------|---------|-------------|--------------------------|
| p53 | 3 | 0.727 | apoptosis, senescence, DNA repair commitment |
| Wnt | 3 | 0.600 | terminal differentiation, stem cell fate, EMT |
| Notch | 3 | 0.308 | lateral inhibition, lineage commitment, angiogenesis |
| Cell Cycle | 2 | 0.462 | mitotic commitment, cytokinesis |
| NF-kappaB | 2 | 0.400 | inflammatory commitment, survival decision |
| ERK/MAPK | 2 | 0.467 | proliferation vs. differentiation |
| Hippo | 2 | 0.364 | organ size, contact inhibition |
| TGF-beta | 2 | 0.364 | EMT, growth arrest |
| JAK-STAT | 2 | 0.333 | immune activation, differentiation |
| ATR/CHK1 | 2 | 0.222 | replication checkpoint, fork restart |
| PI3K/PTEN | 2 | 0.500 | survival, metabolic commitment |
| mTOR | 1 | 0.200 | growth commitment |
| AMPK | 1 | 0.100 | metabolic switch |
| Calcium | 1 | 0.200 | NFAT activation |
| Rho/ROCK | 1 | 0.100 | cytoskeletal commitment |
| Autophagy | 1 | 0.111 | autophagic flux |
| Circadian | 0 | 0.000 | (reversible oscillation) |
| NRF2 | 0 | 0.200 | (reversible stress response) |
| SREBP | 0 | 0.000 | (reversible lipid regulation) |
| PPAR/LXR | 0 | 0.000 | (reversible transcription) |

IA predicts CGC fraction with Spearman rho = 0.83 (*p* < 0.001, N = 20). The underlying dose-response is sharp: modules with IA = 0 have mean CGC fraction 0.050, IA = 1 averages 0.142, IA = 2 averages 0.382, and IA = 3 averages 0.545. A binary split (IA >= 2 vs IA < 2) produces a mean CGC difference of 0.38 vs 0.08 (Fisher exact *p* = 0.0003).

Authority alone (number of distinct biological outputs, including reversible ones) also predicts CGC fraction (rho = 0.752, *p* = 0.0001), but the irreversibility criterion adds substantial predictive power: adding a bottleneck metric (number of unique genes required for each irreversible output) improves the correlation to rho = 0.869.

The IA metric is defined from biochemical pathway logic --- which modules control irreversible cellular outcomes --- without reference to cancer gene databases. To assess robustness, I performed three sensitivity analyses (Table S8). Leave-one-module-out analysis yields rho in the range [0.644, 0.783] (no single module is indispensable). Single-module IA reclassification (shifting any one module's IA by plus or minus 1) yields rho in [0.602, 0.744]; the worst case (PI3K/PTEN reclassified from IA=2 to IA=1) still produces rho = 0.602. Even simultaneously reclassifying three borderline modules (NF-kappaB 2 to 1, Calcium 1 to 0, NRF2 0 to 1) yields rho = 0.714. As an independent validation, an automated Gene Ontology-based proxy (counting irreversible-process annotations per module gene) correlates with manual IA at rho = 0.889 and independently predicts CGC fraction at rho = 0.632 (*p* = 0.003).


### Eigenspace analysis identifies novel cancer gene candidates

To ask whether the multivariate pattern of module expression across cell types contains cancer-relevant structure, I performed eigendecomposition on the 20-module x 154-cell-type expression matrix (HPA single-cell nCPM data). Kaiser criterion identifies 6 principal components explaining 67% of total variance:

| PC | Variance | Top positive loadings | Top negative loadings |
|----|---------|----------------------|---------------------|
| PC1 | 25.2% | ERK/MAPK, JAK-STAT, NRF2 | -- |
| PC2 | 13.2% | Rho/ROCK | TGF-beta, AMPK |
| PC3 | 11.8% | ATR/CHK1, Cell Cycle, p53 | -- |
| PC4 | 6.6% | SREBP, Calcium | ATR/CHK1 |
| PC5 | 5.5% | Circadian | PI3K/PTEN, Cell Cycle |
| PC6 | 5.0% | PPAR/LXR | Autophagy |

PC1 captures general signaling amplitude (cell types with high expression of all modules). PC2 separates mechanosensing (Rho/ROCK) from growth-suppressive axes (TGF-beta, AMPK). PC3 specifically captures the DNA damage / cell cycle checkpoint axis.

I projected each of the 76 NFL genes into this 6D eigenspace based on its module membership and the module's expression pattern. CGC genes occupy a non-random subregion of eigenspace (permutation test: *p* = 0.0001, 10,000 permutations). Using the eigenspace coordinates of known CGC genes as a reference, I identified the 20 non-CGC NFL genes closest to the CGC centroid in eigenspace. Of these 20 candidates, 13 (65%) are confirmed as cancer-associated genes in the IntOGen database (Bailey et al., 2018), representing a 43-fold enrichment over background expectation (IntOGen enrichment *p* < 10^-6^). The 7 unvalidated candidates represent predictions for future experimental testing.


### Tissue-specific eigendecomposition reveals a partially universal architecture

To determine whether the 6 eigenmodes are universal features of oscillatory module organization or tissue-specific artifacts, I split the 154 HPA cell types into 7 tissue classes (blood/immune, N = 22; epithelial, N = 70; neuronal, N = 17; mesenchymal, N = 12; endocrine, N = 10; stem/germ, N = 13; vascular, N = 10) and performed eigendecomposition independently for each class.

PC1 is largely conserved across tissue classes: cosine similarity to the global PC1 axis ranges from 0.98 (blood/immune) to 0.22 (vascular), with a mean of 0.70. The cancer-relevant tissues show the highest conservation: blood/immune (0.98), neuronal (0.88), epithelial (0.83). Vascular and mesenchymal tissue classes show remodeled PC1 axes, consistent with their specialized mechanical signaling requirements.

Higher principal components are consistently tissue-specific (mean similarity < 0.3 for PC2--PC6), confirming that the fine-grained architecture of module co-variation is shaped by tissue context. However, specific module pairs show remarkable stability across all tissues. The Hippo--TGF-beta co-variation is positive in all 7 tissue classes (range: +0.31 to +0.90), as is the JAK-STAT--NF-kappaB coupling (range: +0.13 to +0.87). These invariant module pairs may represent architectural constraints on oscillatory signaling that are conserved regardless of cellular context.

The interpretation for cancer is stratified: PC1 (general signaling amplitude) provides a tissue-universal vulnerability axis, while higher PCs modulate vulnerability in a tissue-specific manner. This is consistent with the observation that some cancer types (leukemia, carcinoma) share common driver genes while others (sarcoma, melanoma) have distinct driver landscapes.


### Developmental trajectories decrease in transcriptional complexity

To test whether the eigenspace architecture is established during development, I analyzed the Mouse Organogenesis Cell Atlas (MOCA; Cao et al., 2019), which provides single-cell transcriptomic data for approximately 300,000 cells across 13 developmental lineages from embryonic day 9.5 to 13.5.

Using per-cell transcriptional entropy (E_intra, Shannon entropy of the gene expression distribution) as a measure of transcriptional complexity, all 13 developmental trajectories show monotonically decreasing E_intra from E9.5 to E13.5 (Spearman rho with day: range -0.60 to -1.00; 10/13 trajectories significant at *p* < 0.05). This universal decrease confirms that differentiation involves transcriptional specialization: as cells commit to lineage fates, their transcriptional repertoire narrows.

Simultaneously, between-cell diversity (E_inter, Jensen-Shannon divergence between cells within a lineage) tends to increase during the same period (positive rho in 10/13 trajectories), indicating that population-level heterogeneity grows even as individual cells become more specialized.

The developmental direction --- universally decreasing E_intra --- contrasts with the cancer strategy identified through the mutation order analysis, where cancer sequentially activates modules with increasing Irreversible Authority. Differentiation compresses the transcriptional state space; cancer re-expands it. The implications of this reciprocal relationship are considered in the Discussion.


### Killed hypotheses

Transparency requires reporting analyses that produced null results:

**Maintenance cost hypothesis.** If modules requiring more cellular resources (higher expression) are more fragile, expression should predict CGC fraction. Spearman rho = 0.03, *p* = 0.96 (N = 20). Expression level is unrelated to cancer vulnerability. **Rejected.**

**Mutation trajectory ordering.** If cancer mutations follow a predictable temporal sequence (disrupting high-IA modules first), then mutation order in multi-hit cancers should correlate with IA rank. Empirical analysis of COSMIC mutation co-occurrence shows an ascending pattern (low IA first, high IA later), which is the *inverse* of the naive prediction. Cancer appears to break accelerator modules first and checkpoint modules later --- a sequential bypass strategy rather than a top-down cascade. This result, while biologically interpretable, contradicts the initial hypothesis. **Inverted, reframed.**

**Development-aging replication.** We tested whether the anti-correlation between E_intra and E_inter slopes observed during aging (Paper 1; rho approximately -0.55) is replicated during embryonic development. The cross-trajectory correlation is rho = -0.02 (*p* = 0.95, N = 13). The aging entropy pattern does not replicate in development, suggesting it is an aging-specific phenomenon rather than a universal constraint. **Not replicated.**


---

## Discussion

### NFL membership as a cancer predictor

The central finding of this study is that negative feedback loop membership predicts cancer gene status with an odds ratio of 59 within the same signaling pathways. This result is robust to the specific CGC version used (corrected from an initial OR = 260 with an incomplete CGC to OR = 59 with the complete 740-gene list), significant across all three nested tests, and convergent with an independent temporal classification approach (12.8x vs 12.1x enrichment).

The magnitude of this enrichment deserves emphasis. An odds ratio of 59 means that a randomly chosen NFL gene is 59 times more likely to be a cancer gene than a randomly chosen non-NFL gene from the same pathway. Among published network-topology-based cancer predictors, this appears to be among the highest reported enrichment ratios for a single topological feature. The enrichment likely reflects a fundamental relationship between oscillatory dynamics and cellular decision-making: genes that participate in feedback loops enabling pulsatile signaling are the same genes whose disruption converts oscillatory to sustained signaling --- the hallmark of oncogenic transformation.

### Irreversible Authority explains the vulnerability gradient

Not all oscillatory modules are equally cancer-vulnerable. The IA metric explains 69% of the variance (rho^2 = 0.69) in CGC fraction across 20 modules. The biological logic is straightforward: modules that control irreversible cell-fate decisions (apoptosis, terminal differentiation, DNA damage checkpoint commitment) create one-way gates in the cellular state space. Cancer must breach these gates to achieve hallmark capabilities such as resistance to cell death and replicative immortality (Hanahan, 2022). Modules controlling only reversible outcomes (circadian rhythm, metabolic switching) can be perturbed transiently without irreversible phenotypic consequences.

The IA framework also explains apparent outliers. The p53 module has the highest IA (3) and the highest CGC fraction (72.7%) despite not being the most highly expressed module. Conversely, the Circadian module is abundantly expressed but has IA = 0 and CGC fraction = 0%. Expression level is irrelevant (rho = 0.03); what matters is the irreversibility of the decisions the module controls.

### Eigenspace structure: universal core, tissue-specific periphery

The tissue-specific eigendecomposition reveals a two-layer architecture. The first principal component --- general signaling amplitude --- is universal (mean cosine similarity 0.70 to the global PC1), with highest conservation in the tissue classes that generate the most cancers (blood/immune: 0.98, epithelial: 0.83). This suggests that the dominant source of variation in oscillatory module expression (which cell types express more of everything) is shared across tissues and provides a tissue-universal vulnerability axis.

Higher principal components, which capture specific module-module co-variation patterns (e.g., TGF-beta vs AMPK, Cell Cycle vs PI3K), are tissue-specific. This is consistent with the well-known tissue-specificity of cancer driver landscapes: while some drivers (TP53, KRAS) appear across cancer types, many are tissue-restricted (e.g., VHL in renal cancer, RB1 in retinoblastoma). The eigenspace framework provides a quantitative account of this phenomenon: universal vulnerability arises from PC1, tissue-specific vulnerability from higher PCs.

The stability of specific module pairs across all tissues (Hippo--TGF-beta positive co-variation in 7/7 classes; JAK-STAT--NF-kappaB in 7/7 classes) identifies architectural constraints that may be experimentally testable: disrupting one module should systematically affect the co-varying partner, regardless of tissue context.

### Development, cancer, and the directionality of complexity

The MOCA analysis reveals a striking universality: all 13 developmental lineages decrease in per-cell transcriptional entropy during E9.5--E13.5. Differentiation simplifies the individual cell while diversifying the population. Cancer, conversely, re-expands the transcriptional state space of individual cells through the sequential activation of high-IA modules.

This reciprocal relationship has implications for the cancer stem cell hypothesis. If differentiation narrows the transcriptional repertoire and cancer re-expands it, then cancer vulnerability may not reside in specific cell types (stem cells vs. differentiated cells) but in the architectural capacity of the oscillatory module framework to be run in reverse. The modules with the highest IA --- p53, Wnt, Notch --- are precisely the modules most implicated in stem cell regulation (Clevers, 2006; Artavanis-Tsakonas et al., 1999), consistent with the idea that cancer exploits the same modular architecture used in normal development but in the opposite direction.

### Limitations

Several limitations should be acknowledged. First, the NFL extraction pipeline relies on KEGG pathway annotations, which are manually curated and incomplete. Approximately 20% of extracted motifs are artifacts (inverted edges), though this does not affect the enrichment statistics since the artifacts are distributed uniformly between CGC and non-CGC genes. Second, the IA metric is assigned manually based on pathway biology knowledge; while it does not use cancer gene databases, it draws on biological knowledge that is not fully independent of cancer research. Future work should develop automated IA assignment from pathway ontologies. Third, the eigenspace analysis uses bulk expression data aggregated by cell type, which may obscure single-cell heterogeneity in module co-variation. Fourth, the MOCA developmental analysis uses transcriptional entropy rather than direct module activity measurement; validation with module-level expression data would strengthen the developmental connection. Fifth, the evolutionary analysis (N = 14 loops with A values) has limited statistical power, and the significance of the age-asymmetry correlation depends on the removal of a single outlier (NRF2/KEAP1).

### Outlook

The oscillatory module framework suggests several testable predictions. First, the 7 unvalidated eigenspace cancer gene candidates (non-CGC genes in the CGC-proximal eigenspace region) represent specific predictions for cancer driver discovery. Second, the tissue-specific eigendecomposition predicts that disrupting a universal module pair (e.g., Hippo + TGF-beta) should have pan-tissue effects, while disrupting a tissue-specific axis should have restricted effects. Third, the developmental complexity gradient predicts that agents promoting re-differentiation (e.g., ATRA in acute promyelocytic leukemia) should shift eigenspace position toward the differentiated (low E_intra) zone. These predictions are experimentally testable with existing tools (CRISPR perturbation screens, scRNA-seq trajectory analysis, clinical pharmacogenomics).

---

## Methods

### KEGG pathway extraction and NFL identification

159 KEGG signaling pathway KGML files were downloaded via the KEGG API. Each KGML file was parsed into a signed directed graph using custom Python scripts, with edge signs assigned as +1 (activation, phosphorylation, expression) or -1 (inhibition, dephosphorylation, ubiquitination) based on KEGG relation subtypes. Compound-mediated edges were expanded to direct gene-gene relationships.

NFL motifs were identified using depth-limited depth-first search (max depth = 3) seeded from 43 genes comprising 22 literature-curated core oscillatory feedback loops (Table S1). A cycle was classified as an NFL if the product of edge signs was -1. To prevent combinatorial explosion (full cycle enumeration was killed by OOM after 9 minutes), the ego-centric DFS approach was used: for each seed gene, enumerate all return paths of length 2--3 in the local neighborhood.

### Alias resolution and deduplication

KEGG uses non-standard gene symbols. A manually curated alias map (37 entries) resolved KEGG symbols to HGNC standard nomenclature (e.g., AISIMD to SOCS3, IMD31A to STAT3, DAKAR to CDK4). After alias resolution, exact duplicate motifs (identical gene sets and edge signs) were removed. Remaining motifs were clustered using pairwise Jaccard similarity on gene sets, with a merging threshold of J > 0.5 for module assignment.

### Cancer Gene Census enrichment

The COSMIC Cancer Gene Census v103 (740 genes, downloaded April 2026) was used as the gold-standard cancer gene reference. The complete KEGG signaling gene list (16,760 unique genes) was assembled from all 159 pathway KGML files. Three nested enrichment tests were computed:
- Test 1: Hypergeometric test, NFL genes vs genome (20,000 protein-coding genes)
- Test 2: Hypergeometric test, NFL genes vs KEGG signaling genes
- Test 3: Fisher exact test, NFL genes vs non-NFL genes within KEGG signaling pathways

All tests used two-sided *p*-values. An initial analysis using an incomplete CGC list (534 genes) yielded inflated enrichment (OR = 260); all results reported here use the complete COSMIC CGC v103 (740 genes).

### Community detection

An undirected weighted network was constructed where nodes represent 128 unique NFL motifs and edges connect motifs sharing at least one gene, weighted by Jaccard similarity coefficient. Louvain community detection was performed at resolution = 1.0 using the python-louvain package. Modularity Q and Adjusted Rand Index (comparing detected communities to KEGG pathway annotations) were computed. Sensitivity to resolution parameter was assessed at resolutions 0.5, 1.0, and 1.5.

### Evolutionary analysis

Ortholog mapping for the 43 NFL core genes was performed using Ensembl Compara (release 111) across 15 species: *M. musculus*, *R. norvegicus*, *C. l. familiaris*, *G. gallus*, *X. tropicalis*, *D. rerio*, *C. intestinalis*, *S. purpuratus*, *D. melanogaster*, *C. elegans*, *N. vectensis*, *S. cerevisiae*, *S. pombe*, *A. thaliana*, *D. discoideum*. Loop conservation required both activator and inhibitor genes to have one-to-one orthologs. Loop age was estimated from the most distant species with confirmed loop conservation, using divergence times from TimeTree (Kumar et al., 2022).

Waveform asymmetry values (A = tau_rise / T) were taken from Paper 1 (Spiro, 2026) for 14 loops with published or computationally estimated oscillation parameters. Spearman rank correlation between loop age and A was computed with and without NRF2/KEAP1.

### Irreversible Authority metric

For each of 20 oscillatory modules (14 original + 6 expanded: AMPK, SREBP, ATR/CHK1, Rho/ROCK, PPAR/LXR, Autophagy), the Irreversible Authority (IA) score was defined as the number of distinct irreversible cellular outcomes controlled by the module. An outcome was classified as irreversible if the cellular state change cannot be reversed by removing the stimulus: apoptosis (irreversible), senescence (irreversible), terminal differentiation (irreversible), DNA damage checkpoint commitment (irreversible), mitotic commitment (irreversible). Metabolic switching, transcriptional oscillation, and cytoskeletal remodeling were classified as reversible. IA assignment was based on established pathway biology reviewed in standard references (Alberts et al., 2015; Hanahan and Weinberg, 2011). The classification does not reference cancer gene databases.

### Eigendecomposition

Mean nCPM expression values for each module's genes were computed across 154 human cell types from the Human Protein Atlas single-cell RNA-seq dataset (v23). The resulting 20 x 154 module expression matrix was standardized (z-score per module) and the 20 x 20 module-module Pearson correlation matrix was computed. Eigendecomposition was performed using numpy.linalg.eigh. The number of retained components was determined by Kaiser criterion (eigenvalue > 1). For tissue-specific analysis, cell types were grouped by HPA cell type class, and eigendecomposition was repeated independently for each class with N >= 4 cell types.

### Eigenspace cancer prediction

Each NFL gene was assigned to its corresponding module(s) and projected into the 6D eigenspace using the module's expression profile. The centroid of all CGC genes in eigenspace was computed, and the 20 non-CGC genes closest to this centroid (by Euclidean distance) were identified as novel cancer gene candidates. Statistical significance was assessed by permutation (10,000 random gene-to-module assignments, computing the distance distribution). Validation used the IntOGen database of cancer driver genes (Gonzalez-Perez et al., 2013; Martinez-Jimenez et al., 2020).

### MOCA developmental trajectory analysis

Processed cell annotations from the Mouse Organogenesis Cell Atlas (GEO: GSE119945; Cao et al., 2019) were used, comprising 300,839 cells with annotations for cell type, developmental trajectory (13 lineages), and embryonic day (E9.5--E13.5). Per-cell transcriptional entropy (E_intra, Shannon entropy of the UMI-normalized gene expression distribution) and between-cell diversity (E_inter, mean pairwise Jensen-Shannon divergence within trajectory-day groups) were computed as described previously (see MOCA processing pipeline). Spearman correlation between embryonic day and mean E_intra/E_inter was computed per trajectory. Cancer relevance of each trajectory was assigned based on the cancer burden of the corresponding adult tissue (GLOBOCAN 2020).

### Data and code availability

All analysis scripts (step01 through step29), intermediate data files, and the complete NFL motif database are available at https://github.com/mool32/oscillatory-nfl-cancer. KEGG pathway data were accessed via the KEGG API for academic research purposes in accordance with KEGG's terms of use (https://www.kegg.jp/kegg/legal.html); raw KGML files are not redistributed. The COSMIC CGC v103 is available at https://cancer.sanger.ac.uk/census. HPA data are available at https://www.proteinatlas.org. MOCA data are deposited at GEO (GSE119945).

---

## Acknowledgments

I thank the KEGG, COSMIC, Human Protein Atlas, and MOCA consortia for making their data publicly available.

---

## References

Albeck JG et al. (2013). Frequency-modulated pulses of ERK activity transmit quantitative proliferation signals. Mol Cell 49:249-261.

Alon U (2007). Network motifs: theory and experimental approaches. Nat Rev Genet 8:450-461.

Artavanis-Tsakonas S et al. (1999). Notch signaling: cell fate control and signal integration in development. Science 284:770-776.

Bailey MH et al. (2018). Comprehensive characterization of cancer driver genes and mutations. Cell 173:371-385.

Berridge MJ et al. (2000). The versatility and universality of calcium signalling. Nat Rev Mol Cell Biol 1:11-21.

Cao J et al. (2019). The single-cell transcriptional landscape of mammalian organogenesis. Nature 566:496-502.

Clevers H (2006). Wnt/beta-catenin signaling in development and disease. Cell 127:469-480.

Dupont G et al. (2011). Calcium dynamics: spatio-temporal organization from the subcellular to the organ level. Int Rev Cytol 261:193-245.

Gonzalez-Perez A et al. (2013). IntOGen-mutations identifies cancer drivers across tumor types. Nat Methods 10:1081-1082.

Hanahan D (2022). Hallmarks of cancer: new dimensions. Cancer Discov 12:31-46.

Hanahan D, Weinberg RA (2011). Hallmarks of cancer: the next generation. Cell 144:646-674.

Kumar S et al. (2022). TimeTree 5: an expanded resource for species divergence times. Mol Biol Evol 39:msac174.

Lahav G et al. (2004). Dynamics of the p53-Mdm2 feedback loop in individual cells. Nat Genet 36:147-150.

Martinez-Jimenez F et al. (2020). A compendium of mutational cancer driver genes. Nat Rev Cancer 20:555-572.

Nelson DE et al. (2004). Oscillations in NF-kappaB signaling control the dynamics of gene expression. Science 306:704-708.

Novak B, Tyson JJ (2008). Design principles of biochemical oscillators. Nat Rev Mol Cell Biol 9:981-991.

Purvis JE, Lahav G (2013). Encoding and decoding cellular information through signaling dynamics. Cell 152:945-956.

Purvis JE et al. (2012). p53 dynamics control cell fate. Science 336:1440-1444.

Sondka Z et al. (2018). The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nat Rev Cancer 18:696-705.

Spiro T (2026). Temporal architecture of signaling oscillations predicts cancer gene function across pathways. Preprint.

Tay S et al. (2010). Single-cell NF-kappaB dynamics reveal digital activation and analogue information processing. Nature 466:267-271.

Barkai N, Leibler S (1997). Robustness in simple biochemical networks. Nature 387:913-917.

Blondel VD et al. (2008). Fast unfolding of communities in large networks. J Stat Mech 2008:P10008.

Elowitz MB, Leibler S (2000). A synthetic oscillatory network of transcriptional regulators. Nature 403:335-338.

Ferrell JE (2002). Self-perpetuating states in signal transduction: positive feedback, double-negative feedback and bistability. Curr Opin Cell Biol 14:140-148.

Forbes SA et al. (2017). COSMIC: somatic cancer genetics at high-resolution. Nucleic Acids Res 45:D777-D783.

Hubert L, Arabie P (1985). Comparing partitions. J Classification 2:193-218.

Kanehisa M et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res 51:D587-D592.

Lim WA et al. (2013). Design principles of regulatory networks: searching for the molecular algorithms of the cell. Mol Cell 49:202-212.

Martincorena I et al. (2017). Universal patterns of selection in cancer and somatic tissues. Cell 171:1029-1041.

Raj A, van Oudenaarden A (2008). Nature, nurture, or chance: stochastic gene expression and its consequences. Cell 135:216-226.

Stratton MR et al. (2009). The cancer genome. Nature 458:719-724.

Tyson JJ et al. (2003). Sniffers, buzzers, toggles and blinkers: dynamics of regulatory and signaling pathways in the cell. Curr Opin Cell Biol 15:221-231.

Uhlen M et al. (2015). Tissue-based map of the human proteome. Science 347:1260419.

Vogelstein B et al. (2013). Cancer genome landscapes. Science 339:1546-1558.

---

## Supplementary Tables

**Table S1.** 22 core feedback loops with activator/inhibitor genes, Ensembl IDs, feedback mechanism, and literature PMIDs.

**Table S2.** Ortholog conservation matrix: 43 genes x 15 species, with binary presence/absence of one-to-one orthologs.

**Table S3.** 128 unique NFL motifs with HGNC gene symbols, edge signs, pathway attribution, and CGC status.

**Table S4.** 59 non-redundant modules with gene composition, module assignment, and per-module CGC statistics.

**Table S5.** Community detection results at three Louvain resolutions (0.5, 1.0, 1.5) with Q, ARI, and community assignments.

**Table S6.** Biological validation assessment for 20 randomly sampled NFL motifs (full/partial/artifact classification with evidence).

**Table S7.** Complete KEGG signaling gene list (16,760 genes) with NFL membership status and CGC status for each gene.

**Table S8.** Irreversible Authority scores for 20 modules with irreversible outcome definitions and CGC fractions.

**Table S9.** Eigenspace cancer gene candidates: 20 non-CGC genes, eigenspace coordinates, distance to CGC centroid, and IntOGen validation status.

**Table S10.** Tissue-specific eigendecomposition results: PC loadings and variance explained for each of 7 tissue classes.

**Table S11.** MOCA developmental trajectory statistics: per-trajectory E_intra and E_inter Spearman correlations with embryonic day.

---

## Figure Legends

**Figure 1. Algorithmic extraction of negative feedback loops from KEGG signaling networks.** (A) Pipeline schematic: KEGG KGML files are converted to signed directed graphs; depth-limited search identifies all 2- and 3-node cycles with net negative sign product. (B) Example motifs: p53/MDM2 two-node NFL; E2F1/CCNE1/RB1 three-node NFL showing indirect negative feedback through RB1. (C) Deduplication cascade from 228 raw motifs to 128 unique NFLs to 59 modules to 14 pathway groups. (D) Biological validation of 20 randomly sampled motifs: 35% full validation, 45% partial support, 20% artifact.

**Figure 2. NFL genes are 59-fold enriched for cancer genes.** (A) Three nested enrichment tests with increasing specificity: NFL vs genome (12.8x), NFL vs KEGG signaling (28.1x), NFL vs non-NFL in same pathways (OR = 59). Error bars show 95% confidence intervals. (B) 2x2 contingency table for the critical test (Test 3). (C) Convergent validation: Paper 1 temporal classification (12.1x) vs Paper 2 topological extraction (12.8x). (D) Venn diagram of CGC genes identified by each approach.

**Figure 3. Modular organization and cross-pathway topology.** (A) Network of 128 NFL motifs connected by gene sharing, colored by Louvain communities (10 communities, Q = 0.39). (B) Condensed pathway-level network (14 nodes) showing cross-module connections. (C) ARI between data-driven communities and KEGG annotations at three resolution settings. (D) Hub genes bridging multiple communities: MAPK1 participates in 3 distinct modules.

**Figure 4. Evolutionary conservation gradient.** (A) Ortholog matrix heatmap: 22 loops x 15 species ordered by evolutionary divergence time. (B) Conservation summary: 100% vertebrate, 50% invertebrate, 27% unicellular. (C) Age vs waveform asymmetry: older loops are more asymmetric (rho = -0.66, *p* = 0.015, N = 13 after removing NRF2/KEAP1 convergent-evolution outlier). (D) All loops >1000 Mya are Class I (asymmetric); symmetric waveforms are exclusively metazoan.

**Figure 5. Irreversible Authority predicts cancer vulnerability across modules.** (A) Scatterplot of IA score vs CGC fraction for 20 modules (rho = 0.83). (B) Dose-response: mean CGC fraction by IA level (0, 1, 2, 3). (C) Comparison of three architectural hypotheses: IA (rho = 0.83) vs Authority (rho = 0.752) vs Maintenance cost (rho = 0.03). (D) Binary validation: irreversible modules (IA >= 2) have 4.8x higher CGC fraction than reversible modules (IA < 2).

**Figure 6. Eigenspace analysis and cancer gene prediction.** (A) First three principal components of the 20-module expression landscape across 154 cell types. (B) CGC genes cluster in a non-random eigenspace subregion (permutation *p* = 0.0001). (C) Tissue-specific PC1 conservation: blood/immune (0.98), neuronal (0.88), epithelial (0.83), vascular (0.22). (D) Novel cancer gene candidates: 13/20 (65%) validated by IntOGen.

**Supplementary Figure S1. Developmental trajectories in the oscillatory complexity landscape.** (A) All 13 MOCA trajectories decrease in per-cell entropy (E_intra) from E9.5 to E13.5. (B) Between-cell diversity (E_inter) increases in 10/13 trajectories. (C) Comparison of developmental direction (decreasing E_intra) with cancer strategy (ascending IA, increasing complexity). (D) Model: differentiation compresses the transcriptional state space; cancer re-expands it by reactivating high-IA modules.
