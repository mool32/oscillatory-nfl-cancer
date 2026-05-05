# Out-of-scope material from Paper 2 → future papers

## Destination: Paper 3 (Perceptual Eigenspace)

### Tissue-specific eigendecomposition (was section 8)
- PC1 universal (mean cosine sim 0.70), blood/immune 0.98, epithelial 0.83, vascular 0.22
- Higher PCs tissue-specific
- Hippo-TGFb invariant pair (positive in 7/7 classes)
- JAK-STAT-NFkB coupling (7/7 classes)
- Data: `tissue_specific_eigen.json`, `step28_tissue_specific_eigen.py`
- Figure: `fig6_eigenspace.pdf` panel C (tissue PC1 bars)

### MOCA developmental trajectory (was section 9)
- All 13 trajectories decrease E_intra (E9.5→E13.5)
- E_inter increases in 10/13
- Dev-aging replication: rho=-0.02 (not replicated)
- Cancer-risk comparison: p=0.61 (not significant)
- Data: `moca_differentiation.json`, `step29_moca_differentiation.py`
- Figure: `fig7_development.pdf`

### Killed hypothesis: Development-aging replication
- rho=-0.02, p=0.95, N=13
- Aging entropy pattern does not replicate in development

## Terminology note
- Paper 2: "module activity eigenspace" or "20-module signaling space"
- Paper 3: "perceptual eigenspace" or "43-module perceptual space"
