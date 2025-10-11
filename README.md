# Optimizing Gastric Organoid Culture Conditions: A Transcriptomic Comparison of 2D vs 3D Systems and Synthetic vs Animal-Derived Growth Factors for Maintaining Primary Cell Fidelity

## Project Description

This study addresses a critical challenge in organoid research: identifying culture conditions that best preserve the biological authenticity of primary gastric cells while minimizing costs. Mouse gastric epithelial cells are cultured under four distinct conditions:

1. Freshly isolated primary cells serving as the gold standard
2. 3D organoid cultures with animal-derived growth factors (AD-GF)
3. 3D organoid cultures with synthetic growth factors (S-GF)
4. 2D monolayer cultures with synthetic growth factors (uncertainty remains regarding purified growth factor usage)

Using bulk RNA sequencing, this project aims to comprehensively compare the transcriptomic profiles across these culture systems to answer three fundamental questions:

1. **Culture dimensionality**: Do 3D organoid cultures better recapitulate the primary cell transcriptome compared to 2D monolayer cultures?
2. **Cell type composition**: Are there differences in cellular heterogeneity and cell type representation between 2D and 3D culture systems that could explain transcriptomic differences?
3. **Growth factor economics**: Can cost-effective synthetic growth factors adequately substitute expensive animal-derived factors without compromising biological fidelity?

The findings will provide evidence-based recommendations for establishing gastric organoid culture protocols that balance biological authenticity, experimental reproducibility, and economic feasibility—a critical consideration for scaling research and potential therapeutic applications.

## Working Plan

### Phase 1: Data Quality Control and Preprocessing

**Objectives:**
- Assess sequencing quality and depth
- Identify potential technical artifacts or batch effects
- Prepare normalized expression matrices

**Deliverables:**
- QC report with sequencing metrics
- Normalized expression data ready for analysis

### Phase 2: Global Transcriptomic Analysis

**Objectives:**
- Compare overall transcriptomic similarity to primary cells
- Identify culture-specific gene signatures
- Determine major sources of variation

**Deliverables:**
- PCA/UMAP plots showing sample relationships
- Correlation heatmaps
- Differential expression results for all comparisons

### Phase 3: Cell Type Deconvolution Analysis

**Objectives:**
- Estimate cell type composition in each culture condition
- Identify cell type-specific markers
- Determine if culture conditions affect cellular heterogeneity

**Deliverables:**
- Cell type proportion estimates
- Cell type marker expression profiles
- Analysis of gastric cell lineages (pit, parietal, chief, enteroendocrine cells, etc.)

### Phase 4: Growth Factor Comparison

**Objectives:**
- Directly compare 3D cultures with animal-derived vs synthetic growth factors
- Identify pathways affected by growth factor source
- Assess cost-benefit trade-offs

**Deliverables:**
- Statistical comparison of AD-GF vs S-GF organoids
- Pathway enrichment analysis
- Recommendations on growth factor selection

### Phase 5: Integration and Reporting

**Objectives:**
- Synthesize findings across all analyses
- Generate high-quality figures
- Provide actionable recommendations

**Deliverables:**
- Comprehensive analysis report
- Manuscript-level figures
- Culture protocol recommendations
