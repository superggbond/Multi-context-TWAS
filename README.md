# Multi-context TWAS
An analytical pipeline for multi-context transcriptome-wide association studies

## 1. Introduction
This GitHub repository contains sample scripts that used for the analyses performed in the manuscript:

Zhang, H., Patrick, M. T., Sarkar, M. K., ... & Tsoi, L. C. (2025). Multi-cytokine based TWAS for seven inflammatory skin disorders identify candidate causal genes in keratinocytes.

We performed a multi-context marginal TWAS analyses on skin disorders. The following sections will go through the workflow of the main analysis step-by-step, using our multi-cytokine based TWAS as an example.

## 2. GReX modeling using BSLMM
For the standard marginal TWAS, the genetically regulated gene expression (GReX) level is reuqired for each gene. To predict GReX for a gene, a fixed window of cis-SNPs (usually +-100kb or +-500kb of a gene's transcription starting site) are used, and the weight for each cis-SNP is required to pair with GWAS summary statistics of a disorder. In our case, we applied [BSLMM](https://github.com/genetics-statistics/GEMMA) for each gene under each cytokine stimulated conditions to generate the cytokine-specific set of cis-SNP weights. Pre-computed SNP weights on different tissues are publicly available [here](https://predictdb.org/) trained by the [PrediXcan](https://github.com/hakyimlab/MetaXcan) family of methods.

Sample scripts we used to apply BSLMM are provided as [run_bslmm.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_bslmm.R)
