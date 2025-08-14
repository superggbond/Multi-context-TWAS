# Multi-context TWAS
An analytical pipeline for multi-context transcriptome-wide association studies

## 1. Introduction
This GitHub repository contains sample scripts that used for the analyses performed in the manuscript:

Zhang, H., Patrick, M. T., Sarkar, M. K., ... & Tsoi, L. C. (2025). Multi-cytokine based TWAS for seven inflammatory skin disorders identify candidate causal genes in keratinocytes.

We performed a multi-context marginal TWAS analyses on skin disorders. The following sections will go through the workflow of the main analysis step-by-step, using our multi-cytokine based TWAS as an example.

## 2. GReX modeling using BSLMM
For the standard marginal TWAS, the genetically regulated gene expression (GReX) level is reuqired for each gene. To predict GReX for a gene, a fixed window of cis-SNPs (usually +-100kb or +-500kb of a gene's transcription starting site) are used, and the weight for each cis-SNP is required to pair with GWAS summary statistics of a disorder. In our case, we applied [BSLMM](https://github.com/genetics-statistics/GEMMA) for each gene under each cytokine stimulated conditions to generate the cytokine-specific set of cis-SNP weights. Pre-computed SNP weights on different tissues are publicly available [here](https://predictdb.org/) trained by the [PrediXcan](https://github.com/hakyimlab/MetaXcan) family of methods.

Sample scripts we used to apply BSLMM are provided as [run_bslmm.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_bslmm.R)

## 3. Single-context TWAS
In order to conduct the multi-context TWAS, we first need to perform the single-context TWAS on each condition using S-PrediXcan. In addition to the GWAS summary statistics, S-PrediXcan also requires the inputs of cis-SNP weights as a .db file and the LD reference file. Sample scripts to generate the cis-SNP weights .db file are provided as [prep_dbfile.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/prep_dbfile.R), and the Sample scripts to generate the LD reference file are provided as [prep_LDref.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/prep_LDref.R).

Then S-PrediXcan can be applied following the sample commands listed in [run_SPrediXcan.sh](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_SPrediXcan.sh).

## 3. Multi-context TWAS
After the S-PrediXcan results are generated for all the conditions, S-MulTiXcan.sh can be applied to perform the multi-context TWAS, following the sample commands listed in [run_SMulTiXcan.sh](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_SMulTiXcan.sh).

