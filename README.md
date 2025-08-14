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

The output will be saved as a .csv file, and here is an example loaded in R:
```r
> head(res)
           gene     gene_name    zscore effect_size       pvalue        var_g
1       KRT8P21       KRT8P21 -7.628920  -10.173711 2.367278e-14 4.013279e-04
2        NMNAT1        NMNAT1  4.696846    3.613510 2.642099e-06 1.392131e-03
3          MCL1          MCL1  4.534705    1.227059 5.768414e-06 1.004760e-02
4        MYBPHL        MYBPHL -4.207041   -2.665227 2.587362e-05 1.789746e-03
5 RP11-296A18.6 RP11-296A18.6 -4.155541   -0.991950 3.245180e-05 1.191158e-02
6          LZIC          LZIC  3.874997   24.363732 1.066262e-04 1.932156e-05
  pred_perf_r2 pred_perf_pval pred_perf_qval n_snps_used n_snps_in_cov
1            0              0              0         601           601
2            0              0              0         160           160
3            0              0              0         299           299
4            0              0              0         321           321
5            0              0              0         133           133
6            0              0              0         187           187
  n_snps_in_model
1             601
2             160
3             299
4             321
5             133
6             187
```
Where each row is a gene's association result:
* `gene`: a gene's id or name used in the GReX model
* `gene_name`: gene name as listed in the GReX model
* `zscore`: S-PrediXcan's association result for the gene
* `effect_size`: S-PrediXcan's association effect size for the gene. Can only be computed when `beta` from the GWAS is used.
* `pvalue`: P-value of the aforementioned statistic
* `pred_perf_r2`: (cross-validated) R2 of tissue model's correlation to gene's measured transcriptome (prediction performance). It is not applicable for our BSLMM model.
* `pred_perf_pval`: pval of tissue model's correlation to gene's measured transcriptome (prediction performance). It is not applicable for our BSLMM model.
* `pred_perf_qval`: qval of tissue model's correlation to gene's measured transcriptome (prediction performance). It is not applicable for our BSLMM model.
* `n_snps_used`: number of snps from GWAS that got used in S-PrediXcan analysis
* `n_snps_in_cov`: number of snps in the LD reference matrix
* `n_snps_in_model`: number of snps in the GReX model
* `var_g`: variance of the gene expression, calculated as `W' * G * W`
(where `W` is the vector of SNP weights in a gene's model,
`W'` is its transpose, and `G` is the LD reference matrix)

## 3. Multi-context TWAS
After the S-PrediXcan results are generated for all the conditions, S-MulTiXcan.sh can be applied to perform the multi-context TWAS, following the sample commands listed in [run_SMulTiXcan.sh](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_SMulTiXcan.sh).

