# Multi-context TWAS
An analytical pipeline for multi-context transcriptome-wide association studies

## 1. Introduction
This GitHub repository contains sample scripts that used for the analyses performed in the manuscript:

Zhang, H., Patrick, M. T., Sarkar, M. K., ... & Tsoi, L. C. (2025). Multi-cytokine based TWAS for seven inflammatory skin disorders identify candidate causal genes in keratinocytes.

We performed a multi-context marginal TWAS analyses on skin disorders. The following sections will go through the workflow of the main analysis step-by-step, using our multi-cytokine based TWAS as an example.

## 2.1 GReX modeling using BSLMM
For the standard marginal TWAS, the genetically regulated gene expression (GReX) level is reuqired for each gene. To predict GReX for a gene, cis-SNPs in a fixed window (usually +-100kb or +-500kb of a gene's transcription starting site) are used, and the weight for each cis-SNP is required to pair with GWAS summary statistics of a disorder. In our case, we applied [BSLMM](https://github.com/genetics-statistics/GEMMA) for each gene under each cytokine stimulated condition to generate the cytokine-specific set of cis-SNP weights. Pre-computed SNP weights on different tissues are publicly available [here](https://predictdb.org/) trained by the [PrediXcan](https://github.com/hakyimlab/MetaXcan) family of methods.

Sample scripts we used to apply BSLMM are provided as [run_bslmm.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_bslmm.R)

BSLMM requires two input files:

(1). A .txt or .txt.gz file containing genotype information. The first column is SNP id, the second and third columns are allele types with minor allele first, and the remaining columns are the posterior/imputed mean genotypes of diﬀerent individuals numbered between 0 and 2. An example mean genotype file with two SNPs and three individuals is as follows:
```txt
rs1, A, T, 0.02, 0.80, 1.50
rs2, G, C, 0.98, 0.04, 1.00
```

(2). A .txt file containing phenotype information. Each line is a number indicating the phenotype value for each individual in turn, in the same order as in the mean genotype file. Notice that only numeric values are allowed and characters will not be recognized by the software. Missing phenotype information is denoted as NA. The number of rows should be equal to the number of individuals in the mean genotype file. An example phenotype file with five individuals and one phenotype is as follows:
```txt
1.2
NA
2.7
-0.2
3.3
```

The BSLMM output file [prefix].param.txt contains the posterior mean estimates for the eﬀect size parameters for each SNP. An example file with a few SNPs is shown below:
```txt
chr rs ps n_miss alpha beta gamma
1 rs3683945 3197400 0 -7.314495e-05 0.000000e+00 0.000000e+00
1 rs3707673 3407393 0 -7.314495e-05 0.000000e+00 0.000000e+00
1 rs6269442 3492195 0 -3.412974e-04 0.000000e+00 0.000000e+00
1 rs6336442 3580634 0 -8.051198e-05 0.000000e+00 0.000000e+00
1 rs13475700 4098402 0 -1.200246e-03 0.000000e+00 0.000000e+00
```
For each SNP, the final estimated effect size can be calculated as: alpha+beta*gamma.

## 2.2 GReX modeling using SuSiE (alternative method)
SuSiE (Sum of Single Effects) is one of the most popular models for sparse multiple regression. It can also be applied for GReX modeling. 

SuSiE model can be implemented by the R package [susieR](https://github.com/stephenslab/susieR). It also requires two input objects in R: (1). A R matrix for genotypes. (2). A R vector for phenotype. And the SuSiE method can be applied in R as:
```R
library(susieR)

# apply SuSiE method for GReX modeling: X is the genotype matrix; y is the phenotype vector.
fit <- susie(X, y)
# calculate posterior mean effect sizes of SNPs
beta_hat <- colSums(fit$alpha * fit$mu)
```

## 3. Perform single-context TWAS
In order to conduct the multi-context TWAS, we first need to perform the single-context TWAS on each condition using S-PrediXcan. In addition to the GWAS summary statistics, S-PrediXcan also requires the inputs of cis-SNP weights as a .db file and the LD reference file. Sample scripts to generate the cis-SNP weights .db file are provided as [prep_dbfile.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/prep_dbfile.R), and the sample scripts to generate the LD reference file are provided as [prep_LDref.R](https://github.com/superggbond/Multi-context-TWAS/blob/main/prep_LDref.R).

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

More details can be found [here](https://github.com/hakyimlab/MetaXcan).

## 4. Perform multi-context TWAS
After the S-PrediXcan results are generated for all the conditions, S-MulTiXcan can be applied to perform the multi-context TWAS, following the sample commands listed in [run_SMulTiXcan.sh](https://github.com/superggbond/Multi-context-TWAS/blob/main/run_SMulTiXcan.sh).

The output will be saved as a tab-separated text file, and here is an example loaded in R:
```r
> head(res)
      gene gene_name       pvalue n n_indep     p_i_best  t_i_best  p_i_worst
1   NMNAT1    NMNAT1 6.436820e-28 8       6 1.248346e-18       IL4 0.01620177
2      FLG       FLG 6.113827e-24 8       5 1.066203e-20      NHEK 0.97930403
3 CTNNBIP1  CTNNBIP1 6.051845e-23 8       5 3.489429e-04 IL17ATNFa 0.24139310
4  KRT8P21   KRT8P21 1.503697e-21 8       4 2.367278e-14      IFNg 0.91871506
5   SDHDP6    SDHDP6 2.921108e-18 8       5 2.936515e-03       IL4 0.69308320
6    UBE4B     UBE4B 6.384423e-17 8       3 9.171576e-08      IFNa 0.77020149
  t_i_worst eigen_max    eigen_min eigen_min_kept     z_min    z_max     z_mean
1      IL13  3.153889 0.0075440227      0.1339378 -5.347993 8.810279  1.0380492
2      IFNg  4.244840 0.0161399116      0.2241544 -9.329252 2.465591 -1.6628579
3      IFNa  3.950087 0.0194258338      0.1373150 -3.569240 3.575962  0.1951729
4      IFNa  5.635723 0.0003576079      0.2489683 -7.628920 4.844059 -1.4727387
5      NHEK  3.042455 0.0005394716      0.1887045 -2.974307 2.623469  0.6251076
6     IL17A  6.278585 0.0035381646      0.2330906 -2.795889 5.342416  1.1213563
      z_sd tmi status
1 5.493796   6      0
2 3.961859   5      0
3 2.532445   5      0
4 3.743329   4      0
5 1.925624   5      0
6 2.497475   3      0
```
Where each row is a gene's association result after leveraging the information across multiple contexts:
* `gene`: a gene's id or name used in the GReX model
* `gene_name`: gene name as listed in the GReX model
* `pvalue`: p-value of S-MultiXcan association
* `n`: number of "tissues"/"conditions" available for this gene
* `n_indep`: number of independent components of variation kept among the tissues' predictions. (Synthetic independent tissues)
* `p_i_best`: best p-value of single-tissue S-PrediXcan association.
* `t_i_best`: name of best single-tissue S-PrediXcan association.
* `p_i_worst`: worst p-value of single-tissue S-PrediXcan association.
* `t_i_worst`: name of worst single-tissue S-PrediXcan association.
* `eigen_max`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the top independent component
* `eigen_min`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the last independent component
* `eigen_min_kept`: In the SVD decomposition of predicted expression correlation: eigenvalue (variance explained) of the smalles independent component that was kept.
* `z_min`: minimum z-score among single-tissue S-Predican associations.
* `z_max`: maximum z-score among single-tissue S-Predican associations.
* `z_mean`: mean z-score among single-tissue S-Predican associations.
* `z_sd`: standard deviation of the mean z-score among single-tissue S-Predican associations.
* `tmi`: trace of `T * T'`, 
where `T`is correlation of predicted expression levels for different tissues 
multiplied by its SVD pseudo-inverse. 
It is an estimate for number of indepent components of variation in predicted expresison across tissues (typically close to `n_indep`)
* `status`: If there was any error in the computation, it is stated here

More details can be found [here](https://github.com/hakyimlab/MetaXcan).

