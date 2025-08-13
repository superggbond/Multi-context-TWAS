# sample bash scripts to run S-PrediXcan
~/MetaXcan/software/SPrediXcan.py \
--keep_non_rsid \
--model_db_path /net/psoriasis/psorgenom/haihan/multi_TWAS/db_files/IFNa.db \
--covariance /net/psoriasis/psorgenom/haihan/multi_TWAS/cov_files/chr$chr.txt.gz \
--gwas_file /net/psoriasis/psorgenom/haihan/multi_TWAS/GWAS_sum/psor/chr$chr \
--snp_column SNP \
--effect_allele_column A2 \
--non_effect_allele_column A1 \
--beta_column BETA \
--pvalue_column P \
--output_file IFNa/chr$chr.csv
