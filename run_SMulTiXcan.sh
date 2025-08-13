# sample scripts to run S-MulTiXcan
~/MetaXcan/software/SMulTiXcan.py \
--models_folder /net/psoriasis/psorgenom/haihan/multi_TWAS/db_files \
--models_name_pattern "(.*).db" \
--snp_covariance /net/psoriasis/psorgenom/haihan/multi_TWAS/cov_files/chr$chr.txt.gz  \
--metaxcan_folder /net/psoriasis/psorgenom/haihan/multi_TWAS/smulti_results/psor/chr$chr \
--metaxcan_file_name_parse_pattern "(.*).csv" \
--gwas_file /net/psoriasis/psorgenom/haihan/multi_TWAS/GWAS_sum/psor/chr$chr \
--snp_column SNP \
--effect_allele_column A2 \
--non_effect_allele_column A1 \
--beta_column BETA \
--pvalue_column P \
--keep_non_rsid \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output chr$chr.csv
