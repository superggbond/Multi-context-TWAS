## scripts used to build up the BSLMM prediction models
library(data.table)
library(snpStats)

# run BSLMM on local snps for each gene +- 100kb range
system(paste0("mkdir -p inter_geno"))  

for (i in all_jobs) {
  gene_name = colnames(gene_tb)
  lower = gene_pos[gene_pos$gene==gene_name,3]
  upper = gene_pos[gene_pos$gene==gene_name,4]
  chr = gene_pos[gene_pos$gene==gene_name,2]
  
  # subset genotype data: used 100kb here
  if (lower < 1e5) {lower = 1e5}
  plink_cmd = paste0("/net/psoriasis/home/hhzhang/plink2 --bfile /net/psoriasis/psorgenom/haihan/multi_TWAS/ref_genotypes/chr",chr," --chr ",chr,
                     " --from-bp ", lower - 1e5," --to-bp ",upper + 1e5," --make-bed --out ./inter_geno/inter_geno 2> plink_out")
  system(plink_cmd)
  # check if there is error in plink step
  system("wc -l plink_out > plink_out2")
  plink.error = as.numeric(read.table("plink_out2"))[1]
  system("rm plink_out*")
  
  if (plink.error == 0) {
    ## compute scaled genotype matrix and write out
    geno_file = read.plink("./inter_geno/inter_geno.bed")
    genotype_mat = t(scale(as(geno_file$genotypes, "numeric")))
    genotype_mat = genotype_mat[complete.cases(genotype_mat), ]
    # subset snps overlapping with GWAS summary data
    gwas_ref = fread(paste0("/net/psoriasis/psorgenom/haihan/multi_TWAS/GWAS_sum/psor/chr",chr))
    genotype_mat = genotype_mat[rownames(genotype_mat)%in%gwas_ref$SNP,]
    
    if (is.null(nrow(genotype_mat))) {
      next
    } else if (nrow(genotype_mat) < 10) {
      next
    } else {
      # add allele info
      allele_info_select = gwas_ref[gwas_ref$SNP%in%rownames(genotype_mat),c(1,3,2)]
      # save genotype and gene expression as BSLMM input files
      genotype_final = cbind(allele_info_select,genotype_mat)
      fwrite(genotype_final,paste0("./geno_",gene_name,".txt"),col.names = F,row.names = F,quote = F, sep = " ")
      write.table(IFNa[,gene_name],paste0("./IFNa_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(IFNg[,gene_name],paste0("./IFNg_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(IL13[,gene_name],paste0("./IL13_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(IL17A[,gene_name],paste0("./IL17A_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(IL17ATNFa[,gene_name],paste0("./IL17ATNFa_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(IL4[,gene_name],paste0("./IL4_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(NHEK[,gene_name],paste0("./NHEK_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      write.table(TNFa[,gene_name],paste0("./TNFa_exp_",gene_name,".txt"),col.names = F,row.names = F,quote = F)
      
      for (smp in smps){
        # run bslmm
        bslmm.cmd = paste0("~/gemma/gemma-0.98.3-linux-static -g ","./geno_",gene_name,".txt ", "-p ","./",smp,"_exp_",gene_name,".txt " ,"-bslmm 1 -w 6000 -s 2000 -rpace 1 -notsnp -o bslmm_output -outdir ./",gene_name)
        system(bslmm.cmd)
        # read and save the output from bslmm file
        paramter.file = read.table(paste0("./",gene_name,"/bslmm_output.param.txt"),header =T)
        beta.hat = paramter.file[,5]+paramter.file[,6]*paramter.file[,7]
        paramter.file$weight = beta.hat
        write.table(paramter.file[,c(2,8)],paste0("../",smp,"/",gene_name),col.names = F,row.names = F,quote = F)
      }
      system(paste0("rm -r *",gene_name,"*"))
    }
  } else {next}
}
