## prep db file for each condition for S-PrediXcan analysis
library(tidyverse)
library(RSQLite)
library(snpStats)
library(data.table)

# check the data structure of the sample database
db <- dbConnect(dbDriver("SQLite"), dbname = "~/MetaXcan/software/data/CD4CLAP_24hst/CD4CLAP_24hst.db")
dbListTables(db)
db.construction <- dbReadTable(db,"construction")
db.extra <- dbReadTable(db,"extra")
db.sample_info <- dbReadTable(db,"sample_info")
db.weights <- dbReadTable(db,"weights")

# generate own customized db files
smps = c("IFNa","IFNg","IL13","IL17A","IL17ATNFa","IL4","NHEK","TNFa")
for (smp in smps){
  setwd(paste0("/net/psoriasis/psorgenom/haihan/multi_TWAS/bslmm_models_filtered/",smp))
  gene.names = list.files()
  
  # generate out.extra table
  out.extra = data.frame(matrix(ncol = 6, nrow = length(gene.names)))
  colnames(out.extra) = c("gene","genename","pred.perf.R2","n.snps.in.model","pred.perf.pval","pred.perf.qval")
  
  out.extra[,c(3,5,6)] = 0
  out.extra$genename = gene.names
  out.extra$gene = gene.names
  
  for (i in 1:length(gene.names)) {
    out.extra[i, "n.snps.in.model"] = length(readLines(gene.names[i]))
  }
  
  # generate out.construction table
  chr = 1:42
  cv.seed = 42
  out.construction = as.data.frame(cbind(chr,cv.seed))
  
  # generate out.sample_info table
  out.sample_info = as.data.frame(50)
  colnames(out.sample_info) = "n.samples"
  
  # generate out.weights table
  out.weights = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(out.weights) = c("rsid","gene","weight","ref_allele","eff_allele")
  for (gene in gene.names){
    weight.file = read.table(gene, header = F)
    weight.add = data.frame(matrix(ncol = 5, nrow = nrow(weight.file)))
    colnames(weight.add) = c("rsid","gene","weight","ref_allele","eff_allele")
    weight.add$gene = gene
    weight.add$rsid = weight.file$V1
    weight.add$weight = weight.file$V2
    
    ref_chr = strsplit(weight.file[1,1],":")[[1]][1]
    ref_file = fread(paste0("/net/psoriasis/psorgenom/haihan/multi_TWAS/GWAS_sum/psor/chr",ref_chr))
    ref_file = ref_file[ref_file$SNP%in%weight.add$rsid,]
    
    weight.add$ref_allele = ref_file$A1
    weight.add$eff_allele = ref_file$A2
    
    out.weights = rbind(out.weights, weight.add)
  }
  
  # output to a new db file
  out <- dbConnect(RSQLite::SQLite(), paste0("/net/psoriasis/psorgenom/haihan/multi_TWAS/db_files/",smp,".db"))
  
  dbWriteTable(out, "extra", out.extra)
  dbWriteTable(out, "weights", out.weights)
  dbWriteTable(out, "construction", out.construction)
  dbWriteTable(out, "sample_info", out.sample_info)
  dbDisconnect(out)
}
