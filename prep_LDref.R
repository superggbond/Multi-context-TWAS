## prep LD ref file for S-PrediXcan analysis
library(tidyverse)
library(RSQLite)
library(snpStats)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
ld_chr = args[1]

# load db file
db <- dbConnect(dbDriver("SQLite"), dbname = "/net/psoriasis/psorgenom/haihan/multi_TWAS/db_files/IFNa.db")
db.weights <- dbReadTable(db,"weights")

pos = strsplit(db.weights$rsid,":")
chr = sapply(pos,"[[",1)
db.weights$chr = chr
db.select = db.weights[db.weights$chr==ld_chr,]

# load 1KG as LD reference panel
ldref = read.plink(paste0("/net/assembly/alextsoi/db/1000g/phase3.v5.reduced.EUR/chr",ld_chr,"_EUR.bed"))
ld.mat = as(ldref$genotypes, "numeric")
snp.names = colnames(ld.mat)
snp.names = str_split(snp.names, "\\:")
snp.names = sapply(snp.names,"[[",2)
snp.names = paste0(ld_chr,":",snp.names)
colnames(ld.mat) = snp.names
ld.mat = ld.mat[,!duplicated(colnames(ld.mat))]

cov.all = data.frame(matrix(ncol = 4, nrow = 0))
colnames(cov.all) = c("GENE","RSID1","RSID2","VALUE")

gene.names = unique(db.select$gene)
for (gene in gene.names){
  db.gene = db.select[db.select$gene==gene,]
  
  ld.mat.cov = ld.mat[,colnames(ld.mat)%in%db.gene$rsid]
  snp.mat.cov = cov(ld.mat.cov)
  snp.cov.uptri = data.frame(GENE=gene, RSID1=rownames(snp.mat.cov)[row(snp.mat.cov)[upper.tri(snp.mat.cov)]], 
                             RSID2=colnames(snp.mat.cov)[col(snp.mat.cov)[upper.tri(snp.mat.cov)]], 
                             VALUE=snp.mat.cov[upper.tri(snp.mat.cov)])
  snp.cov.diag = data.frame(GENE=gene, RSID1=rownames(snp.mat.cov), RSID2=colnames(snp.mat.cov), VALUE=diag(snp.mat.cov))
  cov.all = rbind(cov.all,snp.cov.diag,snp.cov.uptri)
}

fwrite(cov.all,paste0("/net/psoriasis/psorgenom/haihan/multi_TWAS/cov_files/chr",ld_chr,".txt.gz"),row.names = F,col.names = T,quote = F,sep = " ")

