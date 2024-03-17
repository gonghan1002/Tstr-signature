########################  整理数据 ######################## 
#### 

library(data.table)
library(tidyverse)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
getwd()
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/"

checkDir(paste0(res_home,"/code/"))
checkDir(paste0(res_home,"/data/"))
checkDir(paste0(res_home,"/list/"))
checkDir(paste0(res_home,"/res/"))

# read data
data_exprs <- fread("./data/TCGA-STAD.htseq_fpkm.tsv.gz", data.table = F)
head(data_exprs[,1:4])

# data anno
data_gene <- fread("E:/BaiduSyncdisk/r_project/data/TCGA_data/TCGA-mRNA.gene_anno.csv",data.table = F)
head(data_gene[,1:4])

### sub
data_gene_pre <- data_gene[data_gene$gene_type %in% "protein_coding",]
data_gene_pre <- separate(data_gene_pre,V1,into = c("V1"),sep="\\.")

data_exprs <- separate(data_exprs,Ensembl_ID,into = c("Ensembl_ID"),sep="\\.") 
data_exprs <- data_exprs[data_exprs[,1] %in% data_gene_pre$V1,]

### id convert
colnames(data_exprs)[1] <- "ENSEMBL"
gene_id <- id_convert(data_exprs[,1],"ENSEMBL","SYMBOL")
# merge 
data_exprs <- merge(gene_id, data_exprs, "ENSEMBL")
# duplicated
data_exprs <- data_exprs[!duplicated(data_exprs[,2]),]
data_exprs <- data.frame(data_exprs[,-c(1,2)], row.names = data_exprs[,2])
# save
saveRDS(data_exprs, file = "./data/TCGA-STAD.FPKM.RDS")

#### FPMK to TPM
### define function
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# return
data_pre <- 2^data_exprs - 1
range(data_pre)
data_tpm <- apply(data_pre,2,fpkmToTpm)
range(data_tpm)
colSums(data_tpm) # 检查
data_tpm <- as.data.frame(data_tpm)
# save
saveRDS(data_tpm, file = "./data/TCGA-STAD.TPM.RDS")

#### screen tumor
table(substring(colnames(data_tpm),first = 13,last = 16))#观察样本类型
tumor_idx <- grep("01A",substring(colnames(data_tpm),first = 13,last = 16))
normol_idx <- grep("11A",substring(colnames(data_tpm),first = 13,last = 16))
data_pre <- data_tpm[,tumor_idx]
head(data_pre[,1:4])
# save
saveRDS(data_pre, file = "./data/TCGA-STAD.TPM.Tumor.RDS")

