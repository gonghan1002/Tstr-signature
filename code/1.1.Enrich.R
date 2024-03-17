########################   差异分析   ######################## 
#### 先分析tumor与normal的差异

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)
library(GSVA)
library(enrichplot)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "1.degs/"

# read data
res.degs <- read.csv("./res/1.degs/table/degs.tumor_normal.csv")

# interest gene
interest_genes <- read.xlsx("./list/Stress sponse.xlsx")
colnames(interest_genes) <- "gene"
interest_genes$ont <- "Stress sponse"
interest_genes <- interest_genes[,c("ont", "gene")]

##############################  富集分析 -- gsea  ##############################
### 构建list
gene_list <- res.degs$logFC
names(gene_list) <- res.degs[,1]
# 排序  
gene_list <- sort(gene_list,decreasing = T) # 一定要降序排列

### GSEA
kk <- GSEA(gene_list, TERM2GENE = interest_genes, 
           exponent=1, minGSSize=5,  
           pvalueCutoff=0.05, pAdjustMethod="BH",
           seed=F, by="fgsea")
kk.df <- as.data.frame(kk)
# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","GSEA.", "stress", ".csv");outfile_tmp
write.csv(kk.df,file = outfile_tmp)

# plot
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","GSEA.", "stress.2", ".pdf");outfile_tmp
pdf(outfile_tmp,onefile = F, width = 6,height = 6)
gseaplot2(kk, 1, color = "#f89f68",pvalue_table = T, subplots = 1:3, base_size = 10)
dev.off()


