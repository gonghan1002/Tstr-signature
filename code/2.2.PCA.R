########################    验证        ######################## 
####

library(FactoMineR)
library(factoextra)
library(missMDA)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(stringr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "1.degs/"

# read data
data_exprs <- readRDS("./data/TCGA-STAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])
range(data_exprs)
data_exprs <- log2(data_exprs + 1)

# group
data_group <- read.csv("./res/1.degs/table/Group.Tstr.CD8.max.csv")
colnames(data_group)[1] <- "sample"

# 读取list 这里是提取预后基因
interest_genes <- read.xlsx("./list/Tstr.xlsx")
interest_genes <- interest_genes$CD8_Tstr;interest_genes


##############################  PCA  ##############################
data_anno <- data.frame(group = data_group[,4],row.names = data_group$sample)
data_pre <- data_exprs[interest_genes,]
data_pre <- data_pre[,match(rownames(data_anno),colnames(data_pre))]

#### PCA
pca_data <- as.data.frame(t(data_pre))
pca_res <- PCA(pca_data, graph = FALSE) #imputePCA()

# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","PCA.","Tstr_CD8",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = F, width = 6,height = 6)   
fviz_pca_ind(pca_res,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = data_anno$group, # color by groups
             palette = c("#f89f68","#4b84b3"),
             addEllipses = F, # Concentration ellipses
             legend.title = "group",
             title = "")
dev.off()
