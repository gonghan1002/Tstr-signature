########################   差异分析   ######################## 
#### 先分析tumor与normal的差异

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "1.degs/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
data_exprs <- readRDS("./data/TCGA-STAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# interest gene
interest_genes <- read.xlsx("./list/Tstr.xlsx")
interest_genes <- c(unlist(interest_genes),use.names = F)
interest_genes <- interest_genes[!duplicated(interest_genes)]

### 区分tumor normal
### group 
# 观察样本类型
table(substring(colnames(data_exprs),first = 13,last = 16))
# 不要01B B是石蜡样本
tumor_idx <- grep("01A",substring(colnames(data_exprs),first = 13,last = 16))
normol_idx <- grep("11A",substring(colnames(data_exprs),first = 13,last = 16))
# tumor 在前， normal 在后
data_exprs <- data_exprs[,c(tumor_idx, normol_idx)]

### 构建注释组
table(substring(colnames(data_exprs),first = 13,last = 16))
condition <- c(rep("tumor", 373),rep("normal",32))
data_anno <- data.frame(condition, row.names = colnames(data_exprs))

### degs
data_pre <- log2(data_exprs + 1)
res.degs <- diff_limma(inputData = data_pre,colData = data_anno,
                       contrast_name = c("tumor","normal"),
                       log2trans = F,adjust.method = "BH",
                       pval_cutoff = NULL,log2fc_cutoff = NULL)

# 加标签列
res.degs$threshold = factor(ifelse(res.degs$adj.P.Val < 0.05 & abs(res.degs$logFC) >= 0.585, 
                                   ifelse(res.degs$logFC >= 0.585,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
table(res.degs$threshold)

# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","degs.","tumor_normal", ".csv");outfile_tmp
write.csv(res.degs,file = outfile_tmp)

# sub
data_sub <- res.degs[rownames(res.degs) %in% interest_genes,]
table(data_sub$threshold)
# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","degs.Tstr", ".csv");outfile_tmp
write.csv(data_sub,file = outfile_tmp)

###########################     volcano    ###########################  
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Volcano.","tumor_normal",".pdf");outfile_tmp
# plot
volcano_res <- volcano_plot(input=res.degs, output=outfile_tmp, 
                            log2fc_cutoff=1, pval_cutoff=0.05,
                            title=NULL,padj_cutoff=0.05, 
                            genes_highlight=NULL, 
                            color_map=c(Down="#4b84b3", Not="grey", Up="#f89f68"), 
                            width=7, height=7, dpi=300, device="pdf")

