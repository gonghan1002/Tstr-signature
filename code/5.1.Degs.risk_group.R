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
proj_name <- "5.Biology/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
data_exprs <- readRDS("./data/TCGA-STAD.TPM.Tumor.RDS")
data_exprs <- log2(data_exprs + 1)
head(data_exprs[,1:4])

# group
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)

# intesect
commone_col <- intersect(data_group$sample, colnames(data_exprs))
data_exprs <- data_exprs[,commone_col]
data_group <- data_group[data_group$sample %in% commone_col,]

### 构建注释组
data_anno <- data.frame(group = data_group$risk, row.names = colnames(data_exprs))

### degs
res.degs <- diff_limma(inputData = data_exprs,colData = data_anno,
                       contrast_name = c("high","low"),
                       log2trans = F,adjust.method = "BH",
                       pval_cutoff = NULL,log2fc_cutoff = NULL)

# 加标签列
res.degs$threshold = factor(ifelse(res.degs$adj.P.Val < 0.05 & abs(res.degs$logFC) >= 1, 
                                   ifelse(res.degs$logFC >= 1,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
table(res.degs$threshold)

# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","degs.", "fc1", ".csv");outfile_tmp
write.csv(res.degs,file = outfile_tmp)

###########################     volcano    ###########################  
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Volcano.", "fc1",".pdf");outfile_tmp
# plot
volcano_res <- volcano_plot(input=res.degs, output=outfile_tmp, 
                            log2fc_cutoff=1, pval_cutoff=0.05,
                            title=NULL,padj_cutoff=0.05, 
                            genes_highlight=NULL, 
                            color_map=c(Down="#4b84b3", Not="grey", Up="#f89f68"), 
                            width=7, height=7, dpi=300, device="pdf")



