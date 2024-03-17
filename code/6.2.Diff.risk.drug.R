##############################  读入oncoPredict结果 计算相关性 ##############################  

library(tidyverse)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(rstatix)
library(stringr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "7.Drug/"

# load data
data_drug <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
head(data_drug[,1:4])
colnames(data_drug)[1] <- "sample"

### 定义感兴趣的drug
interest_drug <- c("Docetaxel","Oxaliplatin","Paclitaxel","Lapatinib")

# 读取分组数据
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
group <- data_group[,c(1,4,5)]
group_sub <- group

#### 合并数据
data_pre <- data_drug
data_merge <- merge(group_sub, data_pre, by = "sample")
data_merge <- data_merge[!duplicated(data_merge$sample),]
data_merge <- data.frame(data_merge[,-c(1,3)], row.names = data_merge[,1])

##############################  合并数据  ##############################  
data_pre <- data_merge
#### 修改药物名字
colnames(data_pre)[1] <- "RiskScore"
colnames(data_pre) <- str_split(colnames(data_pre), "_",simplify = T)[,1]
data_pre$RiskScore <- log2(data_pre$RiskScore + 1)

data_pre <- data_pre[,c("RiskScore",interest_drug)]
# data_pre <- as.data.frame(scale(data_pre))

# calculate
cor_res_list <- list()
for(drug in interest_drug){
  cor_res_list[[drug]] <- scatter_ggplot(inputArr=data_pre, outfile=NULL,
                                         xcol="RiskScore", ycol= drug,
                                         title = NULL, color="#72BE64",
                                         color_manual=NULL, col_pal_name="lancet",
                                         facet.by=NULL, Marginal_Density=FALSE)
}

# 整合全部癌症的结果到一张图
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Cor.","risk_drug.clinical" , ".pdf");outfile_tmp
multiplot_ggarrange(ggobj_list=cor_res_list, outfile=outfile_tmp, labels=NULL, 
                    ncol=4, nrow=1,
                    legend="bottom", common.legend=TRUE, width=16, height=4)

cor_res_list$IGF1R
dev.off()









