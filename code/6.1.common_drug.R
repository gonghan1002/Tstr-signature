##############################  读入oncoPredict结果 比较auc ##############################  

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

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# load data
data_drug <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
head(data_drug[,1:4])
colnames(data_drug)[1] <- "sample"

# 读取分组数据
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
group <- data_group[,c(1,4,5)]
group_sub <- group

#### 合并数据
data_pre <- data_drug
data_merge <- merge(group_sub, data_pre, by = "sample")
data_merge <- data_merge[!duplicated(data_merge$sample),]
data_merge <- data.frame(data_merge[,-c(1,2)], row.names = data_merge[,1])

##############################  t test  ############################## 
head(data_merge[,1:4])
res_test <- melt(data_merge,id.vars=c("risk"),variable.name = "drug") %>%
  group_by(drug) %>% t_test(value ~ risk) %>% 
  adjust_pvalue(method = "fdr") %>% add_significance("p.adj")

# 统计一下失调基因数目
table(res_test$p<0.05)
# 加threshold列
res_test$threshold = factor(ifelse(res_test$p<0.05 & abs(res_test$statistic)>=0, 
                                   ifelse(res_test$statistic>=0,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
# save
outfile_tmp <-  paste0(res_home,proj_name,"table/","Drug_auc.t_test",".csv");outfile_tmp
write.csv(res_test,file = outfile_tmp,row.names = F)

############################## boxplot ##############################
### 定义感兴趣的drug
interest_drug <- c("Docetaxel","Oxaliplatin","Paclitaxel","Lapatinib")

### 准备绘图数据
data_boxplot <- data_merge
colnames(data_boxplot) <- str_split(colnames(data_boxplot), "_",simplify = T)[,1]
# 转为小写
# colnames(data_boxplot) <- tolower(colnames(data_boxplot))
# extract
common_col <- intersect(interest_drug, colnames(data_boxplot));common_col
data_boxplot <- data_boxplot[,c("risk",common_col)]
head(data_boxplot[,1:4])

# 表达值log2转换
data_boxplot[,common_col] <- log2(data_boxplot[,common_col] + 1)

#### boxplot 分开画
data_pre <- data_boxplot
colnames(data_pre)
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/",
                      "Boxplot.","drug_auc.","Lapatinib",".pdf");outfile_tmp
pdf(outfile_tmp,width = 3,height = 4,onefile = F)
ggplot(data_pre,aes(x=risk,y= Lapatinib )) + 
  geom_boxplot(aes(fill=risk), fill = c("#FFA500", "#72BE64")) + 
  stat_compare_means(method = "wilcox.test",label = "p.format") +
  xlab('risk')+
  ylab('Lapatinib (susceptibility)') + 
  theme_classic()
dev.off()


