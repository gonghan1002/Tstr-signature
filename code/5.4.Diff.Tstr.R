########################   GSVA   ######################## 
#### 

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)
library(ggpubr)


rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "5.Biology/"


# load data
load("./res/1.degs/RData/GSVA.Tstr.RData")

# 读取分组
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
group <- data_group[,c(1,5)]

### merge
data_pre <- rownames_to_column(data_ES_zscore,"sample")
data_merge <- merge(data_pre, group, "sample")

##############################    plot   ##############################
data_plot <- data_merge

# order
data_plot$risk <- factor(data_plot$risk,levels = c("low","high"))
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr.CD8.","risk",".pdf");outfile_tmp
pdf(outfile_tmp,width = 3,height = 4,onefile = F)
ggplot(data_plot,aes(x= risk, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill=risk), fill = c("#4b84b3", "#f89f68")) + 
  stat_compare_means(method = "wilcox.test",label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()





