########################   clinical feature   ######################## 

library(data.table)
library(tidyverse)
library(survival)
library(openxlsx)
library(ggplot2)
library(ggpubr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "2.Prognosis/"

# load data
load("./res/1.degs/RData/GSVA.Tstr.RData")

# read data
data_clinical <- read.csv("./data/TCGA-STAD.clinical_info.csv")

# group 
data_group <- read.csv("./res/1.degs/table/Group.Tstr.CD8.max.csv")

####  merge
data_pre <- rownames_to_column(data_ES_zscore,"sample")
colnames(data_group)[1] <- "sample"
data_merge <- merge(data_pre, data_clinical, "sample")

################   age   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,4)]
data_plot$age <- case_when(data_plot$age >= 67 ~ ">= 65",
                           data_plot$age < 67 ~ "< 65")
colnames(data_plot) <- c("Tstr.CD8", "group")
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Age",".pdf");outfile_tmp
pdf(outfile_tmp,width = 3,height = 4,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#f89f68", "#4b84b3")) + 
  stat_compare_means(method = "t.test",label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()


################   grade   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,5)]
colnames(data_plot) <- c("Tstr.CD8", "group")
data_plot$group <- factor(data_plot$group, levels = c("G1","G2","G3"))
data_plot <- na.omit(data_plot)

# diff
my_lists <- list(c("G1", "G2"), c("G1", "G3"), c("G2", "G3"))

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Grade",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4DA75D", "#4b84b3","#f89f68")) + 
  stat_compare_means(method = "t.test",comparisons = my_lists,
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()

################   new_tumor   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,6)]
colnames(data_plot) <- c("Tstr.CD8", "group")
str(data_plot)

# group
data_plot$group <- factor(data_plot$group, levels = c("NO","YES"))
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Relapse",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4b84b3","#f89f68")) + 
  stat_compare_means(method = "t.test",
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()

################   M   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,7)]
colnames(data_plot) <- c("Tstr.CD8", "group")
str(data_plot)

# group
data_plot$group <- factor(data_plot$group, levels = c("M0","M1"))
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Stage_M",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4b84b3","#f89f68")) + 
  stat_compare_means(method = "t.test",
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()

################   N   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,8)]
colnames(data_plot) <- c("Tstr.CD8", "group")
table(data_plot$group)

# group
data_plot$group <- factor(data_plot$group, levels = c("N0","N1","N2","N3"))
data_plot <- na.omit(data_plot)

# diff
my_lists <- list(c("N0", "N1"), c("N0", "N2"), c("N0", "N3"),c("N1","N2"),c("N1","N3"),c("N2","N3"))

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Stage_N",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#A3E7D3","#4DA75D","#4b84b3","#f89f68")) + 
  stat_compare_means(method = "t.test",comparisons = my_lists,
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()

################   T   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,9)]
colnames(data_plot) <- c("Tstr.CD8", "group")
table(data_plot$group)

# group
data_plot$group <- case_when(data_plot$group == "T1" ~ "T1",
                             data_plot$group == "T2" ~ "T2-4",
                             data_plot$group == "T3" ~ "T2-4",
                             data_plot$group == "T4" ~ "T2-4")
data_plot$group <- factor(data_plot$group, levels = c("T1","T2-4"))
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Stage_T",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4b84b3","#f89f68")) + 
  stat_compare_means(method = "wilcox.test",
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()


################   gender   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,10)]
colnames(data_plot) <- c("Tstr.CD8", "group")
table(data_plot$group)

# group
data_plot$group <- factor(data_plot$group, levels = c("female","male"))
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Gender",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4b84b3","#f89f68")) + 
  stat_compare_means(method = "wilcox.test",
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()

################   stage   ################   
data_plot <- data_merge
colnames(data_plot)
data_plot <- data_plot[,c(2,13)]
colnames(data_plot) <- c("Tstr.CD8", "group")
table(data_plot$group)

# group
data_plot$group <- case_when(data_plot$group == "Stage I" ~ "Stage I",
                             data_plot$group == "Stage II" ~ "Stage II-IV",
                             data_plot$group == "Stage III" ~ "Stage II-IV",
                             data_plot$group == "Stage IV" ~ "Stage II-IV")
data_plot$group <- factor(data_plot$group, levels = c("Stage I","Stage II-IV"))
data_plot <- na.omit(data_plot)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Boxplot.","Tstr_CD8.","Stage",".pdf");outfile_tmp
pdf(outfile_tmp,width = 4,height = 5,onefile = F)
ggplot(data_plot,aes(x= group, y= Tstr.CD8 )) + 
  geom_boxplot(aes(fill = group), fill = c("#4b84b3","#f89f68")) + 
  stat_compare_means(method = "t.test",
                     label = "p.format") +
  xlab('')+
  ylab('Tstr.CD8') + 
  theme_classic()
dev.off()






