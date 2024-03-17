########################   比较突变图谱 + TMB ######################## 
#### 高TMB预后好 与聚类结果相反 
#### 将cluster1和2分开运行 保存结果 之后比较


library(data.table)
library(tidyverse)
library(maftools)
library(openxlsx)
library(survival)
library(survminer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "5.Biology/"

### load data
load("E:/BaiduSyncdisk/r_project/data/TCGA_data/SNV/TCGA-STAD.SNV.Rdata")
head(snv[,1:4])
snv <- snv[,-1]

### 读取分组数据
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
group <- data_group[,c(1,5)]
colnames(group) <- c("Tumor_Sample_Barcode","group")
group$Tumor_Sample_Barcode <- substring(group$Tumor_Sample_Barcode,1,12)
group$Tumor_Sample_Barcode <- gsub("[.]","-",group$Tumor_Sample_Barcode)

# 构建注释矩阵
group <- group[!duplicated(group$Tumor_Sample_Barcode),]
maf_group <- data.frame(group,row.names = group[,1])

###############  根据分组信息读入maf文件   ########################
# 读取maf文件
snv$Tumor_Sample_Barcode <- substring(snv$Tumor_Sample_Barcode,1,12)
maf_data <- read.maf(maf = snv,clinicalData = maf_group)

# set color
my_color <- c("#f89f68","#4b84b3")
names(my_color) = c("high", "low")
group_color = list(group = my_color)

#### plot 
outfile_tmp <- paste0(res_home,proj_name,"pdf/","SNV",".pdf");outfile_tmp
# save
pdf(file = outfile_tmp, onefile = F, height = 6, width = 6)
oncoplot(maf = maf_data,
         clinicalFeatures = "group", sortByAnnotation = T,
         draw_titv = F,
         annotationColor = group_color,
         bgCol = "#F4F4F4"
)
dev.off()

# 提取子集
pres_idx <- group$Tumor_Sample_Barcode[group$group %in% "high"]
maf_sub <-  subsetMaf(maf = maf_data, tsb = pres_idx)

#### plot 
outfile_tmp <- paste0(res_home,proj_name,"pdf/","SNV.","risk_high",".pdf");outfile_tmp
# save
pdf(file = outfile_tmp, onefile = F, height = 6, width = 6)
oncoplot(maf = maf_sub,
         clinicalFeatures = "group", sortByAnnotation = T,
         draw_titv = F,
         annotationColor = group_color,
         bgCol = "#F4F4F4"
)
dev.off()


##########################   TMB比较   ########################## 
tmb_data <- tmb(maf_data)
dev.off()

# 读取risk分组
data_pre <- tmb_data
data_merge <- merge(group,data_pre,by="Tumor_Sample_Barcode")
colnames(data_merge)

### 作图
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","TMB", ".pdf");outfile_tmp
pdf(file = outfile_tmp, width = 3, height = 4)
ggplot(data = data_merge,
       aes(y = total_perMB_log,
           x = group))+
  geom_boxplot(alpha = 1,
               fill = c("#f89f68","#4b84b3"))+
  stat_compare_means(method = "t.test",label = "p.format") + 
  theme_bw()+ theme_classic() +
  ylab('logTMB') +
  xlab('')
dev.off()

##########################   TMB 与 风险得分相关性   ########################## 
library(hrbrthemes)
library(ggplot2)

# merge
data_pre <- data_merge
colnames(data_group)
group <- data_group[,c(1,4)]
colnames(group) <- c("Tumor_Sample_Barcode","score")
group$Tumor_Sample_Barcode <- substring(group$Tumor_Sample_Barcode,1,12)
group$Tumor_Sample_Barcode <- gsub("[.]","-",group$Tumor_Sample_Barcode)
data_merge2 <- merge(data_pre,group,"Tumor_Sample_Barcode")

### plot
cor_res <- scatter_ggplot(inputArr=data_merge2, outfile=NULL,
                          xcol="score", ycol= "total_perMB_log",
                          title = NULL, color="#4b84b3",color_manual=NULL, col_pal_name="lancet",
                          facet.by=NULL, Marginal_Density=FALSE)

# save
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Cor.","TMB_risk", ".pdf");outfile_tmp
pdf(file=outfile_tmp,width=4,height=4)
cor_res + 
  ylab('logTMB') +
  xlab('Risk Score')
dev.off()

