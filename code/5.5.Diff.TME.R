########################   比较TME   ######################## 
#### 

library(tidyverse)
library(data.table)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "5.Biology/"

# 读入estimate res
ESTIMATE_res <- read.table("./data/ESTIMATE.TCGA-STAD.txt",header = T,row.names = 1)


# 读入MCPcounter res
MCPcounter_res <-  read.csv("./data/MCPcounter.TCGA-STAD.csv", row.names = 1)

# 读取分组
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
data_group <- data_group[,c(1,5)]
colnames(data_group) <- c("id","group")

################################   比较ESTIMATE  ################################ 
data_pre <- ESTIMATE_res
data_pre <- rownames_to_column(data_pre,var="id")
### merge
data_pre$id <- gsub("-",".",data_pre$id)
group <- data_group
group$id <- substr(group$id,1,15)
data_pre_merge <- merge(group, data_pre, by="id")

### 转化数据类型
TME.cells <- colnames(data_pre_merge)[3:5];TME.cells
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(data_pre_merge) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = data_pre_merge$id,
    CellType = TME.cells[i],
    group = data_pre_merge$group,
    Composition = data_pre_merge[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

# 调一下risk顺序
plot.info$group <- factor(plot.info$group,levels = c("low","high"))
# 修改一下标签名字
plot.info$CellType <- gsub("Score","",plot.info$CellType)

# boxplot group by risk
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","diff_ESTIMATE", ".pdf")
pdf(file = outfile_tmp,width = 6,height = 4)
ggboxplot(data = plot.info,x = "CellType",y = "Composition",
          palette = c("#4b84b3", "#f89f68"), fill = "group",
          xlab = "",ylab = "scores") +
  stat_compare_means(label = "p.format",method = "t.test",
                     aes(group=group),
                     hide.ns = F) + 
  theme_base() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()


################################   比较MCPcounter  ################################ 
data_pre <- MCPcounter_res
data_pre <- as.data.frame(t(data_pre))
data_pre <- rownames_to_column(data_pre,"id")
colnames(data_pre)

### merge
group <- data_group
data_pre_merge <- merge(group, data_pre, by="id")
colnames(data_pre_merge)

### 转化数据类型
TME.cells <- colnames(data_pre_merge)[3:12];TME.cells
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(data_pre_merge) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = data_pre_merge$id,
    CellType = TME.cells[i],
    group = data_pre_merge$group,
    Composition = data_pre_merge[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

# 调一下risk顺序
plot.info$group <- factor(plot.info$group,levels = c("low","high"))
plot.info$CellType <- gsub("_"," ",plot.info$CellType)

# boxplot group by risk
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","diff_MCPcounter", ".pdf")
pdf(file = outfile_tmp,width = 9,height = 6)
ggboxplot(data = plot.info,x = "CellType",y = "Composition",
          palette = c("#4b84b3", "#f89f68"), fill = "group",
          xlab = "",ylab = "scores") +
  stat_compare_means(label = "p.signif",method = "wilcox.test",
                     aes(group=group),
                     hide.ns = F) + 
  theme_base() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()

