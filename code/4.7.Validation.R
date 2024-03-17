
library(stringr)
library(data.table)
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(survminer)
library(survival)
library(timeROC)
library(GEOmirror)


rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "4.Signature/"


# 读入目的基因
gene_signature <- read.csv("./res/4.Signature/table/LASSO.gene_coef.2.csv")
interest_genes <- gene_signature[,1];interest_genes

########################   下载数据   ########################
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*20) # 不修改会报错
# download
eSet <- geoChina(gse = "GSE15459")
eSet <- eSet[[1]]
data_exprs <- as.data.frame(eSet@assayData$exprs)
data_anno <- eSet@phenoData@data
head(data_exprs[,1:4])

# 数据处理
data_exprs2 <- rownames_to_column(data_exprs,var = "ID_REF")
data_exprs2 <-  getExprs_probes2gene(input = data_exprs2,platform = "hgu133plus2")
head(data_exprs2[,1:4]);range(data_exprs2)
data_exprs2 <- log2(data_exprs2 + 1)

# data anno
data_anno <- read.xlsx("./data/GEO/GSE15459/clinical_info.xlsx")
colnames(data_anno)
data_anno_sub <- data_anno[,c(1,4:8,9,10)]

# save
save(data_exprs2, data_anno, file = "./data/GEO/GSE15459/GSE15459.RData")

### merge
data_pre <- data_exprs2[interest_genes,]
data_pre <- as.data.frame(t(data_pre))
data_pre <- rownames_to_column(data_pre,"GSM.ID")
data_merge <- merge(data_pre, data_anno_sub, "GSM.ID")

########################   risk score   ########################
gene_signature
data_pre <- data_merge
{
  a1 <- 0.06503219*data_pre$SERPINE1
  a2 <- 0.04111661*data_pre$RGS2
  a3 <- 0.01377776*data_pre$PDGFRL
  a4 <- 0.02484820*data_pre$STC1
  a5 <- 0.10362954*data_pre$C5orf46
  a6 <- 0.02717206*data_pre$CST2
  a7 <- 0.03715293*data_pre$GPX3
  a8 <- 0.07365494*data_pre$SNCG
}
data_pre$risk_score <- a1+a2+a3+a4+a5+a6+a7+a8


########################      surv        ########################     
colnames(data_pre)
str(data_pre)
melcli_group <- data_pre
melcli_group <- melcli_group[,c(1,15,16,17)]
colnames(melcli_group)[-1] <- c("months","status","risk_score")

melcli_group$group <- ifelse(melcli_group$risk_score <= median(melcli_group$risk_score),"low","high")

### surv
fit <- survfit(Surv(months, status) ~ group, data = melcli_group)
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Surv.","median",".GSE15459",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, width = 5,height = 6)   
ggsurvplot(fit, 
           data=melcli_group,
           conf.int=F,
           pval=T,
           pval.size=4,
           risk.table=T,
           legend.labs=c("high", "low"), #, "Cluster4"
           legend.title="group",
           xlab="Time(months)",
           ylab="Status",
           break.time.by = 36,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()
