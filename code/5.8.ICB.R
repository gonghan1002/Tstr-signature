########################    ICB      ######################## 

library(openxlsx)
library(tidyverse)
library(survival)
library(data.table)
library(survminer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "5.Biology/"

# 读入目的基因
risk_model <- read.csv("./res/4.Signature/table/LASSO.gene_coef.2.csv")
interest_genes <- risk_model[,1];interest_genes

### 下载数据
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*6)

GSE176307 <- get_Matrix_GEO(GSE_ID = "GSE176307",method = "GEOquery")
GSE176307_anno <- GSE176307$phenoDat
GSE176307_data <- GSE176307$data_matrix

########## 注释矩阵处理 ########## 
colnames(GSE176307_anno)
GSE176307_anno_sub <- GSE176307_anno[,c("title","age:ch1","alive:ch1","d.stage.at.diagnosis:ch1",
                                        "gender:ch1","io.response:ch1","n.stage.at.diagnosis:ch1",
                                        "overall survival:ch1","pfs:ch1","t.stage.at.diagnosis:ch1",
                                        "tmb:ch1","tmb.interpretation:ch1")]
GSE176307_anno_sub$`alive:ch1` <- ifelse(GSE176307_anno_sub$`alive:ch1`=="No",0,1)
GSE176307_anno_sub$`overall survival:ch1` <- as.numeric(GSE176307_anno_sub$`overall survival:ch1`)
colnames(GSE176307_anno_sub)[-1] <- c("age","alive","stage_M","gender","response","stage_N",
                                      "OS_days","PFS_days","stage_T","TMB","TMB_classify")
GSE176307_anno_sub$title <- gsub("Patient sample ","",GSE176307_anno_sub$title)

##########  数据处理 ########## 
data_exprs <- read.csv("./data/GEO/GSE176307/GSE176307_BACI_log_trans_normalized_RNAseq.csv",
                       row.names = 1)
data_exprs_sub <- data_exprs[,interest_genes]
data_exprs_sub <- rownames_to_column(data_exprs_sub,var = "title")

########## 合并 ########## 
data_merge <- merge(GSE176307_anno_sub,data_exprs_sub,by="title")
colnames(data_merge)
data_merge$months <- data_merge$OS_days/30

########################    计算风险评分    ########################
risk_model
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

########################    生存曲线    ########################
colnames(data_pre)
melcli_group <- data_pre[,c(1,3,21,22)]
colnames(melcli_group)[2] <- "OS"

res.cut <- surv_cutpoint(melcli_group,
                         time = "months", 
                         event = "OS", 
                         variables = c("risk_score")
)

res.cat <- surv_categorize(res.cut)
head(res.cat)

### plot 
fit <- survfit(Surv(months, OS) ~ risk_score, data = res.cat)

# save
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Sruv","GSE176307", ".pdf")
pdf(file=outfile_tmp,onefile = FALSE, width = 5,height = 6)   
ggsurvplot(fit, 
           data=res.cat,
           conf.int=F,
           pval=T,
           pval.size=4,
           risk.table=F,
           legend.labs=c("high", "low"), #, "Cluster4"
           legend.title="group",
           xlab="Time(months)",
           ylab="Status",
           break.time.by = 12,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()
