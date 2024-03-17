########################   Surv   ######################## 

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(timeROC)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "1.degs/"

# load data
load("./res/1.degs/RData/GSVA.Tstr.RData")

# surv
data_surv <- fread("./data/TCGA-STAD.survival.tsv",data.table = F)


#### merge
data_pre <- rownames_to_column(data_ES_zscore,"sample")
data_surv_sub <- data_surv[,c(-3)]
data_surv_sub$sample <- gsub("-",".",data_surv_sub$sample)
# merge
data_merge <- merge(data_surv_sub,data_pre,by="sample")

########################      surv        ########################     
colnames(data_merge)
melcli_group <- data.frame(data_merge[,c(2,3,4)], row.names = data_merge[,1])
melcli_group$months <- melcli_group$OS.time/30
melcli_group <- melcli_group[melcli_group$months >= 1,]

res.cut <- surv_cutpoint(melcli_group,
                         time = "months", 
                         event = "OS",
                         variables = c("Tstr.CD8") 
)

res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(months, OS) ~ Tstr.CD8, data = res.cat)
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Sruv.","OS",".","Tstr.CD8",".max",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, width = 5,height = 6)   
ggsurvplot(fit, 
           data=res.cat,
           conf.int=F,
           pval=T,
           pval.size=4,
           risk.table=T,
           legend.labs=c("high CD8Tstr", "low CD8Tstr"), #, "Cluster4"
           legend.title="group",
           xlab="Time(months)",
           ylab="Status",
           break.time.by = 24,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()


### save
outfile_tmp <- paste0(res_home,proj_name,"/table/","Group.","Tstr.CD8.","max",".csv");outfile_tmp
write.csv(res.cat,file = outfile_tmp)




