########################   Surv   ######################## 

library(data.table)
library(tidyverse)
library(survival)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "2.Prognosis/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
data_clinical <- read.xlsx("./data/GEO/GSE26899/Clinical.GSE26899.xlsx")

# 读入目的基因
interest_genes <- read.csv("./res/1.degs/table/degs.Tstr.wilcox.csv")
interest_genes <- interest_genes[interest_genes$padj <= 0.05,]
interest_genes <- interest_genes[abs(interest_genes$FoldChange) >= 1,]
interest_genes <- interest_genes[,1];interest_genes

########################   下载数据   ########################
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*20) # 不修改会报错
# download
eSet <- geoChina(gse = "GSE26899")
eSet <- eSet[[1]]
data_exprs <- as.data.frame(eSet@assayData$exprs)
data_anno <- eSet@phenoData@data
head(data_exprs[,1:4])

### probe anno
data_probe <- fread("./data/GEO/GSE26899/GPL6947-13512.txt",data.table = F)
data_probe_sub <- data_probe[,c("ID","Entrez_Gene_ID")]

data_exprs2 <- rownames_to_column(data_exprs,"ID")
data_merge <- merge(data_probe_sub, data_exprs2, "ID")
head(data_merge[,1:4])

# id convert
gene_id <- id_convert(data_merge[,2],"ENTREZID","SYMBOL")
data_merge <- merge(gene_id, data_merge, by.x= "ENTREZID",by.y ="Entrez_Gene_ID")
data_merge <- data_merge[!duplicated(data_merge$SYMBOL),]
data_merge <- data.frame(data_merge[,-c(1,2,3)], row.names = data_merge[,2])
range(data_merge)
# data_merge <- log2(data_merge + 1)
data_exprs <- data_merge
save(data_exprs, data_anno, file = "./data/GEO/GSE26899/GSE26899.RData")

########################   ssGSEA   ########################
# 转为矩阵
data_pre <- as.matrix(data_exprs)
# ssGSEA
genelist <- list(Tstr.CD8 = interest_genes)
data_ES <- GSVA::gsva(data_pre, genelist, method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
# 计算zscore
data_ES_zscore <- scale(data_ES) %>% data.frame()

#### merge
colnames(data_anno)
data_anno_sub <- data_anno[,c("title","geo_accession","characteristics_ch1.1")]
data_anno_sub <- data_anno_sub[grepl("tumor",data_anno_sub$title),]
data_anno_sub$characteristics_ch1.1 <- gsub("patient: ","",data_anno_sub$characteristics_ch1.1)
data_pre <- rownames_to_column(data_ES, "geo_accession")
# merge
data_merge <- merge(data_pre, data_anno_sub, "geo_accession")
colnames(data_clinical)
colnames(data_merge)[4] <- "Patients_ID"
data_merge <- merge(data_merge, data_clinical, "Patients_ID")

########################      surv        ########################     
colnames(data_merge)
str(data_merge)
melcli_group <- data.frame(data_merge[,c(1,3,12,13)], row.names = data_merge[,1])
colnames(melcli_group) <- c("sample","group","status","months")
melcli_group <- melcli_group[melcli_group$months >= 1,]

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(melcli_group, 
                         time = "months", 
                         event = "status", 
                         variables = c("group")
)

# 2. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 3. KM plot
fit <- survfit(Surv(months, status) ~ group, data = res.cat)
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Surv.",".","Tstr.CD8",".GSE26899",".pdf");outfile_tmp
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
           break.time.by = 36,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()




