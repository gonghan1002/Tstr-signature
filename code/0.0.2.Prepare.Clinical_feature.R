########################   clinical feature   ######################## 

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

# read data
data_clinical <- fread("./data/TCGA-STAD.GDC_phenotype.tsv.gz",data.table = F)

### extract
colnames(data_clinical)
data_clinical_sub <- data_clinical[,c(1,6,34,37,41,42,43,71,83,89,91)]


### 处理临床数据
meta <- data_clinical_sub
colnames(meta) <- c("sample","age","grade","new_tumor","stage_M","stage_N","stage_T",
                    "gender","tumor_type","tumor_origin","stage")
meta$sample <- gsub("-",".",meta$sample)
### grade
idx_tmp <- meta$grade
table(idx_tmp)
idx_tmp <- gsub("GX", NA, idx_tmp)
table(idx_tmp)
meta$grade <- idx_tmp

### grade
idx_tmp <- meta$new_tumor
table(idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$new_tumor <- idx_tmp

### M
idx_tmp <- meta$stage_M
table(idx_tmp)
idx_tmp <- gsub("MX", NA, idx_tmp)
table(idx_tmp)
meta$stage_M <- idx_tmp

### N
idx_tmp <- meta$stage_N
table(idx_tmp)
idx_tmp <- gsub("N3[a-z]", "N3", idx_tmp)
idx_tmp <- gsub("NX", NA, idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$stage_N <- idx_tmp

### T
idx_tmp <- meta$stage_T
table(idx_tmp)
idx_tmp <- gsub("T1[a-z]", "T1", idx_tmp)
idx_tmp <- gsub("T2[a-z]", "T2", idx_tmp)
idx_tmp <- gsub("T4[a-z]", "T4", idx_tmp)
idx_tmp <- gsub("TX", NA, idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$stage_T <- idx_tmp

### stage
idx_tmp <- meta$stage
table(idx_tmp)
idx_tmp <- gsub("stage iii[abc]", "Stage III", 
                gsub("stage ii[abc]", "Stage II",
                     gsub("stage i[abc]", "Stage I", idx_tmp)))
idx_tmp[idx_tmp == "not reported"] = NA
idx_tmp[idx_tmp == "stage iv"] = "Stage IV"
idx_tmp[idx_tmp == "stage iii"] = "Stage III"
idx_tmp[idx_tmp == "stage ii"] = "Stage II"
idx_tmp[idx_tmp == "stage i"] = "Stage I"
table(idx_tmp)
meta$stage <- idx_tmp

# save
write.csv(meta,"./data/TCGA-STAD.clinical_info.csv",row.names = F)


