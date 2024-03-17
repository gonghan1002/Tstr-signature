########################    LASSO回归建模    ######################## 
#### LASSO建模 只剩下6个基因 不需要再cox建模

library(openxlsx)
library(tidyverse)
library(data.table)
library(glmnet)
library(survival)
library(stringr)
library(broom)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "4.Signature/"

# 读取数据
data_exprs <- readRDS("./data/TCGA-STAD.TPM.RDS")
data_exprs <- log2(data_exprs + 1)
head(data_exprs[,1:4])

# 读入临床信息
survival_data <- fread("./data/TCGA-STAD.survival.tsv",data.table = F)

# 读取list
prog_genes <- read.csv("./res/4.Signature/table/Cox_univ.hub_genes.csv",header = T)
interest_genes <- prog_genes[,1];interest_genes

##############################   合并临床数据  ##############################
data_pre <- data_exprs[interest_genes,]
data_pre <- as.data.frame(t(data_pre))
colnames(survival_data)
# merge
data_surv <- survival_data[,c(1,2,4)]
data_surv$sample <- gsub("-",".",data_surv$sample)
data_pre <- data.frame(sample = rownames(data_pre),data_pre)
data_merge <- merge(data_surv, data_pre, by="sample")

# 去掉生存时间少于1个月的患者
data_merge <- data.frame(months = data_merge$OS.time / 30,
                         data_merge)
# data_merge <- data_merge[data_merge$months >= 1,]

##############################   LASSO  ##############################
# 去重
data_pre <- data_merge[!duplicated(data_merge$sample),]
data_pre <- data.frame(data_pre, row.names = data_pre[,2])
data_pre <- na.omit(data_pre)
data_pre <- data_pre[,-c(2,4)]
colnames(data_pre)

## 注意一个大坑，这里的x，y必须为矩阵
colnames(data_pre)
x <- as.matrix(data_pre[,c(3:ncol(data_pre))])
y <- data.matrix(Surv(data_pre$months,data_pre$OS))

# LASSO
fit <- glmnet(x, y, family = "cox", maxit = 1000000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000000)

#################    plot    #################
# extract data
tidy_df <- broom::tidy(fit)
tidy_cvdf <- broom::tidy(cvfit)

# set color
mypalette <- c(brewer.pal(11,"BrBG"),brewer.pal(11,"PiYG"),brewer.pal(11,"PRGn"),
               brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"),brewer.pal(11,"RdGy"),
               brewer.pal(11,"RdYlBu"),brewer.pal(11,"RdYlGn"),brewer.pal(11,"Spectral"),
               brewer.pal(8,"Accent"),brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"),
               brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Set1"),
               brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))
### fit
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","LASSO.","lambda_coef",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 13,height = 6)
ggplot(tidy_df, aes(lambda, estimate, group = term, color = term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+
  scale_color_manual(name="variable",values = mypalette)+
  theme_bw()
dev.off()

### cvfit
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","LASSO.","lambda_cvfit.2",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 8,height = 8)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

#################    计算系数   #################    
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)

# 保留三位有效小数
coef.df <- as.data.frame(geneCoef)
coef.df$Coef <- as.numeric(coef.df$Coef)
coef.df[,2] <- signif(coef.df[,2],3)

# save
geneCoef[,1] <- gsub("[.]","-",geneCoef[,1])
outfile_tmp <- paste0(res_home,proj_name,"/table/","LASSO.","gene_coef",".csv");outfile_tmp
write.csv(geneCoef,file=outfile_tmp,row.names = F)

##############################    计算风险评分  ##############################
# 计算风险得分
riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
# 得到风险矩阵
outCol=c("months","OS",lassoGene)
# 根据中位数区分
risk = as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab <- cbind(riskScore,risk)
colnames(outTab)[1] <- "risk_score"
# 合并临床数据
outTab <- rownames_to_column(as.data.frame(outTab),"sample")
data_merge <- merge(data_surv, outTab, by="sample")

# save
outfile_tmp <- paste0(res_home,proj_name,"/table/","LASSO.","risk_group",".csv");outfile_tmp
write.csv(data_merge,outfile_tmp,row.names = F)

