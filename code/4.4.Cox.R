########################  Cox回归验证risk score  ########################
#### 

library(openxlsx)
library(tidyverse)
library(data.table)
library(glmnet)
library(survival)
library(stringr)
library(broom)
library(ggplot2)
library(forestplot)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "4.Signature/"

#### 读入临床与生存信息
# clinical data
data_clinical <- read.csv("./data/TCGA-STAD.clinical_info.csv")
colnames(data_clinical)

# 读取分组数据
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")

# merge
group_sub <- data_group[,c(-5)]
data_merge <- merge(group_sub,data_clinical,by="sample")

#######################    修改临床参数   ####################### 
data_pre <- data_merge
colnames(data_pre)

{
  ### age
  idx_tmp <- data_pre$age
  idx_tmp <- case_when(idx_tmp < 65 ~ 1,
                       idx_tmp >= 65 ~ 2)
  table(idx_tmp)
  data_pre$age <- idx_tmp
  
  ### gender
  idx_tmp <- data_pre$gender
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "female" ~ 1,
                      idx_tmp %in% "male" ~ 2)
  table(idx_tmp)
  data_pre$gender <- idx_tmp
  
  ### grade
  idx_tmp <- data_pre$grade
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "G1" ~ 1,
                      idx_tmp %in% "G2" ~ 2,
                      idx_tmp %in% "G3" ~ 3)
  table(idx_tmp)
  data_pre$grade <- idx_tmp
  
  ### grade
  idx_tmp <- data_pre$new_tumor
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "NO" ~ 1,
                      idx_tmp %in% "YES" ~ 2)
  table(idx_tmp)
  data_pre$new_tumor <- idx_tmp
  
  ### stage_M
  idx_tmp <- data_pre$stage_M
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "M0" ~ 1,
                      idx_tmp %in% "M1" ~ 2)
  table(idx_tmp)
  data_pre$stage_M <- idx_tmp
  
  ### stage_N
  idx_tmp <- data_pre$stage_N
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "N0" ~ 1,
                      idx_tmp %in% "N1" ~ 2,
                      idx_tmp %in% "N2" ~ 3,
                      idx_tmp %in% "N3" ~ 4)
  table(idx_tmp)
  data_pre$stage_N <- idx_tmp
  
  ### stage_T
  idx_tmp <- data_pre$stage_T
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "T1" ~ 1,
                      idx_tmp %in% "T2" ~ 2,
                      idx_tmp %in% "T3" ~ 3,
                      idx_tmp %in% "T4" ~ 4)
  table(idx_tmp)
  data_pre$stage_T <- idx_tmp
  
  ### stage
  idx_tmp <- data_pre$stage
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "Stage I" ~ 1,
                      idx_tmp %in% "Stage II" ~ 2,
                      idx_tmp %in% "Stage III" ~ 3,
                      idx_tmp %in% "Stage IV" ~ 4,
                      .default = NA)
  table(idx_tmp)
  data_pre$stage <- idx_tmp
  
  ### tumor origin
  idx_tmp <- data_pre$tumor_origin
  table(idx_tmp)
  idx_tmp = case_when(idx_tmp %in% "Cardia, NOS" ~ 1,
                      idx_tmp %in% "Body of stomach" ~ 2,
                      idx_tmp %in% "Gastric antrum" ~ 3,
                      idx_tmp %in% "Fundus of stomach" ~ 4,
                      .default = NA)
  table(idx_tmp)
  data_pre$tumor_origin <- idx_tmp
  
}

####################### cox regression ####################### 
### 选择单因素回归显著的进行多因素回归
data_clinical <- data_pre
str(data_clinical)
colnames(data_clinical)

# Cox
target_col <- colnames(data_clinical)[-c(1,2,3,6,11,12,13)]
HR_res_list <- cox_regression(input=data_clinical,
                              var_cols= target_col,
                              time_col="OS.time", status_col="OS")
hr_res_univ <- HR_res_list$univ_res
hr_res_multiv <- HR_res_list$multiv_res

### plot
# 定义函数
sub_inf <- function(x, sub_method=max){
  sub_value <- sub_method(x[which(is.finite(x))], na.rm=TRUE)
  inf_idx <- which(is.infinite(x))
  if(length(inf_idx)>0)
    x[inf_idx] <- sub_value
  return(x)
}

# 运行时修改
HR_df = hr_res_multiv
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Cox.","multiv.","clinical_factor",".pdf");outfile_tmp
zero_value=1; auto_col=TRUE; width = 5; height = 5; outfile <- outfile_tmp

# 判定CI显示时的区间，避免离群值以及Inf的干扰
HR_df[, 3] <- sub_inf(as.numeric(HR_df[, 3]), sub_method=min)
HR_df[, 4] <- sub_inf(as.numeric(HR_df[, 4]), sub_method=max)

# pvalue取四位有效数字
HR_df[, 5] <- round(HR_df[, 5], digits=4)
# 设定森林图表格，添加标题行
HR_df_table <- as.matrix(HR_df[,c(1, 6, 5)])
HR_df_table <- rbind(colnames(HR_df_table), HR_df_table)
# 设定森林图数值矩阵（添加一行NA，以对应标题行）
HR_df2 <- rbind(NA, HR_df[, c(2:4)])
HR_df2 <- data.frame(HR_df2, box_col='#f89f68') # 第一次univ #71BD63 第二次 multiv #FFA400
if(auto_col){
  box_col <- as.character(HR_df2[, "box_col"])
  box_col[which(HR_df2[, "HR"]<1)] <- "blue"
}
# 设定森林图边界
clip_tmp_min <- min(c(0, HR_df2[, 2]), na.rm=TRUE)
clip_tmp_max <- max(c(2, HR_df2[, 3]), na.rm=TRUE)
clip_tmp <- c(clip_tmp_min, clip_tmp_max)
# clip_tmp <- c(-Inf, Inf)

# order
HR_df_table <- HR_df_table[c(1,4,2,3,7,6,8,5,9),]

pdf(outfile, onefile=FALSE, width = width, height = height)
# 构建森林图对象
forestplot(labeltext = HR_df_table, 
           mean = as.numeric(HR_df2[, 1]), 
           lower = as.numeric(HR_df2[, 2]), 
           upper = as.numeric(HR_df2[, 3]),
           is.summary = c(TRUE, rep(FALSE, nrow(HR_df2)-1)), # 第一行字体加粗
           align=c("l", "c", "c"),
           clip=clip_tmp,  # 设定CI 显示时2的区间
           xlog=FALSE, 
           zero = zero_value,
           boxsize = 0.2, # 设置点估计的方形大小
           new_page = TRUE,
           txt_gp =fpTxtGp(ticks=gpar(cex=0.8)), # x轴字体
           lwd.zero = 1,# 设置参考线的粗细
           lwd.ci = 1,# 设置区间估计线的粗细
           col=fpColors(box=as.character(HR_df2[,"box_col"]), 
                        summary= "#E41A1C",lines = 'lightgray',zero = 'lightgray'))
#使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
dev.off()

