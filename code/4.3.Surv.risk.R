########################     对风险组进行km分析       ######################## 

library(survival)
library(survminer)
library(timeROC)
library(stringr)
library(data.table)
library(tidyverse)
library(openxlsx)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "4.Signature/"

# 读取risk group
group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(group)

################################  KM生存分析   ################################ 
melcli_group <- group
colnames(melcli_group)[c(2,3)] <- c("status","months")
melcli_group$months <- melcli_group$months / 30

# log rank 检验
diff = survdiff(Surv(months, status) ~ risk,data = melcli_group)
# 计算p值，转化为科学计数法
pValue = 1-pchisq(diff$chisq,df=1)
pValue = signif(pValue,2);pValue

# km曲线
fit <- survfit(Surv(months, status) ~ risk, data = melcli_group)

### plot
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Surv.","8_gene",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, width = 5,height = 6)   
ggsurvplot(fit, 
           data=melcli_group,
           conf.int=F,
           pval=T,
           pval.size=4,
           risk.table=T,
           legend.labs=c("high risk", "low risk"), #, "Cluster4"
           legend.title="group",
           xlab="Time(months)",
           ylab="Status",
           break.time.by = 36,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()

################################  ROC   ################################ 
melcli_group <- group
colnames(melcli_group)[c(2,3)] <- c("status","months")
melcli_group$months <- melcli_group$months / 30
colnames(melcli_group)

ROC<- timeROC(T = melcli_group$months, # 结局时间 
              delta = melcli_group$status, # 生存结局 
              marker = melcli_group$risk_score, # 预测变量 ##此处请注意，计算ROC的biomaker，默认是marker值越大，事件越可能发生；反之的话，前面加-号
              cause = 1, #阳性结局赋值，比如死亡与否
              weighting = "marginal", #权重计算方法，marginal是默认值，采用km计算删失分布
              times=c(1*12,5*12,7*12), #时间点，选取5年(60个月)和8年生存率
              ROC = TRUE,
              iid = TRUE)
# 计算AUC
auc_1 = ROC$AUC[[1]]; auc_2 = ROC$AUC[[2]]; auc_3 = ROC$AUC[[3]]
auc_1; auc_2; auc_3

dat = data.frame(tpr_1 = ROC$TP[,1],fpr_1 = ROC$FP[,1],tpr_2 = ROC$TP[,2],
                 fpr_2 = ROC$FP[,2],tpr_3 = ROC$TP[,3],fpr_3 = ROC$FP[,3])

### plot 
my_palet <- c("#7BC67C","#4b84b3","#f89f68")
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","ROC.","8_gene",".years.","3_5_7",".pdf")
pdf(file = outfile_tmp, width=4.5, height=4)
ggplot() + 
  geom_smooth(data = dat,aes(x = fpr_1, y = tpr_1),color = my_palet[1]) + 
  geom_smooth(data = dat,aes(x = fpr_2, y = tpr_2),color = my_palet[2])+
  geom_smooth(data = dat,aes(x = fpr_3, y = tpr_3),color = my_palet[3])+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme(panel.background=element_rect(colour=NA,fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of 5 years = ",round(auc_1,2)),color = my_palet[1])+
  annotate("text",x = .75, y = .15,label = paste("AUC of  8 years = ",round(auc_2,2)),color = my_palet[2])+
  annotate("text",x = .75, y = .05,label = paste("AUC of 10 years = ",round(auc_3,2)),color = my_palet[3])+
  scale_x_continuous(name  = "1-Specificity")+
  scale_y_continuous(name = "Specificity")
dev.off()



