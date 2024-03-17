########################     根据单因素回归分析结果建立nomogram       ######################## 
### ref: https://blog.csdn.net/Ayue0616/article/details/126869563

library(stringr)
library(data.table)
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(rms)
library(foreign)
library(survival)
library(regplot)
library(mstate)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "4.Signature/"

# clinical data
data_clinical <- read.csv("./data/TCGA-STAD.clinical_info.csv")
colnames(data_clinical)

# 读取分组数据
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")

# merge
group_sub <- data_group[,c(-5)]
data_merge <- merge(group_sub,data_clinical,by="sample")
data_merge$months <- data_merge$OS.time/30

####################### 修改临床参数 ####################### 
data_pre <- data_merge
data_pre <- data_pre[,c("OS","months","risk_score","age","new_tumor")]

### age
idx_tmp <- data_pre$age
idx_tmp <- case_when(idx_tmp < 65 ~ "< 65",
                     idx_tmp >= 65 ~ ">= 65")
table(idx_tmp)
data_pre$age <- idx_tmp

###########   nomogram   ###########
ddist <- datadist(data_pre)
options(datadist='ddist')
rt <- data_pre

### 传统Nomogram
cox <- cph(Surv(months,OS) ~ age+risk_score+new_tumor,surv=T,x=T, y=T,data=rt) 
surv <- Survival(cox)
sur_1_year<-function(x)surv(1*12,lp=x)
sur_3_year<-function(x)surv(1*12*5,lp=x)
sur_5_year<-function(x)surv(1*12*7,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),
                    lp= F,funlabel=c('1-Year survival','5-Year survival','7-Year survival'),
                    maxscale=100,
                    fun.at= c('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'))
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Nomogram.","classical",".pdf");outfile_tmp
pdf(outfile_tmp, height = 6, width = 8, onefile = F)
plot(nom_sur,xfrac=0.4)
dev.off()

### 花式Nomogram
rt <- data_pre
rt$id <- 1:nrow(rt)
df.w <- crprep("months", "OS",
               data=rt, trans=c(1,2),
               cens=0, id="id",
               keep=c("age","risk_score","new_tumor"))
df.w$T<- df.w$Tstop - df.w$Tstart

m.cph<-coxph(Surv(months,OS==1)~age+risk_score+new_tumor,
             data=rt)
summary(m.cph)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Nomogram.","new",".pdf");outfile_tmp
regplot(m.cph, observation=rt[rt$id==14,], #可以选择第10个病人，看他的各项得分
        failtime = c(1*12,5*12,7*12), prfail = TRUE,droplines=T)
# ggsave(outfile_tmp,width = 8,height = 6)
dev.off()

###########   Calibration   ###########
### 1-year
cox_1 <- cph(Surv(months,OS) ~ age + risk_score + new_tumor,surv=T,x=T, y=T,time.inc = 1*12,data=rt) 
cal_1 <- calibrate(cox_1, cmethod="KM", method="boot", u=1*12, m= 100, B=1000)
### 2-year
cox_2 <- cph(Surv(months,OS) ~ age+risk_score+new_tumor,surv=T,x=T, y=T,time.inc = 3*12,data=rt) 
cal_2 <- calibrate(cox_2, cmethod="KM", method="boot", u=3*12, m= 100, B=1000)
### 3-year
cox_3 <- cph(Surv(months,OS) ~ age+risk_score+new_tumor,surv=T,x=T, y=T,time.inc = 5*12,data=rt) 
cal_3 <- calibrate(cox_3, cmethod="KM", method="boot", u=5*12, m= 100, B=1000)


### plot
my_color <- c("#7BC67C","#4b84b3","#f89f68")
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Nomogram.","calibration",".pdf");outfile_tmp
pdf(outfile_tmp,onefile = F,width = 8,height = 8)
plot(cal_1, lwd=3, lty=2, errbar.col="black",col= my_color[1],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="Nomogram-predicted OS",ylab="Observed OS")
par(new = T)
plot(cal_2, lwd=3, lty=2, errbar.col="black",col= my_color[2],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="",ylab="")
par(new = T)
plot(cal_3, lwd=3, lty=2, errbar.col="black",col= my_color[3],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="",ylab="")
lines(cal_1[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[1] ,pch = 16)
lines(cal_2[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[2] ,pch = 16)
lines(cal_3[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[3] ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 2,lwd = 2,col = "black")
legend(0.8,0.2,
       c("1-year","3-year","5-year"), 
       lty = c(1,1,1), lwd = c(3,3,3), 
       col = my_color, 
       bty = "n")
dev.off()

################   C index   ################   
f<-coxph(Surv(months,OS==1)~age+new_tumor+risk_score,data=rt)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index  ##


################   ROC   ################  
pret <- predict(cox,data_pre,type = "lp")
ROC <- data.frame(months = data_pre[,"months"],status = data_pre[,"OS"],
                  score = pret)
ROC <- na.omit(ROC)
res_ROC <- timeROC(T = ROC$months,
                   delta = ROC$status,
                   marker = ROC$score, 
                   cause=1, 
                   weighting="marginal", 
                   times=c(1*12,3*12,5*12), 
                   ROC = TRUE,
                   iid = TRUE)

# 计算AUC
ROC <- res_ROC
auc_1 = ROC$AUC[[1]]; auc_2 = ROC$AUC[[2]]; auc_3 = ROC$AUC[[3]]
auc_1; auc_2; auc_3

dat = data.frame(tpr_1 = ROC$TP[,1],fpr_1 = ROC$FP[,1],tpr_2 = ROC$TP[,2],
                 fpr_2 = ROC$FP[,2],tpr_3 = ROC$TP[,3],fpr_3 = ROC$FP[,3])

### plot 
my_palet <- c("#FFA500", "#72BE64", "#97C0E2")
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","ROC.","nomogram.","years.","1_2_3",".pdf")
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
  annotate("text",x = .75, y = .25,label = paste("AUC of 1 years = ",round(auc_1,3)),color = my_palet[1])+
  annotate("text",x = .75, y = .15,label = paste("AUC of 3 years = ",round(auc_2,3)),color = my_palet[2])+
  annotate("text",x = .75, y = .05,label = paste("AUC of 5 years = ",round(auc_3,3)),color = my_palet[3])+
  scale_x_continuous(name  = "1-Specificity")+
  scale_y_continuous(name = "Specificity")
dev.off()









