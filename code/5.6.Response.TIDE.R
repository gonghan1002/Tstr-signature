########################    Immunotherapy   ######################## 
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

# 读取TIDE数据
data_exprs <- read.csv("./data/TIDE_res.TCGA-STAD.SYMBOL.csv")

# 读取分组
data_group <- read.csv("./res/4.Signature/table/LASSO.risk_group.2.csv")
colnames(data_group)
data_group <- data_group[,c(1,4,5)]
colnames(data_group) <- c("sample","score","group")

#### 合并数据
colnames(data_exprs)
colnames(data_exprs)[1] <- "sample"
# merge
data_expr_sub <- data_exprs[,c(1,3,4)]
data_merge <- merge(data_expr_sub, data_group, by = "sample")

#################    卡方检验  #################   
data_pre <- data_merge
table(data_pre$Responder)

# 修改名称
data_pre$Responder <- ifelse(data_pre$Responder=="True","R","NR")

################   比较反应响应情况  ################   
summary(data_pre$score)
cutpoint <- quantile(data_pre[ ,"score"], probs=0.5);cutpoint
# 选择最佳分割点
data_pre$group <- ifelse(data_pre$score < cutpoint,"low","high")

#### plot 
colnames(data_pre)
data_plot <- data_pre[,c("Responder","TIDE","score","group")]
# data_plot <- na.omit(data_plot)
chisq.test(data_plot$group,data_plot$Responder)
pval <- fisher.test(data_plot$group,data_plot$Responder)$p.value; pval

# save
my_platlet <- brewer.pal(4,'Set1')
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Barplot.","response.","TIDE",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, width = 5,height = 6)   
ggplot(data=data_plot, mapping=aes(x=group,fill=Responder))+
  geom_bar(stat="count",width=0.8,position='fill')+
  scale_fill_manual(values=c("#4b84b3","#f89f68"))+ 
  geom_text(stat='count',aes(label=..count..) #scales::percent(..count../sum(..count..)
            , color="white", size=3.5,position=position_fill(0.5))+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(color="black",size=12,
                                   # angle = 60,hjust = 1,vjust = 1
          ),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          strip.background = element_blank())+
  labs(title="Chi-squared test: p=0.002",y="Pecentage(%)")
dev.off()

################   比较m1a  ################   
colnames(data_plot)
data_plot2 <- data_plot

# 调一下risk顺序
data_plot2$group <- factor(data_plot2$group,levels = c("low","high"))

# plot
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Barplot.","TIDE.risk_score",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, height = 4, width = 3)   
ggplot(data = data_plot2,
       aes(y = TIDE,
           x = group))+
  geom_boxplot(alpha = 1,
               fill = c("#4b84b3","#f89f68"))+
  stat_compare_means(method = "t.test",label = "p.format") + 
  theme_bw()+ theme_classic() +
  ylab('TIDE') +
  xlab('group')
dev.off()













