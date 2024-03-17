########################    Enrich     ########################  
#### 

library(data.table)
library(tidyverse)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/Zwy.STAD/res/"
proj_name <- "3.Biology/"

# 读取差异基因
degs_list <- read.csv("./res/3.Biology/table/degs.Tstr.fc0.585.csv")
degs <- degs_list[degs_list$threshold != "NoSig",]

## id convert 
gene_id <- id_convert(input = degs[,1],from = "SYMBOL",to = "ENTREZID")

########################    富集分析    -- go   ######################## 
# go
kk <- enrichGO(gene_id$ENTREZID, OrgDb=org.Hs.eg.db,
               pvalueCutoff=0.05, qvalueCutoff=0.05, 
               ont="all",readable =T)
kk.df <- as.data.frame(GO_res)
# save
outfile_tmp <- paste0(res_home,proj_name, "/table/", "GO.","Tstr_CD8",".csv")
write.csv(kk.df, file=outfile_tmp)

### plot
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","GO.","Tstr_CD8", ".pdf")
pdf(file = outfile_tmp, width = 8, height = 6)
dotplot(kk, showCategory =5,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
dev.off()

########################    富集分析    -- kegg   ########################
# kegg
kk <- enrichKEGG(gene_id$ENTREZID,organism = "hsa",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
kk.df <- as.data.frame(KEGG_res)
# save
outfile_tmp <- paste0(res_home,proj_name, "/table/", "KEGG.","Tstr_CD8",".csv");outfile_tmp
write.csv(KEGG_res.df,file = outfile_tmp)

### plot
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","KEGG.","Tstr_CD8", ".pdf")
pdf(file = outfile_tmp, width = 8, height = 6)
dotplot(kk, showCategory =10)
dev.off()


########################    富集分析    -- GSEA   ######################## 
gene_id <- id_convert(input = degs_list[,1],from = "SYMBOL",to = "ENTREZID")
data_pre <- merge(gene_id, degs_list, by.x = "SYMBOL", by.y = "X")
colnames(data_pre)
# 提取用于GSEA的基因集
logFC.list <- data_pre$logFC
names(logFC.list) <- data_pre$ENTREZID
logFC.list = sort(logFC.list,decreasing = T)

####### GSEA
### 0.001可以保留比较少的通路 方便画一起
# read file
TERM2GENE_set <- read.gmt("E:/r_project/data/Msigdb/h.all.v2022.1.Hs.entrez.gmt")
colnames(TERM2GENE_set) <- c("ont", "gene")
# gsea
kk <- GSEA(logFC.list, TERM2GENE = TERM2GENE_set, 
           exponent=1, minGSSize=5,  
           pvalueCutoff=0.05, pAdjustMethod="BH",
           seed=F, by="fgsea")
kk.df <- as.data.frame(kk)

# save
outfile_tmp <- paste0(res_home, proj_name,"/table/","Hallmark.","Tstr_CD8",".csv");outfile_tmp
write.csv(kk.df, file = outfile_tmp)

#### plot
my_palet <- brewer.pal(5,"Set1")
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","GSEA.","Hallmark_top5",".pdf");outfile_tmp
# save
pdf(file = outfile_tmp, onefile = F,width = 7,height = 6)
gseaplot2(kk,c(1:5),color = my_palet,
          subplots = 1:2, pvalue_table = F)
dev.off()




