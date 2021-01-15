#加载包
library(DESeq2)
#导入数据（有重复值，所以平时我们都是用geneid，但这里给的基因名，放在excel去重下再导入）
mRNA_count<- read.table("GSE122709_count.txt",header = T, stringsAsFactors = F,row.names = 1) 
#新增样本信息表
sample_info<-data.frame(ID=colnames(mRNA_count))

#样本信息表新增分组
#第一种方法：字符串截取
sample_info$group<-substr(sample_info$ID,1,2)
sample_info$group2<-substr(sample_info$ID,1,2)
sample_info[,2]<-c("NC","NC","NC","NC","NC","D","D","D","D","D","D","D","D","D","D")
#转化成因子类型，方便分组
sample_info$group<-factor(sample_info$group)
#查看表格的类型
str(sample_info)
#分组因子
group_list<-sample_info$group

#数据格式处理
dds <- DESeqDataSetFromMatrix(countData = mRNA_count,
                              colData = data.frame(group_list=sample_info$group),
                              design= ~ group_list)
#DESeq2运行,标准化
dds2<- DESeq(dds) 
#查看结果(对照组和实验组进行对比,注意实验组在前)
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","D","NC"))
#标准化后的矩阵
TMM<-counts(dds2,normalized=T)

#转化为data.frame形式
NC_D_deseq2<-data.frame(res)
#标准化的矩阵
rld<-rlogTransformation(dds2)

##筛选差异表达基因
#首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序

NC_D_deseq2<- NC_D_deseq2[order(NC_D_deseq2$padj, NC_D_deseq2$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.05标识up,代表显著上调的基因
#log2FC<1$padj<0.05标识down,代表显著下调的基因
#其余标识 none，代表非差异的基因

NC_D_deseq2[which(NC_D_deseq2$log2FoldChange >= 1 & NC_D_deseq2$padj < 0.05),'direction'] <- 'up'
NC_D_deseq2[which(NC_D_deseq2$log2FoldChange <= -1 & NC_D_deseq2$padj < 0.05),'direction'] <- 'down'
NC_D_deseq2[which(abs(NC_D_deseq2$log2FoldChange) <= 1 | NC_D_deseq2$padj >= 0.05),'direction'] <- 'none'

#输出选择的差异基因总表  subset求子集
library(tidyverse)
#方法一用filter
#NC_D_deseq2_select1<-NC_D_deseq2 %>% filter(direction=="up" )
NC_D_deseq2_select2<-NC_D_deseq2 %>% filter(direction %in% c('up', 'down'))
#方法二用subset
NC_D_deseq2_select <- subset(NC_D_deseq2, direction %in% c('up', 'down'))
#保存数据
write.table(NC_D_deseq2_select, file = 'D_NC.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#画火山图看差异基因
library(tidyverse)
NC_D_deseq2<-rownames_to_column(NC_D_deseq2,var = "ID")
top_de <- filter(NC_D_deseq2,abs(log2FoldChange) > 4)
#
library(ggrepel)
my_palette <- c('#E64B35FF','#999999','#4DBBD5FF')
#创建一个画布
library(ggplot2)
#添加几何对象  geom_point散点图，将direction映射给点颜色  aes映射颜色
ggplot(data = NC_D_deseq2,aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = direction,size = abs(log2FoldChange))) + 
  geom_label_repel(data = top_de,aes(label = ID))+
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed') +
  scale_color_manual(values = my_palette) +
  scale_size(range = c(0.1,2)) +
  labs(x = 'log2 fold change',
       y = '-log10(pvalue)',
       title = 'Vocano plot',
       size = 'log2 fold change') +
  ylim(c(0,20)) +
  guides(size = FALSE) +
  theme_classic() + 
  theme(plot.title = element_text(size = 18,hjust = 0.5),
        legend.background = element_blank(),
        legend.position = c(0.93,0.85)) 
