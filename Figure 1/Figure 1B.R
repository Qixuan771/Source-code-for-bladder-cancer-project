library(DESeq2)
library(ggrepel)
counts = read.table("RNAseq_4celllines_counts.dedup.txt",header = F, row.names = 1)
colnames(counts) = c("HT1376.rep1",	"HT1376.rep2",	"SCABER_rep1",	"SCABER_rep2",	"RT4_rep1",	"RT4_rep2",	"SW780_rep1",	"SW780_rep2")
countdata=sapply(round(counts), as.integer)
condition = c(rep('basal',4), rep('luminal',4))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
rownames(dds) = rownames(counts)
## Pre-filter
keep = rowSums(counts(dds)) >1
dds = dds[keep,]
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition", "basal", "luminal"))
res <- res[order(res$pvalue),]
threshold_OE <- res$padj < 0.01 & res$log2FoldChange > 2
res = res[! is.na(res$padj),]
res2 = as.data.frame(res)
res2$sig <- -log10(res2$padj)
res2$condition <- ifelse((res2$log2FoldChange > 2) & (res2$padj < 0.01), "Up",
                         ifelse((res2$log2FoldChange < -2) & (res2$padj < 0.01), "Down",
                                ifelse((res2$log2FoldChange > 2) & (res2$padj > 0.01), "NS",
                                       ifelse((res2$log2FoldChange < -2) & (res2$padj > 0.01), "NS",
                                              ifelse((res2$log2FoldChange < 2) & (res2$padj > 0.01), "NS", "NS")))))

dim(res2[res2$condition == 'Up',]) #427 8
dim(res2[res2$condition == 'Down',]) # 524 8
# get up-regulated and down-regulated genes
write.table(res2[res2$condition == 'Up',],file = "basal.specific.genes.txt",sep = "\t",
            quote=F, col.names = T, row.names = T)

write.table(res2[res2$condition == 'Down',],file = "luminal.specific.genes.txt",sep = "\t",
            quote=F, col.names = T, row.names = T)

#To generate plot#

res2$genelabels <- ""
res2[res2$condition=='Up',]$genelabels[1:10] <- rownames(res2[res2$condition=='Up',])[1:10]
res2[res2$condition=='Down',]$genelabels[1:10] <- rownames(res2[res2$condition=='Down',])[1:10]

##volcano plot
pdf("figure 1B.pdf")
ggplot(res2) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = condition)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj),label = genelabels), 
                  box.padding  = 0.5) +
  #ggtitle("Mov10 overexpression") 
  xlab("log2 fold change")+ggtitle("Basal vs Luminal")+ theme(plot.title = element_text(hjust = 0.5))+ 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + 
  geom_vline(aes(xintercept=2), linetype = "dashed") + 
  geom_vline(aes(xintercept= -2), linetype = "dashed") +
  scale_color_manual(values=c("blue","grey", "red" ))+
  xlim(c(-15,15))+
  ylim(c(0,30))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_hline(aes(yintercept= -log10(0.01)), linetype = "dashed") #To add lines#

dev.off()
