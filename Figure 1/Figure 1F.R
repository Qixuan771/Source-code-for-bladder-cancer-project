library(BiocParallel)
numCores <- 8
register(MulticoreParam(workers = numCores), default = TRUE)

args = commandArgs(trailingOnly=TRUE)
fl=args[1]

library(ggplot2)
library(ggrepel)
library('plyr')
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)

#enh_genes=read.csv("Allmerged_ChIP_K27ac_matrix.txt",
#                         header=F,sep="\t",numerals = c("allow.loss", "warn.loss", "no.loss"))
#fl="RSClus_enh_all_ID_cor_sign_poscorr.txt"

enh_genes=read.csv(fl,
                         header=F,sep="\t",numerals = c("allow.loss", "warn.loss", "no.loss"))

colnames(enh_genes)=c("chr","start","stop","ID","gene","CS_RT4_1","CS_RT4_2","CS_SW780_1","CS_SW780_2","CS_SCaBER_1","CS_SCaBER_2","CS_HT1376_1","CS_HT1376_2","cor","rand_cor_mean","rand_cor_sd","Log10pval","Log10BH")
head(enh_genes)

enh_genes_l2=enh_genes
enh_genes_l2[,6:13]=log2(enh_genes[,6:13]+1)

var=apply(enh_genes_l2[,6:13],1,var)
enh_genes_f1=enh_genes_l2[order(var,decreasing=T)[1:10000],]

enh_genes_l2$CS_cell=""

enh_genes_l2$CS_cell[grep("RT4",enh_genes_l2$ID)]="RT4"
enh_genes_l2$CS_cell[grep("SC",enh_genes_l2$ID)]="SC"
enh_genes_l2$CS_cell[grep("HT",enh_genes_l2$ID)]="HT"
enh_genes_l2$CS_cell[grep("SW",enh_genes_l2$ID)]="SW"

ID_freq=as.data.frame(table(as.vector(enh_genes_f1$ID)))
ID_freq=ID_freq[order(ID_freq$Freq,decreasing=T),]
gene_freq=as.data.frame(table(as.vector(enh_genes_f1$gene)))
gene_freq=gene_freq[order(gene_freq$Freq,decreasing=T),]

dim(ID_freq) #2229 2
dim(gene_freq) #2781 2

AS_mat=as.matrix(enh_genes_l2[,6:13])
rownames(AS_mat)=enh_genes_l2$gene
AS_mat_nm=AS_mat-apply(AS_mat,1,mean)
AS_mat_nm_f1=AS_mat_nm[order(var,decreasing=T)[1:10000],]
AS_mat_nm_f1=AS_mat_nm_f1[(enh_genes_f1$cor>=0.5),] #5805 8

enh_genes_f2=enh_genes[order(var,decreasing=T)[1:10000],]
enh_genes_f2=enh_genes_f2[(enh_genes_f2$cor>=0.5),]

write.table(enh_genes_f2,
file=paste("top5805_",fl,sep=""),
sep="\t",row.names=F,col.names=F,quote = FALSE)

################# Deal with our RNASeq now ###################
TJ_4cells_genes=read.csv("191104_RNASeq_aggregate_TPM.txt",header=T,sep="\t",numerals = c("allow.loss", "warn.loss", "no.loss"))
colnames(TJ_4cells_genes)=c("Genes","RS_RT4_1","RS_RT4_2","RS_SW780_1","RS_SW780_2","RS_SCaBER_1","RS_SCaBER_2","RS_HT1376_1","RS_HT1376_2")
head(TJ_4cells_genes)
GSM_ENSG=read.csv("GSM_ENSG_file_GeneID.txt",
header=T,sep="\t",numerals = c("allow.loss", "warn.loss", "no.loss"))
ENSEMBL_annotations=read.csv("/projects/b1100/TJ/Bladder_Cancer/invitro/RNA_SEQ/RNASeq_Matrices/Homo_sapiens.GRCh37.87_Annotations.txt",header=FALSE,sep=" ",
                             ,stringsAsFactors = F)
ENSEMBL_annotations_prot=unique(ENSEMBL_annotations[ENSEMBL_annotations$V3=="protein_coding",])

TJ_4cells_genes$Genes=unlist(lapply(as.vector(TJ_4cells_genes$Genes),function(x) unlist(strsplit(x,"[.]"))[1]))
TJ_4cells_genes1=subset(TJ_4cells_genes,TJ_4cells_genes$Genes %in% as.vector(GSM_ENSG$initial_alias))

m=match(TJ_4cells_genes1$Genes,as.vector(GSM_ENSG$initial_alias))

rmat=as.matrix(TJ_4cells_genes1[,2:ncol(TJ_4cells_genes)])
row.names(rmat)=as.vector(GSM_ENSG$name)[m]
head(rmat)
dim(rmat)

#Filter for protein coding

rmat_protcod=rmat[which(rownames(rmat) %in% ENSEMBL_annotations_prot$V2),]
dim(rmat_protcod) #18084 8

#Remove <2 FPKM genes
rmat_protcod_f1=rmat_protcod[(rowSums(rmat_protcod>2)>=2),]
dim(rmat_protcod_f1) #11675     8

rmat_f1=rmat[(rowSums(rmat>2)>=2),]
dim(rmat_f1) #13484     8

#Quantile Normalize data
require("preprocessCore")
rmat_protcod_f1_QN=normalize.quantiles(rmat_protcod_f1)
row.names(rmat_protcod_f1_QN)=row.names(rmat_protcod_f1)
colnames(rmat_protcod_f1_QN)=colnames(rmat_protcod_f1)

### ComBat BATCH removal
#https://github.com/kundajelab/mesoderm/blob/master/RNA-seq/bulkDataViz.Rmd
library(sva)
rmat_protcod_f1_QN_batch=ComBat(log2(rmat_protcod_f1_QN+1),
                                  batch=c("TJ_2016","TJ_2017","TJ_2019","TJ_2019","TJ_2016","TJ_2017","TJ_2019","TJ_2019"),
                                  mod=NULL, par.prior = TRUE, prior.plots = FALSE)

rmat_protcod_f1_QN_batch[rmat_protcod_f1_QN_batch<0]=0
rmat_protcod_f1_QN_batch_m2=rmat_protcod_f1_QN_batch-apply(rmat_protcod_f1_QN_batch,1,mean)

################### Plot ranked plot ################
options(scipen=999) # turn-off scientific notation like 1e+48
theme_set(theme_bw())  # pre-set the bw theme.

#Ranked plot
ID_freq$rank=1:nrow(ID_freq)
gene_freq$rank=1:nrow(gene_freq)

ggplot(ID_freq, aes(y=Freq,x=rank)) + geom_point(aes(colour = "black",size=1))
#+scale_colour_manual(values = c("A"="purple","B"= "black"))+
#geom_text_repel(data=t[(t$threshold=="A"),], aes(label=Gene),
#direction="y",hjust = 1,xlim  = c(NA, 0.66*max(t$rank)),ylim=c(0.2*max(t$`Log10(P)`), NA),nudge_x=10)+
#labs(subtitle=unique(dff_str$Group)[i], x="Rank", y="-Log10(Pvalue)", title="MOTIF Analysis Results", caption = "Source: 190405")
#Save plot as PDF
ggsave(paste("TJ_",format(Sys.Date(), format="%Y%m%d"),"_ChIPSeq_corr_RankedIDs",".pdf",sep=""),width = 6, height = 6, units = 'in')

ggplot(gene_freq, aes(y=Freq,x=rank)) + geom_point(aes(colour = "black",size=1))
ggsave(paste("TJ_",format(Sys.Date(), format="%Y%m%d"),"_ChIPSeq_corr_RankedGenes",".pdf",sep=""),width = 6, height = 6, units = 'in')
################## HEATMAP #################

MDACC_Lum=c("FGFR3","FOXA1","GPX2","ERBB2","ERBB3","CYP2J2","GATA3","PPARG","KRT19","KRT7","KRT8","KRT20","FABP4","CD24","KRT18","XBP1")
MDACC_Bas=c("CD44","CDH3","KRT14","KRT16","KRT5","KRT6A","KRT6B","KRT6C","KRT1")

mL=which(enh_genes_l2$gene %in% MDACC_Lum)
mB=which(enh_genes_l2$gene %in% MDACC_Bas)

mL1=which(rownames(AS_mat_nm_f1) %in% MDACC_Lum)
mB1=which(rownames(AS_mat_nm_f1) %in% MDACC_Bas)

pdf(file = paste(strsplit(fl,split="[.]")[[1]][1],"_hmap_Lum_Bas_only_withintop5K.pdf",sep=""),
    width = 4, height = 8) # 1.46

ht1=Heatmap(AS_mat_nm_f1[mL1,],
col = colorRamp2(c(-1.5, 0, 1.5),c("blue", "#EEEEEE", "red")),
            show_column_names = T, show_row_names = T,
            cluster_columns = T,cluster_rows = T,
            show_column_dend = T, show_row_dend = T,
            clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
            show_heatmap_legend =T)

ht2=Heatmap(AS_mat_nm_f1[mB1,],
col = colorRamp2(c(-1.5, 0, 1.5),c("blue", "#EEEEEE", "red")),
show_column_names = T, show_row_names = T,
cluster_columns = T,cluster_rows = T,
show_column_dend = T, show_row_dend = T,
clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
show_heatmap_legend =T)

ht1 %v% ht2
dev.off()


################ Connect ENH correlation values to TPM values of that Luminal/Basal gene

eng_Lum_TPM=rmat_protcod_f1_QN_batch_m2[match(rownames(AS_mat_nm_f1[mL1,]),rownames(rmat_protcod_f1_QN_batch_m2)),]
eng_Bas_TPM=rmat_protcod_f1_QN_batch_m2[match(rownames(AS_mat_nm_f1[mB1,]),rownames(rmat_protcod_f1_QN_batch_m2)),]

dim(eng_Lum_TPM) #12   8
dim(eng_Bas_TPM) #2   8

#### Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ha1_L = HeatmapAnnotation(which ="row",
                                Lum_TPM = eng_Lum_TPM[,c(1,2,7,8,4,3,6,5)],
                                col = list(Lum_TPM = col_fun)
)

ha1_B = HeatmapAnnotation(which ="row",
                                Bas_TPM = eng_Bas_TPM[,c(1,2,7,8,4,3,6,5)],
                                col = list(Bas_TPM = col_fun)
)

ht1=Heatmap(AS_mat_nm_f1[mL1,],
col = colorRamp2(c(-1.5, 0, 1.5),c("blue", "#EEEEEE", "red")),
            show_column_names = F, show_row_names = T,
            cluster_columns = T,cluster_rows = T,
            show_column_dend = T, show_row_dend = T,
            clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
            show_heatmap_legend =T,right_annotation=ha1_L)

ht2=Heatmap(AS_mat_nm_f1[mB1,],
col = colorRamp2(c(-1.5, 0, 1.5),c("blue", "#EEEEEE", "red")),
show_column_names = T, show_row_names = T,
cluster_columns = T,cluster_rows = T,
show_column_dend = T, show_row_dend = F,
clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
show_heatmap_legend =T,right_annotation=ha1_B)

pdf(file = paste(strsplit(fl,split="[.]")[[1]][1],"_hmap_Lum_Bas_only_wTPM_1.pdf",sep=""),
width = 6, height = 8) # 1.46
ht1 %v% ht2
dev.off()

######## Plot genome-wide

eng_TPM=rmat_protcod_f1_QN_batch_m2[match(rownames(AS_mat_nm_f1),rownames(rmat_protcod_f1_QN_batch_m2)),]


dim(eng_TPM) #5805   8

#### Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

mLL=match(MDACC_Lum,rownames(eng_TPM))
mBB=match(MDACC_Bas,rownames(eng_TPM))

ha1 = HeatmapAnnotation(which ="row",
                                TPM = eng_TPM[,c(4,3,1,2,6,5,7,8)],
                                col = list(TPM = col_fun),
                                LumBas = anno_mark(at = which(rownames(eng_TPM) %in% c(MDACC_Lum,MDACC_Bas)),
                                                    labels = rownames(eng_TPM)[which(rownames(eng_TPM) %in% c(MDACC_Lum,MDACC_Bas))])
)

ht=Heatmap(AS_mat_nm_f1,
            col = colorRamp2(c(-1.5, 0, 1.5),c("blue", "#EEEEEE", "red")),
            show_column_names = T, show_row_names = F,
            cluster_columns = T,cluster_rows = T,
            show_column_dend = T, show_row_dend = T,
            clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",
            show_heatmap_legend =T,right_annotation=ha1)

pdf(file = paste(strsplit(fl,split="[.]")[[1]][1],"_hmap_top5K_wTPM_1.pdf",sep=""),
width = 6, height = 8) # 1.46
ht
dev.off()
