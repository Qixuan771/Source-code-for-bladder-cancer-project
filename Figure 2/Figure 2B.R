library(ggrepel)
luminal = read.table("luminal.pvalue.0.01.txt",header = T)
basal = read.table("basal.0.0.1.pvalue.txt",header=T)
luminal.dat = data.frame(luminal)
basal.dat = data.frame(basal)
######luminal specific motif
luminal.dat$label = ""
luminal.dat$label[1:20] = as.character(luminal.dat$Motif)[1:20]
pdf("luminal.specific.motifrank.0.01.pdf")
ggplot(luminal.dat, aes(x=rank, y=adjusted.p.value)) + geom_point(size=3)+ 
  geom_text_repel(aes(x=rank, y =adjusted.p.value, label=label, size=10), box.padding = 1)+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.5)))+
  theme(axis.title = element_text(size = 20)) +    
  xlab("Rank")+ylab("-Log10(Pvalue)")
dev.off()
######basal specific motif
basal.dat$label = ""
basal.dat$label[1:20] = as.character(basal.dat$Motif)[1:20]
pdf("basal.specific.motifrank.0.01.pdf")
ggplot(basal.dat, aes(x=rank, y=adjusted.p.value)) + geom_point(size=3)+ 
  geom_text_repel(aes(x=rank, y =adjusted.p.value, label=label, size=10), box.padding = 1)+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.5)))+
  theme(axis.title = element_text(size = 20)) +    
  xlab("Rank")+ylab("-Log10(Pvalue)")

dev.off()