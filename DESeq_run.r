library(DESeq2)
library("RColorBrewer")
library("gplots")

sampleTable <- read.csv("design_table.csv",stringsAsFactors = FALSE)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "./", 
                                  design= ~ Genotype)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
rld <- varianceStabilizingTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(Genotype, rownames(colData(dds)), sep=" : "))

pdf("MA_plot_before.pdf",10,10)
plotMA(dds,ylim=c(-1,1), main = "MA plot ( before remove outlers)")
dev.off()






pdf("cluster_plot.pdf",10,10)
heatmap.2(mat,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

pdf("PCA_plot.pdf",6,6)
print(plotPCA(rld, intgroup=c("Genotype")))
dev.off()


tiff("cluster_plot.tiff",10,10)
heatmap.2(mat,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

tiff("PCA_plot.tiff",6,6)
print(plotPCA(rld, intgroup=c("Genotype")))
dev.off()


withoutTtn <- assay(rld[-26687,]) 
distsRLWithouTtn <- dist(t(withoutTtn))
matWithoutTtn <- as.matrix(distsRLWithouTtn)
rownames(matWithoutTtn) <- colnames(matWithoutTtn) <- with(colData(dds),
                                       paste(Genotype, rownames(colData(dds)), sep=" : "))


pdf("cluster_plot_withoutTTN.pdf",10,10)
heatmap.2(matWithoutTtn,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 




## remove the outlers
sampleTableReduce <- sampleTable[sampleTable$SampleID != "543-1-1",]
ddsR <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableReduce,
                                   directory = "./", 
                                   design= ~ Genotype)
rldR <- varianceStabilizingTransformation(ddsR, blind=TRUE)
distsRLR <- dist(t(assay(rldR)))
matR <- as.matrix(distsRLR)
rownames(matR) <- colnames(matR) <- with(colData(ddsR),
                                       paste(Genotype, rownames(colData(ddsR)), sep=" : "))



pdf("cluster_plot_R.pdf",10,10)
heatmap.2(matR,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

pdf("PCA_plot_R.pdf",6,6)
print(plotPCA(rldR, intgroup=c("Genotype")))
dev.off()


pdf("cluster_plot_R.pdf",10,10)
heatmap.2(matR,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

pdf("PCA_plot_R.pdf",6,6)
print(plotPCA(rldR, intgroup=c("Genotype")))
dev.off()





tiff("cluster_plot_R.tiff")
heatmap.2(matR,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

tiff("PCA_plot_R.tiff")
print(plotPCA(rldR, intgroup=c("Genotype")))
dev.off()


tiff("cluster_plot.tiff")
heatmap.2(mat,trace="none", col = rev(hmcol),margin=c(13, 13))
dev.off() 

tiff("PCA_plot.tiff")
print(plotPCA(rld, intgroup=c("Genotype")))
dev.off()


tiff("MA_plot_before.tiff")
plotMA(dds,ylim=c(-1,1), main = "MA plot ( before remove outlers)")
dev.off()

tiff("MA_plot_after.tiff")
plotMA(ddsMT_WT,ylim=c(-1,1), main = "MA plot ( after remove outlers)")
dev.off()


## Analysis using sample after remove outlers 
library(rtracklayer)
annot <- import("Rattus_norvegicus.Rnor_5.0.79_withSTMtitin.gtf")
geneInfo <- as.data.frame(annot@elementMetadata@listData)[,c("gene_id","gene_name")]
geneInfo <- unique(geneInfo)

ddsR <- DESeq(ddsR)
resR <- results(ddsR)
resR <- resR[order(resR$padj),]

pdf("MA_plot_after.pdf",10,10)
plotMA(ddsR,ylim=c(-1,1), main = "MA plot ( after remove outlers)")
dev.off()

##  compare up with wt 
resR_up_wt <- results(ddsR, contrast=c("Genotype","TTN_up","WT"))
resR_up_wt <- resR_up_wt[order(resR_up_wt$pvalue),]
up_wt <- as.data.frame(resR_up_wt)

up_wt <- up_wt[!is.na(up_wt$padj),]
up_wt$gene_id <- rownames(up_wt)
up_wt <- merge(geneInfo,up_wt,by="gene_id")
up_wt <- up_wt[order(up_wt$padj),]
up_wt_sig <- up_wt[up_wt$padj <0.05,]
write.csv(up_wt,"result_UP_WT.csv",row.names = FALSE)
write.csv(up_wt_sig,"result_UP_WT_significant.csv", row.names = FALSE)

## compare down with wt 
resR_down_wt <- results(ddsR, contrast=c("Genotype","TTN_down","WT"))
resR_down_wt <- resR_down_wt[order(resR_down_wt$pvalue),]
down_wt <- as.data.frame(resR_down_wt)

down_wt$gene_id <- rownames(down_wt)
down_wt <- merge(geneInfo, down_wt, by = "gene_id")
down_wt <- down_wt[!is.na(down_wt$padj),]
down_wt <- down_wt[order(down_wt$padj),]
down_wt_sig <- down_wt[down_wt$padj < 0.05, ]
write.csv(down_wt, "result_DOWN_WT.csv", row.names  = FALSE )
write.csv(down_wt_sig, "result_DOWN_WT_significant.csv", row.names = FALSE)



## compare MT with WT
colData(ddsR)$Condition <- as.factor(colData(ddsR)$Condition)
design(ddsR) <-formula(~Condition)
ddsMT_WT <- DESeq(ddsR)
resMT_WT <- results(ddsMT_WT,contrast = c("Condition","MT","WT"))
resMT_WT <- resMT_WT[order(resMT_WT$pvalue),]
MT_WT <- as.data.frame(resMT_WT)

MT_WT$gene_id <- rownames(MT_WT)
MT_WT <- merge(geneInfo, MT_WT, by = "gene_id")
MT_WT <- MT_WT[order(MT_WT$padj),]
MT_WT <- MT_WT[!is.na(MT_WT$padj),]
MT_WT_sig <- MT_WT[MT_WT$padj < 0.05,]
write.csv(MT_WT,"result_MT_WT.csv")
write.csv(MT_WT_sig, "result_WT_MT_significant.csv") 



resR_up_down <- results(ddsR, contrast=c("Genotype","TTN_up","TTN_down"))
resR_up_down <- resR_up_down[order(resR_up_down$pvalue),]
up_down <- as.data.frame(resR_up_down)

up_down$gene_id <- rownames(up_down)
up_down <- merge(geneInfo, up_down, by = "gene_id")
up_down <- up_down[!is.na(up_down$padj),]
up_down <- up_down[order(up_down$padj),]
up_down_sig <- up_down[up_down$padj < 0.05, ][1,]
write.csv(up_down, "result_UP_DOWN.csv", row.names  = FALSE )
write.csv(up_down_sig, "result_UP_DOWN_significant.csv", row.names = FALSE)





TTN_up_wt <- as.data.frame(resR_up_wt["ENSRNOG99999999999",])
TTN_down_wt <- as.data.frame(resR_down_wt["ENSRNOG99999999999",]) 
TTN_MT_wt <- as.data.frame(resMT_WT["ENSRNOG99999999999",]) 
TTN_up_down <-as.data.frame(resR_up_down["ENSRNOG99999999999",]) 

TTN_stat <- rbind(TTN_up_wt,TTN_down_wt, TTN_MT_wt, TTN_up_down) 
write.csv(TTN_stat, "TTN_statistic.csv")
rownames(TTN_stat) <- c("UP_WT","DOWN_WT", "MT_WT", "UP_DOWN")
#TTN_expression <- counts(ddsR, normalized = TRUE)["ENSRNOG99999999999",]
TTN_expression <- counts(dds, normalized = TRUE)["ENSRNOG99999999999",]
# col <- as.data.frame(colData(ddsR)) 

TTN_df <- data.frame(sampleID = names(TTN_expression), value = TTN_expression, group = col[names(TTN_expression),"Genotype"] )




library(ggplot2)
library(rtracklayer)
library(plyr)

raw_count <- counts(ddsR)
#raw_count <- counts(dds)
col <- as.data.frame(colData(ddsR))
annot_flat <- import("Rnor_5.0.79_withSTMtitin_merged.gtf")
annot_exon <- as.data.frame(annot_flat[annot_flat$type == "exonic_part",]) 
gene_length <- ddply(annot_exon, .(gene_id),summarize, length = sum(width))
rownames(gene_length) <- gene_length$gene_id

gene_overlap <- gene_length[rownames(raw_count),] 
library_size <- apply(raw_count,2,sum)
fpkm <- raw_count *1000000000/library_size/gene_overlap$length
TTN_fpkm <- fpkm["ENSRNOG99999999999",]
TTN_fpkm_df <- data.frame(sampleID = names(TTN_fpkm), value = TTN_fpkm, GenoType = col[names(TTN_fpkm),"Genotype"] )


ggplot(data = TTN_fpkm_df, aes(x = GenoType, y = value,colour = GenoType)) + 
        scale_x_discrete(label = c("A+/- (TTN down)","Z+/- (TTN up)","WT")) + 
        geom_point(size = 5,position = position_jitter(width = 0.1)) + 
        ylab("RNA level(FPKM)") + 
        xlab("Genotype") + 
        ylim(1,50)+
        theme(axis.text = element_text(size = 15, colour = "black"), 
              axis.title = element_text(size = 20, colour = "black"),
              panel.background = element_blank(),
              panel.grid.major = element_line(linetype= "dashed",colour = "gray"),
              legend.position="none")

ggsave("pTTN_scaled.pdf")



## without removing outliers 
raw_count <- counts(dds)
#raw_count <- counts(dds)
annot_flat <- import("Rnor_5.0.79_withSTMtitin_merged.gtf")
annot_exon <- as.data.frame(annot_flat[annot_flat$type == "exonic_part",]) 
gene_length <- ddply(annot_exon, .(gene_id),summarize, length = sum(width))
rownames(gene_length) <- gene_length$gene_id
col <- as.data.frame(colData(dds))
gene_overlap <- gene_length[rownames(raw_count),] 
library_size <- apply(raw_count,2,sum)
fpkm <- t(t(raw_count *1000000000)/library_size)/gene_overlap$length
TTN_fpkm <- fpkm["ENSRNOG99999999999",]
TTN_fpkm_df <- data.frame(sampleID = names(TTN_fpkm), value = TTN_fpkm, GenoType = col[names(TTN_fpkm),"Genotype"] )




ggplot(data = TTN_fpkm_df, aes(x = GenoType, y = value,colour = GenoType)) + 
        scale_x_discrete(label = c("A+/- (TTN down)","Z+/- (TTN up)","WT")) + 
        geom_point(size = 5,position = position_jitter(width = 0.1)) + 
        ylab("RNA level(FPKM)") + 
        xlab("Genotype") + 
        ylim(1,50)+
        theme(axis.text = element_text(size = 15, colour = "black"), 
              axis.title = element_text(size = 20, colour = "black"),
              panel.background = element_blank(),
              panel.grid.major = element_line(linetype= "dashed",colour = "gray"),
              legend.position="none")



ggsave("pTTN_scaled.pdf")
ggsave("pTTN_scaled.tiff")




md <- aov(value ~ GenoType, data = TTN_fpkm_df )
TTN_aov <- summary(md)
TTN_tukey <- TukeyHSD(md)
write.csv(TTN_tukey$GenoType, "TTN_tukey.csv")
write.csv(TTN_aov[[1]], "TTN_anov.csv")







MT_WT_sig$gene_id 
down_wt_sig$gene_id 
up_wt_sig$gene_id 


library(rtracklayer)
library(plyr)
annot_flat <- import("gene.merged.gtf")
annot_exon <- as.data.frame(annot_flat[annot_flat$type == "exonic_part",]) 
gene_length <- ddply(annot_exon, .(gene_id),summarize, length = sum(width))
gene_overlap <- gene_length[gene_length$gene_id %in% rownames(raw_count),] 


raw_count <- counts(ddsMT_WT)[gene_overlap$gene_id,]
library_size <- apply(raw_count,2,sum)
fpkm <- t(t(raw_count *1000000000)/library_size)/gene_overlap$length

write.csv(fpkm,"fpkm.csv")
write.csv(counts(ddsMT_WT,normalized = TRUE), "counts.csv") 









