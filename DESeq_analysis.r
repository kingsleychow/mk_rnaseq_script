library(DESeq2)
library(plyr)
library(gplots)
## Only for extimating and nurmalize the counts  
## creating a heatmap between cases and control
## median ratio method

sampleFiles <- list.files(pattern = "*.count") 
sampleNames <- sub(".count","",sampleFiles)
sampleNames <- sub("^Sample_","",sampleNames)
sampleNames <- sub("B$","",sampleNames)

pheno <- read.csv("../TOF.csv",stringsAsFactors = FALSE)[1:27,]
phenoCtr <- read.csv("../RV_control.csv",stringsAsFactors = FALSE)
phenoUse <- pheno[,c("No",
                     "Age_at_Operation",
                     "Sex",
                     "Remol_RVEF",
                     "Remodel_RVvols",
                     "Restrictive.physiology")]
phenoUse$QRS.duration <- pheno$QRS.duration > 160 
names(phenoUse)[1:3] <- c("sampleName","Age","Gender")
names(phenoCtr)[1] <- "sampleName"


condition <- c(rep("case",27),rep("control",11))
sampleTable <- data.frame(sampleName = sampleNames, 
                          fileName = sampleFiles,
                          condition = condition) 

## selecting the all case sample 
sampleCase <- merge(sampleTable[sampleTable$condition == "case",],phenoUse,all = TRUE,by = "sampleName")

## selecting all contral sample 
sampleCtr <- merge(sampleTable[sampleTable$condition == "control",], phenoCtr, all = TRUE, by = "sampleName")

## all the sample 
sampleTable <- rbind.fill(sampleCase, sampleCtr)

## reading the pheno type data 
sampleCase$Remol_RVEF <- factor(sampleCase$Remol_RVEF)
sampleCase$Remodel_RVvols <- factor(sampleCase$Remodel_RVvols)
sampleCase$Gender <- factor(sampleCase$Gender)
sampleCase$Restrictive.physiology <- factor(sampleCase$Restrictive.physiology)
sampleCase$QRS.duration <- factor(sampleCase$QRS.duration)

# ddsCase <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase[-5,],directory=".",
#                                       design=~Remol_RVEF + Age + Gender) 




## All the case data 
ddsCase <- estimateSizeFactors(ddsCase )
df.counts <- counts(ddsCase, normalized = TRUE)
d <- dist(t(df.counts))  
mat <- as.matrix(d)



library("RColorBrewer")
rld <- varianceStabilizingTransformation(ddsCase, blind=TRUE)

d <- dist(t(assay(rld)))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
mat <- as.matrix(d)
rownames(mat) <- colnames(mat) <- paste(rownames(colData(rld)),paste(":",colData(rld)$Remodel_RVvols))

jpeg("heatmap_case_cluster.jpeg",600,600)
heatmap.2(mat, trace="none", col =rev(hmcol),margin=c(13, 13))
dev.off()

pdf("heatmap_cluster.pdf",12,12)
heatmap.2(mat, trace="none", col =rev(hmcol),margin=c(13, 13))
dev.off()

pdf("pca_cluster.pdf",12,12)
print(plotPCA(rld, intgroup = ""))
dev.off()

countMatrix <- counts(dds, normalized = TRUE)
write.csv(countMatrix,"countMatrix.csv")






##  comparision base on Remol_RVEL phenotype, 0 is the base control 
ddsCase <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase[-5,],directory=".",
                                      design=~Remol_RVEF + Age + Gender) 
ddsCase <- DESeq(ddsCase)
res_RVEF <- results(ddsCase,contrast = c("Remol_RVEF","1","0"))
res_RVEF <- res_RVEF[order(res_RVEF$pvalue),]

RVEF_sig <- data.frame(res_RVEF[res_RVEF$pvalue < 0.05 & !is.na(res_RVEF$pvalue),]) 
RVEF_Gsig <- data.frame(res_RVEF[res_RVEF$padj < 0.05 & !is.na(res_RVEF$padj),]) 
write.csv(RVEF_sig,"RVEF_sig.csv",row.names = TRUE)
write.csv(RVEF_Gsig, "RVEF_Gsig.csv",row.names = TRUE)



## comparison based on Remodel_RVvols, 0 is the base control 
ddsCase2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase[-5,],directory=".",
                                       design=~Remodel_RVvols + Age + Gender) 
ddsCase2 <- DESeq(ddsCase2)
res_RVvols <- results(ddsCase2,contrast = c("Remodel_RVvols","1","0"))
res_RVvols <- res_RVvols[order(res_RVvols$pvalue),]


RVvols_sig <- data.frame(res_RVvols[res_RVvols$pvalue < 0.05 & !is.na(res_RVvols$pvalue),]) 
RVvols_Gsig <- data.frame(res_RVvols[res_RVvols$padj < 0.05 & !is.na(res_RVvols$padj),])
write.csv(RVvols_sig, "RVvols_sig.csv",row.names = TRUE)
write.csv(RVvols_Gsig, "RVvols_Gsig.csv", row.names = TRUE)







## comparion based on QRS.duration , FALSE is the based control 
ddsCase4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase,directory=".",
                                       design=~QRS.duration + Age + Gender)
ddsCase4 <- DESeq(ddsCase4)
res_QRS <- results(ddsCase4,contrast = c("QRS.duration","TRUE","FALSE"))
res_QRS <- res_QRS[order(res_QRS$pvalue),]
QRS_sig <- data.frame(res_QRS[res_QRS$pvalue < 0.05 & !is.na(res_QRS$pvalue),])
QRS_Gsig <- data.frame(res_QRS[res_QRS$padj < 0.05 & !is.na(res_QRS$padj),])
write.csv(QRS_sig, "QRS_sig.csv", row.names = TRUE)
write.csv(QRS_Gsig, "QRS_Gsig.csv", row.names =TRUE)


ddsCase3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase[-20,],directory=".",
                                       design=~Restrictive.physiology + Age + Gender)





##  comparsion based patient and healthy control 
condition <- c(rep("case",27),rep("control",11))
sampleTable <- data.frame(sampleName = sampleNames, 
                          fileName = sampleFiles,
                          condition = condition) 
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory=".",
                                  design=~condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$pvalue),]



dataCountsNormalized <- counts(dds, normalized = TRUE)
write.csv(dataCountsNormalized, "CountsNormalized.csv")

##  loading phenotype data 
mainPhe <- read.csv("main_pheno.csv",stringsAsFactors = FALSE)[1:27,]
cmr <- read.csv("CMR.csv",stringsAsFactors = FALSE)
general <- read.csv("general_check.csv",stringsAsFactors = FALSE)
cmrSelect <- read.csv("cmr_regress.csv",stringsAsFactors = FALSE)
meta <- read.csv("meta_data.csv",stringsAsFactors = FALSE)[1:27,]
meta$Sex <- as.factor(meta$Sex) 

datCmr <- cmr[,c("No","testing_bench",cmrSelect$phe[cmrSelect$sig])]
names(datCmr) <- gsub("_A","",names(datCmr))

binarPhe <- mainPhe[,c("No","Remodel_RVvols",
                       "Remol_RVEF",
                       "Restrictive.physiology")]

numericPhe <- mainPhe[,c("No","QRS.duration",
                         "Pre_RV_LGE_score",
                         "Pre_LV_LGE_score",
                         "RVOT_scar_sag_Pre",
                         "Post_RV_LGE_score",
                         "Post_LV_LGE_score",
                         "RVOT_scar_sag_post",
                         "RV_fibrosis")]


datCmr <- merge(datCmr, meta[,c(1,2,3)],by="No")
binarPhe <- merge(binarPhe, meta[,c(1,2,3)],by="No")
numericPhe <- merge(numericPhe, meta[,c(1,2,3)],by="No")
datCt <- t(dataCountsNormalized)
rownames(datCt) <- gsub("Sample_","",rownames(datCt))
rownames(datCt) <- gsub("B$","",rownames(datCt))
datCase <- datCt[as.character(numericPhe$No),]

datCaseMean <- apply(datCase, 2, mean )
datCaseNoZero <- datCase[,datCaseMean != 0]

resBinr <- list() 
for(gene in colnames(datCaseNoZero)){
        tmpRes <- list()
        for(phen in names(binarPhe)[2:4]){
                tmp <- summary(glm(binarPhe[,phen] ~ datCaseNoZero[,gene] + 
                                           binarPhe[,"Sex"] + 
                                           binarPhe[,"Age_at_Operation"],family=binomial(link="logit")))
                print(nrow(tmp$coefficients)) 
                if(nrow(tmp$coefficients) < 4) { 
                        tmpRes[[phen]] <- NA
                } else tmpRes[[phen]] <- tmp$coefficients["datCaseNoZero[, gene]",4]
        }
        resBinr[[gene]] <- do.call(c,tmpRes)
}

resBinDf <- do.call(rbind,resBinr)
write.csv(resBinDf,"result_binary_data.csv")

resNumeric <- list()
for(gene in colnames(datCaseNoZero)){
        tmpRes <- list()
        for(phen in names(numericPhe)[2:9]){
                tmp <- summary(glm(numericPhe[,phen] ~ datCaseNoZero[,gene] + 
                                           numericPhe[,"Sex"] + 
                                           numericPhe[,"Age_at_Operation"],
                                   family = gaussian("identity")))
                print(nrow(tmp$coefficients)) 
                if(nrow(tmp$coefficients) < 4) { 
                        tmpRes[[phen]] <- NA
                }else tmpRes[[phen]] <- tmp$coefficients["datCaseNoZero[, gene]",4]
        }
        resNumeric[[gene]] <- do.call(c,tmpRes)
}
resNumeric <- do.call(rbind,resNumeric)
write.csv(resNumeric, "result_numeric_data.csv")



datCaseBatch <- data.frame(datCaseNoZero) 
datCaseBatch$No <- rownames(datCaseNoZero)
datCaseBatch <- merge(datCmr, datCaseBatch, by = "No")





resCmr <- list()
for(gene in colnames(datCaseNoZero)){
        tmpRes <- list()
        for(phen in names(datCmr)[3:13]){
                tmp <- summary(glm(datCaseBatch[,phen] ~ datCaseBatch[,gene] + 
                                           datCaseBatch[,"testing_bench"] + 
                                           datCaseBatch[,"Sex"] + 
                                           datCaseBatch[,"Age_at_Operation"],
                                   family = gaussian("identity")))
                print(nrow(tmp$coefficients))
                print(phen)
                print(gene) 
                if(nrow(tmp$coefficients) < 6) { 
                        tmpRes[[phen]] <- NA 
                }else tmpRes[[phen]] <- tmp$coefficients["datCaseBatch[, gene]",4]
        }
        resCmr[[gene]] <- do.call(c,tmpRes)
}
resCmr <- do.call(rbind,resCmr)
write.csv(resCmr,"result_Cmr.csv")



res <- res[order(res$padj),]
res.df <- data.frame(res)
res.df.Gsig <- res.df[res.df$padj < 0.05 & !is.na(res.df$padj), ] 
res.df.nona <- res.df[!is.na(res.df$padj),]

write.csv(res.df,"case_control_Gsig.csv")
write.csv(res.df, "case_control.csv")


gene_discription <- read.table("./gene_type/gene_type.csv", sep = "\t",stringsAsFactors = FALSE)





rldcase_control <- varianceStabilizingTransformation(dds, blind=TRUE)

dcase_control <- dist(t(assay(rldcase_control)))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
matcase_control <- as.matrix(dcase_control)
rownames(matcase_control) <- colnames(matcase_control) <- paste(rownames(colData(rld)),paste(":",colData(rldcase_control)$condition  ))

jpeg("heatmap_case_control_cluster.jpeg",600,600)
heatmap.2(matcase_control, trace="none", col =rev(hmcol),margin=c(13, 13))
dev.off()





## comparsion based on restrictive 1 and none restrictive 0, 0 is the based control 
ddsCase3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleCase[-20,],directory=".",
                                        design=~Restrictive.physiology + Age + Gender)
ddsCase3 <- DESeq(ddsCase3)
res_Restrictive <- results(ddsCase3,contrast = c("Restrictive.physiology","1","0"))
res_Restrictive <- res_Restrictive[order(res_Restrictive$pvalue),]
Restrictive_sig <- data.frame(res_Restrictive[res_Restrictive$pvalue < 0.05 & !is.na(res_Restrictive$pvalue),])
Restrictive_Gsig <- data.frame(res_Restrictive[res_Restrictive$padj < 0.05 & !is.na(res_Restrictive$padj),])
write.csv(Restrictive_sig, "Restrictive_sig.csv", row.names = TRUE)
write.csv(Restrictive_Gsig, "Restrictive_Gsig.csv", row.names =TRUE)





### pathway analysis  between restrictive and none restrictive 

count_restric <- counts(ddsCase3,normalized = TRUE )

library(gage)
kg.human <- kegg.gsets("human", id.type = "entrez")
kegg.gs <- kg.human$kg.sets 


library(biomaRt)
hsMart <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
hsId <- getBM(mart = hsMart,attributes = c("ensembl_gene_id", "entrezgene"))



ensID <- hsId$ensembl_gene_id
names(ensID) <- hsId$entrezgene 
kegg.gs.ensID <- lapply(kegg.gs, function(x) ensID[x])


                        
                        
tmp <-  colData(ddsCase3)



res.gs <- gage(count_restric, 
               gsets = kegg.gs.ensID, 
               ref = which( tmp$Restrictive.physiology != 1), 
               samp = which( tmp$Restrictive.physiology == 1), 
               compare = "unpaired",
               same.dir = FALSE
)


res.kegg <- data.frame(res.gs$greater ) 

res.kegg.Gsig <- res.kegg[ res.kegg$q.val < 0.05 & !is.na(res.kegg$q.val), ]
write.csv(res.kegg.Gsig, "restrictive.kegg.Gsig.csv ")


go.hm <- go.gsets("human",id.type = "EG")
goset.hm <- go.hm$go.sets 
goset.hm.ensID <-  lapply(goset.hm, function(x) ensID[x])

res.go <- gage(count_restric, 
               gsets = goset.hm.ensID, 
               ref = which( tmp$Restrictive.physiology != 1), 
               samp = which( tmp$Restrictive.physiology == 1), 
               compare = "unpaired",
               same.dir = FALSE
)
res.go <- data.frame(res.go$greater ) 
res.go <- res.go[order(res.go$q.val),]

res.go.Gsig <- res.go[ res.go$q.val < 0.05 & !is.na(res.go$q.val), ]
write.csv(res.go.Gsig, "restrictive.go.Gsig.csv")


countsMatrix <- as.matrix(ratIdCounts[,c(-1,-2)])
rownames(countsMatrix) <- ratIdCounts$entrezgene



# > sessionInfo()
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# 
# locale:
#         [1] LC_COLLATE=English_Singapore.1252  LC_CTYPE=English_Singapore.1252    LC_MONETARY=English_Singapore.1252
# [4] LC_NUMERIC=C                       LC_TIME=English_Singapore.1252    
# 
# attached base packages:
#         [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#         [1] gplots_2.16.0             plyr_1.8.1                DESeq2_1.6.3              RcppArmadillo_0.4.600.4.0
# [5] Rcpp_0.11.5               GenomicRanges_1.18.4      GenomeInfoDb_1.2.4        IRanges_2.0.1            
# [9] S4Vectors_0.4.0           BiocGenerics_0.12.1       biomaRt_2.22.0            gage_2.16.0              
# [13] reshape2_1.4.1           
# 
# loaded via a namespace (and not attached):
#         [1] acepack_1.3-3.3      annotate_1.44.0      AnnotationDbi_1.28.1 base64enc_0.1-2      BatchJobs_1.5        BBmisc_1.9          
# [7] Biobase_2.26.0       BiocParallel_1.0.3   Biostrings_2.34.1    bitops_1.0-6         brew_1.0-6           caTools_1.17.1      
# [13] checkmate_1.5.1      cluster_1.15.2       codetools_0.2-8      colorspace_1.2-4     DBI_0.3.1            digest_0.6.8        
# [19] fail_1.2             foreach_1.4.2        foreign_0.8-61       Formula_1.2-0        gdata_2.13.3         genefilter_1.48.1   
# [25] geneplotter_1.44.0   ggplot2_1.0.0        graph_1.44.1         grid_3.1.1           gtable_0.1.2         gtools_3.4.1        
# [31] Hmisc_3.15-0         httr_0.6.1           iterators_1.0.7      KEGGREST_1.6.4       KernSmooth_2.23-12   lattice_0.20-29     
# [37] latticeExtra_0.6-26  locfit_1.5-9.1       MASS_7.3-33          munsell_0.4.2        nnet_7.3-8           png_0.1-7           
# [43] proto_0.3-10         RColorBrewer_1.1-2   RCurl_1.95-4.5       rpart_4.1-8          RSQLite_1.0.0        scales_0.2.4        
# [49] sendmailR_1.2-1      splines_3.1.1        stringr_0.6.2        survival_2.37-7      tools_3.1.1          XML_3.98-1.1        
# [55] xtable_1.7-4         XVector_0.6.0        zlibbioc_1.12.0     






