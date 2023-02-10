# A3A cut tag H3k27ac downstream analysis #
library(DiffBind)

sampleinfo<-read.csv("/public/home/nieyg/project/A3A/cuttag/20220826-ARID3A-CUTTAG/5_analysis/A3A_input_sampleinfo.csv")
A3A <- dba(sampleSheet=sampleinfo)

A3A <- dba.count(A3A, bUseSummarizeOverlaps=TRUE);
pdf("H3K27ac-PCA.pdf")
dba.plotPCA(A3A,  attributes=DBA_CONDITION, label=DBA_ID)
plot(A3A)
dev.off()
#A3A <- dba.normalize(A3A)

#Establishing a model design and contrast
A3A <- dba.contrast(A3A,categories=DBA_CONDITION,minMembers = 2);
dbObj <- dba.analyze(A3A, method=DBA_ALL_METHODS)
pdf("./A3A-H3k27ac-DiffBind-results.pdf")
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()

res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=0.05,fold=1)

# Write to file

out <- as.data.frame(res_deseq)
write.table(out, file="./WT_vs_A3AKO_deseq2_FDR05FC1.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
library(dplyr)
WT_enrich <- out %>% filter(FDR < 0.05 & Fold > 1) %>%   select(seqnames, start, end)
write.table(WT_enrich, file="WT_enrich.bed", sep="\t", quote=F, row.names=F, col.names=F)
A3AKO_enrich <- out %>%   filter(FDR < 0.05 & Fold < -1) %>%   select(seqnames, start, end)
write.table(A3AKO_enrich, file="A3AKO_enrich.bed", sep="\t", quote=F, row.names=F, col.names=F)

#get the sample normlizecount for heatmap 
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=0.05,fold=1,bCounts=T)
out <- as.data.frame(res_deseq)
data<-out[,12:15];
library(pheatmap);

filtermat=t(scale(t(data),scale = T,center = F))
pdf("./A3A-H3k27ac-DEseq-heatmap.pdf",width=5,heigh=8)
pheatmap(filtermat, cluster_rows=TRUE,
         show_rownames=F,border_color = NA,clustering_method="centroid",
         color = colorRampPalette(colors = c("white","blue"))(100),
         cluster_cols=T,cutree_rows =1,main = "WT vs A3AKO",
         show_colnames = T,fontsize_row = 7.5)
dev.off()




# annotate the peak and gene functional analysis 

# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
# Load data
samplefiles <- list.files("/public/home/nieyg/project/A3A/cuttag/20220826-ARID3A-CUTTAG/5_analysis", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("A3AKO", "WT")

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
pdf("./Anno_peak_distribution.pdf",width=6,heigh=6)
plotAnnoPie(peakAnnoList[["WT"]]);
plotAnnoPie(peakAnnoList[["A3AKO"]]);
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList, title="Distribution of DEP relative to TSS")
dev.off()
# Get annotation data frame
WT_peakAnno <- annotatePeak(samplefiles[["WT"]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
A3AKO_peakAnno <- annotatePeak(samplefiles[["A3AKO"]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
WT_annot <- as.data.frame(WT_peakAnno@anno)
A3AKO_annot <- as.data.frame(A3AKO_peakAnno@anno)
write.csv(WT_annot,"WT-peak-anno.csv")
write.csv(A3AKO_annot,"A3AKO-peak-anno.csv")


pdf("WT-H3K27ac-GO.pdf")

ego <- enrichGO(WT_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"WT-H3K27ac-BP.csv")

ego <- enrichGO(WT_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"WT-H3K27ac-MF.csv")

ego <- enrichGO(WT_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"WT-H3K27ac-CC.csv")
dev.off()
######KEGG##########
pdf("WT-H3K27ac_KEGG.pdf")
ego <- enrichKEGG(
  gene = WT_annot$geneId,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"WT-H3K27ac-KEGG.csv")
dev.off()

pdf("A3AKO-H3K27ac-GO.pdf")
ego <- enrichGO(A3AKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"A3AKO-H3K27ac-BP.csv")

ego <- enrichGO(A3AKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"A3AKO-H3K27ac-MF.csv")

ego <- enrichGO(A3AKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"A3AKO-H3K27ac-CC.csv")
dev.off()
######KEGG##########
pdf("A3AKO-H3K27ac_KEGG.pdf")

ego <- enrichKEGG(
  gene = A3AKO_annot$geneId,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"A3AKO-H3K27ac-KEGG.csv")

dev.off();

# Motif enrichment 

findMotifsGenome.pl A3AKO_enrich.bed  hg19 A3AKO-KO_motifDir_Homer -len 6,8,10,12  
findMotifsGenome.pl WT_enrich.bed  hg19 WT_motifDir_Homer -len 6,8,10,12  

