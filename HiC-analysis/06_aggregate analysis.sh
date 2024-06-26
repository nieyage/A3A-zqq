# The APA plot in Loops:

# Hi-C aggregate analysis
usage: fanc aggregate [-h] [-m MATRIX_FILE] [-p PLOT_FILE] [--tads]
                      [--tads-imakaev] [--loops]
                      [--loop-strength LOOP_STRENGTH_FILE]
                      [--tad-strength TAD_STRENGTH_FILE] [-w WINDOW]
                      [--pixels PIXELS] [-v REGION_VIEWPOINT]
                      [-b BOUNDARY_MODE] [-i INTERPOLATION] [-r RELATIVE]
                      [-a ABSOLUTE] [-e] [-l] [--rescale]
                      [--colormap COLORMAP] [--vmin VMIN] [--vmax VMAX] [-tmp]
                      [-C] [--keep-submatrices] [-s] [--lower-triangular-plot]
                      [--labels LABELS] [--label-locations LABEL_LOCATIONS]
                      input [regions] [output]

cd /md01/nieyg/project/A3A-Hic/06_aggregate_analysis
conda activate fanc
A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/A3A.hic
A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/A3A3.hic
KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/KO.hic
WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/WT.hic
# region file:
#global_loop:
global=/md01/nieyg/project/A3A-Hic/05_Loop_calling/WT_merged.loops.bedpe

# get the A3A binding region
scp -r /public/home/chenxy/ChIP-seq-cuttag_data/Mr.zhao/A3A/A3A-IgG_peaks.narrowPeak nieyg@202.116.90.56:/md01/nieyg/project/A3A-Hic/05_Loop_calling/
conda activate fanc

# annotation the loop file 
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
#library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
loops <- read.table("WT_10kb_merged.loops.bedpe")   
peaks <- read.table("A3A-IgG_peaks.narrowPeak")  
peaks_gr<- GRanges(seqnames=peaks$V1, ranges=IRanges(start=peaks$V2, end=peaks$V3))
loops_end1 <- GRanges(seqnames=loops$V1, ranges=IRanges(start=loops$V2, end=loops$V3))  
loops_end2 <- GRanges(seqnames=loops$V4, ranges=IRanges(start=loops$V5, end=loops$V6))  
overlaps_end1 <- findOverlaps(loops_end1, peaks_gr, type="any")  
overlaps_end2 <- findOverlaps(loops_end2, peaks_gr, type="any")  
overlapping_loops_indices <- intersect(queryHits(overlaps_end1),queryHits(overlaps_end2))
overlapping_loops <- loops[overlapping_loops_indices, ]  
A3A_mediated_loop<- as.data.frame(overlapping_loops)
write.table(A3A_mediated_loop,"A3A_binding_WT_loops.bedpe",sep="\t",row.names=F,col.names=F)

# get the cardiac develop related gene locations 
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="gene",]
gene_related_heart_morphogenesis<-getGO("GO:0003007")
gene_related_cardiac_development<-getGO("GO:0048738")
heart_genes <- unique(unlist(c(gene_related_heart_morphogenesis,gene_related_cardiac_development)))
heart_gene_ranges <- gtf[gtf$gene_name %in% heart_genes]
overlaps_end1 <- findOverlaps(loops_end1, heart_gene_ranges, type="any")  
overlaps_end2 <- findOverlaps(loops_end2, heart_gene_ranges, type="any")  
Cardiac_develop_loops_indices <- unique(c(queryHits(overlaps_end1), queryHits(overlaps_end2)))  
Cardiac_develop_loops <- loops[Cardiac_develop_loops_indices, ]  
Cardiac_develop_loops<- as.data.frame(Cardiac_develop_loops)
write.table(Cardiac_develop_loops,"Cardiac_develop_loops.bedpe",sep="\t",row.names=F,col.names=F)

# Get the SE-P loops 

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=3000, downstream=3000)
SE <- read.table("WTac-1_peaks_Gateway_SuperEnhancers.bed")   
SE_gr<- GRanges(seqnames=SE$V1, ranges=IRanges(start=SE$V2, end=SE$V3))

overlaps_end1 <- findOverlaps(loops_end1, SE_gr, type="any")  
overlaps_end2 <- findOverlaps(loops_end2, promoter, type="any")  
SE_loops_indices1 <- intersect(queryHits(overlaps_end1),queryHits(overlaps_end2))
overlaps_end1 <- findOverlaps(loops_end1, promoter, type="any")  
overlaps_end2 <- findOverlaps(loops_end2, SE_gr, type="any")  
SE_loops_indices2 <- intersect(queryHits(overlaps_end1),queryHits(overlaps_end2))
SE_loops_indices<- c(SE_loops_indices1,SE_loops_indices2)
SE_P_loops <- loops[SE_loops_indices, ]  
write.table(SE_P_loops,"SE_P_loops.bedpe",sep="\t",row.names=F,col.names=F)

# the upsetR for the four indices 
library(UpSetR)

listInput <- list(
        All_loops = 1:2571, 
        A3A_mediated_loops = overlapping_loops_indices, 
        Cardiac_development_loops = Cardiac_develop_loops_indices, 
        SE_loops = SE_loops_indices)
data<- fromList(listInput)

pdf("./loops-upsetR.pdf", width=8, height=4)
UpSetR::upset(data)
dev.off()

A3A_SE_P_loop<- intersect(overlapping_loops_indices,SE_loops_indices)
A3A_SE_P_loop <- loops[A3A_SE_P_loop, ]  
write.table(A3A_SE_P_loop,"A3A_SE_P_loops.bedpe",sep="\t",row.names=F,col.names=F)

A3A_SE_P_loop<- intersect(overlapping_loops_indices,Cardiac_develop_loops_indices)

intersect(A3A_SE_P_loop,Cardiac_develop_loops_indices)


sed -i 's/"//g' A3A_binding_WT_loops.bedpe

# Loops and other pairwise genomic regions
fanc aggregate /md01/nieyg/project/A3A-Hic/05_Loop_calling/A3A_100kb_merged.loops \
               WT_100kb_merged.loops.bedpe \
               A3A.loops_no_singlets.agg \
               -p A3A.loops_no_singlets.agg.pdf \
               --loops -e -l -r 1.0

fanc aggregate /md01/nieyg/project/A3A-Hic/05_Loop_calling/A3A3_100kb_merged.loops \
               WT_100kb_merged.loops.bedpe \
               A3A.loops_no_singlets.agg \
               -p A3A3.loops_no_singlets.agg.pdf \
               --loops -e -l -r 1.0

# # Aggregate TADs
# fanc aggregate $A3A \
#                architecture/domains/gm12878_tads.bed \
#                architecture/aggregate/fanc_example_100kb.agg \
#                -p architecture/aggregate/fanc_example_100kb_oe.agg.png \
#                -m architecture/aggregate/fanc_example_100kb_oe.agg.txt \
#                -e -l -r 1.0
sed -i '/chrUn_gl000220/d;' A3A_10kb_merged.loops.bedpe

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic \
               WT_10kb_merged.loops.bedpe \
               KO_10kb_merged_no_singlets.agg \
               -p KO_KOloop_10kb_merged_no_singlets.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic \
               WT_10kb_merged.loops.bedpe \
               WT_10kb_merged_no_singlets.agg \
               -p WT_10kb_merged_no_singlets.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic \
               WT_10kb_merged.loops.bedpe \
               A3A_A3A_10kb_merged_no_singlets.agg \
               -p A3A_A3A_10kb_merged_no_singlets.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic \
               WT_10kb_merged.loops.bedpe \
               A3A3_10kb_merged_no_singlets.agg \
               -p A3A3_10kb_merged_no_singlets.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

# The A3A binding region 

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic \
               A3A_binding_WT_loops.bedpe \
               KO_10kb_merged_no_singlets.agg \
               -p KO_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic \
               A3A_binding_WT_loops.bedpe \
               WT_10kb_merged_no_singlets.agg \
               -p WT_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic \
               A3A_binding_WT_loops.bedpe \
               A3A_A3A_10kb_merged_no_singlets.agg \
               -p A3A_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic \
               A3A_binding_WT_loops.bedpe \
               A3A3_10kb_merged_no_singlets.agg \
               -p A3A3_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2


fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic \
               A3A_SE_P_loops.bedpe \
               KO_10kb_merged_no_singlets.agg \
               -p KO_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic \
               A3A_SE_P_loops.bedpe \
               WT_10kb_merged_no_singlets.agg \
               -p WT_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic \
               A3A_SE_P_loops.bedpe \
               A3A_A3A_10kb_merged_no_singlets.agg \
               -p A3A_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2

fanc aggregate /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic \
               A3A_SE_P_loops.bedpe \
               A3A3_10kb_merged_no_singlets.agg \
               -p A3A3_A3A_loops.agg.pdf \
               --loops  -e -l  --rescale -v center -i 1 --vmin 0 --vmax 2


fancplot chr11:1mb-6mb -p triangular -vmax 0.05 /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic -o A3A3.png
fancplot chr11:1mb-6mb -p triangular -c Reds -vmax 0.05 /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic
fancplot chr11:1mb-6mb -p triangular -c Reds -vmax 0.05 /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic
fancplot chr11:1mb-6mb -p triangular -c Reds -vmax 0.05 /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic


fancplot chr11:2.2mb-3mb -o A3A3.png -p triangular -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic
fancplot chr11:2.2mb-3mb -o A3A.png -p triangular -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic
fancplot chr11:2.2mb-3mb -o KO.png -p triangular -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic
fancplot chr11:2.2mb-3mb -o WT.png -p triangular -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic


# TH
chr11:2.2mb-3mb

# MYH7 chr14 23433876 23476124 chr14 23503876 23546124  . 3.063005e-11
chr14:23433876-23546124
fancplot chr14:23.3MB-23.6MB -o MYH7_A3A3.png -p triangular -c Reds -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic
fancplot chr14:23.3MB-23.6MB  -o MYH7_A3A.png -p triangular -c Reds -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic
fancplot chr14:23.3MB-23.6MB  -o MYH7_KO.png -p triangular -c Reds -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic
fancplot chr14:23.3MB-23.6MB  -o MYH7_WT.png -p triangular -c Reds -vmin 0 -vmax 0.05  /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10
# MYL3 chr3 46722501 46767500 chr3 46872501 46917500  . 0.07323615
fancplot chr3:46.6MB-47MB -o MYL3_A3A3.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic
fancplot chr3:46.6MB-47MB -o MYL3_A3A.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic
fancplot chr3:46.6MB-47MB -o MYL3_KO.png -p triangular -vmin 0 -vmax 0.05  -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic
fancplot chr3:46.6MB-47MB -o MYL3_WT.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic

# RBM15 chr1 110310001 110340000 chr1 110510001 110540000  . 3.530003e-05

fancplot chr1:110.30MB-110.6MB -o RBM15_A3A3.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic
fancplot chr1:110.30MB-110.6MB -o RBM15_A3A.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic
fancplot chr1:110.30MB-110.6MB -o RBM15_KO.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic
fancplot chr1:110.30MB-110.6MB -o RBM15_WT.png -p triangular -vmin 0 -vmax 0.05 -c Reds /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic


