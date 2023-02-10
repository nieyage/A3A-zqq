#Step1: compare with input 
macs2 callpeak -t ../2_bam/A3A_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n A3A-D5-IG 

macs2 callpeak -t ../2_bam/KOac-1_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n KOac-1-D5-IG 
macs2 callpeak -t ../2_bam/KOac-2_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n KOac-2-D5-IG 

nohup macs2 callpeak -t ../2_bam/WTac-1_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n WTac-1-D5-IG &
nohup macs2 callpeak -t ../2_bam/WTac-2_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n WTac-2-D5-IG &

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./A3A-D5-IG_peaks.narrowPeak > A3A-D5-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./KOac-1-D5-IG_peaks.narrowPeak > KOac-1-D5-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./KOac-2-D5-IG_peaks.narrowPeak > KOac-2-D5-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./WTac-1-D5-IG_peaks.narrowPeak > WTac-1-D5-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./WTac-2-D5-IG_peaks.narrowPeak > WTac-2-D5-IG_peaks_homer.tmp

nohup findMotifsGenome.pl A3A-D5-IG_peaks_homer.tmp  hg19 A3A-D5-IG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl KOac-1-D5-IG_peaks_homer.tmp  hg19 KOac-1-D5-IG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl KOac-2-D5-IG_peaks_homer.tmp  hg19 KOac-2-D5-IG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl WTac-1-D5-IG_peaks_homer.tmp  hg19 WTac-1-D5-IG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl WTac-2-D5-IG_peaks_homer.tmp  hg19 WTac-2-D5-IG_motifDir_Homer -len 8,10,12  &

#Step2: Basical Figure



computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/chip/SA-Chip/rawdata/3.align/hg19_RefSeq.bed  \
-S A3A_CPM_normalized.bw D5-IG_CPM_normalized.bw KOac-1_CPM_normalized.bw KOac-2_CPM_normalized.bw WTac-1_CPM_normalized.bw WTac-2_CPM_normalized.bw \
--skipZeros  -o A3A-CPM.gz  \
--outFileSortedRegions A3A-bedtools-out.bed
plotHeatmap -m A3A-CPM.gz  -out A3A-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m A3A-CPM.gz  -out A3A-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)

peakAnno <- annotatePeak("A3A-D5-IG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("A3A-D5-IG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of A3A-binding loci\nrelative to TSS")

dev.off()

write.csv(as.data.frame(peakAnno),"A3A-cuttag-peak.annotation.csv",row.names = F)

ascp -v -QT -l 300m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR130/023/SRR13079923/SRR13079923_1.fastq.gz ./

#!/bin/bash
for i in $(cat SraAccList.txt)
do 
x=$(echo $i | cut -b1-6)
y=`echo ${i: -2}`
echo "vol1/fastq/${x}/0${y}/${i}/${i}_1.fastq.gz" >>fastqid_trim.txt
echo "vol1/fastq/${x}/0${y}/${i}/${i}_2.fastq.gz" >>fastqid_trim.txt
done

ascp -v -QT -l 300m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp --file-list fastqid_trim.txt ./
rm fastqid_trim.txt

#2022.12.06
#A3A TF cut tag No.3
macs2 callpeak -t ../2_bam/A3A_sorted_rmDup_mapped_rmbl.bam \
               -c ../2_bam/IgG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n A3A-IgG
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' A3A-IgG_peaks.narrowPeak > A3A-IgG_homer.tmp
nohup findMotifsGenome.pl A3A-IgG_homer.tmp  hg19 A3A-IgG_motifDir_Homer -len 6,8,10,12  &


#macs2 callpeak -t ../2_bam/IgG_sorted_rmDup_mapped_rmbl.bam \
#               -c ../2_bam/A3A_sorted_rmDup_mapped_rmbl.bam  \
#               -f BAMPE \
#               -g hs \
#               -n IgG-A3A
#awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' IgG-A3A_peaks.narrowPeak > IgG-A3A_homer.tmp
#nohup findMotifsGenome.pl IgG-A3A_homer.tmp  hg19 IgG-A3A_motifDir_Homer -len 6,8,10,12  &

#awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' A3A_peaks.narrowPeak > A3A_homer.tmp
#nohup findMotifsGenome.pl A3A_homer.tmp  hg19 A3A_motifDir_Homer -len 6  &


computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/reference/ref_gene/hg19_RefSeq.bed  \
-S A3A_CPM_normalized.bw IgG_CPM_normalized.bw \
--skipZeros  -o A3A-IgG-CPM.gz  \
--outFileSortedRegions A3A-bedtools-out.bed
plotHeatmap -m A3A-IgG-CPM.gz  -out A3A-IgG-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m A3A-IgG-CPM.gz  -out A3A-IgG-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)

peakAnno <- annotatePeak("A3A-IgG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("A3A-IgG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of A3A-binding loci\nrelative to TSS")

dev.off()
write.csv(as.data.frame(peakAnno),"A3A-cuttag-peak.annotation.csv",row.names = F)



##MEME 

bedtools getfasta -fi /public/home/nieyg/reference/genome/hg19/hg19.fa \
-bed A3A-IgG_peaks.narrowPeak -fo ARID3A.fa
conda activate meme
# meme denovo motif 
nohup meme  ARID3A.fa -oc meme_denovo -dna -p 12 -maxw 10 -minw 6 -mod zoops -nmotifs 50 -revcomp &

nohup meme-chip \
-meme-p 12 \
-oc meme-chip-A3A \
-maxw 10 \
-db /public/home/nieyg/database/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme \
ARID3A.fa > meme-chip.out &

#conda activate python27
#awk '{print $2"\t"$3"\t"$4}' A3A-IgG_homer.tmp > A3A-IgG_peaks_BETA.bed
#nohup BETA plus \
#    -p A3A-IgG_peaks_BETA.bed -e A3A-IgG_peaks_BETA.txt #DEG file    -k O \
#    --info 1,3,6\
#    -g hg19 \
#    --gs /public/home/nieyg/reference/genome/hg19/hg19.fa \
#    --pn 4000 \
#    --gname2 \
#    -n A3A-IgG_peaks_BETA \
#    --df 0.05 \
#    --da 1 \
#    -o BETA_A3A_rmbg_pvalue > BETA_plus.out &


