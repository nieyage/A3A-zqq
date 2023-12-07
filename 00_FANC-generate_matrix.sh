
## PBS configure
#PBS -N hic_fanc
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=60G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/A3A-Hic
conda activate fanc
    fanc auto \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A-1/A3A-1-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A-1/A3A-1-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A-2/A3A-2-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A-2/A3A-2-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A \
    -g /md01/nieyg/ref/genome/hg19/hg19.fa \
    -i /md01/nieyg/ref/bowtie2/hg19/hg19 \
    -n A3A -t 20 -r MboI --iterative --split-ligation-junction -q 30

    fanc auto \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A3-1/A3A3-1-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A3-1/A3A3-1-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A3-2/A3A3-2-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/A3A3-2/A3A3-2-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3 \
    -g /md01/nieyg/ref/genome/hg19/hg19.fa \
    -i /md01/nieyg/ref/bowtie2/hg19/hg19 \
    -n A3A3 -t 20 -r MboI --iterative --split-ligation-junction -q 30

    fanc auto \
    /md01/nieyg/project/A3A-Hic/merged_data/KO-1/KO-1-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/KO-1/KO-1-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/KO-2/KO-2-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/KO-2/KO-2-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/02.fanc_out/KO \
    -g /md01/nieyg/ref/genome/hg19/hg19.fa \
    -i /md01/nieyg/ref/bowtie2/hg19/hg19 \
    -n KO -t 20 -r MboI --iterative --split-ligation-junction -q 30

    fanc auto \
    /md01/nieyg/project/A3A-Hic/merged_data/WT-1/WT-1-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/WT-1/WT-1-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/WT-2/WT-2-merged_R1.fastq.gz \
    /md01/nieyg/project/A3A-Hic/merged_data/WT-2/WT-2-merged_R2.fastq.gz \
    /md01/nieyg/project/A3A-Hic/02.fanc_out/WT \
    -g /md01/nieyg/ref/genome/hg19/hg19.fa \
    -i /md01/nieyg/ref/bowtie2/hg19/hg19 \
    -n WT -t 20 -r MboI --iterative --split-ligation-junction -q 30

echo "All samples processed"