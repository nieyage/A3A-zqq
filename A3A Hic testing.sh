# A3A Hic testing 

# random select some reads for testing 
seqtk sample -s100 /md01/nieyg/project/A3A/01.RawData/A76/A76_clean/A76/A76clean_R1.fastq.gz 10000 > /md01/nieyg/project/A3A/01.RawData/A76/test_clean/test/test_R1.fq
seqtk sample -s100 /md01/nieyg/project/A3A/01.RawData/A76/A76_clean/A76/A76clean_R2.fastq.gz 10000 > /md01/nieyg/project/A3A/01.RawData/A76/test_clean/test/test_R2.fq

# Hic Pro testing 
conda activate /data/R02/nieyg/software/hic-pro
time /data/R02/nieyg/software/HiC-Pro_3.1.0/bin/HiC-Pro -c config_test_latest.txt -i /md01/nieyg/project/A3A/01.RawData/A76/test_clean/ -o hicpro_latest_test



fastp -i A76_1.fq.gz \
        -I A76_2.fq.gz \
        -o A76_clean_R1.fastq.gz \
        -O A76_clean_R2.fastq.gz \
        -h A76_clean_fastp.html \
        -j A76_clean_fastp.json


## PBS configure 
#PBS -N hic_test
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=30
#PBS -l mem=20G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/A3A-Hic
conda activate /data/R02/nieyg/software/hic-pro
/data/R02/nieyg/software/HiC-Pro_3.1.0/bin/HiC-Pro -i /md01/nieyg/project/A3A-Hic/00.CleanData/  -o hiC_out_A3A  -c config_A3A_testdata.txt 

#/public/home/yangjw28/software/HiC-Pro-master/bin/HiC-Pro -i /public/home/yangjw28/data/OSN/HIC-F1-19d/cleandata -o out -c config-hicpro.txt -s quality_checks  1> run.log 2>&1
#/public/home/yangjw28/software/HiC-Pro-master/bin/HiC-Pro -i /public/home/yangjw28/data/OSN/HIC-F1-19d/cleandata -o out -c config-hicpro.txt -s build_contact_maps 1> run.log 2>&1
#/public/home/yangjw28/software/HiC-Pro-master/bin/HiC-Pro -i /public/home/yangjw28/data/OSN/HIC-F1-19d/cleandata -o out -c config-hicpro.txt -s ice_norm  1> run.log 2>&1

perl hic_summary.pl hiC_out_A3A/ >  HiC-Pro-split.Efficiency.xls
