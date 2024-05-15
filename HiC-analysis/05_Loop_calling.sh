
# PART1: Loop calling
# We can use fanc loops to call loops in Hi-C matrices using the HICCUPS algorithm
## PBS configure
#PBS -N Loop_calling
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=100G
source /public/home/nieyg/.bash_profile
conda activate fanc

cd /md01/nieyg/project/A3A-Hic/05_Loop_calling
# A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/A3A.hic
# A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/A3A3.hic
# KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/KO.hic
# WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/WT.hic
A3A_100kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_100kb.hic
A3A3_100kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_100kb.hic
KO_100kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_100kb.hic
WT_100kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_100kb.hic

# res:100kb 
# 1.Annotating pixels for loop calling
fanc loops $A3A_100kb \
           A3A_100kb.loops \
           -t 10

# 2.Filtering annotated pixels
fanc loops A3A_100kb.loops \
           A3A_100kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops A3A_100kb_filtered.loops \
           A3A_100kb_merged.loops \
           -j --remove-singlets

# 4.Exporting to BEDPE
fanc loops A3A_100kb_merged.loops \
           -b A3A_100kb_merged.loops.bedpe


fanc loops $A3A3_100kb \
           A3A3_100kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops A3A3_100kb.loops \
           A3A3_100kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops A3A3_100kb_filtered.loops \
           A3A3_100kb_merged.loops \
           -j --remove-singlets
# 4.Exporting to BEDPE
fanc loops A3A3_100kb_merged.loops \
           -b A3A3_100kb_merged.loops.bedpe


fanc loops $WT_100kb \
           WT_100kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops WT_100kb.loops \
           WT_100kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops WT_100kb_filtered.loops \
           WT_100kb_merged.loops \
           -j --remove-singlets

# 4.Exporting to BEDPE
fanc loops WT_100kb_merged.loops \
           -b WT_100kb_merged.loops.bedpe



fanc loops $KO_100kb \
           KO_100kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops KO_100kb.loops \
           KO_100kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops KO_100kb_filtered.loops \
           KO_100kb_merged.loops \
           -j --remove-singlets
# 4.Exporting to BEDPE
fanc loops KO_100kb_merged.loops \
           -b KO_100kb_merged.loops.bedpe

# 5kb 

## PBS configure
#PBS -N Loop_calling
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10
#PBS -l mem=100G
source /public/home/nieyg/.bash_profile
conda activate fanc
cd /md01/nieyg/project/A3A-Hic/05_Loop_calling
A3A_5kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_5kb.hic
A3A3_5kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_5kb.hic
KO_5kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_5kb.hic
WT_5kb=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_5kb.hic
# res:5kb 
# 1.Annotating pixels for loop calling
fanc loops $A3A_5kb \
           A3A_5kb.loops \
           -t 10

# 2.Filtering annotated pixels
fanc loops A3A_5kb.loops \
           A3A_5kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops A3A_5kb_filtered.loops \
           A3A_5kb_merged.loops \
           -j --remove-singlets

# 4.Exporting to BEDPE
fanc loops A3A_5kb_merged.loops \
           -b A3A_5kb_merged.loops.bedpe


fanc loops $A3A3_5kb \
           A3A3_5kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops A3A3_5kb.loops \
           A3A3_5kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops A3A3_5kb_filtered.loops \
           A3A3_5kb_merged.loops \
           -j --remove-singlets
# 4.Exporting to BEDPE
fanc loops A3A3_5kb_merged.loops \
           -b A3A3_5kb_merged.loops.bedpe


fanc loops $WT_5kb \
           WT_5kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops WT_5kb.loops \
           WT_5kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops WT_5kb_filtered.loops \
           WT_5kb_merged.loops \
           -j --remove-singlets

# 4.Exporting to BEDPE
fanc loops WT_5kb_merged.loops \
           -b WT_5kb_merged.loops.bedpe



fanc loops $KO_5kb \
           KO_5kb.loops \
           -t 5

# 2.Filtering annotated pixels
fanc loops KO_5kb.loops \
           KO_5kb_filtered.loops \
           --rh-filter -d 5 -o 5

# 3.Merging unfiltered pixels into loops
fanc loops KO_5kb_filtered.loops \
           KO_5kb_merged.loops \
           -j --remove-singlets
# 4.Exporting to BEDPE
fanc loops KO_5kb_merged.loops \
           -b KO_5kb_merged.loops.bedpe
