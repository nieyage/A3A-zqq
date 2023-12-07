# Matrix analysis: chromatin compartments FAN-C includes implementations of the most established analyses and measures for the characterisation of Hi-C matrix properties. 
# Contact strength and the preference of contacts between certain genomic regions are particularly useful measures for gaining a global view of chromatin organisation. 

# 1. Expected values

# Contact distance decay plots: the average contact strength between loci separated by a certain distance, also called “expected contacts”, is typically shown in a log-log plot of expected contacts vs distance (Fig. 3c). 
# The slope and shape of the curve can inform about compaction of chromatin at various distance scales
# Hi-C矩阵中两个基因座之间的接触强度随着它们之间的距离增加而逐渐减弱。对于Hi-C矩阵，预期值随距离的变化呈现出独特的轮廓，可以近似为一个幂律，并在对数-对数图中形成几乎直线。

usage: fanc expected [-h] [-p PLOT_FILE] [-l LABELS [LABELS ...]]
                     [-c CHROMOSOME] [-tmp] [--recalculate] [-N]
                     input [input ...] output

## PBS configure
#PBS -N fanc_expected_10kb_chr1
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=5
#PBS -l mem=60G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/A3A-Hic/03_chromatin_compartments
conda activate fanc
output=/md01/nieyg/project/A3A-Hic/03_chromatin_compartments/01_expected
fanc expected -p $output/A3A_10kb_expected_values.pdf -l A3A  -c chr1 /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_10kb.hic $output/A3A_10kb_expected_values.txt
fanc expected -p $output/A3A3_10kb_expected_values.pdf -l A3A3  -c chr1 /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_10kb.hic $output/A3A3_10kb_expected_values.txt
fanc expected -p $output/KO_10kb_expected_values.pdf -l KO  -c chr1 /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_10kb.hic $output/KO_10kb_expected_values.txt
fanc expected -p $output/WT_10kb_expected_values.pdf -l WT  -c chr1 /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_10kb.hic $output/WT_10kb_expected_values.txt

echo "All samples processed"

# 2. Comparing expected values

samples=("A3A" "A3A-2" "A3A3-1" "A3A3-2" "KO-1" "KO-2" "WT-1" "WT-2")
dir=/md01/nieyg/project/A3A-Hic/02.fanc_out
fanc expected -p $output/A3A_10kb_expected_values_compare.pdf \
-l "A3A" "A3A3"  "KO" "WT"  -c chr1\
$dir/A3A/hic/binned/A3A_10kb.hic \
$dir/A3A3/hic/binned/A3A3_10kb.hic \
$dir/KO/hic/binned/KO_10kb.hic \
$dir/WT/hic/binned/WT_10kb.hic \
$output/A3A_10kb_expected_values_compare.txt

# 3. O/E matrices
# Observed/expected (O/E) transformation: a central transformation used by many analyses in which each pixel represents the (log2-)fold-change enrichment over the expected contact intensity for a region at that distance (Fig. 3d). 
# Expected values are stored by FAN-C inside each matrix, allowing a fast, dynamic conversion of normalised into O/E contacts for various applications;

fancplot -o architecture/expected/fanc_example_500kb_chr18_oe.png \
     chr18:1-78mb -p triangular -e output/hic/binned/fanc_example_500kb.hic \
     -vmin -2 -vmax 2



















