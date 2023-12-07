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
dir=/md01/nieyg/project/A3A-Hic/02.fanc_out
output=/md01/nieyg/project/A3A-Hic/03_chromatin_compartments/02_PCA
mkdir $output

fanc pca -n "A3A" "A3A3"  "KO" "WT" \
         -Z -s 100000 -p $output/A3A_10kb_pca.pdf \
$dir/A3A/hic/binned/A3A_10kb.hic \
$dir/A3A3/hic/binned/A3A3_10kb.hic \
$dir/KO/hic/binned/KO_10kb.hic \
$dir/WT/hic/binned/WT_10kb.hic \
$output/A3A_10kb_pca_lowc.pca

