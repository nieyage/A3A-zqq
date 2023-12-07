# Saddle plots showing contacts relative to expectation among regions with different compartment eigenvector values (binned by 2% percentiles).

#
# Matrix plots
#
## PBS configure
#PBS -N fanc_Matrixplots
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=5
#PBS -l mem=60G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/A3A-Hic/04_comparison_workflow
conda activate fanc
#samples=("A3A" "A3A-2" "A3A3-1" "A3A3" "KO-1" "KO" "WT-1" "WT")
output=/md01/nieyg/project/A3A-Hic/04_comparison_workflow/01_Matrixplots
mkdir $output
fancplot -o $output/A3A_1mb.chr1.pdf chr1:1-195mb -p square /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_1mb.hic  -vmax 0.01 -vmin 0 -f
fancplot -o $output/A3A3_1mb.chr1.pdf chr1:1-195mb -p square /md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_1mb.hic  -vmax 0.01 -vmin 0 -f
fancplot -o $output/KO_1mb.chr1.pdf chr1:1-195mb -p square /md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_1mb.hic  -vmax 0.01 -vmin 0 -f
fancplot -o $output/WT_1mb.chr1.pdf chr1:1-195mb -p square /md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_1mb.hic  -vmax 0.01 -vmin 0 -f
echo "All samples processed"

#
# Compartments
#

output=/md01/nieyg/project/A3A-Hic/04_comparison_workflow/02_Compartments
#mkdir $output
GENOME=/md01/nieyg/ref/genome/hg19/hg19.fa
A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/binned/A3A_1mb.hic
A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/binned/A3A3_1mb.hic
KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/binned/KO_1mb.hic
WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/binned/WT_1mb.hic

fanc compartments $A3A A3A_1mb.ab -e A3A_1mb.ab.pdf -g $GENOME -tmp -p 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 --enrichment-colormap RdBu_r
fanc compartments $A3A3 A3A3_1mb.ab -e A3A3_1mb.ab.pdf -g $GENOME -tmp -p 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 --enrichment-colormap RdBu_r
fanc compartments $KO KO_1mb.ab -e KO_1mb.ab.pdf -g $GENOME -tmp -p 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 --enrichment-colormap RdBu_r
fanc compartments $WT WT_1mb.ab -e WT_1mb.ab.pdf -g $GENOME -tmp -p 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 --enrichment-colormap RdBu_r


#
# Compartment strength
#
output=/md01/nieyg/project/A3A-Hic/04_comparison_workflow/03_Compartments_strength
#mkdir $output
A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/A3A.hic
A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/A3A3.hic
KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/KO.hic
WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/WT.hic

nohup fanc compartments $A3A $output/A3A_all_cs.ab --compartment-strength $output/A3A_all.compartment_strength.txt -g $GENOME -tmp &
nohup fanc compartments $A3A3 $output/A3A3_all_cs.ab --compartment-strength $output/A3A3_all.compartment_strength.txt -g $GENOME -tmp &
nohup fanc compartments $KO  $output/KO_all_cs.ab --compartment-strength $output/KO_all.compartment_strength.txt -g $GENOME -tmp &
nohup fanc compartments $WT  $output/WT_all_cs.ab --compartment-strength $output/WT_all.compartment_strength.txt -g $GENOME -tmp &

#
# Insulation scores
#
A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/A3A.hic
A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/A3A3.hic
KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/KO.hic
WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/WT.hic

nohup fanc insulation  -g $A3A  A3A_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation  -g $A3A3 A3A3_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation  -g $KO   KO_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation  -g $WT   WT_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb &

# --impute？

nohup fanc insulation --impute -g $A3A  A3A_all_imputed.ins -tmp -w  100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation --impute -g $A3A3 A3A3_all_imputed.ins -tmp -w  100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation --impute -g $KO   KO_all_imputed.ins -tmp -w  100kb 250kb 500kb 1mb 2mb 5mb &
nohup fanc insulation --impute -g $WT   WT_all_imputed.ins -tmp -w  100kb 250kb 500kb 1mb 2mb 5mb &

#
# Boundaries
#

fanc boundaries A3A_all.ins A3A.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb 
fanc boundaries A3A3_all.ins A3A3.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb 
fanc boundaries WT_all.ins WT.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb 
fanc boundaries KO_all.ins KO.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb 

# Let’s plot the boundaries from the 1mb scores:
nohup fanc insulation  -g $A3A  A3A_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb -o bed &
nohup fanc insulation  -g $A3A3 A3A3_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb -o bed &
nohup fanc insulation  -g $KO   KO_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb -o bed &
nohup fanc insulation  -g $WT   WT_all.ins -tmp -w 100kb 250kb 500kb 1mb 2mb 5mb -o bed &

fancplot --width 6 -o architecture/domains/fanc_example_100kb_tads_insulation_1mb_boundaries.png \
     chr18:18mb-28mb \
     -p triangular output/hic/binned/fanc_example_100kb.hic -m 4000000 -vmin 0 -vmax 0.05 \
     -p line architecture/domains/fanc_example_100kb.insulation_1mb.bed -l "1mb" \
     -p bar architecture/domains/fanc_example_100kb.insulation_boundaries_1mb.bed





#
# Insulation score differences
#
# A3A3 vs A3A
fanc compare -tmp -c difference A3A3_all.ins A3A_all.ins A3A3_minus_A3A_all.ins
# A3A vs WT
fanc compare -tmp -c difference A3A_all.ins WT_all.ins A3A_minus_WT_all.ins
# A3A3 vs WT
fanc compare -tmp -c difference A3A3_all.ins WT_all.ins A3A3_minus_WT_all.ins

# A3A vs KO
fanc compare -tmp -c difference A3A_all.ins KO_all.ins A3A_minus_KO_all.ins
# A3A3 vs KO
fanc compare -tmp -c difference A3A3_all.ins KO_all.ins A3A3_minus_KO_all.ins
# KO vs WT
fanc compare -tmp -c difference KO_all.ins WT_all.ins KO_minus_WT_all.ins


# Boundaries for different region 
fanc boundaries KO_minus_WT_all.ins KO_minus_WT_all.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb
fanc boundaries A3A_minus_WT_all.ins A3A_minus_WT_all.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb
fanc boundaries A3A3_minus_WT_all.ins A3A3_minus_WT_all.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb
fanc boundaries A3A_minus_KO_all.ins A3A_minus_KO_all.ins.boundaries.bed  -w 100kb,250kb,500kb,1mb,2mb,5mb
fanc boundaries A3A3_minus_KO_all.ins A3A3_minus_KO_all.ins.boundaries.bed  -w 100kb,250kb,500kb,1mb,2mb,5mb
fanc boundaries A3A3_minus_A3A_all.ins A3A3_minus_A3A_all.ins.boundaries.bed -w 100kb,250kb,500kb,1mb,2mb,5mb


# show some region TAD 
fancplot -o architecture/domains/fanc_example_100kb_tads.png chr18:18mb-28mb \
     -p triangular output/hic/binned/fanc_example_100kb.hic -m 4000000 \
     -vmin 0 -vmax 0.05


# aggregate analysis

#
# Sample boundary plots
#
A3A=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A/hic/A3A.hic
A3A3=/md01/nieyg/project/A3A-Hic/02.fanc_out/A3A3/hic/A3A3.hic
KO=/md01/nieyg/project/A3A-Hic/02.fanc_out/KO/hic/KO.hic
WT=/md01/nieyg/project/A3A-Hic/02.fanc_out/WT/hic/WT.hic
GENES=/md01/nieyg/ref/10X/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz
conda activate fanc

fancplot --width 4 -w 100000 \
  -o MEIS1.pdf \
  chr11:18786657-19113663 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A3 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $KO \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $WT \
  -p gene $GENES --aspect-ratio 0.3 --group-by "gene_id" -sq --label-field gene_name -a 0

fancplot --width 4 -w 100000 \
  -o ACTN2.pdf \
  chr13:12224440-123857183 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A3 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $KO \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $WT \
  -p gene $GENES --aspect-ratio 0.3 --group-by "gene_id" -sq --label-field gene_name -a 0


fancplot --width 4 -w 100000 \
  -o ACTN2.pdf \
  chr18:16456109-16942015 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $A3A3 \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $KO \
  -p triangular -vmin 0 -vmax 0.02 -m 500kb $WT \
  -p gene $GENES --aspect-ratio 0.3 --group-by "gene_id" -sq --label-field gene_name -a 0












