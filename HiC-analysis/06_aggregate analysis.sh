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
conda activate r4-base



A3A_mediated_loop=

scp -r /public/home/chenxy/ChIP-seq_data/Mr.zhao/A3A/SE/WTac-1/WTac-1_peaks_Gateway_SuperEnhancers.bed nieyg@202.116.90.56:/md01/nieyg/project/A3A-Hic/05_Loop_calling/
SE_P_loop=

# get the heart_development gene
heart_development_region=

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


fanc aggregate /md01/nieyg/project/A3A-Hic/05_Loop_calling/KO_10kb_merged.loops \
               WT_10kb_merged.loops.bedpe \
               KO_10kb_merged_no_singlets.agg \
               -p KO_10kb_merged_no_singlets.agg.pdf \
               --loops -e -l -r 1.0 --rescale

fanc aggregate /md01/nieyg/project/A3A-Hic/05_Loop_calling/WT_10kb_merged.loops \
               WT_10kb_merged.loops.bedpe \
               WT_10kb_merged_no_singlets.agg \
               -p WT_10kb_merged_no_singlets.agg.pdf \
               --loops -e -l -r 1.0