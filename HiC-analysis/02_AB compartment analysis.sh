# AB compartment analysis
# Regions in a Hi-C matrix can generally be assigned to either the active or the inactive compartment, also called ‘A’ and ‘B’ compartments, respectively.


usage: fanc compartments [-h] [-d DOMAINS] [-v EIGENVECTOR]
                         [-e ENRICHMENT_FILE] [-m MATRIX_FILE] [-g GENOME]
                         [-w] [-r REGION] [-i EIGENVECTOR_INDEX]
                         [-p PERCENTILES [PERCENTILES ...]] [-c COLORMAP]
                         [-s SYMMETRIC_AT] [--enrichment-min VMIN]
                         [--enrichment-max VMAX] [-G]
                         [-x EXCLUDE [EXCLUDE ...]]
                         [--compartment-strength COMPARTMENT_STRENGTH_FILE]
                         [-tmp] [-f] [--recalculate]
                         matrix [ab_compartments]

# 1. Correlation matrix: Compartments are derived from a correlation matrix, in which each entry i, j corresponds to the Pearson correlation between row i and column j of the (Hi-C) matrix.
# include AB Eigenvector,AB domains :

fanc -v compartments output/hic/binned/fanc_example_1mb.hic architecture/compartments/fanc_example_1mb.ab

fancplot -o architecture/compartments/fanc_example_1mb.ab_and_ev.png chr18 \
     -p square architecture/compartments/fanc_example_1mb.ab \
     -vmin -0.75 -vmax 0.75 -c RdBu_r \
     -p line architecture/compartments/fanc_example_1mb.ev.txt


# 2. AB domains




