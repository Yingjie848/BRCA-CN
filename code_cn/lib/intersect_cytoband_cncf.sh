# Intersect the segments covering 50% of the gene
# -- keep any intersected region
bedtools intersect -a analysis/cellular_fraction_cn/input/all_sample_cncf.bed -b input/cytoBand.bed -wa -wb -F 0.5 > analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband_50pct_ovlp.txt
