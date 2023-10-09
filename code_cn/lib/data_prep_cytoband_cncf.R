# Prep bed files with all samples cncf and cytoband labels

rm(list = ls())
library(data.table)
library(ipfun)

save_all_sample_cncf_bed <- function() {
    # Save all the monoallelic sample's cncf file into bed format for bedtools
    cncf <- fread("input/final_all_tumor_cncf_combined.csv")
    codedcncf <- cncf
    codedcncf$chromosome <- ifelse(codedcncf$chrom == 23, "X",
        as.character(codedcncf$chrom)
    )
    m <- codedcncf[
        ,
        .(chromosome, start, end, sample, tcn.em, lcn.em, ploidy, cf.em)
    ]
    fwrite(m, "analysis/cellular_fraction_cn/input/all_sample_cncf.bed", sep = "\t", col.names = F)
}
save_cytoband_bed <- function() {
    cyto <- fread("input/cytoBand.txt")
    cyto <- cyto[V1 != "chrY"]
    cyto$V1 <- sub("chr", "", cyto$V1)
    fwrite(cyto, "input/cytoBand.bed", sep = "\t", col.names = F)
}

save_cncf_cytoband_dataset <- function() {
    d <- fread("analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband_50pct_ovlp.txt")
    names(d) <- c("chromosome", "start", "end", "sample", "tcn.em", "lcn.em", "ploidy", "cf.em", "chrom_txt", "cyto_start", "cyto_end", "cytoband", "gieStain")
    fwrite(d, "analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband.csv")
}
############################################################### 3
# main execution logic
# save_all_sample_cncf_bed()
# save_cytoband_bed()
# save_cncf_cytoband_dataset()