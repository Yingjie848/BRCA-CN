# 2021-01-13 This script uses different copy number loss definition:
#                 as:   <= ploidy -1
#
# Find enriched cn change cytobands for biallelic samples
# -- this script starts as to check cellular fraction
# -- then later evolve into testing the cytobands signifcance 2020-06
rm(list = ls())
library(data.table)
library(ipfun)
library(ggpubr)
library(ggplot2)
source("lib/lib_get_interested_biallelic_sample.R")

GAIN_CUTOFF <- 1 # Recent DDR group discussion - thinking of trying ploidy+1
# Final biallelic sample category classification
# sample_cat <- get_biallelic_samples_control_w_VUS_noHiLST_domSig3()

## Further tweaking different positive and negative control
# - trying with BRCA1 vs control with VUS
# sample_cat <- get_samples_basedon_gene_and_function(
#     c("BRCA1"),
#     get_biallelic_samples_control_w_VUS_noHiLST_domSig3
# )
sample_cat <- get_samples_basedon_gene_and_function(
    c("BRCA2", "PALB2"),
    get_biallelic_samples_control_w_VUS_noHiLST_domSig3
)

## start data loading and processing
load_processed_data <- function() {
    # Loading copy number and cytoband info
    d <- fread("analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband.csv")

    # Loading master info (sample categories)
    m <- fread("input/master_score_MC3_PanImmun.csv")
    m$cancer <- m$Cancer
    m$cat <- "all"
    m[Id %in% sample_cat$biallelic, cat := "biallelicBRCA12"]
    m[Id %in% sample_cat$control, cat := "control"]

    # m <- fread("analysis/explore_monoallele_aneuploidy_HR/input/master_novus_brca12_cn_category.csv")
    m_merge_columns <- m[, .(Id, cat, cancer)]
    dm <- merge(d, m_merge_columns, by.x = "sample", by.y = "Id")
    dm$gen_loc <- paste(dm$chromosome, dm$cytoband, sep = "_")

    # Simplify comparison between biallelic samples and controls only
    # spd <- dm[cat == "biallelicBRCA12" | cat == "control"]
    spd <- dm
    return(spd)
}

spd <- load_processed_data()

clonal <- spd
# clonal$ploidy <- round(clonal$ploidy)
clonal$ploidy <- ifelse(clonal$ploidy < 2, round(clonal$ploidy), clonal$ploidy)


## selecting OV clonal bands
ov_data <- clonal[cancer == "OV"]
brca_data <- clonal[cancer == "BRCA"]
erpos_data <- clonal[sample %in% sample_cat$er_positive]
tri_data <- clonal[sample %in% sample_cat$triple_negative]

save_tri_data <- function() {
    # Save the tcga triple negative data for combined analysis with Serena's dataset
    tri_data$group <- "TCGA_triple_negative"
    fwrite(tri_data, "output/copy_number_cytoband_group/triple_negative_cytoband_cn_TCGA.csv")
}

tcn_compare_cytoband <- function(inputband) {
    # function to get summary about cytoband, cloncal, subclone diff and plot
    s_clonal <- clonal[gen_loc == inputband]
    s_subclonal <- subclonal[gen_loc == inputband]
    s_all <- spd[gen_loc == inputband]

    writeLines(paste("There are", nrow(s_all), "all cytobands."))
    writeLines(paste("- There are", nrow(s_clonal), "clonal cytobands."))
    writeLines(paste("- There are", nrow(s_subclonal), "subclonal cytobands."))

    flist <- list()
    flist$clonal <- plot_tcn_cat_boxplot(clonal, inputband)
    flist$subclonal <- plot_tcn_cat_boxplot(subclonal, inputband)
    flist$all <- plot_tcn_cat_boxplot(spd, inputband)
    return(flist)
}

# Exploring cncf segments number count vs cytoband
# -- There are about 860 cytoband
# cncf=fread("input/final_all_tumor_cncf_combined.csv")
# -- in comparision, segments are much fewer, each sample median 45 segments, so only about 2 segments per chromosome

test_clonal_band <- function(inputband, testdata = clonal) {
    # main test to return the stats using testdata on inputband(cytoband)
    # generating the stats for a single inputband
    writeLines(paste("Testing cytoband:", inputband))
    band_clonal <- testdata[gen_loc == inputband]

    # Look at gain
    bi_gain <- nrow(band_clonal[cat == "biallelicBRCA12" &
        tcn.em >= (ploidy + GAIN_CUTOFF)])
    bi_nogain <- nrow(band_clonal[cat == "biallelicBRCA12" &
        tcn.em < (ploidy + GAIN_CUTOFF)])

    control_gain <- nrow(band_clonal[cat == "control" &
        tcn.em >= (ploidy + GAIN_CUTOFF)])
    control_nogain <- nrow(band_clonal[cat == "control" &
        tcn.em < (ploidy + GAIN_CUTOFF)])

    others_gain <- nrow(band_clonal[cat == "all" &
        tcn.em >= (ploidy + GAIN_CUTOFF)])
    others_nogain <- nrow(band_clonal[cat == "all" &
        tcn.em < (ploidy + GAIN_CUTOFF)])

    # Look at loss
    bi_loss <- nrow(band_clonal[cat == "biallelicBRCA12" &
        tcn.em <= (ploidy - 1)])
    bi_noloss <- nrow(band_clonal[cat == "biallelicBRCA12" & tcn.em > (ploidy - 1)])
    control_loss <- nrow(band_clonal[cat == "control" & tcn.em <= (ploidy - 1)])
    control_noloss <- nrow(band_clonal[cat == "control" & tcn.em > (ploidy - 1)])
    others_loss <- nrow(band_clonal[cat == "all" & tcn.em <= (ploidy - 1)])
    others_noloss <- nrow(band_clonal[cat == "all" & tcn.em > (ploidy - 1)])

    # Look at loh
    bi_loh <- nrow(band_clonal[cat == "biallelicBRCA12" & lcn.em == 0])
    bi_noloh <- nrow(band_clonal[cat == "biallelicBRCA12" & lcn.em > 0])
    control_loh <- nrow(band_clonal[cat == "control" & lcn.em == 0])
    control_noloh <- nrow(band_clonal[cat == "control" & lcn.em > 0])

    # Look at homdel
    bi_homdel <- nrow(band_clonal[cat == "biallelicBRCA12" & tcn.em == 0])
    bi_nohomdel <- nrow(band_clonal[cat == "biallelicBRCA12" & tcn.em > 0])
    control_homdel <- nrow(band_clonal[cat == "control" & tcn.em == 0])
    control_nohomdel <- nrow(band_clonal[cat == "control" & tcn.em > 0])

    # calc p values
    pval_gain <- fisher.test(matrix(c(
        bi_gain, bi_nogain,
        control_gain, control_nogain
    ),
    nrow = 2
    ))$p.value
    pval_loss <- fisher.test(matrix(c(
        bi_loss, bi_noloss,
        control_loss, control_noloss
    ),
    nrow = 2
    ))$p.value

    pval_loh <- fisher.test(matrix(c(
        bi_loh, bi_noloh,
        control_loh, control_noloh
    ),
    nrow = 2
    ))$p.value

    pval_homdel <- fisher.test(matrix(c(
        bi_homdel, bi_nohomdel,
        control_homdel, control_nohomdel
    ),
    nrow = 2
    ))$p.value

    row <- data.table(
        band = inputband,
        p_gain = pval_gain, p_loss = pval_loss, p_loh = pval_loh, p_homdel = pval_homdel,
        n_biallelic_gain = bi_gain, n_biallelic_nogain = bi_nogain, n_control_gain = control_gain, n_control_nogain = control_nogain, n_biallelic_loss = bi_loss, n_biallelic_noloss = bi_noloss, n_control_loss = control_loss, n_control_noloss = control_noloss,
        n_biallelic_loh = bi_loh, n_biallelic_noloh = bi_noloh, n_control_loh = control_loh, n_control_noloh = control_noloh,
        n_biallelic_homdel = bi_homdel, n_biallelic_nohomdel = bi_nohomdel, n_control_homdel = control_homdel, n_control_nohomdel = control_nohomdel,
        n_others_gain = others_gain, n_others_nogain = others_nogain,
        n_others_loss = others_loss, n_others_noloss = others_noloss
    )

    return(row)
}

test_sample_cytobands <- function(inputdata, resultfile) {
    # Return a object and save results in the resultfile
    # -- result calculated all the pval and number of cases for all cytobands for clonal copy number changes
    result_file_path <- paste0(TEST_OUTPUT_FOLDER, resultfile, ".csv")
    all_band <- unique(inputdata$gen_loc)

    # Define a function to use the inputdata to test the band
    test_specific_sample_clonal_band <- function(x) test_clonal_band(x, testdata = inputdata)

    all_result <- rbindlist(lapply(all_band, test_specific_sample_clonal_band))
    fwrite(all_result, result_file_path)
    return(all_result)
}

save_fdr_adjustment <- function(inputfile, outputfile) {
    # Save the FDR ajudst p values for the cloncal cytoband cn change tests
    # ds <- fread("analysis/cellular_fraction_cn/output/all_clonal_cn_change_test_stats.csv")
    input_file_path <- paste0(TEST_OUTPUT_FOLDER, inputfile)
    output_file_path <- paste0(TEST_OUTPUT_FOLDER, outputfile)
    ds <- fread(input_file_path)

    # FDR adjustment pval
    ds$p_gain_FDR <- p.adjust(ds$p_gain, , method = "fdr")
    ds$p_loss_FDR <- p.adjust(ds$p_loss, , method = "fdr")
    ds$p_loh_FDR <- p.adjust(ds$p_loh, , method = "fdr")
    ds$p_homdel_FDR <- p.adjust(ds$p_homdel, , method = "fdr")

    # Calculate percentage among biallelic samples vs. control samples
    ds$biallelic_gain_pct <- with(ds, n_biallelic_gain / (n_biallelic_gain + n_biallelic_nogain))
    ds$control_gain_pct <- with(ds, n_control_gain / (n_control_gain + n_control_nogain))
    ds$biallelic_loss_pct <- with(ds, n_biallelic_loss / (n_biallelic_loss + n_biallelic_noloss))
    ds$control_loss_pct <- with(ds, n_control_loss / (n_control_loss + n_control_noloss))
    ds$biallelic_loh_pct <- with(ds, n_biallelic_loh / (n_biallelic_loh + n_biallelic_noloh))
    ds$control_loh_pct <- with(ds, n_control_loh / (n_control_loh + n_control_noloh))
    ds$biallelic_homdel_pct <- with(ds, n_biallelic_homdel / (n_biallelic_homdel + n_biallelic_nohomdel))
    ds$control_homdel_pct <- with(ds, n_control_homdel / (n_control_homdel + n_control_nohomdel))

    # Add the rest of the samples (not biallelic nor control) pct of gain/loss
    ds$others_gain_pct <- with(ds, n_others_gain / (n_others_gain + n_others_nogain))
    ds$others_loss_pct <- with(ds, n_others_loss / (n_others_loss + n_others_noloss))

    # fwrite(ds, "analysis/cellular_fraction_cn/output/all_clonal_cn_change_test_stats_FDR.csv")
    fwrite(ds, output_file_path)
    return(ds)
}

save_result_tumortypes <- function() {
    ## Testing specific samples (OV/BRCA ER+) samples
    test_sample_cytobands(clonal, "all_clonal_cn_test_stats")
    save_fdr_adjustment("all_clonal_cn_test_stats.csv", "FDR_all_clonal_cn_test_stats.csv")

    test_sample_cytobands(ov_data, "ov_clonal_cn_test_stats")
    save_fdr_adjustment("ov_clonal_cn_test_stats.csv", "FDR_ov_clonal_cn_test_stats.csv")

    test_sample_cytobands(brca_data, "brca_clonal_cn_test_stats")
    save_fdr_adjustment("brca_clonal_cn_test_stats.csv", "FDR_brca_clonal_cn_test_stats.csv")

    test_sample_cytobands(erpos_data, "erpos_clonal_cn_test_stats")
    save_fdr_adjustment("erpos_clonal_cn_test_stats.csv", "FDR_erpos_clonal_cn_test_stats.csv")

    test_sample_cytobands(tri_data, "tri_clonal_cn_test_stats")
    save_fdr_adjustment("tri_clonal_cn_test_stats.csv", "FDR_tri_clonal_cn_test_stats.csv")
}

load_and_process_fdr_result <- function(result_folder) {
    # Loading tht FDR result from given result_folder
    # To return a list with objects and overlapping cytobands
    folder <- paste0(TEST_OUTPUT_FOLDER, result_folder)
    all_path <- paste0(folder, "/FDR_all_clonal_cn_test_stats.csv")
    ov_path <- paste0(folder, "/FDR_ov_clonal_cn_test_stats.csv")
    brca_path <- paste0(folder, "/FDR_brca_clonal_cn_test_stats.csv")
    er_path <- paste0(folder, "/FDR_erpos_clonal_cn_test_stats.csv")
    tri_path <- paste0(folder, "/FDR_tri_clonal_cn_test_stats.csv")
    fr <- list()
    fr$fdr <- fread(all_path)
    fr$fdr_ov <- fread(ov_path)
    fr$fdr_brca <- fread(brca_path)
    fr$fdr_er <- fread(er_path)
    fr$fdr_tri <- fread(tri_path)
    fr$fdr_ov <- fr$fdr_ov[order(band), ]
    fr$fdr_brca <- fr$fdr_brca[order(band), ]
    fr$fdr_er <- fr$fdr_er[order(band), ]
    fr$fdr_tri <- fr$fdr_tri[order(band), ]

    # Check the overlapping loss armband between brca and ov
    fr$loss <- list()
    fr$loss$all <- fr$fdr[(p_loss_FDR < 0.05) &
        (biallelic_loss_pct > control_loss_pct)]
    fr$loss$brca <- fr$fdr_brca[(p_loss_FDR < 0.05) &
        (biallelic_loss_pct > control_loss_pct)]
    fr$loss$ov <- fr$fdr_ov[(p_loss_FDR < 0.05) &
        (biallelic_loss_pct > control_loss_pct)]
    fr$loss$er <- fr$fdr_er[(p_loss_FDR < 0.05) &
        (biallelic_loss_pct > control_loss_pct)]
    fr$loss$tri <- fr$fdr_tri[(p_loss_FDR < 0.05) &
        (biallelic_loss_pct > control_loss_pct)]

    # Finding overlapping cytoband
    fr$loss$bo <- fr$loss$brca[band %in% fr$loss$ov$band]
    fr$loss$ov_er <- fr$loss$ov[band %in% fr$loss$er$band]
    # fr$loss$bo_in_yj <- fr$loss$bo[band %in% yj_loss$gen_loc]

    # Check fdr gain region
    fr$gain <- list()
    fr$gain$all <- fr$fdr[(p_gain_FDR < 0.15) &
        (biallelic_gain_pct > control_gain_pct)]
    fr$gain$brca <- fr$fdr_brca[(p_gain_FDR < 0.15) &
        (biallelic_gain_pct > control_gain_pct)]
    fr$gain$ov <- fr$fdr_ov[(p_gain_FDR < 0.15) &
        (biallelic_gain_pct > control_gain_pct)]
    fr$gain$er <- fr$fdr_er[(p_gain_FDR < 0.15) &
        (biallelic_gain_pct > control_gain_pct)]
    fr$gain$tri <- fr$fdr_tri[(p_gain_FDR < 0.15) &
        (biallelic_gain_pct > control_gain_pct)]

    fr$gain$bo <- fr$gain$brca[band %in% fr$gain$ov$band]
    fr$gain$ov_er <- fr$gain$ov[band %in% fr$gain$er$band]
    # fr$gain$bo_in_yj <- fr$gain$bo[band %in% yj_gain$gen_loc]


    fr$combined_gain_pval <- data.table(
        band = fr$fdr_ov$band,
        # band = gsub("^X_", "23_", fr$fdr_ov$band),
        # -- needed to be consistant with other data
        OV = fr$fdr_ov$p_gain_FDR,
        BRCA = fr$fdr_brca$p_gain_FDR,
        ER = fr$fdr_er$p_gain_FDR,
        TRI = fr$fdr_tri$p_gain_FDR
    )

    fr$combined_gain_pct_change <- data.table(
        band = fr$fdr_ov$band,
        # band = gsub("^X_", "23_", fr$fdr_ov$band),
        OV = fr$fdr_ov$biallelic_gain_pct - fr$fdr_ov$control_gain_pct,
        BRCA = fr$fdr_brca$biallelic_gain_pct - fr$fdr_brca$control_gain_pct,
        ER = fr$fdr_er$biallelic_gain_pct - fr$fdr_er$control_gain_pct,
        TRI = fr$fdr_tri$biallelic_gain_pct - fr$fdr_tri$control_gain_pct
    )

    fr$combined_loss_pval <- data.table(
        band = fr$fdr_ov$band,
        # band = gsub("^X_", "23_", fr$fdr_ov$band),
        OV = fr$fdr_ov$p_loss_FDR,
        BRCA = fr$fdr_brca$p_loss_FDR,
        ER = fr$fdr_er$p_loss_FDR,
        TRI = fr$fdr_tri$p_loss_FDR
    )

    fr$combined_loss_pct_change <- data.table(
        band = fr$fdr_ov$band,
        # band = gsub("^X_", "23_", fr$fdr_ov$band),
        OV = fr$fdr_ov$biallelic_loss_pct - fr$fdr_ov$control_loss_pct,
        BRCA = fr$fdr_brca$biallelic_loss_pct - fr$fdr_brca$control_loss_pct,
        ER = fr$fdr_er$biallelic_loss_pct - fr$fdr_er$control_loss_pct,
        TRI = fr$fdr_tri$biallelic_loss_pct - fr$fdr_tri$control_loss_pct
    )

    return(fr)
}

sho_selected_columns <- function(x) {
    # interactive inspecting function to show the selected columns in x
    a <- x[, .(band, p_loss_FDR, biallelic_loss_pct, control_loss_pct, n_biallelic_loss, n_control_loss)]
    a[order(p_loss_FDR), ]
}

save_processed_cytoband_result <- function() {
    # Using the loaded & processed cytoband results to save combined pval for further plotting
    fwrite(r0$combined_gain_pval, "analysis/cellular_fraction_cn/input/combined_cytobands_result/combined_gain_pval.csv")
    fwrite(r0$combined_gain_pct_change, "analysis/cellular_fraction_cn/input/combined_cytobands_result/combined_gain_pct_change.csv")
    fwrite(r0$combined_loss_pval, "analysis/cellular_fraction_cn/input/combined_cytobands_result/combined_loss_pval.csv")
    fwrite(r0$combined_loss_pct_change, "analysis/cellular_fraction_cn/input/combined_cytobands_result/combined_loss_pct_change.csv")
}

load_process_save_cytoband_result <- function(test_folder = "clonal_all_clones") {
    # Performing loading from the test_folder result
    # and then process, FDR adjust, save the result in the input/sig_cytobands_ ... test_folder data

    r0 <- load_and_process_fdr_result(test_folder)
    result_folder <- paste0("analysis/cellular_fraction_cn/input/sig_cytobands_", test_folder)
    dir.create(result_folder)
    result_f <- paste0(result_folder, "/cytobands_level_fdr_tests.rda")
    save(r0, file = result_f)
}

get_cohort_number <- function() {
    # Get the study cohort number for the manuscript
    d <- fread("analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband.csv")

    # Loading master info (sample categories)
    m <- fread("input/master_score_MC3_PanImmun.csv")
    m$cancer <- m$Cancer
    tcga.ov <- m[cancer == "OV"]
    fr <- list()
    fr$ov.number <- nrow(tcga.ov[Id %in% d$sample])
    fr$er.number <- length(sample_cat$er_positive[sample_cat$er_positive %in% d$sample])
    fr$tri.number <- length(sample_cat$triple_negative[sample_cat$triple_negative %in% d$sample])
    return(fr)
}

#####################################################
main <- function() {
    ## main section
    # -- there are mainly 2 logic steps:
    ## 1. testing cytoband level cn changes for all samples
    # save_result_tumortypes() # this saves to default location, later to organize

    # 2. further combinding and processings tested pvals
    # load_process_save_cytoband_result("cnloss_newRoundedPloidy1")
    # load_process_save_cytoband_result("cnloss_newRoundedPloidy1_BRCA1")
    load_process_save_cytoband_result("cnloss_newRoundedPloidy1_BRCA2")
}

# r <- main()
a <- get_cohort_number()
