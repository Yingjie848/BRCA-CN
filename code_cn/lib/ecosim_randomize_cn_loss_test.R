# Randomize the cytoband using the same ratio/number of loss per sample to assesss the background level cn loss for biallelic/control samples
rm(list = ls())
library(data.table)
library(ipfun)
library(EcoSimR)
source("lib/fun_get_interested_biallelic_sample.R")

# Giving the runNumber argument to make runs
runNumber <- commandArgs(trailingOnly = TRUE)
sample_cat <- get_biallelic_samples_control_w_VUS_noHiLST_domSig3() # control remove high LST/domSig3

## start data loading and processing
load_processed_data <- function() {
    # Loading copy number and cytoband info
    d <- fread("analysis/cellular_fraction_cn/input/all_sample_cncf_cytoband.csv")
    d$ploidy <- ifelse(d$ploidy < 2, round(d$ploidy), d$ploidy)

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
    spd <- dm[cat == "biallelicBRCA12" | cat == "control"]
    return(spd)
}

## Previous code testing further of 'clonal' copy number changes based on CCF, here we just use the same variable/but include all CN
spd <- load_processed_data()

# -------------------------------------------------
## selecting OV clonal bands
clonal <- spd
ov_data <- clonal[cancer == "OV"]
brca_data <- clonal[cancer == "BRCA"]
erpos_data <- clonal[sample %in% sample_cat$er_positive]
tri_data <- clonal[sample %in% sample_cat$triple_negative]

# -------------------------------------------------
## Use ecosimR to simutation/randomize the data, need to sim each tumor type
start_burn_in_sim_n <- function(input, n = 50) { # n could be 50 as in the Nils' code
    # Convert the input data.table object to matrix and simulate n times
    # TODO: convert from tall to wide format for each cytoband, and 0, 1(for loss), sample by cytoband
    # EcoSimR doesn't allow input to have negative value or NA, so recode all NA as 0 here
    input$cnloss <- fifelse(input$tcn.em <= input$ploidy - 1, 1, 0, na = 0)
    scm <- dcast(input, sample ~ gen_loc, value.var = "cnloss", fill = 0)
    scm_m <- as.matrix(scm, rownames = "sample")
    ## This is the key step of generating the simutation matrix
    for (i in 1:n) {
        scm_m <- sim9_single(scm_m)
    }
    return(scm_m)
}

get_dt_from_sim <- function(scm_m) {
    # From the simuation scm_m to obtain the stat test object
    scm_dt <- as.data.table(scm_m, keep.rownames = T)
    scm_dt_m <- melt(scm_dt, id.vars = c("rn"))
    names(scm_dt_m) <- c("sample", "gen_loc", "cnloss")
    scm_dt_m$gen_loc <- as.character(scm_dt_m$gen_loc)
    scm_dt_m$cat <- "all"
    scm_dt_m[sample %in% sample_cat$biallelic, cat := "biallelicBRCA12"]
    scm_dt_m[sample %in% sample_cat$control, cat := "control"]
    return(scm_dt_m)
}

# Exploring cncf segments number count vs cytoband
# -- There are about 860 cytoband
# cncf=fread("input/final_all_tumor_cncf_combined.csv")
# -- in comparision, segments are much fewer, each sample median 45 segments, so only about 2 segments per chromosome

test_clonal_band <- function(inputband, testdata = clonal) {
    # main test to return the stats using testdata on inputband(cytoband)
    # generating the stats for a single inputband

    # writeLines(paste("Testing cytoband:", inputband))
    # -- disable the line above to make the system run faster

    band_clonal <- testdata[gen_loc == inputband]

    # Look at loss
    bi_loss <- nrow(band_clonal[cat == "biallelicBRCA12" & cnloss == 1])
    bi_noloss <- nrow(band_clonal[cat == "biallelicBRCA12" & cnloss == 0])
    control_loss <- nrow(band_clonal[cat == "control" & cnloss == 1])
    control_noloss <- nrow(band_clonal[cat == "control" & cnloss == 0])

    pval_loss <- fisher.test(matrix(c(
        bi_loss, bi_noloss,
        control_loss, control_noloss
    ),
    nrow = 2
    ))$p.value

    row <- data.table(
        band = inputband,
        p_loss = pval_loss
    )
    return(row)
}

test_sample_cytobands <- function(inputdata, run_n = 0) {
    # Return a object and save results in the resultfile
    # -- result calculated all the pval and number of cases for all cytobands for clonal copy number changes
    all_band <- unique(inputdata$gen_loc)

    # Define a function to use the inputdata to test the band
    test_specific_sample_clonal_band <- function(x) test_clonal_band(x, testdata = inputdata)
    all_result <- rbindlist(lapply(all_band, test_specific_sample_clonal_band))
    names(all_result) <- c("band", paste0("p_loss_run_", as.character(run_n)))
    return(all_result)
}

## Need to run these result a 100 times and combined the reulst to save
# setup initial result - burn_in: run 500 times to generate some randomness
# clonal_result <- test_sample_cytobands(sim_n(clonal))
# brca_result <- test_sample_cytobands(sim_n(brca_data))
ov_sim <- start_burn_in_sim_n(ov_data, 500)
erpos_sim <- start_burn_in_sim_n(erpos_data, 500)
ov_result <- test_sample_cytobands(get_dt_from_sim(ov_sim))
erpos_result <- test_sample_cytobands(get_dt_from_sim(erpos_sim))
# tri_result <- test_sample_cytobands(sim_n(tri_data))

for (localrun in 2:100) {
    # for (localrun in 2:3) { # This is for debugging to see how long it takes to run
    # Converting the different tumor data to simuation data

    # This is the local jump to enable some randomness
    for (j in 1:50) {
        ov_sim <- sim9_single(ov_sim)
        erpos_sim <- sim9_single(erpos_sim)
    }
    ov_result <- merge(ov_result,
        test_sample_cytobands(get_dt_from_sim(ov_sim), localrun),
        by = "band"
    )
    erpos_result <- merge(erpos_result,
        test_sample_cytobands(get_dt_from_sim(erpos_sim), localrun),
        by = "band"
    )

    # progress and debug code
    writeLines(paste(" -- Working on localrun: ", localrun))
    writeLines(paste("    ncol =", ncol(ov_result)))
}

run_folder <- paste0("analysis/cellular_fraction_cn/output/randomize_background/runs/run_", runNumber, "/")
mkdirp(run_folder)

save_tumor_result <- function(result, tumor) {
    result_file_path <- paste0(run_folder, tumor, ".csv")
    fwrite(result, result_file_path)
}

## Save the calculation result
save_tumor_result(ov_result, "ov_result")
save_tumor_result(erpos_result, "erpos_result")