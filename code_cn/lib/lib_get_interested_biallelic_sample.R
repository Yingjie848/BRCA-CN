library(data.table)
#setwd("..")

get_interested_biallelic_samples <- function() {
    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)

    # hrd = fread("input/mc3_hrd_homdels_newSig3.csv")
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat != "Wild_type", ]
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    sample_cat <- list()
    # -- these are in the output
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples

    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_biallelic_gene_samples <- function(interested_gene_set) {
    # With input: interested_gene_set,
    # - return the biallelic samples of the geneset

    # Loading and prepare biallelic samples
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)

    # hrd = fread("input/mc3_hrd_homdels_newSig3.csv")
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    return(interested_biallelic_samples)
}

get_monoallelic_gene_samples <- function(interested_gene_set) {
    # With input: interested_gene_set,
    # - return the monoallelic samples (LoF) of the geneset
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd
    # Select the mono allelic samples
    hrd_mono <- hrd[gene %in% interested_gene_set & mc3_cat == "Mono_allelic_path"]
    mono_samples <- unique(hrd_mono$sample)
    return(mono_samples)
}

get_biallelic_samples_control_w_VUS <- function() {
    # Include VUS in the returned controls

    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat == "Biallelic_path" | mc3_cat == "Mono_allelic_path"] # This is different
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    # control samples in this case include VUS cases
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))
    sample_cat <- list()
    # -- these are in the output
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples


    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_biallelic_samples_control_w_VUS_highLST_domSig3 <- function() {
    # Remove samples of highLST, domSig3 from control to biallelic samples
    # Include VUS in the returned controls

    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat == "Biallelic_path" | mc3_cat == "Mono_allelic_path"] # This is different
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    # control samples in this case include VUS cases
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    # Identify high LST and domSig3 control samples
    m <- fread("input/master_score_MC3_PanImmun.csv")
    m$sample <- m$Id
    m <- m[sample %in% control_samples]
    mdom <- m[LST >= 15 & dominant == "Signature.3"]
    control_samples <- control_samples[!control_samples %in% mdom$sample]
    interested_biallelic_samples <- unique(c(interested_biallelic_samples, mdom$sample))


    # -- these are in the output
    sample_cat <- list()
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples


    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_biallelic_samples_control_w_VUS_noHiLST_domSig3 <- function() {
    # Remove samples of highLST, domSig3 from control to biallelic samples
    # Include VUS in the returned controls

    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat == "Biallelic_path" | mc3_cat == "Mono_allelic_path"] # This is different
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    # control samples in this case include VUS cases
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    # Identify high LST and domSig3 control samples
    m <- fread("input/master_score_MC3_PanImmun.csv")
    m$sample <- m$Id
    m <- m[sample %in% control_samples]
    mdom <- m[LST >= 15 & dominant == "Signature.3"]
    control_samples <- control_samples[!control_samples %in% mdom$sample]

    # -- these are in the output
    sample_cat <- list()
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples


    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    ## Adding a other category for samples which are not biallelic nor control
    others <- all_samples[!(all_samples %in% interested_biallelic_samples)
    & !(all_samples %in% control_samples)]
    sample_cat$others <- others
    return(sample_cat)
}

#a= get_biallelic_samples_control_w_VUS_noHiLST_domSig3()
#saveRDS(a, file="~/biallelic_tcga_sample_lists_biallelic-cn-project.rds")

get_core_biallelic_samples <- function() {
    # Return sample category of biallelic core genes
    # Loading and prepare biallelic samples
    # -- this could be further simplified
    core <- fread("input/gene_lists/genelist_core.txt", header = F)
    interested_gene_set <- core$V1

    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)

    # hrd = fread("input/mc3_hrd_homdels_newSig3.csv")
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat != "Wild_type", ]
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    sample_cat <- list()
    # -- these are in the output
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples

    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_biallelic_samples_control_w_monoLoF <- function() {
    # Include VUS in the returned controls

    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat == "Biallelic_path"]
    hrd_not_control_samples <- unique(hrd_not_control$sample)
    # control samples in this case include VUS cases
    control_samples <- all_samples[!all_samples %in% hrd_not_control_samples]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    sample_cat <- list()
    # -- these are in the output
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples

    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]
    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_biallelic_samples_control_w_monoVUS <- function() {
    # Include VUS in the returned controls

    # Loading and prepare biallelic samples
    interested_gene_set <- c("BRCA1", "BRCA2", "PALB2")
    bi <- fread("input/Biallelic_Pathogenic_hits_MC3.txt")
    homdel <- bi[gene %in% interested_gene_set & status == "homdel"]
    homdel_samples <- gsub(".", "-", unique(homdel$Id), fixed = T)
    hrd <- fread("input/hrd_mc3_master102gene.csv") # update using the fixed hrd

    # calculate control samples
    all_samples <- unique(hrd$sample)
    hrd_not_control <- hrd[mc3_cat == "Biallelic_path" | mc3_cat == "Mono_allelic_path"] # This is different
    hrd_not_control_samples <- unique(hrd_not_control$sample)

    # Calcuating for biVUS, monoVUS
    monoLoF <- hrd[mc3_cat == "Mono_allelic_path"] # all mono allelic path
    biVUS <- hrd[mc3_cat == "Biallelic_VUS" &
        !(sample %in% monoLoF$sample)] # biVUS but not LoF
    monoVUS <- hrd[mc3_cat == "Mono_allelic_VUS" &
        !(sample %in% monoLoF$sample) &
        !(sample %in% biVUS$sample)] # monoVUS but not Lof

    # control samples in this case include VUS cases
    control_samples <- all_samples[!(all_samples %in% hrd_not_control_samples |
        all_samples %in% biVUS$sample)]

    # calculate biallelic samples
    hrd_bi <- hrd[gene %in% interested_gene_set & mc3_cat == "Biallelic_path"]
    biallelic_samples <- unique(hrd_bi$sample)
    interested_biallelic_samples <- unique(c(homdel_samples, biallelic_samples))

    sample_cat <- list()
    # -- these are in the output
    sample_cat$biallelic <- interested_biallelic_samples
    sample_cat$control <- control_samples
    sample_cat$all <- all_samples

    # Adding triple negative sample and ER pos sample in this object
    brca <- fread("input/brca_tcga_clinical_data.tsv")
    brca$sample <- substr(brca[[3]], 1, 12)
    triple_negative <- brca[Receptor_Status_Patient == "Triple Negative"]
    er_positive <- brca[grepl("ER+", Receptor_Status_Patient, fixed = T)]

    # -- these are in the output
    sample_cat$er_positive <- unique(er_positive$sample)
    sample_cat$triple_negative <- unique(triple_negative$sample)

    ## Further identify monoallelic VUS/monoLoF
    sample_cat$monoLoF <- unique(monoLoF$sample)
    sample_cat$biVUS <- unique(biVUS$sample)
    sample_cat$monoVUS <- unique(monoVUS$sample)

    return(sample_cat)
}

get_samples_basedon_gene_and_function <- function(interested_genes, base_function) {
    # return similar object as get_interested_biallelic_samples
    # However, for the controls, it will be the same as in base_function
    sample_cat <- base_function()
    sample_cat$biallelic <- get_biallelic_gene_samples(interested_genes)
    return(sample_cat)
}

get_icgc_sig3 <- function() {
    # Return the signature 3 calculation for ICGC samples
    s <- fread("/data/share/icgc/pcwag/mutational_signatures_in_samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")

    s[, total := rowSums(.SD), .SDcols = 4:68]
    s$sig3 <- s$SBS3 / s$total
    a <- s[, 4:48]
    a$domSig <- names(a)[apply(a, 1, which.max)]
    s$domSig <- a$domSig
    snames <- names(s)
    snames[1] <- "cancer_type"
    snames[2] <- "sample_name"
    names(s) <- snames
    ss <- s[, .(cancer_type, sample_name, sig3, total, SBS3, domSig)]
    return(ss)
}

get_serena_signatures <- function() {
    ## Signature and other info from Serena's paper
    serena_sigs <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table21_signature_contribution.csv")
    serena_lst <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table17_HRD_intermediary.csv")
    serena_sigs[, total := rowSums(.SD), .SDcols = 2:13]
    serena_sigs$sig3 <- serena_sigs[[4]] / serena_sigs$total
    # -- For these signature cutoff, for now just using sig3 >=0.3, LST>=20
    serena_sigs$sample <- gsub("(a|b).*$", "", serena_sigs[[1]])
    serena_lst$sample <- gsub("(a|b).*$", "", serena_lst$Sample)
    signatures <- merge(serena_sigs, serena_lst[, .(sample, LST)], by = "sample")
    return(signatures)
}

get_davies_biallelic_samples <- function() {
    sample <- fread("analysis/icgc_dataset/analysis/brca_eu_biallelic_calc/input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)
    # This is the cohort we are using for final testing 320 ER+/HER2-
    allSamples <- unique(sample$sampleid)

    # Return the biallelic/control samples list based on Davies HRDetect paper
    ## Loading mutation information to determine biallelic vs. control samples
    brca_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4a_mutation_details.csv")
    brca_mut$sample <- gsub("(a|b).*$", "", brca_mut$Sample)
    # -- use this to define the 77 cases with biallelic
    ddr_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4b_otherHRgene_mutation_details.csv")
    ddr_mut$sample <- gsub("(a|b).*$", "", ddr_mut$Sample)

    # -- use this to remote the other DDR gene (PALB2) mutations
    davies_sample <- fread("input/Davies_HRDetect_2017/Davies_suppltable_1_sampleInfo.csv")

    # Organize breast and OV samples list
    biallelic <- brca_mut[isBrcaMonoallelic == F]
    monoallelic <- ddr_mut[Gene == "PALB2"]
    signatures <- get_serena_signatures()
    dom_signatures <- signatures[signatures$LST >= 15 & signatures$sig3 >= 0.2]

    # Finalize sample categories
    control_samples <- allSamples[!(allSamples %in% biallelic$sample) &
        !(allSamples %in% monoallelic$sample) &
        !(allSamples %in% dom_signatures$sample)]

    biallelic_samples <- unique(biallelic$sample)

    fr <- list()
    fr$control <- control_samples
    fr$biallelic <- biallelic_samples
    return(fr)
}
