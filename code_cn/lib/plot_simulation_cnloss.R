# Plot the histogram of the simulation result from ecoSimR and naive sampleing method

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(ggpubr)

save_sig_fdrpval_counts <- function(type='loss') {
    get_sig_pval_count <- function(data_filepath) {
        # Given our simutation reuslt data file path, return a count object
        # d <- fread("analysis/cellular_fraction_cn/output/randomize_background/runs/run_1/ov_result.csv")
        d <- fread(data_filepath)
        d2 <- d[, 2:ncol(d)]

        count_fdrpval_num <- function(x) {
            # with the input x, do fdr adjust, and return the number of significant counts
            ax <- p.adjust(x, method = "fdr")
            nsig <- sum(ax <= 0.05)
            return(nsig)
        }

        pcount <- data.table(pval_count = apply(d2, 2, count_fdrpval_num))
        return(pcount)
    }

    # list all the run folders and obtain each tumor type
    get_run_result_tumor <- function(rN, tumor) {
        fpath <- paste0("analysis/cellular_fraction_cn/output/randomize_background_", type,"/runs/run_", rN, "/", tumor, "_result.csv")
        pval <- get_sig_pval_count(fpath)
        return(pval)
    }

    get_run_ov <- function(x) get_run_result_tumor(x, "ov")
    get_run_erpos <- function(x) get_run_result_tumor(x, "erpos")

    # Run through the 100x runs to get ov, and erpos values
    ov <- rbindlist(lapply(1:100, get_run_ov))
    erpos <- rbindlist(lapply(1:100, get_run_erpos))
    # Saving the result for plotting next step
    path = paste0("analysis/cellular_fraction_cn/output/randomize_background_", type, "/sig_fdrpval_distribution/")
    mkdirp(path)
    fwrite(ov, paste0( path, "ov_sig_fdrpval_counts.csv"))
    fwrite(erpos, paste0(path, "erpos_sig_fdrpval_counts.csv"))
}

save_ecosim_figure <- function(type='loss', ovband=41, erband=25, xlimmax=100) {

    # Loading saved data for plotting
    path = paste0("analysis/cellular_fraction_cn/output/randomize_background_", type, "/sig_fdrpval_distribution/")
    ov <- fread(paste0(path, "ov_sig_fdrpval_counts.csv"))
    erpos <- fread(paste0(path, "erpos_sig_fdrpval_counts.csv"))

    xlabel = paste("Number of significant CN", type,  "cytobands")
    p.ov <- gghistogram(ov,
        x = "pval_count", add = "mean", rug = T,
        xlab = xlabel,
        color = "lightblue",
        fill = "lightblue",
        bins = 100,
        xlim = c(0, xlimmax),
        ylim = c(0, 3000),
        title = "Histogram of significant FDR cytoband counts in OV"
    ) +
        geom_vline(xintercept = ovband, color = "green")

    xlabel = paste("Number of significant CN", type,  "cytobands")
    p.er <- gghistogram(erpos,
        x = "pval_count", add = "mean", rug = T,
        xlab = xlabel,
        color = "orange",
        fill = "orange",
        bins = 100,
        xlim = c(0, xlimmax),
        ylim = c(0, 3000),
        title = "Histogram of significant FDR cytoband counts in ER+ breast tumors"
    ) +
        geom_vline(xintercept = erband, color = "green")

    fig <- ggarrange(p.ov, p.er,
        labels = c("OV", "ER+"),
        ncol = 1, nrow = 2
    )

    ggsave(paste0(path, "ecosim_sig_fdrpval_distribution.png"), fig)
}

save_naivesim_figure <- function() {

    # Loading saved data for plotting
    ov <- rbindlist(list(
        fread("analysis/cellular_fraction_cn/output/naive_sampling_simutation/batch1/ov_sig_fdrp_loss_counts.csv"),
        fread("analysis/cellular_fraction_cn/output/naive_sampling_simutation/batch2/ov_sig_fdrp_loss_counts.csv")
    ))

    erpos <- rbindlist(list(
        fread("analysis/cellular_fraction_cn/output/naive_sampling_simutation/batch1/er_sig_fdrp_loss_counts.csv"),
        fread("analysis/cellular_fraction_cn/output/naive_sampling_simutation/batch2/er_sig_fdrp_loss_counts.csv")
    ))

    # Remove unfinished runs
    ov <- ov[n_row == 829]
    erpos <- erpos[n_row == 829]
    ov <- ov[1:10000]
    erpos <- erpos[1:10000]

    p.ov <- gghistogram(ov,
        x = "nsig", add = "mean", rug = T,
        xlab = "Number of significant CN loss cytobands",
        color = "lightblue",
        fill = "lightblue",
        bins = 100,
        xlim = c(0, 100),
        ylim = c(0, 3000),
        title = "Histogram of significant FDR cytoband counts in OV"
    ) +
        geom_vline(xintercept = 91, color = "green")

    p.er <- gghistogram(erpos,
        x = "nsig", add = "mean", rug = T,
        xlab = "Number of significant CN loss cytobands",
        color = "orange",
        fill = "orange",
        bins = 100,
        xlim = c(0, 100),
        ylim = c(0, 3000),
        title = "Histogram of significant FDR cytoband counts in ER+ breast tumors"
    ) +
        geom_vline(xintercept = 32, color = "green")

    fig <- ggarrange(p.ov, p.er,
        labels = c("OV", "ER+"),
        ncol = 1, nrow = 2
    )

    ggsave("analysis/cellular_fraction_cn/output/randomize_background/sig_fdrpval_distribution/naivesim_sig_fdrpval_distribution.png", fig)
}

######################################################################
main <- function() {
    # Parse through the simulation result
    # save_sig_fdrpval_counts(type="loss")
    # save_sig_fdrpval_counts(type="gain")

    # Plot figures
    save_ecosim_figure(type="loss")
    save_ecosim_figure(type="gain", ovband=9, erband=145,xlimmax=150)
    # save_naivesim_figure()
}
main()