# Steps in how to calculate significant changed bands
# -- These are pesudo make code, to note the steps in generating the final results
# -- the path folder are all relevant to this file 


# Generating the calculation result
# -- the folder "cnloss_rounded_ploidy_minus_1" is the particular version definition for the script
output/cnloss_rounded_ploidy_minus_1/*.csv:
	analysis/cellular_fraction_cn/lib/find_enriched_cn_change_ploidy_cytobands.R::save_result_tumortypes()
	mv analysis/cellular_fraction_cn/output/*.csv \
		output/cnloss_rounded_ploidy_minus_1/

# Compile the results to a R object
analysis/cellular_fraction_cn/input/sig_cytobands_cnloss_rounded_ploidy_minus_1/cytobands_level_fdr_tests.rda:
	analysis/cellular_fraction_cn/lib/find_enriched_cn_change_ploidy_cytobands.R:: load_process_save_cytoband_result("cnloss_rounded_ploidy_minus_1")

# Generating result for ICGC samples (OV) 
analysis/icgc_dataset/output/icgc_cn_biallelic_tests/*.csv:
	analysis/icgc_dataset/lib/s6_calc_cn_change_cytobands.R::save_result_tumortypes()

# Generating Serena ICGC ER+ breast tumor result
analysis/icgc_dataset/output/icgc_cn_biallelic_tests/*.csv:
	analysis/icgc_dataset/lib/s8_calc_serena_breast560_cn_change_cytobands.R::	test_sample_cytobands(erpos_data, "serena_ERpos_breastTumor_cn_test_stats")
