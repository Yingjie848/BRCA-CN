## -----------------------------------------------------------------------------
# identify potential genes those enhance proliferation in BRCA1-null and BRCA2-null cells

## -----------------------------------------------------------------------------
library(dplyr)

## -----------------------------------------------------------------------------
define_BF_signif <- function(x,x_cutoff=c(0,0),use_quantile_cutoff=FALSE,quartile=0.25){

    if(isTRUE(use_quantile_cutoff)){
      x_negative = x[x<0]
      x_positive = x[x>0]
      x_cutoff=c(quantile(x_negative,probs=1-quartile),quantile(x_positive,probs=quartile))
    }

    print(x_cutoff)

    q <- ifelse(x > 0 & x > x_cutoff[2],'signif +',
         ifelse(x < 0 & x < x_cutoff[1],'signif -',
                'ns'))
    q
}

## -----------------------------------------------------------------------------
load("../crispr_results.Rdata")


## -----------------------------------------------------------------------------
# define extreme BRCA - WT
quartile_cutoff = 0.5

crispr <- crispr %>% mutate(

	DLD1_MSK_Brca2_minus_wt.BF.BFextreme       = define_BF_signif(DLD1_MSK_Brca2_minus_wt.BF,use_quantile_cutoff = T,quartile = quartile_cutoff),

	DLD1_REPARE_Brca2_minus_wt.BF.BFextreme    = define_BF_signif(DLD1_REPARE_Brca2_minus_wt.BF,use_quantile_cutoff = T,quartile = quartile_cutoff),

	RPE1_REPARE_Brca1_minus_wt.BF.BFextreme    = define_BF_signif(RPE1_REPARE_Brca1_minus_wt.BF,use_quantile_cutoff = T,quartile = quartile_cutoff)

)

## -----------------------------------------------------------------------------
# genes lost in OV BRCA1 tumors, and had down-regulated expression
crispr.cnloss.brca1 = crispr[crispr$cLLL.OV_BRCA1=='loss,low (CNloss cytoband)',]; nrow(crispr.cnloss.brca1)
# genes lost in OV or ER+ BRCA2 tumors, and had down-regulated expression
crispr.cnloss.brca2 = crispr[(crispr$cLLL.OV_BRCA2=='loss,low (CNloss cytoband)' | crispr$cLLL.ER_BRCA2=='loss,low (CNloss cytoband)'),]; nrow(crispr.cnloss.brca2)

## -----------------------------------------------------------------------------
# BRCA1, only based on method 2
crispr.cnloss.brca1_m2 = crispr.cnloss.brca1[crispr.cnloss.brca1$RPE1_REPARE_Brca1_minus_wt.BF.BFextreme == 'signif -',]; nrow(crispr.cnloss.brca1_m2)
crispr.cnloss.brca1_m2$Gene

## -----------------------------------------------------------------------------
# BRCA2, only based on method 2
crispr.cnloss.brca2_m2 = crispr.cnloss.brca2[crispr.cnloss.brca2$DLD1_REPARE_Brca2_minus_wt.BF.BFextreme == 'signif -',]; nrow(crispr.cnloss.brca2_m2)
crispr.cnloss.brca2_m2$Gene

## -----------------------------------------------------------------------------
save(crispr,crispr.cnloss.brca1,crispr.cnloss.brca2,crispr.cnloss.brca1_m2,crispr.cnloss.brca2_m2,file="crispr_results.20230501.Rdata")

write.csv(crispr.cnloss.brca1_m2      ,"crispr_results.cnloss.brca1_potential_277genes.csv",row.names=F)
write.csv(crispr.cnloss.brca2_m2      ,"crispr_results.cnloss.brca2_potential_218genes.csv",row.names=F)
