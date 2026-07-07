############################################################################
#       --- results : COMPARISON ---
############################################################################

# upload data
setwd("~/Documents/PROJECTS/social_mobility/simulations_data")

load("scenario_1.RData")
load("scenario_2.RData")
load("scenario_3.RData")
load("scenario_4.RData")

load("BPCF_1.RData")
load("BPCF_2.RData")
load("BPCF_3.RData")
load("BPCF_4.RData")

load("results_doubleBNP_s1.RData")
load("results_doubleBNP_s2.RData")
load("results_doubleBNP_s3.RData")
load("results_doubleBNP_s4.RData")

load("SLM_X_1.RData")
load("SLM_X_2.RData")


############################################################################

# see 'table' functions in 2_stan_estimation.R and 2_competitor_estimation.R


total_table<- function(models){
  
  table_all <- do.call(cbind, models)
  names_col <- paste0(rep(names(models),each =2), "-", colnames(table_all))
  
  colnames(table_all) <- names_col
  
  return(round(table_all,4))
}

total_table(models=list(BNP = table_doubleBNP_s1, BPCF = table_BPCF_s1, SLM = table_SLM_s1))
total_table(models=list(BNP = table_doubleBNP_s2, BPCF = table_BPCF_s2, SLM = table_SLM_s2))
total_table(models=list(BNP = table_doubleBNP_s3, BPCF = table_BPCF_s3))
total_table(models=list(BNP = table_doubleBNP_s4, BPCF = table_BPCF_s4))
