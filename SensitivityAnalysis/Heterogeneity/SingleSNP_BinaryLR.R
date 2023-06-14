# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                             Logistic regression                              #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #      
# Generalized MR on binary data using UK Biobank patient level data. Here, it  #
# is used a logistic regression to calculate the effect of the instruments on  #
# the outcomes.                                                                #
# Also, sensitivity analysis "leave-one-out" is performed.
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- read_delim(here("SensitivityAnalysis","SingleSNP","exposure_data.csv")) %>%
  select(-"...1")

# Outcome data -----------------------------------------------------------------
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

unlink(here("SensitivityAnalysis","SingleSNP","MR_res_BinaryLR.xlsx"))
unlink(here("SensitivityAnalysis","SingleSNP","MR_dat_BinaryLR.xlsx"))

t <- read_delim(paste0(pathData,"UKB\\cohort.csv"))

for (i in 1:length(outc)){
  t_outcome <- t %>%
    left_join(
      read_delim(paste0(pathData,"MR_UKBiobank\\BinaryData\\Logistic\\Phenotype_",outc[i],".csv")) %>%
        select('eid',
               'outcome' = 'state'),
      by = 'eid') %>%
    filter(!is.na(outcome))
  
  n <- nrow(t_outcome)
  
  # LOGISTIC REGRESSION ----------------------------------------------------------
  # Eliminate those rows with NaN
  # Eliminate those rows with NaN
  tab <- t_outcome %>%
    select("eid","sex","outcome", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
           "PC6","PC7","PC8","PC9","PC10")
  
  # Substitute the 1 for 2
  eaf <- sum(tab$snp)/(2*length(tab$snp))
  tab <- tab %>%
    mutate(snp = if_else(snp == 2, 1, snp))
  
  #Logistic regression - Adjusted model for sex and the first 10 principal components
  logistic_model <- glm(outcome ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab,
                        family = "binomial")
  logistic_model <- coefficients(summary(logistic_model))
  
  # Coefficients
  beta <- logistic_model[2,1]
  se <- logistic_model[2,2]
  p  <- logistic_model[2,4]
  
  # OUTCOME DATA -----------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = exposure_dat$SNP, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION
  dat <- harmonise_data(exposure_dat, outcome_dat)
  res <- mr(dat) %>%
    mutate('OR' = exp(-b),
           'CI_LOW_OR'  = exp(-b-1.96*se),
           'CI_HIGH_OR' = exp(-b+1.96*se))
  
  write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_BinaryLR.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_BinaryLR.xlsx"), sheetName = outc[i], append = TRUE)
}







