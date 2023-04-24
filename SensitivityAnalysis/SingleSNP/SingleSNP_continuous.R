# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for continuous data                            #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
#  MR for continuous data using UK Biobank patient-level data.                 #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- read_delim(here("SensitivityAnalysis","SingleSNP","exposure_data.csv")) %>%
  select(-"...1")

# Outcome data -----------------------------------------------------------------
out  <- c('eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A','Apo-B','CRP','Lipoprotein','HbA1c','Glucose')
outc <- c('3148-0.0','30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0','30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0')

unlink(here("SensitivityAnalysis","SingleSNP","MR_res_continuous.xlsx"))
unlink(here("SensitivityAnalysis","SingleSNP","MR_dat_continuous.xlsx"))
t <- read_delim(paste0(pathData,"UKB\\cohort.csv"))

for (i in 1:11){
  t_outcome <- t %>%
    left_join(
      read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>%
        select('eid',
               'outcome' = outc[i]),
      by = 'eid') %>%
    filter(!is.na(outcome)) %>%
    mutate(outcome = (outcome - mean(outcome,na.rm = TRUE))/sd(outcome, na.rm = TRUE)) %>% # SD Units
    select('eid','outcome',exposure_dat$SNP,'sex',contains('PC'))
  
    n <- nrow(t_outcome)
  # LINEAR REGRESSION ----------------------------------------------------------
  # Eliminate those rows with NaN
  tab <- t_outcome %>%
    select("eid","sex","outcome", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
           "PC6","PC7","PC8","PC9","PC10")
  
  # Substitute the 1 for 2
  eaf <- sum(tab$snp)/(2*length(tab$snp))
  tab <- tab %>%
    mutate(snp = if_else(snp == 2, 1, snp))
  
  # Linear regression - Adjusted model for sex and the first 10 principal components
  mlr <- lm(outcome ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab)
  mlr <- coefficients(summary(mlr))
  
  # Coefficients:
  beta  <- mlr[2,1]
  se    <- mlr[2,2]
  p     <- mlr[2,4]

  # OUTCOME DATA -----------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = exposure_dat$SNP, 
                            "beta.outcome" = beta, 
                            "se.outcome" = se,
                            "p.outcome" = p, 
                            "eaf.outcome" = eaf, 
                            "n.outcome" = n) %>%
    mutate(id.outcome = out[i],
           outcome    = out[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  dat <- harmonise_data(exposure_dat, outcome_dat)
  res <- mr(dat) %>%
    rename('Estimate' = 'b') %>%
    mutate(Estimate   = -Estimate, # Lowering sclerostin
           CI_LOW   = Estimate-1.96*se,
           CI_HIGH  = Estimate+1.96*se)
  
  write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_continuous.xlsx"), sheetName = out[i], append = TRUE)
  write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_continuous.xlsx"), sheetName = out[i], append = TRUE)
}




