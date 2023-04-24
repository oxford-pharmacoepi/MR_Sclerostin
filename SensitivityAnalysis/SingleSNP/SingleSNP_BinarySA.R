# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                               Survival analysis                              #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok")))

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- read_delim(here("SensitivityAnalysis","SingleSNP","exposure_data.csv")) %>%
  select(-"...1")

# Outcome data -----------------------------------------------------------------
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

unlink(here("SensitivityAnalysis","SingleSNP","MR_res_BinarySA.xlsx"))
unlink(here("SensitivityAnalysis","SingleSNP","MR_dat_BinarySA.xlsx"))

t <- read_delim(paste0(pathData,"UKB\\cohort.csv"))

for (i in 1:length(outc)){
  t_outcome <- t %>%
    inner_join(
      read_delim(paste0(pathData,"MR_UKBiobank\\BinaryData\\SA_Enrolment\\Phenotype_",outc[i],".csv")) %>%
        select("eid",
               "state",
               "age"),
      by = 'eid')
  
  n        <- nrow(t_outcome)
  
  # SURVIVAL ANALYSIS ------------------------------------------------------------
  # Eliminate those rows with NaN
  tab <- t_outcome %>%
    select("eid","sex","state","age", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
           "PC6","PC7","PC8","PC9","PC10")
  
  # Substitute the 1 for 2
  eaf <- sum(tab$snp)/(2*length(tab$snp))
  tab <- tab %>%
    mutate(snp = if_else(snp == 2, 1, snp))
  
  # Survival analysis - Adjusted model for sex and the first 10 principal components
  cox <- coxph(Surv(age,state) ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab)
  cox <- coefficients(summary(cox))
  
  # Coefficients
  beta <- cox[1,1]
  hr   <- cox[1,2]
  se   <- cox[1,3]
  p    <- cox[1,5]
  
  # OUTCOME DATA ---------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = exposure_dat$SNP, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION ----------------------------------------------------
  dat <- harmonise_data(exposure_dat,outcome_dat)
  res <- mr(dat) %>%
    mutate(OR = exp(-b),
           CI_LOW_OR = exp(-b-1.96*se),
           CI_HIGH_OR = exp(-b+1.96*se))
    
  
  write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_BinarySA.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_BinarySA.xlsx"),  sheetName = outc[i], append = TRUE)
}



















