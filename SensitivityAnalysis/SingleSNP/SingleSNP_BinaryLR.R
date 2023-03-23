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

for (i in 1:length(outc)){
  t <- read_delim(here("MR_UKBiobank","BinaryData","Logistic",paste0("Phenotype_",outc[i],".csv"))) %>%
    select("eid",
           "outcome" = "state") %>%
    inner_join(read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>%
                 select("eid",
                        "Caucasian" = "22006-0.0"), # Caucasian
               by = "eid") %>%
    filter(!is.na(outcome)) %>%
    filter(Caucasian == 1) %>%
    inner_join(read_delim(paste0(pathData,"genetics.csv")) %>% # SNPs expression
                 select(-"...1"),
               by = "eid") %>%
    inner_join(read_delim(paste0(pathData,"UKB\\PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS
  
  # LOGISTIC REGRESSION ----------------------------------------------------------
  # Eliminate those rows with NaN
  tab <- t %>%
    select("eid","sex","outcome", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
           "PC6","PC7","PC8","PC9","PC10") %>%
    filter(!is.na(sex),
           !is.na(snp),!is.na(PC1),!is.na(PC2),!is.na(PC3),!is.na(PC4),!is.na(PC5),
           !is.na(PC6),!is.na(PC7),!is.na(PC8),!is.na(PC9),!is.na(PC10))
    
  # Substitute the 1 for 2
  tab1 <- tab
  tab1$snp[tab$snp == 2] <- 1
    
  #Logistic regression - Adjusted model for sex and the first 10 principal components
  logistic_model <- glm(outcome ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab1,
                        family = "binomial")
  logistic_model <- summary(logistic_model)
  
  # Coefficients
  n    <- nrow(tab1)
  eaf  <- (sum(tab1$snp)/(2*length(tab1$snp)))
  beta <- logistic_model$coefficients[2,1]
  se   <- logistic_model$coefficients[2,2]
  p  <- logistic_model$coefficients[2,4]

  
  # OUTCOME DATA -----------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = exposure_dat$SNP, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION
  dat <- harmonise_data(exposure_dat, outcome_dat)
  res <- mr(dat)
  
  write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_BinaryLR.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_BinaryLR.xlsx"), sheetName = outc[i], append = TRUE)
}







