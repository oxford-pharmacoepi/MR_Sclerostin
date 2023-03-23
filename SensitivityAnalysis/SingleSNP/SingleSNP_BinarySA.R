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

for (i in 1:length(outc)){
  t    <- read_delim(here("MR_UKBiobank","BinaryData","SA_Birth",paste0("Phenotype_",outc[i],".csv"))) %>%
    select("eid",
           "state",
           "age") %>%
    inner_join(read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>%
                 select("eid",
                        "Caucasian" = "22006-0.0"), # Caucasian
               by = "eid") %>%
    filter(Caucasian == 1) %>%
    inner_join(read_delim(paste0(pathData,"genetics.csv")) %>% # SNPs expression
                 select(-"...1"),
               by = "eid") %>%
    inner_join(read_delim(paste0(pathData,"UKB\\PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS
  
  # SURVIVAL ANALYSIS ------------------------------------------------------------

  # Eliminate those rows with NaN
  tab <- t %>%
    select("eid","sex","state","age", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
           "PC6","PC7","PC8","PC9","PC10") %>%
    filter(!is.na(age),!is.na(sex),
           !is.na(snp),!is.na(PC1),!is.na(PC2),!is.na(PC3),!is.na(PC4),!is.na(PC5),
           !is.na(PC6),!is.na(PC7),!is.na(PC8),!is.na(PC9),!is.na(PC10))
  
  # Substitute the 1 for 2
  tab1 <- tab
  tab1$snp[tab$snp == 2] <- 1
  
  #Logistic regression - Adjusted model for sex and the first 10 principal components
  cox <- coxph(Surv(age,state) ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab1)
  cox <- summary(cox)
  
  # Coefficients
  n    <- nrow(tab1)
  eaf  <- (sum(tab1$snp)/(2*length(tab1$snp)))
  beta <- cox$coefficients[1,1]
  hr   <- cox$coefficients[1,2]
  se <- cox$coefficients[1,3]
  p  <- cox$coefficients[1,5]
  
  
  # OUTCOME DATA ---------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = exposure_dat$SNP, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION ----------------------------------------------------
  dat <- harmonise_data(exposure_dat,outcome_dat)
  res <- mr(dat)
  
  write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_BinarySA.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_BinarySA.xlsx"),  sheetName = outc[i], append = TRUE)
}



















