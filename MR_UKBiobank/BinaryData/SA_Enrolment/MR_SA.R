# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                               Survival analysis                              #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok")))

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv"))
snps <- ldrho %>% select("SNP"= RS_number) # SNPs

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1")

# Outcome data -----------------------------------------------------------------
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

unlink(here("SensitivityAnalysis","Enrolment","gMR_res.xlsx"))
unlink(here("SensitivityAnalysis","Enrolment","gMR_dat.xlsx"))
for (i in 1:length(outc)){
  t    <- read_delim(here("SensitivityAnalysis","Enrolment","SA_Enrolment",paste0("Phenotype_",outc[i],".csv"))) %>%
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
  beta  <- matrix(0,nrow(snps),1)
  se    <- matrix(0,nrow(snps),1)
  p     <- matrix(0,nrow(snps),1)
  eaf   <- matrix(0,nrow(snps),1)
  n     <- matrix(0,nrow(snps),1)
  hr    <- matrix(0,nrow(snps),1)
  
  for (j in(c(1:nrow(snps)))){
    # Eliminate those rows with NaN
    tab <- t %>%
      select("eid","sex","state","age", "snp" = snps[j,1],"PC1","PC2","PC3","PC4","PC5",
             "PC6","PC7","PC8","PC9","PC10") %>%
      filter(!is.na(age),!is.na(sex),
             !is.na(snp),!is.na(PC1),!is.na(PC2),!is.na(PC3),!is.na(PC4),!is.na(PC5),
             !is.na(PC6),!is.na(PC7),!is.na(PC8),!is.na(PC9),!is.na(PC10))
    
    # Substitute the 1 for 2
    tab1 <- tab
    tab1$snp[tab$snp == 2] <- 1
    
    # Survival analysis - Adjusted model for sex and the first 10 principal components
    cox <- coxph(Surv(age,state) ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab1)
    cox <- summary(cox)
    
    # Coefficients
    n[j]    <- nrow(tab1)
    eaf[j]  <- (sum(tab1$snp)/(2*length(tab1$snp)))
    beta[j] <- cox$coefficients[1,1]
    hr[j]   <- cox$coefficients[1,2]
    se[j] <- cox$coefficients[1,3]
    p[j]  <- cox$coefficients[1,5]
  }
  
  # OUTCOME DATA ---------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = snps, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION ----------------------------------------------------
  source(here("Functions","gMR.R"))
  a <- gMR(exposure_dat,outcome_dat,ldrho)
  tab <- a[[2]] %>% mutate(OR = exp(Estimate),
                           CI_LOW_OR = exp(Estimate - 1.96*SE),
                           CI_HIGH_OR = exp(Estimate + 1.96*SE))
  
  write.xlsx(tab,here("SensitivityAnalysis","Enrolment","gMR_res.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(a[[1]],here("SensitivityAnalysis","Enrolment","gMR_dat.xlsx"),  sheetName = outc[i], append = TRUE)
}



















