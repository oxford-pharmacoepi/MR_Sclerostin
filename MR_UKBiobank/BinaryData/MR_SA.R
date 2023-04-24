# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                               Survival analysis                              #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv"))
snps <- ldrho %>% select("SNP"= RS_number) # SNPs

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1")

# Outcome data -----------------------------------------------------------------
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

unlink(here("MR_UKBiobank","BinaryData","SA_Birth","gMR_res.xlsx"))
unlink(here("MR_UKBiobank","BinaryData","SA_Birth","gMR_dat.xlsx"))
unlink(here("SensitivityAnalysis",'LeaveOneOut',"BinaryData_SA","loo.xlsx"))

t <- read_delim(paste0(pathData,"UKB\\cohort.csv"))

for (i in 1:length(outc)){
  t_outcome <- t %>%
    inner_join(
      read_delim(paste0(pathData,"MR_UKBiobank\\BinaryData\\SA_Birth\\Phenotype_",outc[i],".csv")) %>%
        select("eid",
               "state",
               "age"),
      by = 'eid') 
  
  cases    <- sum(t_outcome$state == 1)
  controls <- sum(t_outcome$state == 0)
  n        <- nrow(t_outcome)
  
  # SURVIVAL ANALYSIS ------------------------------------------------------------
  beta  <- matrix(0,nrow(snps),1)
  se    <- matrix(0,nrow(snps),1)
  p     <- matrix(0,nrow(snps),1)
  eaf   <- matrix(0,nrow(snps),1)
  hr    <- matrix(0,nrow(snps),1)
 
  for (j in(c(1:nrow(snps)))){
    # Eliminate those rows with NaN
    tab <- t_outcome %>%
      select("eid","sex","state","age", "snp" = snps[j,1],"PC1","PC2","PC3","PC4","PC5",
             "PC6","PC7","PC8","PC9","PC10")
    
    # Substitute the 1 for 2
    eaf[j] <- sum(tab$snp)/(2*length(tab$snp))
    tab <- tab %>%
      mutate(snp = if_else(snp == 2, 1, snp))
    
    # Survival analysis - Adjusted model for sex and the first 10 principal components
    cox <- coxph(Surv(age,state) ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab)
    cox <- coefficients(summary(cox))
    
    # Coefficients
    beta[j] <- cox[1,1]
    hr[j]   <- cox[1,2]
    se[j] <- cox[1,3]
    p[j]  <- cox[1,5]
  }
  
  # OUTCOME DATA ---------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = snps, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome  = exposure_dat$other_allele.exposure,
           n.outcome  = n,
           cases.outcome    = cases,
           controls.outcome = controls)
  
  # MENDELIAN RANDOMIZATION ----------------------------------------------------
  source(here("Functions","gMR.R"))
  a <- gMR(exposure_dat,outcome_dat,ldrho)
  tab <- a[[2]] %>% mutate(OR = exp(Estimate),
                           CI_LOW_OR = exp(Estimate - 1.96*SE),
                           CI_HIGH_OR = exp(Estimate + 1.96*SE))
  
  write.xlsx(tab,here("MR_UKBiobank","BinaryData","SA_Birth","gMR_res.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(a[[1]],here("MR_UKBiobank","BinaryData","SA_Birth","gMR_dat.xlsx"),  sheetName = outc[i], append = TRUE)
  
  # Leave-one-out-analysis
  loo <- matrix(0,7,3)
  loo[1,1:3] <- c(tab$Estimate,tab$SE,tab$Pval)
  
  source(here("Functions","leaveoneout.R"))
  leaveoneout(loo, a[[1]], ldrho, 'BinaryData_SA')
}



















