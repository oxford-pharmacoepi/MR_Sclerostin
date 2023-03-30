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
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok")))

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv"))
snps  <- ldrho %>% select("SNP" = RS_number) # SNPs

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1") %>% 
  right_join(data.frame("SNP" = snps), by = "SNP")

# Outcome data -----------------------------------------------------------------
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

unlink(here("MR_UKBiobank","BinaryData","Logistic","MR_res.xlsx"))
unlink(here("MR_UKBiobank","BinaryData","Logistic","MR_dat.xlsx"))
for (i in 1:length(outc)){
  t <- read_delim(paste0(pathData,"MR_UKBiobank\\BinaryData\\Logistic\\Phenotype_",outc[i],".csv")) %>%
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
  beta  <- matrix(0,nrow(snps),1)
  se    <- matrix(0,nrow(snps),1)
  p     <- matrix(0,nrow(snps),1)
  eaf   <- matrix(0,nrow(snps),1)
  n     <- matrix(0,nrow(snps),1)
  
  for (j in(c(1:nrow(snps)))){
    # Eliminate those rows with NaN
    tab <- t %>%
      select("eid","sex","outcome", "snp" = snps[j,1],"PC1","PC2","PC3","PC4","PC5",
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
    n[j]    <- nrow(tab1)
    eaf[j]  <- (sum(tab1$snp)/(2*length(tab1$snp)))
    beta[j] <- logistic_model$coefficients[2,1]
    se[j] <- logistic_model$coefficients[2,2]
    p[j]  <- logistic_model$coefficients[2,4]
  }
  
  # OUTCOME DATA -----------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = snps, "beta.outcome" = beta, "se.outcome" = se,
                            "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
    mutate(id.outcome = outc[i],
           outcome    = outc[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  # MENDELIAN RANDOMIZATION
  source(here("Functions","gMR.R"))
  a <- gMR(exposure_data = exposure_dat,outcome_data = outcome_dat,ldrho)
  tab <- a[[2]] %>% mutate(OR = exp(Estimate),
                           CI_LOW_OR = exp(Estimate - 1.96*SE),
                           CI_HIGH_OR = exp(Estimate + 1.96*SE))
  
  write.xlsx(tab,here("MR_UKBiobank","BinaryData","Logistic","MR_res.xlsx"), sheetName = outc[i], append = TRUE)
  write.xlsx(a[[1]],here("MR_UKBiobank","BinaryData","Logistic","MR_dat.xlsx"),  sheetName = outc[i], append = TRUE)
  
  # Leave-one-out analysis
  loo <- matrix(0,7,3)
  loo[1,1:3] <- c(tab$Estimate,tab$SE,tab$Pval)
  source(here("Functions","leaveoneout.R"))
  loo <- leaveoneout(loo, a[[1]], ldrho,'BinaryData_LR')
}







