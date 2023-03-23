# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for continuous data                            #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Generalized MR for continuous data using UK Biobank patient-level data.      #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv")) 
snps  <- ldrho %>% select("SNP" = RS_number) 

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1") %>% 
  right_join(data.frame("SNP" = snps), by = "SNP")

# Outcome data -----------------------------------------------------------------
out  <- c('eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A','Apo-B','CRP','Lipoprotein','HbA1c','Glucose')
outc <- c('3148-0.0','30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0','30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0')

unlink(here("MR_UKBiobank","ContinuousData","gMR_res.xlsx"))
unlink(here("MR_UKBiobank","ContinuousData","gMR_dat.xlsx"))

for (i in 1:11){
  t    <-  read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>% # Outcome
    select("eid",
           "outcome" = outc[i],
           "Caucasian" = "22006-0.0") %>%
    filter(!is.na(outcome)) %>%
    filter(Caucasian == 1) %>%
    mutate(outcome = (outcome - mean(outcome,na.rm = TRUE))/sd(outcome, na.rm = TRUE)) %>% # SD units
    inner_join(read_delim(paste0(pathData,"genetics.csv")) %>% # SNPs expression
                 select(-"...1"),
               by = "eid") %>%
    inner_join(read_delim(paste0(pathData,"UKB\\PC.csv")), by = "eid")  # 10 PRINCIPAL COMPONENTS
  
  # LINEAR REGRESSION ----------------------------------------------------------
  beta  <- matrix(0,nrow(snps),1)
  se    <- matrix(0,nrow(snps),1)
  p     <- matrix(0,nrow(snps),1)
  eaf   <- matrix(0,nrow(snps),1)
  n     <- matrix(0,nrow(snps),1)
  
  for (j in 1:nrow(snps)){
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
    
    # Linear regression - Adjusted model for sex and the first 10 principal components
    mlr <- lm(outcome ~ snp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = tab1)
    
    # Coefficients:
    beta[j]  <- summary(mlr)$coefficients[2,1]
    se [j]   <- summary(mlr)$coefficients[2,2]
    p[j]     <- summary(mlr)$coefficients[2,4]
    eaf[j]   <- (sum(tab1$snp)/(2*length(tab1$snp)))
    n[j]     <- length(tab1$snp)
  }
  
  # OUTCOME DATA -----------------------------------------------------------------
  outcome_dat <- data.frame("SNP" = snps, 
                            "beta.outcome" = beta, 
                            "se.outcome" = se,
                            "p.outcome" = p, 
                            "eaf.outcome" = eaf, 
                            "n.outcome" = n) %>%
    mutate(id.outcome = out[i],
           outcome    = out[i],
           effect_allele.outcome = exposure_dat$effect_allele.exposure,
           other_allele.outcome = exposure_dat$other_allele.exposure)
  
  source(here("Functions","gMR.R"))
  a <- gMR(exposure_data = exposure_dat,outcome_data = outcome_dat,ldrho)
  
  write.xlsx(a[[2]],here("MR_UKBiobank","ContinuousData","gMR_res.xlsx"), sheetName = out[i],append = TRUE)
  write.xlsx(a[[1]],here("MR_UKBiobank","ContinuousData","gMR_dat.xlsx"), sheetName = out[i], append = TRUE)
  
  # Leave-one-out-analysis
  loo <- matrix(0,7,3)
  loo[1,1:3] <- c(a[[2]]$Estimate[1],a[[2]]$SE[1],a[[2]]$Pval)
  source(here("Functions","leaveoneout.R"))
  loo <- leaveoneout(loo, a[[1]], ldrho,'ContinuousData')
}




