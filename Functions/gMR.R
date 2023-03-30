gMR <- function(exposure_data,outcome_data,correl_matrix){
  snps <- correl_matrix %>% select(SNP = RS_number) 
  
  # Harmonize data
  dat <- snps %>%
    left_join(harmonise_data(exposure_data,outcome_data), by = "SNP") %>%
    filter(remove == FALSE) %>%
    distinct() 
  
  # Correlation matrix
  corr <- snps %>% 
    left_join(correl_matrix %>%
                rename("SNP" = "RS_number"), by = "SNP") %>%
    select(-SNP) %>%
    data.matrix() 
  rownames(corr) = colnames(corr)
  
  # Mendelian randomization
  mr <- mr_input(
    bx = dat$beta.exposure,
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome,
    corr = corr,
    snps = snps
  )
  
  ivw <- MendelianRandomization::mr_ivw(mr,correl = TRUE)
  
  tab <- data.frame("Method" = c("ivw"),
                    "Estimate" = -c(ivw$Estimate)) %>%
    mutate(CI_LOW  = Estimate-1.96*ivw$StdError, 
           CI_HIGH = Estimate+1.96*ivw$StdError,
           SE      = ivw$StdError,
           Pval    = ivw$Pvalue)
  
  return(list(dat,tab))
  
}