leaveoneout <- function(loo, dat, ldrho,type){
  ldrho <- ldrho %>% 
    select(-RS_number) %>%
    data.matrix() 
  rownames(ldrho) <- colnames(ldrho)
  
  for (j in c(1:length(dat$SNP))){
    dat1  <- dat %>% filter(SNP != dat$SNP[j])
    corr1 <- ldrho[-which(rownames(ldrho) == dat$SNP[j]),-which(colnames(ldrho) == dat$SNP[j])]
    
    mr <- mr_input(bx = dat1$beta.exposure, bxse = dat1$se.exposure, by = dat1$beta.outcome,
                   byse = dat1$se.outcome, corr = corr1, snps = dat1$SNP
                   )
    
    ivw <- mr_ivw(mr,correl = TRUE)
    loo[j+1,] <- c(-ivw$Estimate[1],ivw$StdError[1],ivw$Pvalue[1])
  }
  
  write.xlsx(data.frame(SNPS = c("All",dat$SNP),
                        Estimate = loo[,1],
                        CI_LOW   = loo[,1]-1.96*loo[,2],
                        CI_HIGH = loo[,1]+1.96*loo[,2],
                        SE   = loo[,2],
                        pval = loo[,3]),
             here("SensitivityAnalysis","LeaveOneOut",type,"loo.xlsx"),  sheetName = dat$outcome[1], append = TRUE)
}