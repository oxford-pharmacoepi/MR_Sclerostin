table <- function(ivw,egger){
  table <- data.frame(c("ivw",        "egger",           "eggerIntercept"),
                      c(ivw$Estimate,  egger$Estimate,    egger$Intercept),
                      c(ivw$StdError,  egger$StdError.Est,egger$StdError.Int),
                      c(ivw$Pvalue,    egger$Causal.pval,  egger$Pleio.pval),
                      c(ivw$SNPs,      egger$SNPs,        egger$SNPs))
  names(table) <- c("Method","Estimate","SE","Pval","SNPs")
  table <- table %>%
    mutate(OR = exp(Estimate), C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))
}