# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                             Logistic regression                              #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","ukb")))
library(pacman)
pacman::p_load('readr','trio','tibble','dplyr','readxl','lubridate',
               'MendelianRandomization','TwoSampleMR','xlsx')

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv")) %>%
  select(-X)
snps <- data.frame("SNP" = ldrho$RS_number) # SNPs

# Outcome data -----------------------------------------------------------------
outc <- "T2DM"
type <- "Logistic"

t    <- read_delim(here("MR_UKBiobank","BinaryData",type,paste0("Phenotype_",outc,".csv"))) %>%
  select("eid",
         "outcome" = "state") %>%
  inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
               select("eid",
                      "Caucasian" = "22006-0.0"), # Caucasian
             by = "eid") %>%
  filter(!is.na(outcome)) %>%
  filter(Caucasian == 1) %>%
  inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
               select(-"...1"),
             by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS
t
sum(t$outcome == 1)
sum(t$outcome == 1)/nrow(t)*100
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
    filter(!is.na(outcome),
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

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1") %>% 
  right_join(data.frame("SNP" = snps), by = "SNP")

# OUTCOME DATA -----------------------------------------------------------------
outcome_dat <- data.frame("SNP" = snps, "beta.outcome" = beta, "se.outcome" = se,
                          "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
  mutate(id.outcome = outc,
         outcome    = outc,
         effect_allele.outcome = exposure_dat$effect_allele.exposure,
         other_allele.outcome = exposure_dat$other_allele.exposure)

# HARMONIZE DATA AND CORRELATION MATRIX ----------------------------------------
dat <- snps %>%
  left_join(harmonise_data(exposure_dat,outcome_dat), by = "SNP") %>%
  filter(remove == FALSE) %>%
  distinct()

corr <- dat %>%
  select("SNP") %>%
  left_join(ldrho %>%
              rename("SNP" = "RS_number"),
            by = "SNP") %>%
  select(dat$SNP) %>%
  data.matrix() 
rownames(corr) = colnames(corr)

mr <- mr_input(
  bx = dat$beta.exposure,
  bxse = dat$se.exposure,
  by = dat$beta.outcome,
  byse = dat$se.outcome,
  corr = corr,
  snps = dat$SNP
)

# GENERALIZED MENDELIAN RANDOMIZATION ------------------------------------------
ivw   <- mr_ivw(mr,correl = TRUE)
egger <- mr_egger(mr,correl = TRUE)

table <- data.frame(c("ivw",        "egger",           "eggerIntercept"),
                    c(ivw$Estimate,  egger$Estimate,    egger$Intercept),
                    c(ivw$StdError,  egger$StdError.Est,egger$StdError.Int),
                    c(ivw$Pvalue,    egger$Causal.pval,  egger$Pleio.pval),
                    c(ivw$SNPs,      egger$SNPs,        egger$SNPs))
names(table) <- c("Method","Estimate","SE","Pval","SNPs")
table <- table %>%
  mutate(OR = exp(Estimate), C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))

table

# write.xlsx(table,here("MR_UKBiobank","BinaryData","Logistic","MR_res.xlsx"),sheetName = outc, append = TRUE)
# write.xlsx(dat,here("MR_UKBiobank","BinaryData","Logistic","MR_dat.xlsx"),  sheetName = outc, append = TRUE)

# Leave-one-out-analysis
loo <- matrix(0,7,3)
loo[1,1:3] <- c(ivw$Estimate[1],ivw$StdError[1],ivw$Pvalue)

for (i in c(1:length(dat$SNP))){
  dat1  <- dat %>% filter(SNP != dat$SNP[i])
  corr1 <- corr[,-which(colnames(corr) == dat$SNP[i])]
  corr1 <- corr1[-which(rownames(corr) == dat$SNP[i]),]
  
  mr <- mr_input(bx = dat1$beta.exposure, bxse = dat1$se.exposure, by = dat1$beta.outcome,
                 byse = dat1$se.outcome, corr = corr1, snps = dat1$SNP
  )
  
  ivw <- mr_ivw(mr,correl = TRUE)
  loo[i+1,] <- c(ivw$Estimate[1],ivw$StdError[1],ivw$Pvalue[1])
}

write.xlsx(data.frame(SNPS = c("All",dat$SNP),
                      Estimate = loo[,1],
                      C1 = loo[,1]-1.96*loo[,2],
                      C2 = loo[,1]+1.96*loo[,2],
                      SE = loo[,2],
                      pval = loo[,3]),
           here("SensitivityAnalysis","LeaveOneOut","BinaryData_LR","loo.xlsx"),  sheetName = outc, append = TRUE)
