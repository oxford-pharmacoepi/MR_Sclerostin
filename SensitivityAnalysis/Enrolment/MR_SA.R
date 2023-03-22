# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                               Survival analysis                              #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('readr','trio','tibble','dplyr','readxl','lubridate',
               'MendelianRandomization','TwoSampleMR','xlsx','survival')

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv")) %>%
  select(-X)
snps <- data.frame("SNP" = ldrho$RS_number) # SNPs

# Outcome data -----------------------------------------------------------------
outc <- "Fracture"
type <- "SA_Enrolment" # /SA_Enrolment"

t    <- read_delim(here("MR_UKBiobank","BinaryData",type,paste0("Phenotype_",outc,".csv"))) %>%
  select("eid",
         "state",
         "age") %>%
  inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
               select("eid",
                      "Caucasian" = "22006-0.0"), # Caucasian
             by = "eid") %>%
  filter(Caucasian == 1) %>%
  inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
               select(-"...1"),
             by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS
325713-nrow(t)
nrow(t)

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
    filter(!is.na(state), !is.na(age),
           !is.na(snp),!is.na(PC1),!is.na(PC2),!is.na(PC3),!is.na(PC4),!is.na(PC5),
           !is.na(PC6),!is.na(PC7),!is.na(PC8),!is.na(PC9),!is.na(PC10))
  
  # Substitute the 1 for 2
  tab1 <- tab
  tab1$snp[tab$snp == 2] <- 1
  
  #Logistic regression - Adjusted model for sex and the first 10 principal components
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


write.xlsx(table,here("SensitivityAnalysis","Enrolment",type,"MR_res.xlsx"),sheetName = outc, append = TRUE)
write.xlsx(dat,here("SensitivityAnalysis","Enrolment",type,"MR_dat.xlsx"),  sheetName = outc, append = TRUE)





















# 
# tab1 <- as_tibble(read.csv(here("Figures","MI.csv")))
# 
# fit <- survfit(Surv(age,state)  ~ pr, data = tab1)
# eth.death.surv <- ggsurvplot(fit, linetype = "strata", # Change line type by groups
#                              # xlim = c(0,413),
#                               ylim = c(0.88,1),
#                              # pval = FALSE, conf.int = T,
#                              # risk.table = TRUE, # Add risk table
#                              # risk.table.col = "strata", # Change risk table color by groups
#                              
#                              # ggtheme = theme_bw(), # Change ggplot2 theme)
#                              legend.labs = c("Lower sclerostin","Normal sclerostin","High sclerostin")
# )
# eth.death.surv





