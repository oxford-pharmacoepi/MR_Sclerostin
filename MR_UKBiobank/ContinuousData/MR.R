# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for continuous data                            #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx','LDlinkR')

# Correlation matrix -----------------------------------------------------------
ldrho <- read.csv(here("Pruning","correlation.csv")) %>%
  select(-X)
snps <- data.frame("SNP" = ldrho$RS_number) # SNPs\

# Outcome data -----------------------------------------------------------------
# DICTIONARY OF THE OUTCOMES
  # "Cholesterol" = "30690-0.0", 
  # "LDL" = "30780-0.0", 
  # "HDL" = "30760-0.0",
  # "Triglycerides" = "30870-0.0",
  # "Apo-A" = "30630-0.0",
  # "Apo-B" = "30640-0.0",
  # "CRP" = "30710-0.0",
  # "Lipoprotein" = "30790-0.0",
  # "HbA1c" = "30750-0.0",
  # "Glucose" = "30740-0.0",
  # "eBMD" = "3148-0.0",
  # "Height" = "50-0.0"

outc <- "30760-0.0"
t    <-  read_delim(here("Data","UKB","ukb669864_1.csv")) %>% # Outcome
  select("eid",
         "outcome" = outc,
         "Caucasian" = "22006-0.0") %>%
  filter(!is.na(outcome)) %>%
  filter(Caucasian == 1) %>%
  mutate(outcome = (outcome - mean(outcome,na.rm = TRUE))/sd(outcome, na.rm = TRUE)) %>% # SD units
  inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
               select(-"...1"),
             by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid")  # 10 PRINCIPAL COMPONENTS


# LINEAR REGRESSION ------------------------------------------------------------
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
    filter(!is.na(outcome), !is.na(sex),
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

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- snps %>%
  left_join(read_delim(here("Pruning","exposure_data.csv"))) %>%
  select(-"...1") %>% 
  right_join(data.frame("SNP" = snps), by = "SNP")

# OUTCOME DATA -----------------------------------------------------------------
outcome_dat <- data.frame("SNP" = snps, 
                          "beta.outcome" = beta, 
                          "se.outcome" = se,
                          "p.outcome" = p, 
                          "eaf.outcome" = eaf, 
                          "n.outcome" = n) %>%
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


source(here("Functions","table.R"))
table <- table(ivw,egger)
table <- table %>% mutate(Beta = Estimate, C1 = Estimate-1.96*SE, C2 = Estimate+1.96*SE)

write.xlsx(table,here("MR_UKBiobank","ContinuousData","gMR_res.xlsx"),sheetName = outc,append = TRUE)
write.xlsx(dat,here("MR_UKBiobank","ContinuousData","gMR_dat.xlsx"),  sheetName = outc, append = TRUE)

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
           here("SensitivityAnalysis","LeaveOneOut","ContinuousData","loo.xlsx"),  sheetName = outc, append = TRUE)
