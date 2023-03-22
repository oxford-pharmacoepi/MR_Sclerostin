# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for continuous data                            #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx','LDlinkR')

# EXPOSURE DATA ----------------------------------------------------------------
exposure_dat <- read_delim(here("SensitivityAnalysis","SingleSNP","exposure_data.csv")) %>%
  select(-"...1")

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

outc <- "30740-0.0"
t    <-  read_delim(here("Data","UKB","ukb669864_1.csv")) %>% # Outcome
  select("eid",
         "outcome" = outc,
         "Caucasian" = "22006-0.0") %>%
  filter(!is.na(outcome)) %>%
  filter(Caucasian == 1) %>%
  mutate(outcome = (outcome - mean(outcome,na.rm = TRUE))/sd(outcome, na.rm = TRUE)) %>% # SD units
  inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
               select("eid", exposure_dat$SNP),
             by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid")  # 10 PRINCIPAL COMPONENTS


# LINEAR REGRESSION ------------------------------------------------------------
# Eliminate those rows with NaN
tab <- t %>%
  select("eid","sex","outcome", "snp" = exposure_dat$SNP,"PC1","PC2","PC3","PC4","PC5",
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
beta  <- coefficients(mlr)[2]
se    <- summary(mlr)$coefficients[2,2]
p     <- summary(mlr)$coefficients[2,4]
eaf   <- (sum(tab1$snp)/(2*length(tab1$snp)))
n     <- length(tab1$snp)


# OUTCOME DATA -----------------------------------------------------------------
outcome_dat <- data.frame("SNP" = exposure_dat$SNP, 
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
dat <- harmonise_data(exposure_dat,outcome_dat) %>%
  filter(remove == FALSE) %>%
  distinct()

res <- mr(dat)
res <- res %>% mutate(Beta = b, C1_LOW =b-1.96*se, C2_LOW = b+1.96*b)

write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_continuous.xlsx"),sheetName = outc,append = TRUE)
write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_continuous.xlsx"),  sheetName = outc, append = TRUE)
