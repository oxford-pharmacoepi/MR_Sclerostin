# ============================================================================ #
#                       Generalized Mendelian Randomization                    #
#                               for binary data                                #
#                             Logistic regression                              #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','tibble','dplyr','tidyverse',
               'readr','here','xlsx')

# Exposure data ----------------------------------------------------------------
exposure_data <- read.csv(here("SensitivityAnalysis","SingleSNP","exposure_data.csv"))

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
               select("eid", exposure_data$SNP),
             by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS
t
sum(t$outcome == 1)
sum(t$outcome == 1)/nrow(t)*100
# LOGISTIC REGRESSION ----------------------------------------------------------
# Eliminate those rows with NaN
tab <- t %>%
  select("eid","sex","outcome","snp" = exposure_data$SNP,"PC1","PC2","PC3","PC4","PC5",
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
n    <- nrow(tab1)
eaf  <- (sum(tab1$snp)/(2*length(tab1$snp)))
beta <- logistic_model$coefficients[2,1]
se <- logistic_model$coefficients[2,2]
p <- logistic_model$coefficients[2,4]


# OUTCOME DATA -----------------------------------------------------------------
outcome_dat <- data.frame("SNP" = exposure_data$SNP, "beta.outcome" = beta, "se.outcome" = se,
                          "p.outcome" = p, "eaf.outcome" = eaf, "n.outcome" = n) %>%
  mutate(id.outcome = outc,
         outcome    = outc,
         effect_allele.outcome = exposure_data$effect_allele.exposure,
         other_allele.outcome = exposure_data$other_allele.exposure)

# HARMONIZE DATA AND CORRELATION MATRIX ----------------------------------------
dat <- exposure_data %>% select(SNP) %>%
  left_join(harmonise_data(exposure_data,outcome_dat), by = "SNP") %>%
  filter(remove == FALSE) %>%
  distinct()

res <- mr(dat)
res <- res %>% mutate(OR = exp(b), C1_L = exp(b-1.96*se), C2_U = exp(b+1.96*se))

# MENDELIAN RANDOMIZATION ------------------------------------------

write.xlsx(res,here("SensitivityAnalysis","SingleSNP","MR_res_binaryLR.xlsx"),sheetName = outc,append = TRUE)
write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_binaryLR.xlsx"),  sheetName = outc, append = TRUE)
