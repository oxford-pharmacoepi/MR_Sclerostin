# ============================================================================ #
#                                 PRUNING                                      #
# This file is to select the instruments for the Mendelian randomization. It   #
# performs clumping and the mendelian randomization with eBMD.                 #
# It generates the file "exposure_data" containing the information for the     #
# instruments.                                                                 #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx','LDlinkR')

# Read meta-analysis results
t <- read_delim(here("Metaanalysis","Fixed_results","Fixed.csv"), delim = ",", show_col_types = FALSE) %>%
  select(-"...1") %>%
  left_join(read.delim(here("Metaanalysis","sc.txt")) %>% select("Pos","rs_number" = "MARKERNAME"),
            by = "rs_number") %>%
  filter(Pos <= 41836156+500000, Pos >= 41831099-500000, # SNPs in the SOST region
         p.value <= 1e-6) %>% # Select those SNPs significantly correlated with sclerostin
  rename("SNP" = "rs_number",
         "effect_allele.exposure" = "reference_allele",
         "other_allele.exposure" = "other_allele",
         "beta.exposure" = "beta",
         "se.exposure" = "se",
         "pval.exposure" = "p.value",
         "N.exposure" = "n_samples") %>%
  inner_join(read_delim(here("Metaanalysis","Fixed_results","Fixed_metal.csv"), delim = ",", show_col_types = FALSE) %>%
               select("SNP" = "MarkerName",
                      "eaf.exposure" = "Freq1"),
             by = "SNP") %>%
  mutate(id.exposure = "Sclerostin",
         exposure = "Sclerostin")

# Pruning
exposure_dat <- clump_data(t,clump_r2 = 0.8, clump_kb = 500, pop = "EUR")
  
# LD correlation matrix
# Sometimes there is no connection with the server, so the correlation matrix has
# has been stored in the following file:
ldrho <- LDmatrix(exposure_dat$SNP,"EUR",token = "93873e604a3e")
write_delim(ldrho,here("Pruning","ldrho.txt"))
ldrho <- read_delim(here("Pruning","ldrho.txt"))
snps  <- data.frame("SNP" = ldrho$RS_number)

# Outcome
# out <- "ebi-a-GCST006979" # eBMD
# outcome_dat <- extract_outcome_data(snps = snps$SNP, outcomes = out)
out <- "HF"
outcome_dat <- read_delim(here("Data","GWAS_HipFracture","GCST90161240_buildGRCh37.tsv")) %>%
  rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
         "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
         "se.outcome" = "standard_error") %>%
  filter(chromosome == 17, SNP %in% exposure_dat$SNP) %>%
  mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")

# Harmonize data
dat <- snps %>%
  left_join(harmonise_data(exposure_dat,outcome_dat), by = "SNP") %>%
  filter(remove == FALSE) %>%
  distinct() 

# Correlation matrix
corr <- dat %>%
  select("SNP") %>%
  left_join(ldrho %>%
              rename("SNP" = "RS_number"), by = "SNP") %>%
  select(dat$SNP) %>%
  data.matrix() 
rownames(corr) = colnames(corr)

# Mendelian randomization
mr <- mr_input(
  bx = dat$beta.exposure,
  bxse = dat$se.exposure,
  by = dat$beta.outcome,
  byse = dat$se.outcome,
  corr = corr,
  snps = dat$SNP
)

ivw   <- MendelianRandomization::mr_ivw(mr,correl = TRUE)
egger <- MendelianRandomization::mr_egger(mr,correl = TRUE)

table <- data.frame("Method" = c("ivw",        "egger",           "eggerIntercept"),
                    "Estimate" = c(ivw$Estimate,  egger$Estimate,    egger$Intercept),
                    "SE" = c(ivw$StdError,  egger$StdError.Est,egger$StdError.Int),
                    "Pval" = c(ivw$Pvalue,    egger$Causal.pval,  egger$Pleio.pval)) %>%
  mutate(OR = exp(Estimate), 
         C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))

write.xlsx(table,here("Pruning","gMR_res.xlsx"),sheetName = paste0("results_",out),append = TRUE)
write.xlsx(dat,here("Pruning","gMR_res.xlsx"),sheetName = paste0("data_",out),append = TRUE)

# Selection of the instruments
dat <- dat %>%
  mutate(beta.exposure = -beta.exposure) %>% # Lowering sclerostin effect
  select("SNP", contains("exposure"),"pos")
write.csv(dat,here("Pruning","exposure_data.csv"))

ldrho <- ldrho %>%
  filter(RS_number %in% dat$SNP) %>%
  select(RS_number,contains(dat$SNP))
write.csv(ldrho, file=here("Pruning","correlation.csv"))

