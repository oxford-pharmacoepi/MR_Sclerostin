rm(list = ls())

library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx','LDlinkR')

# SELECT INSTRUMENT ------------------------------------------------------------
# Read meta-analysis results
t <- read_delim(here("Metaanalysis","Fixed_results","Fixed.csv"), delim = ",", show_col_types = FALSE) %>%
  select(-"...1") %>%
  filter(p.value <= 1e-6) %>% # Select those SNPs significantly correlated with sclerostin
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
exposure_dat <- clump_data(t,clump_r2 = 0.1, clump_kb = 500, pop = "EUR")

# LD correlation matrix
# Sometimes there is no connection with the server, so the correlation matrix has
# has been stored in the following file:
# ldrho <- LDmatrix(exposure_dat$SNP,"EUR",token = "93873e604a3e")
# write_delim(ldrho,here("Pruning","ldrho.txt"))
ldrho <- read_delim(here("Pruning","ldrho.txt"))
snps  <- data.frame("SNP" = ldrho$RS_number)


# Outcome
out <- "ebi-a-GCST006979" # eBMD
outcome_dat <- extract_outcome_data(snps = snps$SNP, outcomes = out)

# Harmonize data
dat <- snps %>%
  left_join(harmonise_data(exposure_dat,outcome_dat), by = "SNP") %>%
  filter(remove == FALSE) %>%
  distinct()

dat <- dat %>%
  filter(sign(beta.exposure) != sign(beta.outcome)) %>% # Select those snps with inverse direction between exposure and outcome
  select("SNP", contains("exposure"),"pos") %>%
  mutate(beta.exposure = -beta.exposure) # Lowering sclerostin effect
write.csv(dat,here("SensitivityAnalysis","SingleSNP","exposure_data.csv"))
