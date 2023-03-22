# ============================================================================ #
#                        MENDELIAN RANDOMIZATION                               #
# MR analysis using GWAS for the outcomes                                      #
# ============================================================================ #


rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','tibble','dplyr','tidyverse',
               'readr','here','xlsx')

# Exposure data ----------------------------------------------------------------
exposure_data <- read.csv(here("SensitivityAnalysis","SingleSNP","exposure_data.csv"))


out <- "ebi-a-GCST006867"
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = out)

# Outcome data -----------------------------------------------------------------
# Copy the id of the outcome of interest into the "out" variable:
# - diabetes: ebi-a-GCST006867
# - hypertension: ukb-b-14057
# - coronary artery disease: ebi-a-GCST005194
# - myocardial infarction: ebi-a-GCST011365
# - eBMD: ebi-a-GCST006979



# Hip fracture (uncomment lines from 24 to 30)
# out <- "HF"
# outcome_data <- read_delim(here("Data","GWAS_HipFracture","GCST90161240_buildGRCh37.tsv")) %>%
#   rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
#          "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
#          "se.outcome" = "standard_error") %>%
#   filter(chromosome == 17, SNP %in% exposure_data$SNP) %>%
#   mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")

# Ischaemic stroke (uncomment lines from 33 to 40)
out <- "IS"
outcome_data <- read_delim(here("Data","GWAS_IschaemicStroke","Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
  rename("Pos" = "base_pair_location") %>%
  left_join(exposure_data %>%
              select("SNP","Pos" = "pos"), by = "Pos") %>%
  select("SNP","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta","se.outcome" = "standard_error","pval.outcome" = "p_value",
         "effect_allele.outcome" = "effect_allele","other_allele.outcome" = "other_allele") %>%
  filter(SNP %in% exposure_data$SNP) %>%
  mutate(outcome = "IS", id.outcome = "IS")

# MR ANALYSIS ------------------------------------------------------------------
dat <-harmonise_data(exposure_data,outcome_data) %>%
  filter(remove == FALSE) %>%
  distinct() 

res <- mr(dat)

table <- data.frame("Method" = c("Wald ratio"),
                    "Estimate" = c(res$b),
                    "SE" = c(res$se),
                    "Pval" = c(res$pval),
                    "SNPs" = c(res$nsnp)) %>%
  mutate(C1_beta = Estimate-1.96*SE, C2_beta = Estimate + 1.96*SE) %>%
  mutate(OR = exp(Estimate), C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))

write.xlsx(table,here("SensitivityAnalysis","SingleSNP","MR_res_metaanalysis.xlsx"),sheetName = out,append = TRUE)
write.xlsx(dat,here("SensitivityAnalysis","SingleSNP","MR_dat_metaanalysis.xlsx"),  sheetName = out, append = TRUE)
