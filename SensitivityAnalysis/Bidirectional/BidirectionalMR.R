# ============================================================================ # 
#                             Bidirectional MR                                 #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx','LDlinkR')

# Exposure data ----------------------------------------------------------------
# Copy the id of the outcome of interest into the "out" variable:
# - diabetes: ebi-a-GCST006867
# - hypertension: ukb-b-14057
# - coronary artery disease: ebi-a-GCST005194
# - myocardial infarction: ebi-a-GCST011365
# - eBMD: ebi-a-GCST006979

out <- "ebi-a-GCST006867"
exposure_data <- extract_instruments(outcomes = out) %>%
  filter(chr.exposure == 17,
         pval.exposure <= 1e-6)

# Hip fracture
# out <- "HF"
# exposure_data <- read_delim(here("Data","GWAS_HipFracture","GCST90161240_buildGRCh37.tsv")) %>%
#   rename("SNP" = "variant_id","pval.exposure" = "p_value","effect_allele.exposure" = "effect_allele",
#          "other_allele.exposure" = "other_allele","eaf.exposure" = "effect_allele_frequency","beta.exposure" = "beta",
#          "se.exposure" = "standard_error",'pos' = 'base_pair_location') %>%
#   filter(chromosome == 17, pval.exposure <= 1e-6) %>%
#   mutate(exposure = "Hip fracture", id.exposure = "Hip fracture")

# Ischaemic stroke
# out <- "IS"
# exposure_data <- read_delim(here("Data","GWAS_IschaemicStroke","Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
#   rename("Pos" = "base_pair_location") %>%
#   left_join(read_delim(here("SensitivityAnalysis","Bidirectional","sc.txt")) %>% select("Pos","SNP" = "MARKERNAME"), by = "Pos") %>%
#   select("SNP","eaf.exposure" = "effect_allele_frequency","beta.exposure" = "beta","se.exposure" = "standard_error","pval.exposure" = "p_value",
#          "effect_allele.exposure" = "effect_allele","other_allele.exposure" = "other_allele", "chromosome") %>%
#   filter(chromosome == 17, pval.exposure <= 1e-6) %>%
#   mutate("id.exposure" = "IS","exposure" = "IS")


exposure_data <- clump_data(exposure_data,clump_r2 = 0.8, clump_kb = 500, pop = "EUR")
ldrho <- LDmatrix(exposure_data$SNP,"EUR",token = "93873e604a3e")

# Outcome data -----------------------------------------------------------------
# Read meta-analysis results
outcome_data <- read_delim(here("SensitivityAnalysis","Bidirectional","Fixed_results","Fixed.csv"), delim = ",", show_col_types = FALSE) %>%
  select(-"...1") %>%
  filter(rs_number %in% exposure_data$SNP) %>% # Select those SNPs significantly correlated with sclerostin
  rename("SNP" = "rs_number",
         "effect_allele.outcome" = "reference_allele",
         "other_allele.outcome" = "other_allele",
         "beta.outcome" = "beta",
         "se.outcome" = "se",
         "pval.outcome" = "p.value",
         "N.outcome" = "n_samples") %>%
  inner_join(read_delim(here("SensitivityAnalysis","Bidirectional","Fixed_results","Fixed_metal.csv"), delim = ",", show_col_types = FALSE) %>%
               select("SNP" = "MarkerName",
                      "eaf.outcome" = "Freq1"),
             by = "SNP") %>%
  mutate(id.outcome = "Sclerostin",
         outcome = "Sclerostin")

# MR ANALYSIS ------------------------------------------------------------------
dat <- ldrho %>% select("SNP" = "RS_number") %>%
  left_join(harmonise_data(exposure_data,outcome_data),
            by = "SNP") %>%
  filter(remove == FALSE) %>%
  distinct() 

corr <- dat %>%
  select("SNP") %>%
  left_join(ldrho %>%
              rename(SNP = RS_number),
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

ivw   <- mr_ivw(mr,correl = TRUE)
egger <- mr_egger(mr,correl = TRUE)

table <- data.frame("Method" = c("ivw",        "egger",           "eggerIntercept"),
                    "Estimate" = c(ivw$Estimate,  egger$Estimate,    egger$Intercept),
                    "SE" = c(ivw$StdError,  egger$StdError.Est,egger$StdError.Int),
                    "Pval" = c(ivw$Pvalue,    egger$Causal.pval,  egger$Pleio.pval),
                    "SNPs" = c(ivw$SNPs,      egger$SNPs,        egger$SNPs)) %>%
  mutate(beta = Estimate, C1_Beta = Estimate - 1.96*SE, C2_Beta = Estimate + 1.96*SE) %>%
  mutate(OR = exp(Estimate), C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))

write.xlsx(table,here("SensitivityAnalysis","Bidirectional","gMR_res.xlsx"),sheetName = out,append = TRUE)
write.xlsx(dat,here("SensitivityAnalysis","Bidirectional","gMR_dat.xlsx"),  sheetName = out, append = TRUE)

