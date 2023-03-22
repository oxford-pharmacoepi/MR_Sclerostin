# ============================================================================ #
#                        MENDELIAN RANDOMIZATION                               #
# MR analysis using GWAS for the outcomes                                      #
# ============================================================================ #
rm(list = ls())
library(pacman)
pacman::p_load('TwoSampleMR','MendelianRandomization','tibble','dplyr','tidyverse',
               'readr','here','xlsx')

# Exposure data ----------------------------------------------------------------
exposure_data <- read.csv(here("Pruning","exposure_data.csv"))
ldrho         <- read.csv(here("Pruning","correlation.csv")) %>%
  select(-X)
snps <- data.frame("SNP" = ldrho$RS_number)

# Outcome data -----------------------------------------------------------------
# Copy the id of the outcome of interest into the "out" variable:
# - diabetes: ebi-a-GCST006867
# - hypertension: ukb-b-14057
# - coronary artery disease: ebi-a-GCST005194
# - myocardial infarction: ebi-a-GCST011365
# - eBMD: ebi-a-GCST006979

out <- "ebi-a-GCST005194"
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = out)

# Hip fracture (uncomment lines from 24 to 30)
# out <- "HF"
# outcome_data <- read_delim(here("Data","GWAS_HipFracture","GCST90161240_buildGRCh37.tsv")) %>%
#   rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
#          "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
#          "se.outcome" = "standard_error") %>%
#   filter(chromosome == 17, SNP %in% exposure_data$SNP) %>%
#   mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")

# Ischaemic stroke (uncomment lines from 33 to 40)
# out <- "IS"
# outcome_data <- read_delim(here("Data","GWAS_IschaemicStroke","Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
#   rename("Pos" = "base_pair_location") %>%
#   left_join(exposure_data %>%
#               select("SNP","Pos" = "pos"), by = "Pos") %>%
#   select("SNP","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta","se.outcome" = "standard_error","pval.outcome" = "p_value",
#          "effect_allele.outcome" = "effect_allele","other_allele.outcome" = "other_allele") %>%
#   filter(SNP %in% exposure_data$SNP) %>%
#   mutate(outcome = "IS", id.outcome = "IS")

# MR ANALYSIS ------------------------------------------------------------------
dat <- snps %>%
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
  mutate(OR = exp(Estimate), C1 = exp(Estimate - 1.96*SE), C2 = exp(Estimate + 1.96*SE))

write.xlsx(table,here("MR_metaanalysis","gMR_res.xlsx"),sheetName = out,append = TRUE)
write.xlsx(dat,here("MR_metaanalysis","gMR_dat.xlsx"),  sheetName = out, append = TRUE)

# Leave-one-out-analysis
loo <- matrix(0,7,3)
loo[1,1:3] <- c(ivw$Estimate[1],ivw$StdError[1],ivw$Pvalue[1])

for (i in c(1:length(dat$SNP))){
  dat1  <- dat %>% filter(SNP != dat$SNP[i])
  corr1 <- corr[-which(rownames(corr) == dat$SNP[i]),-which(colnames(corr) == dat$SNP[i])]

  mr <- mr_input(bx = dat1$beta.exposure, bxse = dat1$se.exposure, by = dat1$beta.outcome,
    byse = dat1$se.outcome, corr = corr1, snps = dat1$SNP
  )

  ivw <- mr_ivw(mr,correl = TRUE)
  loo[i+1,] <- c(ivw$Estimate[1],ivw$StdError[1],ivw$Pvalue[1])
}

write.xlsx(data.frame(SNPS = c("All",dat$SNP),
                      loo),
           here("SensitivityAnalysis","LeaveOneOut","Metaanalysis","loo.xlsx"),  sheetName = out, append = TRUE)

