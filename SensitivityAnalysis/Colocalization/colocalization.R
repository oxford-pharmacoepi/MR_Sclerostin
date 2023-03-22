# COLOCALIZATION
rm(list = ls())
library(pacman)
pacman::p_load('coloc','ggplot2','dplyr','readr','xlsx')

genstart <- 41831099
genend   <- 41836156
window   <- 20

# EXPOSURE ---------------------------------------------------------------------
exposure_data <- read_delim(here("Metaanalysis","Fixed_results","Fixed.csv")) %>% # Results from the Meta-analysis
  select(-"...1") %>%
  select("SNP" = "rs_number","pval.exposure" = "p.value","beta.exposure" = "beta","se.exposure" = "se",
         "effect_allele.exposure" = "reference_allele","other_allele.exposure" = "other_allele", "n.exposure" = "n_samples") %>%
  mutate(CHR = 17,
         id.exposure = "exposure",
         exposure = "exposure") %>%
  inner_join(read_delim(here("Metaanalysis","Fixed_results","Fixed_metal.csv")) %>%
             select(SNP = MarkerName,
                    eaf.exposure = "Freq1"), by = "SNP") %>%
  inner_join(read_delim(here("Metaanalysis","sc.txt")) %>% select(SNP = MARKERNAME,Pos), by = "SNP") %>%
  filter(Pos <= genend+window*1000, Pos >= genstart-window*1000)

# OUTCOME ----------------------------------------------------------------------
# diabetes: 'ebi-a-GCST006867'
# hypertension: ukb-b-14057
# coronary artery disease: ebi-a-GCST005194
# myocardial infarction: ebi-a-GCST011365
# eBMD: ebi-a-GCST006979

out <- "ebi-a-GCST006979"
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = out)

# Hip fracture
# out <- "HF"
# outcome_data <- read_delim(here("Data","GWAS_HipFracture","GCST90161240_buildGRCh37.tsv")) %>%
#   rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
#          "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
#          "se.outcome" = "standard_error",'pos' = 'base_pair_location') %>%
#   filter(SNP %in% exposure_data$SNP) %>%
#   mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")

# Ischaemic stroke (uncomment lines from 33 to 40)
# out <- "IS"
# outcome_data <- read_delim(here("Data","GWAS_IschaemicStroke","Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
#   rename("Pos" = "base_pair_location") %>%
#   left_join(read_delim(here("SensitivityAnalysis","Bidirectional","sc.txt")) %>% select("Pos","SNP" = "MARKERNAME"), by = "Pos") %>%
#   select("SNP","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta","se.outcome" = "standard_error","pval.outcome" = "p_value",
#          "effect_allele.outcome" = "effect_allele","other_allele.outcome" = "other_allele", "chromosome") %>%
#   filter(SNP %in% exposure_data$SNP) %>%
#   mutate("id.outcome" = "IS","outcome" = "IS")
# write.csv(outcome_data,here("SensitivityAnalysis","Colocalization","outcomes",paste0("outcome_",out,".csv")))

outcome_data <- read_delim(here("SensitivityAnalysis","Colocalization","outcomes",paste0("outcome_",out,".csv"))) %>%
  filter(eaf.outcome > 0,
         eaf.outcome < 1) %>%
  mutate(id.outcome = "outcome",
          outcome    = "outcome")

dat <- harmonise_data(exposure_data, outcome_data) %>%
  filter(remove == FALSE) %>%
  distinct(SNP, .keep_all = TRUE)

write.xlsx(dat,here("SensitivityAnalysis","Colocalization","data.xlsx"),sheetName = out,append = TRUE)

# Continuous variable (quantitative)
res <- coloc.abf(dataset1 = list(snp = dat$SNP, beta = dat$beta.outcome, varbeta = dat$se.outcome^2, pvalues = dat$pval.outcome, MAF = dat$eaf.outcome, N = dat$samplesize.outcome, type = "quant"),
                 dataset2 = list(snp = dat$SNP, beta = dat$beta.exposure, varbeta = dat$se.exposure^2, pvalues = dat$pval.exposure, MAF = dat$eaf.exposure, N = dat$n.exposure, type = "quant"),
                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

p1 <- ggplot(dat, aes(x = pos, y = -log10(pval.exposure))) +
  geom_point() +
  labs(x = "", y = bquote(-log[10](italic(p))),
       title = "Sclerostin") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(74500000, 74600000, 74700000, 74800000),
                     labels = c(74500, 74600, 74700, 74800))

p2 <- ggplot(dat, aes(x = pos, y = -log10(pval.outcome))) +
  geom_point() +
  labs(x = "", y = bquote(-log[10](italic(p))),
       title = "eBMD") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(name = paste0("Chromosome ", 17, " Position (kb)"),
                     breaks = c(74500000, 74600000, 74700000, 74800000),
                     labels = c(74500, 74600, 74700, 74800))
ggpubr::ggarrange(p1, p2,
                  heights = c(1, 1), nrow = 2, 
                  ncol = 1, align = "hv")


# Binary variable
res <- coloc.abf(dataset1 = list(snp = dat$SNP, beta = dat$beta.exposure, varbeta = dat$se.exposure^2, pvalues = dat$pval.exposure, MAF = dat$eaf.exposure, N = dat$n.exposure, type = "quant"),
                 dataset2 = list(snp = dat$SNP, beta = dat$beta.outcome,  varbeta = dat$se.outcome^2,  pvalues = dat$pval.outcome,  MAF = dat$eaf.outcome, N = dat$samplesize.outcome, s = 0.3, type = "cc"),
                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
res










