# ============================================================================ #
#                               Exposure data for                              #
#                      non correlated SNPs sensitivy analysis                  #
#                        2023 - Marta Alcalde-Herraiz                          #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

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


# Outcome
out <- c('eBMD','HF')
for (i in 1:2){
  switch(i,
         "1" = {
           outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes = "ebi-a-GCST006979")
         },
         "2" = {
           outcome_dat <- read_delim(paste0(pathData,"GWAS_HipFracture\\GCST90161240_buildGRCh37.tsv")) %>%
             rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
                    "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
                    "se.outcome" = "standard_error") %>%
             filter(SNP %in% exposure_dat$SNP) %>%
             mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")
         }
  )
  
  # Harmonize data
  dat <- harmonise_data(exposure_dat,outcome_dat) %>%
    filter(remove == FALSE) %>%
    distinct()
  
  res <- TwoSampleMR::mr(dat)
}


exposure_dat$beta.exposure <- -exposure_dat$beta.exposure
write.csv(exposure_dat,here("SensitivityAnalysis","SingleSNP","exposure_data.csv"))
