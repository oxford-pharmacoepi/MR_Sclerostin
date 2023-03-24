# ============================================================================ #
#                                 PRUNING                                      #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This file selects the instruments for the Mendelian randomization. It        #
# performs clumping and the mendelian randomization with eBMD and hip fracture.#
# It generates the file "exposure_data" containing the information for the     #
# instruments.                                                                 #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# Read meta-analysis results
t <- read_delim(here("Metaanalysis","Fixed_results","Fixed.csv"), delim = ",", show_col_types = FALSE) %>%
  select(-"...1") %>%
  left_join(read.delim(here("Metaanalysis","sc.txt")) %>% select("Pos","rs_number" = "MARKERNAME"),
            by = "rs_number") %>%
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
exposure_dat <- clump_data(t,clump_r2 = 0.8, clump_kb = 500000, pop = "EUR")
  
# LD correlation matrix
# Sometimes there is no connection with the server, so the correlation matrix has
# has been stored in the following file:
ldrho <- LDmatrix(exposure_dat$SNP,"EUR",token = tok)
write.csv(ldrho, file=here("Pruning","correlation.csv"), row.names = FALSE)
ldrho <- read_delim(here("Pruning","correlation.csv"))
unlink(here('Pruning','gMR_res.xlsx'))

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
  
  source(here('Functions','gMR.R'))
  a <- gMR(exposure_data = exposure_dat, outcome_data = outcome_dat, correl_matrix = ldrho)
  
  write.xlsx(a[[1]],here("Pruning","gMR_res.xlsx"),sheetName = paste0("data_",out[i]),append = TRUE)
  write.xlsx(a[[2]],here("Pruning","gMR_res.xlsx"),sheetName = paste0("results_",out[i]),append = TRUE)
}

# Selection of the instruments
dat <- ldrho %>% select(SNP = RS_number) %>%
    left_join(exposure_dat, by = "SNP") %>%
  mutate(beta.exposure = -beta.exposure) %>% # Lowering sclerostin effect
  select("SNP", contains("exposure"))
write.csv(dat,here("Pruning","exposure_data.csv"))



