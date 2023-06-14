# ============================================================================ #
#                           MENDELIAN RANDOMIZATION                            #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# MRanalysis using GWAS for the outcomes and sensitivity analysis              #
# "leave-one-out"                                                              #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# Exposure data ----------------------------------------------------------------
exposure_data <- read.csv(here("SensitivityAnalysis","Heterogeneity","exposure_data.csv"))
ldrho         <- read.csv(here("SensitivityAnalysis","Heterogeneity","correlation.csv")) 
snps <- ldrho %>% select(SNP = RS_number)

# Outcome data -----------------------------------------------------------------
out  <- c('eBMD','HF','CAD','MI','IS','Hypertension','T2DM')
outc <- c('ebi-a-GCST006979','','ebi-a-GCST005194','ebi-a-GCST011365','','ukb-b-14057','')

unlink(here('MR_gwas','gMR_res.xlsx'))
unlink(here('MR_gwas','gMR_dat.xlsx'))

for (i in 1:length(out)){
  if(i %in% c(1,3,4,6)){
    outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = outc[i]) %>%
      select(SNP,pos,beta.outcome, se.outcome, samplesize.outcome, pval.outcome, eaf.outcome, effect_allele.outcome,
             other_allele.outcome) %>%
      mutate(outcome = out[i], id.outcome = out[i])
  }else if(i == 2){
    outcome_data <- read_delim(paste0(pathData,"GWAS_HipFracture\\GCST90161240_buildGRCh37.tsv")) %>%
      rename("SNP" = "variant_id","pval.outcome" = "p_value","effect_allele.outcome" = "effect_allele",
             "other_allele.outcome" = "other_allele","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta",
             "se.outcome" = "standard_error") %>%
      filter(SNP %in% exposure_data$SNP) %>%
      mutate(outcome = "Hip fracture", id.outcome = "Hip fracture")
  }else if(i == 5){
    outcome_data <- read_delim(paste0(pathData,"GWAS_IschaemicStroke\\Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
      rename("Pos" = "base_pair_location") %>%
      left_join(
        exposure_data %>% select("SNP") %>%
                  left_join(read_delim(here("Metaanalysis","sc.txt")) %>% 
                              select("SNP" = "MARKERNAME","Pos"), by = "SNP"),
        by = "Pos") %>%
      select("SNP","eaf.outcome" = "effect_allele_frequency","beta.outcome" = "beta","se.outcome" = "standard_error","pval.outcome" = "p_value",
             "effect_allele.outcome" = "effect_allele","other_allele.outcome" = "other_allele") %>%
      filter(SNP %in% exposure_data$SNP) %>%
      mutate(outcome = "IS", id.outcome = "IS")
  }else if(i == 7){
    outcome_data <- read_delim(paste0(pathData,"GWAS_T2DM\\DIAMANTE-EUR.sumstat.txt")) %>%
      select('SNP' = 'rsID','beta.outcome' = 'Fixed-effects_beta','effect_allele.outcome' = 'effect_allele',
             'other_allele.outcome' = 'other_allele','eaf.outcome' = 'effect_allele_frequency',
             'se.outcome' = 'Fixed-effects_SE','pval.outcome' = 'Fixed-effects_p-value') %>%
      right_join(
        exposure_data %>% select('SNP'),
        by = 'SNP'
      ) %>%
      mutate(outcome = "T2DM", id.outcome = "T2DM")
  }
  
  source(here('Functions','gMR.R'))
  
  loo <- matrix(0,7,3) # leave-one-out analysis
  
  a <- gMR(exposure_data = exposure_data, outcome_data = outcome_data, correl_matrix = ldrho)

  if(i %in% c(2:7)){
    tab <- a[[2]] %>%
      mutate(OR = exp(Estimate), CI_LOW_OR = exp(CI_LOW), CI_HIGH_OR = exp(CI_HIGH))
  }else{
    tab <- a[[2]]
  }
  
  write.xlsx(a[[1]],here("SensitivityAnalysis","Heterogeneity","gMR_dat.xlsx"),sheetName = out[i],append = TRUE)
  write.xlsx(tab,here("SensitivityAnalysis","Heterogeneity","gMR_res.xlsx"),sheetName = out[i],append = TRUE)
}






