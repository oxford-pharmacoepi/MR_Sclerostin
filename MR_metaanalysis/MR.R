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
exposure_data <- read.csv(here("Pruning","exposure_data.csv"))
ldrho         <- read.csv(here("Pruning","correlation.csv")) 
snps <- ldrho %>% select(SNP = RS_number)

# Outcome data -----------------------------------------------------------------
out  <- c('eBMD','HF','CAD','MI','IS','Hypertension','T2DM')
outc <- c('ebi-a-GCST006979','','ebi-a-GCST005194','ebi-a-GCST011365','','ukb-b-14057','ebi-a-GCST006867')
unlink(here('MR_metaanalysis','gMR_res.xlsx'))
unlink(here('MR_metaanalysis','gMR_dat.xlsx'))
unlink(here("SensitivityAnalysis","LeaveOneOut","Metaanalysis","loo.xlsx"))

for (i in 1:length(out)){
  if(i %in% c(1,3,4,6,7)){
    outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = outc[i])
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
  }
  
  source(here('Functions','gMR.R'))
  if (i == length(out)){
    outcome_data <- snps %>% inner_join(outcome_data, by = "SNP")
    ldrho        <- ldrho %>% select(-"rs66838809") %>% filter(RS_number != "rs66838809")
    snps         <- ldrho %>% select(SNP = RS_number)
    loo <- matrix(0,6,3) # leave-one-out analysis
  }else{
    loo <- matrix(0,7,3) # leave-one-out analysis
  }
  
  a <- gMR(exposure_data = exposure_data, outcome_data = outcome_data, correl_matrix = ldrho)
  
  if(i %in% c(2:7)){
    tab <- a[[2]] %>%
      mutate(OR = exp(Estimate), CI_LOW_OR = exp(Estimate - 1.96*SE), CI_HIGH_OR = exp(Estimate + 1.96*SE))
  }else{
    tab <- a[[2]]
  }
  
  write.xlsx(a[[1]],here("MR_metaanalysis","gMR_dat.xlsx"),sheetName = out[i],append = TRUE)
  write.xlsx(tab,here("MR_metaanalysis","gMR_res.xlsx"),sheetName = out[i],append = TRUE)
  
  # Leave-one-out-analysis
  loo[1,1:3] <- c(a[[2]]$Estimate,a[[2]]$SE,a[[2]]$Pval)
  
  source(here('Functions','leaveoneout.R'))
  loo <- leaveoneout(loo, a[[1]], ldrho,"Metaanalysis")

  
}






