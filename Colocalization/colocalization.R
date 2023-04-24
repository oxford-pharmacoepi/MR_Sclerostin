# ============================================================================ #
#                               COLOCALIZATION                                 #
#                        2023 - Marta Alcalde-Herraiz                          #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

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
out  <- c('eBMD','HF','CAD','MI','IS','Hypertension','T2DM')
outc <- c('ebi-a-GCST006979','','ebi-a-GCST005194','ebi-a-GCST011365','','ukb-b-14057','')

unlink(here("Colocalization","coloc.xlsx"))
unlink(here("Colocalization","data.xlsx"))

ss <- c(0,0.0157,0.2243,0.097,0.3317,0.2586,0.0954)

for (i in 1:7){
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
  
  outcome_data <- outcome_data %>%
    filter(eaf.outcome > 0,
           eaf.outcome < 1) %>%
    mutate(id.outcome = "outcome",
           outcome    = "outcome")
  
  dat <- harmonise_data(exposure_data, outcome_data) %>%
    filter(remove == FALSE) %>%
    distinct(SNP, .keep_all = TRUE)
  write.xlsx(dat,here("Colocalization","data.xlsx"),sheetName = out[i],append = TRUE)
  
  if (i == 1){
    # Continuous variable (quantitative)
    res <- coloc.abf(dataset1 = list(snp = dat$SNP, 
                                         beta = dat$beta.outcome, 
                                         varbeta = dat$se.outcome^2, 
                                         pvalues = dat$pval.outcome,
                                         MAF = dat$eaf.outcome,
                                         N = dat$samplesize.outcome,
                                         type = "quant"),
                     dataset2 = list(snp = dat$SNP, 
                                     beta = dat$beta.exposure, 
                                     varbeta = dat$se.exposure^2, 
                                     pvalues = dat$pval.exposure,
                                     MAF = dat$eaf.exposure, 
                                     N = dat$n.exposure, 
                                     type = "quant"),
                     p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  }else{
    # Binary variable
    res <- coloc.abf(dataset1 = list(snp = dat$SNP,
                                     beta = dat$beta.exposure,
                                     varbeta = dat$se.exposure^2,
                                     pvalues = dat$pval.exposure,
                                     MAF = dat$eaf.exposure,
                                     N = dat$n.exposure,
                                     type = "quant"),
                     dataset2 = list(snp = dat$SNP,
                                     beta = dat$beta.outcome,
                                     varbeta = dat$se.outcome^2,
                                     pvalues = dat$pval.outcome,
                                     MAF = dat$eaf.outcome,
                                     s = ss[i], type = "cc"),
                     p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  }
  tab <- data.frame(nsnps = res$summary[1],
                    H0 = res$summary[2],
                    H1 = res$summary[3],
                    H2 = res$summary[4],
                    H3 = res$summary[5],
                    H4 = res$summary[6])
  write.xlsx(tab,here("Colocalization","coloc.xlsx"),sheetName = out[i],append = TRUE)
}








