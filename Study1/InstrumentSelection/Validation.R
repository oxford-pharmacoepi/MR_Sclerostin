# ============================================================================ #
#                                  Validation                                  #
#                             Marta Alcalde-Herraiz                            #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This scripts validates the instruments based on the positive control         #
# association.                                                                 #
# ============================================================================ #

source(here::here('Functions','harmonisingData.R'))
dir.create(paste0(pathResults,"InstrumentSelection/Validation"))
metaMethod <- c("Fixed","Random")

for(metaMethod_i in metaMethod){
  list_of_files <- list.files(path = paste0(pathResults,"InstrumentSelection/Instruments"),
                              recursive  = TRUE,
                              pattern    = paste0("iv_",metaMethod_i,"_gene"),
                              full.names = TRUE)
  list_of_files <- list_of_files[list_of_files != "C:/Users/martaa/Desktop/Projects/MR_Sclerostin_NEW/Results/Study1/InstrumentSelection/Instruments/iv_Fixed_gene100000kb_r20.8_clumpw500000.txt"]

  exposure <- readr::read_delim(list_of_files) %>% dplyr::distinct()

  # Step 1 of validation: Check the association of this SNPs with BMD
  # Article: https://www.nature.com/articles/s41588-018-0302-x
  # Download GWAS: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006979/
  bmd <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/bmd/30598549-GCST006979-EFO_0009270.h.tsv"))
  bmd <- bmd %>%
    dplyr::select(
      SNP = hm_rsid,
      effect_allele = hm_effect_allele,
      other_allele  = hm_other_allele,
      beta = hm_beta,
      eaf  = hm_effect_allele_frequency,
      se   = standard_error,
      pval = p_value) %>%
    dplyr::mutate(
      samplesize = 426824,
      id = 'outcome'
    ) %>%
    dplyr::filter(SNP %in% exposure$SNP)
  bmd <- TwoSampleMR::format_data(dat = bmd, type = "outcome")

  dat_harmonised <- harmonisingData(exposure = exposure,
                                    outcome = bmd,
                                    gwasID = "ebi-a-GCST006979")

  readr::write_delim(dat_harmonised,paste0(pathResults,"InstrumentSelection/Validation/Validation_",metaMethod_i,"step1.txt"))

  # Step 2 of validation: Check the association of this SNPs with hip fracture risk
  # Article: https://www.sciencedirect.com/science/article/pii/S2666379122003317?via%3Dihub
  # Download GWAS: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90161001-GCST90162000/GCST90161240/
  hf <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/hipFracture/GCST90161240_buildGRCh37.tsv"))
  hf <- hf %>%
    dplyr::select(
      SNP = variant_id,
      pval = p_value,
      effect_allele, other_allele,
      eaf = effect_allele_frequency,
      beta,
      se = standard_error) %>%
    dplyr::mutate(
      samplesize = 735354,
      id = "outcome") %>%
    dplyr::filter(SNP %in% exposure$SNP)
  hf <- TwoSampleMR::format_data(dat = hf, type = "outcome")

  dat_harmonised <- harmonisingData(exposure = exposure,
                                    outcome = hf,
                                    gwasID = "ebi-a-GCST90161240")

  readr::write_delim(dat_harmonised,paste0(pathResults,"InstrumentSelection/Validation/Validation_",metaMethod_i,"step2.txt"))
}
# Step 3
# https://gtexportal.org/home/

