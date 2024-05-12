# ============================================================================ #
#                                  Processing                                  #
#                            Marta Alcalde-Herraiz                             #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This file reads the files created after using GWAS and METAL and converts    #
# the results into .txt                                                        #
# ============================================================================ #
fixed_gwama  <- tibble::as_tibble(readr::read_delim(here::here('Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Fixed.out')))
random_gwama <- tibble::as_tibble(readr::read_delim(here::here('Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Random.out')))
fixed_metal  <- tibble::as_tibble(readr::read_delim(here::here('Study1/MetaAnalysis/ResultsMetaAnalysis/METAL_Fixed.TBL'), delim = "\t"))

fixed <- fixed_gwama %>%
  dplyr::mutate(se = (beta_95U - beta)/1.95) %>%
  dplyr::select(SNP = rs_number,
                effect_allele = reference_allele,
                other_allele  = other_allele,
                beta = beta,
                se   = se,
                pval = `p-value`,
                i2   = i2,
                samplesize  = n_samples,
                q_statistic = q_statistic,
                q_pvalue    = `q_p-value`,
                zscore = z) %>%
  dplyr::inner_join(
    fixed_metal %>%
      dplyr::select(SNP = MarkerName,
                    effect_allele1 = Allele1,
                    other_allele1  = Allele2,
                    eaf = Freq1) %>%
      dplyr::mutate(effect_allele1 = toupper(effect_allele1),
                    other_allele1  = toupper(other_allele1)),
    by = 'SNP') %>%
  dplyr::mutate(eaf = dplyr::if_else(effect_allele == effect_allele1, eaf, 1-eaf)) %>%
  dplyr::select(-effect_allele1, -other_allele1)

random <- random_gwama %>%
  dplyr::mutate(se = (beta_95U - beta)/1.95) %>%
  dplyr::select(SNP = rs_number,
                effect_allele = reference_allele,
                other_allele  = other_allele,
                beta = beta,
                se   = se,
                pval = `p-value`,
                i2   = i2,
                samplesize = n_samples,
                q_statistic = q_statistic,
                q_pvalue    = `q_p-value`,
                zscore = z) %>%
  dplyr::inner_join(
    fixed_metal %>%
      dplyr::select(SNP = MarkerName,
                    effect_allele1 = Allele1,
                    other_allele1  = Allele2,
                    eaf = Freq1) %>%
      dplyr::mutate(effect_allele1 = toupper(effect_allele1),
                    other_allele1  = toupper(other_allele1)),
    by = 'SNP') %>%
  dplyr::mutate(eaf = dplyr::if_else(effect_allele == effect_allele1, eaf, 1-eaf)) %>%
  dplyr::select(-effect_allele1, -other_allele1)

readr::write_tsv(fixed, paste0(pathResults,'MetaAnalysis/Fixed.txt'))
readr::write_tsv(random,paste0(pathResults,'MetaAnalysis/Random.txt'))

