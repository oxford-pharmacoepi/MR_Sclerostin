# ============================================================================ #
#                                 findProxies                                  #
#                           Marta Alcalde-Herraiz                              #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This script finds proxies for the SNPs that cannot be found in an outcome    #
# gwas.                                                                        #
# ============================================================================ #
findProxies <- function(exposure,outcome, gwasID = NULL){
  # Extract proxies
  exposure_proxies <- exposure %>% dplyr::anti_join(outcome %>%
                                                      dplyr::select(SNP),
                                                    by = 'SNP')

  if(is.null(gwasID)){
    missage("The SNP ", exposure_proxies$SNP, " needs a proxy. Please provide the gwasID")
  }else{
    proxies <- TwoSampleMR::extract_outcome_data(snps = exposure_proxies$SNP,
                                                 outcomes = gwasID,
                                                 proxies  = TRUE,
                                                 rsq = 0.8,
                                                 align_alleles = 1,
                                                 palindromes   = 1,
                                                 maf_threshold = 0.01,
                                                 proxy_splitsize = 500,
                                                 splitsize = 10000)
    # Add columns of interest into the outcome table
    outcome <- outcome %>%
      dplyr::mutate(target_snp.outcome = NA,
                    proxy_snp.outcome  = NA)

    # Add missing columns of outcome to proxies
    proxies <- proxies %>%
      dplyr::mutate(pval_origin.outcome = NA)

    outcome <- outcome %>%
      dplyr::union_all(proxies %>%
                         dplyr::select(colnames(outcome)))

    return(outcome)
  }
}


