# ============================================================================ #
#                            mendelianRandomisation                            #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #


mendelianRandomisation <- function(exposure, outcome, gwasID = NULL, binary = TRUE, name = NULL){
  # Single variant or more than one variant
  if(nrow(exposure) == 1){
    dat <- harmonisingData(exposure = exposure,
                           outcome  = outcome,
                           gwasID   = gwasID)

    mrObject <- MendelianRandomization::mr_input(
      bx = dat$beta.exposure,
      bxse = dat$se.exposure,
      by = dat$beta.outcome,
      byse = dat$se.outcome,
      snps = dat$SNP)

    mrResults <- MendelianRandomization::mr_ivw(
      mrObject,
      correl = FALSE,
    )

  }else{

    ld_matrix <- ld_matrix_local(
      exposure$SNP,
      bfile    = paste0(pathData,'LD_ReferencePannel/EUR'),
      plink_bin = paste0(pathData,'Plink/plink.exe'))

    pattern <- ld_matrix %>%
      tibble::as_tibble() %>%
      colnames() %>%
      tibble::tibble() %>%
      tidyr::separate(col = ".", into = c("SNP", "effect_allele","other_allele"), sep = "_")

    dat_harmonised <- pattern %>%
      dplyr::left_join(
        harmonisingData(exposure = exposure,
                        outcome  = outcome,
                        gwasID   = gwasID),
        by = "SNP"
      ) %>%
      dplyr::mutate(
        beta.exposure = dplyr::if_else(effect_allele.exposure != effect_allele, -beta.exposure, beta.exposure),
        eaf.exposure  = dplyr::if_else(effect_allele.exposure != effect_allele, 1-eaf.exposure, eaf.exposure),
        effect_allele.exposure  = dplyr::if_else(effect_allele.exposure != effect_allele, effect_allele, effect_allele.exposure),
        other_allele.exposure   = dplyr::if_else(other_allele.exposure  != other_allele,  other_allele,  other_allele.exposure),
      ) %>%
      dplyr::mutate(
        beta.outcome = dplyr::if_else(effect_allele.outcome != effect_allele, -beta.outcome, beta.outcome),
        eaf.outcome  = dplyr::if_else(effect_allele.outcome != effect_allele, 1-eaf.outcome, eaf.outcome),
        effect_allele.outcome  = dplyr::if_else(effect_allele.outcome != effect_allele, effect_allele, effect_allele.outcome),
        other_allele.outcome   = dplyr::if_else(other_allele.outcome  != other_allele,  other_allele,  other_allele.outcome),
      )

    dat <- TwoSampleMR::harmonise_ld_dat(dat_harmonised,ld_matrix)

    mrObject <- MendelianRandomization::mr_input(
      bx = dat$x$beta.exposure,
      bxse = dat$x$se.exposure,
      by = dat$x$beta.outcome,
      byse = dat$x$se.outcome,
      corr = dat$ld,
      snps = dat$x$SNP)

    mrResults <- MendelianRandomization::mr_ivw(
      mrObject,
      correl = TRUE,
    )

    dat <- dat$x
   }

  tableResults <- list()

  tableResults$Dat <- dat
  if(binary){
    tableResults$Results <- tibble::tibble(
      "Instruments" = nrow(exposure),
      "Method"      = "IVW",
      "Outcome"     = .env$name,
      "OR" = exp(-mrResults$Estimate),
      "U95CI" = exp(-mrResults$Estimate+1.96*mrResults$StdError),
      "L95CI" = exp(-mrResults$Estimate-1.96*mrResults$StdError),
      "SE" = mrResults$StdError,
      "Pval" = mrResults$Pvalue)
  }else{
    tableResults$Results <- tibble::tibble(
      "Instruments" = nrow(exposure),
      "Method"      = "IVW",
      "Outcome"     = .env$name,
      "Beta" = -mrResults$Estimate,
      "U95CI" = -mrResults$Estimate+1.96*mrResults$StdError,
      "L95CI" = -mrResults$Estimate-1.96*mrResults$StdError,
      "SE" = mrResults$StdError,
      "Pval" = mrResults$Pvalue)
  }

  return(tableResults)
}
