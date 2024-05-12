# ============================================================================ #
#                               COLOCALIZATION                                 #
#                            Marta Alcalde-Herraiz                             #
#                                                                              #
# This script runs the colocalisation analysis betwheen sclerostin encoding    #
# gene and the outcomes of interest.                                           #
# ============================================================================ #
source(here::here("Functions/loadOutcomeGwas.R"))

genstart <- 43753738
genend   <- 43758791
window   <- 20e3

# EXPOSURE ---------------------------------------------------------------------
exposure_fixed  <- readr::read_delim(paste0(pathData,"Results/Study1/InstrumentSelection/FixedMapped.txt")) %>%
  dplyr::filter(chr == 17, pos >= genstart - window, pos <= genend + window) %>%
  dplyr::mutate(exposure = "Fixed", id.exposure = "Fixed")  %>%
  dplyr::rename("samplesize" = "sample_size")
# exposure_random <- readr::read_delim(paste0(pathData,"Results/Study1/InstrumentSelection/RandomMapped.txt")) %>%
#   dplyr::filter(chr == 17, pos >= genstart - window, pos <= genend + window) %>%
#   dplyr::mutate(exposure = "Random", id.exposure = "Random") %>%
#   dplyr::rename("samplesize" = "sample_size")

# COLOCALISATION ---------------------------------------------------------------
exposure_i <- TwoSampleMR::format_data(exposure_fixed,"exposure", phenotype_col = "exposure")

for(outcome_i  in c("bmd","hf","cad","mi","is","hypertension","t2dm","ldl","hdl","glucose","hba1c")){

  if(outcome_i == "hypertension"){
    outcome <- TwoSampleMR::extract_outcome_data(snps = exposure_i$SNP, outcomes = 'ukb-b-14057') %>%
      dplyr::select(SNP,pos,beta.outcome, se.outcome, samplesize.outcome, pval.outcome, eaf.outcome, effect_allele.outcome,
                    other_allele.outcome) %>%
      dplyr::mutate(outcome = "outcome", id.outcome = "outcome")

  }else{
    outcome    <- loadOutcomeGwas(outcome_i, exposure_i) %>% tibble::as_tibble() %>% dplyr::distinct() %>%
      dplyr::inner_join(exposure_i %>% dplyr::select(SNP), by = "SNP")
  }

  outcome <- outcome %>%
    dplyr::filter(SNP != "rs77697917")

  data <- TwoSampleMR::harmonise_data(
    exposure_i %>%
      dplyr::filter(!is.na(eaf.exposure), eaf.exposure > 0, eaf.exposure < 1),
    outcome %>%
      dplyr::filter(!is.na(eaf.outcome), eaf.outcome > 0, eaf.outcome < 1)) %>%
    tibble::as_tibble() %>%
    dplyr::filter(remove == FALSE)

  readr::write_delim(data,paste0(pathData,"Results/Colocalization/",outcome_i,".txt"))

  if(outcome_i %in% c("bmd","hdl","ldl","glucose","hba1c")){
    res <- coloc::coloc.abf(
      dataset2 = list("snp"  = data$SNP,
                      "beta" = data$beta.outcome,
                      "varbeta" = data$se.outcome^2,
                      "pvalues" = data$pval.outcome,
                      "N"    = data$samplesize.outcome,
                      "MAF"  = data$eaf.outcome,
                      "type" = "quant"),
      dataset1 = list("snp"  = data$SNP,
                      "beta" = data$beta.exposure,
                      "varbeta" = data$se.exposure^2,
                      "pvalues" = data$pval.exposure,
                      "MAF"  = data$eaf.exposure,
                      "N"    = data$samplesize.exposure,
                      "type" = "quant"),
      p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  }else{
    res <- coloc::coloc.abf(
      dataset2 = list("snp"  = data$SNP,
                      "beta" = data$beta.outcome,
                      "varbeta" = data$se.outcome^2,
                      "pvalues" = data$pval.outcome,
                      "N"    = data$samplesize.outcome,
                      "MAF"  = data$eaf.outcome,
                      "type" = "cc"),
      dataset1 = list("snp"  = data$SNP,
                      "beta" = data$beta.exposure,
                      "varbeta" = data$se.exposure^2,
                      "pvalues" = data$pval.exposure,
                      "MAF"  = data$eaf.exposure,
                      "N"    = data$samplesize.exposure,
                      "type" = "quant"),
      p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  }

  tab <- data.frame(nsnps = res$summary[1],
                    H0 = res$summary[2],
                    H1 = res$summary[3],
                    H2 = res$summary[4],
                    H3 = res$summary[5],
                    H4 = res$summary[6],
                    exposure_method = exposure_i$exposure %>% unique(),
                    outcome_method  = outcome_i)

  if("coloc_table" %in% ls()){
    coloc_table <- coloc_table %>%
      dplyr::union_all(tab)
  }else{
    coloc_table <- tab
  }
 # png(file = paste0(pathResults,"/SensitivityAnalysis/",outcome_i,".png"), width = 3000, height = 2750, units = "px",res = 300)
 # coloc::sensitivity(res, "H4>H3",   preserve.par = FALSE)
 # dev.off()
}

readr::write_delim(coloc_table, paste0(pathResults,"/colocalisation.txt"))




