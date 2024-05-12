# ============================================================================ #
#                               COLOCALIZATION                                 #
#                        2024 - Marta Alcalde-Herraiz                          #
#                                                                              #
# This script uses SuSiE to perform the colocalization analysis.               #
# ============================================================================ #
library(coloc)
source(here::here("Functions/ld_matrix_local.R"))
source(here::here("Functions/loadOutcomeGwas.R"))

genstart <- 43753738
genend   <- 43758791
window   <- 20e3

# EXPOSURE ---------------------------------------------------------------------
exposure_fixed  <- readr::read_delim(paste0(pathData,"Results/Study1/InstrumentSelection/FixedMapped.txt")) %>%
  dplyr::filter(chr == 17, pos >= genstart - window, pos <= genend + window) %>%
  dplyr::mutate(exposure = "Fixed", id.exposure = "Fixed")  %>%
  dplyr::rename("samplesize" = "sample_size")

ld_matrix <- ld_matrix_local(exposure_fixed$SNP,
                      bfile    = paste0(pathData,'LD_ReferencePannel/EUR'),
                      plink_bin = paste0(pathData,'Plink/plink.exe'))

D1 <- prepare_data(ld_matrix, exposure_fixed, type = "quant")

check_dataset(D1, req="LD")

# OUTCOME ----------------------------------------------------------------------
outcome_names <- list("bmd" = "Heel bone mineral density",
                      "hf"  = "Hip fracture",
                      "cad" = "Coronary artery disease",
                      "mi"  = "Myocardial infarction",
                      "is"  = "Ischaemic stroke",
                      "hypertension" = "Hypertension",
                      "t2dm" = "Type 2 diabetes mellitus",
                      "hdl" = "HDL cholesterol",
                      "ldl" = "LDL cholesterol",
                      "glucose" = "Fasting glucose",
                      "hba1c" = "HbA1c")

outcome_list <- list()
for(i in names(outcome_names)){
  outcome <- loadOutcomeGwas(i,exposure_fixed %>%
                               dplyr::inner_join(tibble::tibble("SNP" = D1$snp)) %>%
                               TwoSampleMR::format_data(type = "exposure")) %>%
    dplyr::tibble() %>%
    dplyr::distinct() %>%
    dplyr::rename_with(~stringr::str_remove(.,".outcome"), everything())

 outcome_list[[i]] <- prepare_data(ld_matrix, outcome, type = "cc")

 check_dataset(outcome_list[[i]] , req="LD")
}

# SuSiE ------------------------------------------------------------------------
D1$N <- max(D1$N)
S1 <- runsusie(D1)



prepare_data <- function(ld_matrix, data, type = "quant"){
  # Harmonise alleles to the ld_matrix alleles
  pattern <- ld_matrix %>%
    tibble::as_tibble() %>%
    colnames() %>%
    tibble::tibble() %>%
    tidyr::separate(col = ".", into = c("SNP", "effect_allele.ld","other_allele.ld"), sep = "_")

  data1 <- pattern %>%
    dplyr::left_join(
      data, by = "SNP"
    ) %>%
    dplyr::mutate(
      beta = dplyr::if_else(effect_allele != effect_allele, -beta, beta),
      eaf  = dplyr::if_else(effect_allele != effect_allele, 1-eaf, eaf),
      effect_allele  = dplyr::if_else(effect_allele != effect_allele, effect_allele, effect_allele),
      other_allele   = dplyr::if_else(other_allele  != other_allele,  other_allele,  other_allele),
    ) %>%
    dplyr::select(-"effect_allele.ld", -"other_allele.ld")

  data1_no_passing_the_filters <- data1 %>%
    dplyr::filter(eaf <= 0 | eaf >= 1 | is.na(eaf)) %>%
    dplyr::summarise("SNP" = paste0(SNP,"_", effect_allele,"_", other_allele))

  data1 <- data1 %>%
    dplyr::filter(eaf > 0 & eaf < 1 & !is.na(eaf))

  # Remove non-valid snps from ld_matrix
  ld_matrix <- ld_matrix[row.names(ld_matrix) != data1_no_passing_the_filters$SNP,]
  ld_matrix <- ld_matrix[,colnames(ld_matrix) != data1_no_passing_the_filters$SNP]

  colnames(ld_matrix)  <- gsub("_.*","",colnames(ld_matrix))
  row.names(ld_matrix) <- gsub("_.*","",row.names(ld_matrix))

  D1 <- list("snp"  = data1$SNP,
             "beta" = data1$beta,
             "varbeta" = data1$se^2,
             "pvalues" = data1$pval,
             "N"    = data1$samplesize,
             "MAF"  = data1$eaf,
             "type" = "quant",
             "LD"   = ld_matrix)
  return(D1)
}
