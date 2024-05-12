# ============================================================================ #
#                             Survival package                                 #
#                           Marta Alcalde-Herraiz                              #
#                                                                              #
# This script performs the sensitivity analysysis "Survival package", where we #
# repeat the survival analysis using the "gwassurvivr" package".               #
# ============================================================================ #
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("gwasurvivr")
# BiocManager::install("GWASdata")
# BiocManager::install("SNPRelate")

library(plinkFile)
library(gwasurvivr)
library(SNPRelate)
source(here::here("Functions/mendelianRandomisation.R"))
source(here::here("Functions/harmonisingData.R"))
source(here::here("Functions/ld_matrix_local.R"))
# Create covariate file --------------------------------------------------------
genetics <- tibble::as_tibble(read.table(paste0(pathData,"Results/Study2/Cohort/genetic_data.txt"), header = TRUE)) %>%
  dplyr::mutate("rs9910625_0"  = "T", "rs80107551_0" = "T")

covariate_file_all <- genetics %>%
  dplyr::select("eid", "sex", "batch", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                "PC7", "PC8", "PC9", "PC10")

outcomes_names <- list("Fracture" = "Fracture",
                       "CAD"   = "Coronary artery disease",
                       "MI"    = "Myocardial infarction",
                       "IS"    = "Ischaemic stroke",
                       "Hypertension" = "Hypertension",
                       "T2DM"  = "Type 2 diabetes")


rm("survival_results")
for(on_i in names(outcomes_names)){
  covariate_file <- covariate_file_all %>%
    dplyr::inner_join(
      read.table(paste0(pathData,"Results/Study2/Phenotyping/Phenotype_Birth_",on_i,".txt"), header = TRUE),
      by = "eid"
    ) %>%
  dplyr::mutate(
    "eid" = as.character(eid)
  )
  # Afterwards, convert them into plink files using PLINK. Read the generated files:

  bed <- readBED(paste0(pathData,"UKBiobank/selected_snps_c17_2.bed"),
                 iid = 1, vid = 1, vfr = NULL, vto = NULL, quiet = TRUE)

  plinkCoxSurv(
    bed.file = paste0(pathData,"UKBiobank/selected_snps_c17_2.bed"),
    covariate.file = covariate_file,
    id.column = "eid",
    sample.ids = covariate_file$eid,
    time.to.event = "age",
    event = "state",
    covariates = c("sex", "batch", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
    print.covs = "only",
    out.file   = paste0(pathData,"Results/SensitivityAnalyses/survival_results"),
    maf.filter = NULL
  )

  surv_results <- read.table(paste0(pathData,"Results/SensitivityAnalyses/survival_results.coxph"),
                             sep="\t", header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select("SNP" = "RSID",
                  "CHR" = "CHR",
                  "effect_allele.outcome" = "A1",
                  "other_allele.outcome"  = "A0",
                  "eaf.outcome"  = "exp_freq_A1",
                  "pval.outcome" = "PVALUE",
                  "beta.outcome" = "COEF",
                  "se.outcome"   = "SE.COEF",
                  "N" = "N",
                  "Cases" = "N.EVENT") %>%
    dplyr::mutate(
      "Outcome" = outcomes_names[[on_i]]
    )

  if("survival_results" %in% ls()){
    survival_results <- survival_results %>%
      dplyr::union_all(surv_results)
  }else{
    survival_results <- surv_results
  }
}


instruments <- read.table(paste0(pathData,"Results/Study1/InstrumentSelection/Instruments/iv_Fixed_gene500000kb_r20.3_clumpw250000.txt"), header = TRUE)

rm("mr_results")
for(on_i in names(outcomes_names)){
  mr_res <- mendelianRandomisation(instruments, survival_results %>%
                           dplyr::filter(Outcome == outcomes_names[[on_i]]) %>%
                           dplyr::mutate(id.outcome = Outcome, outcome = Outcome),
                         gwasID = NULL, binary = TRUE, name = NULL)

  if("mr_results" %in% ls()){
    mr_results <- mr_results %>% dplyr::union_all(mr_res$Results %>%
                                                    dplyr::mutate(Outcome = outcomes_names[[on_i]]))
  }else{
    mr_results <- mr_res$Results %>%
      dplyr::mutate(Outcome = outcomes_names[[on_i]])
  }
}

readr::write_delim(mr_results,paste0(pathData,"Results/SensitivityAnalyses/SurvivalPackage_results.txt"))
