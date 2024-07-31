# ============================================================================ #
#                             MendelianRandomisation                           #
#                             Marta Alcalde-Herraiz                            #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This script performs generalised IVW when there is more than one instrument, #
# and Wald Ratio when there is only one instrument.                            #
# ============================================================================ #
source(here::here("Functions","loadOutcomeGwas.R"))
source(here::here("Functions","harmonisingData.R"))
source(here::here("Functions","mendelianRandomisation.R"))
source(here::here("Functions","ld_matrix_local.R"))

instruments <- list.files(path = paste0(pathData,"Results/Study1/InstrumentSelection/Instruments/"),
                          recursive  = TRUE,
                          pattern    = "_gene",
                          full.names = TRUE)

outcomes    <- list.files(path = paste0(pathResults,"Phenotyping/"),
                          full.names = TRUE,
                          pattern = "Phenotype")

genetics    <- tibble::as_tibble(read.table(paste0(pathResults,"Cohort/genetic_data.txt"), header = TRUE)) %>%
  dplyr::mutate(rs9910625_0  = "T",
                rs80107551_0 = "T")

resultsMR <- list()
datMR     <- list()

for(i in instruments){
  # Select instruments ---------------------------------------------------------
  file <- i

  # Extract meta-analysis method (fixed or random)
  metaMethod <- gsub(".*iv_","",file)
  metaMethod <- gsub("_gene.*","",metaMethod)

  # Extract r2 value from the file name
  r2 <- gsub(".*r2","",file)
  r2 <- gsub("_.*","",r2)

  # Read exposure file
  exposure_table <- readr::read_delim(file)

  # Select only snps of interest
  genetics_table <- genetics %>%
    dplyr::select("eid", exposure_table$SNP, paste0(exposure_table$SNP,"_0"), "sex", "age_at_assessment", dplyr::starts_with("PC"), "batch")

  for(outcome_i in outcomes){
    # Extract regression method
    regressionMethod <- gsub(".*Phenotype_","",outcome_i)
    regressionMethod <- gsub("_.*","",regressionMethod)

    # Extract outcome
    outcomeName <- gsub(paste0(".*",regressionMethod,"_"),"",outcome_i)
    outcomeName <- gsub(".txt.*","",outcomeName)

    outcome_table <- tibble::as_tibble(read.table(outcome_i, header = TRUE))

    # Build outcome dataset
    if(regressionMethod == "Birth" | regressionMethod == "Enrol"){
      for(j in 1:nrow(exposure_table)){
        # Table with all the necessary variables
        table <- genetics_table %>%
          dplyr::inner_join(outcome_table, by = "eid") %>%
          dplyr::select("eid", "state", "age", "snp" = exposure_table$SNP[j], "sex", "age_at_assessment", "batch", dplyr::starts_with("PC")) %>%
          dplyr::filter(!is.na(snp))

        # Survival analysis - Adjusted model for sex and the first 10 principal components
        regression <- survival::coxph(survival::Surv(age,state) ~ snp + sex + age_at_assessment + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                                      data = table)

        regression <- coefficients(summary(regression))

        tableResults_i <- tibble::tibble(
          SNP      = exposure_table$SNP[j],
          samplesize.outcome = length(table$snp),
          cases    = sum(table$state == 1),
          controls = sum(table$state == 0),
          beta.outcome  = regression[1,1],
          se.outcome    = regression[1,3],
          pval.outcome  = regression[1,5],
          eaf.outcome   = sum(table$snp)/(2*length(table$snp)),
          effect_allele.outcome = genetics_table %>% dplyr::select("ea" = paste0(exposure_table$SNP[j],"_0")) %>% dplyr::pull() %>% unique()
        )

        if(j == 1){
          tableResults <- tableResults_i
        }else{
          tableResults <- tableResults %>%
            dplyr::union_all(tableResults_i)
        }
      }
    }else if(regressionMethod == "Logaritmic"){
      for(j in 1:nrow(exposure_table)){
        # Table with all the necessary variables
        table <- genetics_table %>%
          dplyr::inner_join(outcome_table, by = "eid") %>%
          dplyr::select("eid", "state", "snp" = exposure_table$SNP[j], "sex", "age_at_assessment", "batch", dplyr::starts_with("PC")) %>%
          dplyr::filter(!is.na(snp))

        # Logaritmic regression - Adjusted model for sex and the first 10 principal components
        regression <- stats::glm(state ~ snp + sex + age_at_assessment + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                                 data = table,
                                 family = "binomial")
        regression <- coefficients(summary(regression))

        tableResults_i <- tibble::tibble(
          SNP = exposure_table$SNP[j],
          samplesize.outcome  = length(table$snp),
          cases    = sum(table$state == 1),
          controls = sum(table$state == 0),
          beta.outcome = regression[2,1],
          se.outcome   = regression[2,2],
          pval.outcome    = regression[2,4],
          eaf.outcome  = sum(table$snp)/(2*length(table$snp)),
          effect_allele.outcome = genetics_table %>% dplyr::select("ea" = paste0(exposure_table$SNP[j],"_0")) %>% dplyr::pull() %>% unique()
        )

        if(j == 1){
          tableResults <- tableResults_i
        }else{
          tableResults <- tableResults %>%
            dplyr::union_all(tableResults_i)
        }
      }
    }else if(regressionMethod == "Logistic"){
      for(j in 1:nrow(exposure_table)){
        # Table with all the necessary variables
        table <- genetics_table %>%
          dplyr::inner_join(outcome_table, by = "eid") %>%
          dplyr::select("eid", "outcome", "snp" = exposure_table$SNP[j], "age_at_assessment", "sex", "batch", dplyr::starts_with("PC")) %>%
          dplyr::filter(!is.na(snp))

        # Logistic regression - Adjusted model for sex and the first 10 principal components
        regression <- stats::lm(outcome ~ snp + sex + age_at_assessment + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                                data = table)
        regression <- coefficients(summary(regression))

        tableResults_i <- tibble::tibble(
          SNP = exposure_table$SNP[j],
          samplesize.outcome  = length(table$snp),
          cases    = NA,
          controls = NA,
          beta.outcome = regression[2,1],
          se.outcome   = regression[2,2],
          pval.outcome    = regression[2,4],
          eaf.outcome     = sum(table$snp)/(2*length(table$snp)),
          effect_allele.outcome = genetics_table %>% dplyr::select("ea" = paste0(exposure_table$SNP[j],"_0")) %>% dplyr::pull() %>% unique()
        )

        if(j == 1){
          tableResults <- tableResults_i
        }else{
          tableResults <- tableResults %>%
            dplyr::union_all(tableResults_i)
        }
      }
    }

    # Format output data
    tableResults <- tableResults %>%
      dplyr::inner_join(exposure_table %>%
                          dplyr::select("SNP","effect_allele.exposure","other_allele.exposure")) %>%
      dplyr::mutate(other_allele.outcome = dplyr::if_else(effect_allele.outcome == effect_allele.exposure, other_allele.exposure, effect_allele.exposure)) %>%
      dplyr::select(-dplyr::ends_with(".exposure")) %>%
      dplyr::mutate(id.outcome = paste0(outcomeName,"_",regressionMethod),
                    outcome    = paste0(outcomeName,"_",regressionMethod))

    # Mendelian randomisation
    mr_res <- mendelianRandomisation(exposure = exposure_table,
                                     outcome  = tableResults,
                                     gwasID   = "" ,
                                     binary   = (regressionMethod != "Logistic"), name = paste0(outcomeName,"_",regressionMethod))

    varname <- paste0(metaMethod,"_",r2,"_",regressionMethod,"_",outcomeName)
    resultsMR[[varname]] <- mr_res$Results
    datMR[[varname]]     <- mr_res$Dat
  }
}


dir.create(paste0(pathResults,"MendelianRandomisation/"))
for(list_name_i in 1:length(resultsMR)){
  name_i <- names(resultsMR)
  readr::write_tsv(resultsMR[[list_name_i]],paste0(pathResults,"MendelianRandomisation/",name_i[list_name_i],".txt"))
}

dir.create(paste0(pathResults,"MendelianRandomisation/Data"))
for(list_name_i in 1:length(datMR)){
  name_i <- names(datMR)
  readr::write_tsv(datMR[[list_name_i]],paste0(pathResults,"MendelianRandomisation/Data/",name_i[list_name_i],".txt"))
}

