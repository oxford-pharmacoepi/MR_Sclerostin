# ============================================================================ #
#                                 Phenotyping                                  #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

pop <- readr::read_delim(paste0(pathResults,"Cohort/ukb.txt")) %>% dplyr::select("eid")

# Binary outcomes --------------------------------------------------------------
outcomes <- c("Fracture", "CAD", "MI", "IS", "Hypertension", "T2DM")

ukb        <- readr::read_delim(paste0(pathResults,"Cohort/ukb.txt"))
gp         <- readr::read_delim(paste0(pathResults,"Cohort/gp.txt"))
hesin      <- readr::read_delim(paste0(pathResults,"Cohort/hesin.txt"))
hesin_diag <- readr::read_delim(paste0(pathResults,"Cohort/hesin_diag.txt"))

source(here::here("Functions/getCodes.R"))
source(here::here("Functions/phenotypingHes.R"))
source(here::here("Functions/phenotypingGP.R"))
source(here::here("Functions/phenotypingUKB.R"))

for(i in 1:length(outcomes)){
  hes_table <- phenotypingHes(outcome = outcomes[[i]],
                              ukb     = ukb,
                              hesin   = hesin,
                              hesin_diag = hesin_diag,
                              pop     = pop)

  gp_table  <- phenotypingGP(outcome = outcomes[[i]],
                             ukb     = ukb,
                             gp      = gp,
                             pop     = pop)

  ukb_table <- phenotypingUKB(outcome = outcomes[[i]],
                              ukb = ukb,
                              pop = pop)

  # Merging tables
  logistic_table <- hes_table %>%
    dplyr::left_join(gp_table,  by = "eid") %>%
    dplyr::left_join(ukb_table, by = "eid") %>%
    dplyr::distinct() %>%
    dplyr::filter(state_ukb != -1, state_hes != -1, state_gp != -1) %>%
    dplyr::mutate(state = dplyr::if_else(state_hes == 1 | state_gp == 1 | state_ukb == 1,1,0)) %>%
    dplyr::select("eid", "state_hes", "state_gp", "state_ukb", "state")

  birth_table <- hes_table %>%
    dplyr::left_join(gp_table,  by = "eid") %>%
    dplyr::left_join(ukb_table, by = "eid") %>%
    dplyr::distinct() %>%
    dplyr::filter(state_ukb_birth != -1, state_hes_birth != -1, state_gp_birth != -1) %>%
    dplyr::mutate(state = dplyr::if_else(state_hes_birth == 1 | state_gp_birth == 1 | state_ukb_birth == 1,1,0),
                  age   = dplyr::if_else(state == 1, pmin(age_hes_birth, age_gp_birth, age_ukb_birth, na.rm = TRUE), age_hes_birth)) %>%
    dplyr::select("eid", "state_hes_birth", "state_gp_birth", "state_ukb_birth", "state", "age")

  enrol_table <- hes_table %>%
    dplyr::left_join(gp_table,  by = "eid") %>%
    dplyr::left_join(ukb_table, by = "eid") %>%
    dplyr::distinct() %>%
    dplyr::filter(state_ukb_enrol != -1, state_hes_enrol != -1, state_gp_enrol != -1) %>%
    dplyr::mutate(state = dplyr::if_else(state_hes_enrol == 1 | state_gp_enrol == 1 | state_ukb_enrol == 1,1,0),
                  age   = dplyr::if_else(state == 1, pmin(age_hes_enrol, age_gp_enrol, na.rm = TRUE), age_hes_enrol)) %>%
    dplyr::select("eid", "state_hes_enrol", "state_gp_enrol", "state_ukb_enrol", "state", "age")

  # Statistics tables
  stats_table <- logistic_table %>%
    dplyr::group_by(state_hes, state_gp, state_ukb) %>%
    dplyr::tally() %>%
    dplyr::mutate(Phenotye = paste0(outcomes[[i]],"_logistic")) %>%
    dplyr::union_all(
      birth_table %>%
        dplyr::rename("state_hes" = "state_hes_birth", "state_gp" = "state_gp_birth", "state_ukb" = "state_ukb_birth") %>%
        dplyr::group_by(state_hes, state_gp, state_ukb) %>%
        dplyr::tally() %>%
        dplyr::mutate(Phenotye = paste0(outcomes[[i]],"_birth"))
    ) %>%
    dplyr::union_all(
      enrol_table %>%
        dplyr::rename("state_hes" = "state_hes_enrol", "state_gp" = "state_gp_enrol", "state_ukb" = "state_ukb_enrol") %>%
        dplyr::group_by(state_hes, state_gp, state_ukb) %>%
        dplyr::tally() %>%
        dplyr::mutate(Phenotye = paste0(outcomes[[i]],"_enrol"))
    )

  if(i == 1){
    st <- stats_table
  }else{
    st <- st %>%
      dplyr::union_all(stats_table)
  }
  readr::write_delim(logistic_table %>% dplyr::select("eid","state"), paste0(pathResults,"Phenotyping/Phenotype_Logaritmic_",outcomes[[i]],".txt"))
  readr::write_delim(birth_table %>% dplyr::select("eid","state","age"), paste0(pathResults,"Phenotyping/Phenotype_Birth_",outcomes[[i]],".txt"))
  readr::write_delim(enrol_table %>% dplyr::select("eid","state","age"), paste0(pathResults,"Phenotyping/Phenotype_Enrol_",outcomes[[i]],".txt"))
}

readr::write_delim(st, paste0(pathResults,"Phenotyping/Statistics_Categorical.txt"))

# Continuous outcomes ----------------------------------------------------------
outcomes <- c("ebmd", "cholesterol", "ldl", "hdl", "triglycerides", "apoA",
              "apoB", "crp", "lipoprotein", "hba1c", "glucose")

names <- c("Heel bone mineral density [g/cm2]",
           "Cholesterol [mmol/L]",
           "LDL-Cholesterol [mmol/L]",
           "HDL-Cholesterol [mmol/L]",
           "Triglycerides [mmol/L]",
           "Apolipoprotein A [g/L]",
           "Apolipoprotein B [g/L]",
           "C-Reactive protein [mg/L]",
           "Lipoprotein A [nmol/L]",
           "Glycated haemoglobin (HbA1c) [mmol/mol]",
           "Glucose [mmol/L]")

for(i in 1:length(outcomes)){
  ukb_table <- ukb %>%
    dplyr::select("eid", "outcome" = getCodes(outcomes[[i]])) %>%
    dplyr::filter(!is.na(outcome))

  # Statistics tables
  stats_table <- tibble::tibble(outcome = outcomes[[i]],
                                N = ukb_table %>% dplyr::tally() %>% dplyr::pull(),
                                Mean = mean(ukb_table$outcome),
                                SD   = sd(ukb_table$outcome),
                                q0.05 = quantile(ukb_table$outcome, probs = 0.05),
                                q0.25 = quantile(ukb_table$outcome, probs = 0.25),
                                q0.75 = quantile(ukb_table$outcome, probs = 0.75),
                                q0.95 = quantile(ukb_table$outcome, probs = 0.95)
  )

  # Histogram to explore outliers
  ggplot2::ggplot(data = ukb_table, ggplot2::aes(x = outcome)) +
    ggplot2::geom_histogram(color = "black", fill = "grey") +
    ggplot2::ylab(label = "Counts") +
    ggplot2::xlab(label = names[[i]])

  ukb_table <- ukb_table %>%
    dplyr::mutate(outcome = (outcome-mean(outcome))/sd(outcome))

  stats_table <- stats_table %>%
    dplyr::mutate(N_2 = ukb_table %>% dplyr::tally() %>% dplyr::pull(),
                  Mean_2 = mean(ukb_table$outcome),
                  SD_2   = sd(ukb_table$outcome),
                  q0.05_2 = quantile(ukb_table$outcome, probs = 0.05),
                  q0.25_2 = quantile(ukb_table$outcome, probs = 0.25),
                  q0.75_2 = quantile(ukb_table$outcome, probs = 0.75),
                  q0.95_2 = quantile(ukb_table$outcome, probs = 0.95))


  if(i == 1){
    st <- stats_table
  }else{
    st <- st %>%
      dplyr::union_all(stats_table)
  }
  readr::write_delim(ukb_table %>% dplyr::select("eid","outcome"), paste0(pathResults,"Phenotyping/Phenotype_Logistic_",outcomes[[i]],".txt"))
}

readr::write_delim(st, paste0(pathResults,"Phenotyping/Statistics_Continuous.txt"))
