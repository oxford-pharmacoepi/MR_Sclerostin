# ============================================================================ #
#                                phenotypingHes                                #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

phenotypingHes <- function(outcome, ukb, hesin, hesin_diag, pop){
  hes_codes <- getCodes(outcome)$HES

  hes_table <- hesin_diag %>%
    dplyr::select("eid","ins_index","diag_icd10") %>%
    # Filter patients having a code
    dplyr::filter(diag_icd10 %in% hes_codes) %>%
    dplyr::inner_join(hesin %>%
                        dplyr::select("eid","ins_index","epistart","admidate"),
                      by = c("eid","ins_index")) %>%
    # Minimum date between epistart and admidate
    dplyr::mutate(epistart = as.Date(epistart, "%d/%m/%Y")) %>%
    dplyr::mutate(admidate = as.Date(admidate, "%d/%m/%Y")) %>%
    dplyr::mutate(event_date = pmin(epistart,admidate, na.rm = TRUE)) %>%
    dplyr::group_by(eid) %>%
    # Year of first event
    dplyr::summarise(first_event = min(event_date, na.rm = TRUE)) %>%
    dplyr::left_join(ukb %>% dplyr::select("eid",
                                           "year_of_birth"   = "f.34.0.0",
                                           "date_assessment" = "f.53.0.0"),
                     by = "eid") %>%
    # Time of first event since birth
    dplyr::mutate(age_hes_birth = lubridate::year(first_event) - year_of_birth) %>%
    # Time of first event since first assessment
    dplyr::mutate(age_hes_enrol = lubridate::year(first_event) - lubridate::year(date_assessment)) %>%
    # Logistic state
    dplyr::mutate(state_hes     = 1) %>%
    # Remove participants with missing dates or negative ages - Build survival states
    dplyr::mutate(state_hes_birth     = dplyr::if_else(is.na(age_hes_birth)     | age_hes_birth < 0, -1, state_hes)) %>%
    dplyr::mutate(state_hes_enrol = dplyr::if_else(is.na(age_hes_enrol) | age_hes_enrol < 0, -1, state_hes)) %>%
    dplyr::select("eid","state_hes","state_hes_birth","age_hes_birth","state_hes_enrol","age_hes_enrol") %>%
    # If death, censure it
    dplyr::right_join(pop, by = "eid") %>%
    dplyr::left_join(
      ukb %>% dplyr::select("eid",
                            "year_of_birth" = "f.34.0.0",
                            "age_of_death"  = "f.40007.0.0"),
      by = "eid") %>%
    dplyr::mutate(age_hes_birth = dplyr::if_else(is.na(age_hes_birth) & !is.na(age_of_death), round(age_of_death), age_hes_birth),
                  age_hes_enrol = dplyr::if_else(is.na(age_hes_enrol) & !is.na(age_of_death), round(age_of_death), age_hes_enrol),
                  age_hes_birth = dplyr::if_else(is.na(age_hes_birth), 2020 - year_of_birth, age_hes_birth),
                  age_hes_enrol = dplyr::if_else(is.na(age_hes_enrol), 2020 - year_of_birth, age_hes_enrol),
                  state_hes       = dplyr::if_else(is.na(state_hes),0,state_hes),
                  state_hes_birth = dplyr::if_else(is.na(state_hes_birth),0,state_hes_birth),
                  state_hes_enrol = dplyr::if_else(is.na(state_hes_enrol),0,state_hes_enrol)) %>%

    dplyr::select("eid", "state_hes", "state_hes_birth", "age_hes_birth", "state_hes_enrol", "age_hes_enrol")

  return(hes_table)
}
