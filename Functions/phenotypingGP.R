# ============================================================================ #
#                                phenotypingGP                                 #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

phenotypingGP <- function(outcome, ukb, gp, pop){
  gp_codes <- getCodes(outcome)$GP

  gp_table <- gp %>%
    # Filter patients with the code
    dplyr::filter(read_2 %in% gp_codes | read_3 %in% gp_codes) %>%
    dplyr::group_by(eid) %>%
    # Event dates
    dplyr::reframe(event_date = as.Date(event_dt,"%d/%m/%Y")) %>%
    # First event date
    dplyr::group_by(eid) %>%
    dplyr::mutate(first_event = min(event_date, na.rm = TRUE)) %>%
    dplyr::left_join(ukb %>%
                       dplyr::select("eid",
                                     "year_of_birth"   = "f.34.0.0",
                                     "date_assessment" = "f.53.0.0"),
                     by = "eid") %>%
    # Time since first event since birth
    dplyr::mutate(age_gp_birth = lubridate::year(first_event) - year_of_birth) %>%
    # Time since first event after enrolment
    dplyr::mutate(age_gp_enrol = lubridate::year(first_event) - lubridate::year(date_assessment)) %>%
    dplyr::mutate(state_gp = 1) %>%
    # Participant with missing dates or negative age
    dplyr::mutate(state_gp_birth = dplyr::if_else(is.na(age_gp_birth) | age_gp_birth < 0, -1, state_gp)) %>%
    dplyr::mutate(state_gp_enrol = dplyr::if_else(is.na(age_gp_enrol) | age_gp_enrol < 0, -1, state_gp)) %>%
    # All participants
    dplyr::right_join(pop, by = "eid") %>%
    dplyr::mutate(state_gp       = dplyr::if_else(is.na(state_gp),0,state_gp)) %>%
    dplyr::mutate(state_gp_birth = dplyr::if_else(is.na(state_gp_birth),0,state_gp_birth)) %>%
    dplyr::mutate(state_gp_enrol = dplyr::if_else(is.na(state_gp_enrol),0,state_gp_enrol)) %>%
    dplyr::ungroup() %>%
    dplyr::select("eid", "state_gp", "state_gp_birth", "age_gp_birth", "state_gp_enrol", "age_gp_enrol")

  return(gp_table)

}
