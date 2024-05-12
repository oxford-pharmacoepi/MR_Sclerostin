# ============================================================================ #
#                               phenotypingUKB                                 #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

phenotypingUKB <- function(outcome, ukb, pop){
  ukb_table = switch(
    outcome,

    "CAD" = {
      ukb_table <- data.frame("eid" = pop$eid,
                              "state_ukb" = 0,
                              "state_ukb_birth" = 0,
                              "age_ukb_birth"   = 1000,
                              "state_ukb_enrol" = 0)
      return(ukb_table)
    },

    "Fracture" = {
      # Field 3005 ------------------------------------------------------------
      ukb1 <- ukb %>%
        dplyr::select("eid",
                      "date_assesment" = "f.53.0.0",
                      "year_of_birth"  = "f.34.0.0",
                      "state_ukb_1"    = "f.3005.0.0") %>%
        dplyr::filter(state_ukb_1 == 1) %>%
        dplyr::mutate(state_ukb_1 = as.double(state_ukb_1)) %>%
                      # For birth cohort, calculate age
        dplyr::mutate(age_ukb_1_birth = lubridate::year(date_assesment) - year_of_birth,
                      # For enrolment cohort, remove people with previous fractures
                      age_ukb_1_enrol = -1) %>%
        # Remove participants with negative ages or missing dates
        dplyr::mutate(state_ukb_1_birth = dplyr::if_else(is.na(age_ukb_1_birth) | age_ukb_1_birth < 0, -1, state_ukb_1),
                      state_ukb_1_enrol = dplyr::if_else(is.na(age_ukb_1_enrol) | age_ukb_1_enrol < 0, -1, state_ukb_1))

      # Field 20002.0.0 --------------------------------------------------------
      codes_ukb <- getCodes(outcome)$UKB_CODE2

      ukb2 <- searchUKBcode(ukb, codes_ukb) %>%
        dplyr::rename("state_ukb_2" = "state_ukb",
                      "state_ukb_2_birth" = "state_ukb_birth",
                      "age_ukb_2_birth"   = "age_ukb_birth",
                      "state_ukb_2_enrol" = "state_ukb_enrol")

      # Merge the two fields
      ukb_table <- ukb1 %>%
        dplyr::full_join(ukb2,
                         by = "eid") %>%
        # States
        dplyr::mutate(state_ukb_1 = dplyr::if_else(is.na(state_ukb_1),0,state_ukb_1),
                      state_ukb_1_birth = dplyr::if_else(is.na(state_ukb_1_birth), 0, state_ukb_1_birth),
                      state_ukb_1_enrol = dplyr::if_else(is.na(state_ukb_1_enrol), 0, state_ukb_1_enrol),
                      state_ukb_2 = dplyr::if_else(is.na(state_ukb_2),0,state_ukb_2),
                      state_ukb_2_birth = dplyr::if_else(is.na(state_ukb_2_birth), 0, state_ukb_2_birth),
                      state_ukb_2_enrol = dplyr::if_else(is.na(state_ukb_2_enrol), 0, state_ukb_2_enrol)) %>%
        # Remove those with missing dates or negative age
        dplyr::mutate(state_ukb       = dplyr::if_else(state_ukb_1 == -1 | state_ukb_2 == -1, -1, 1),
                      state_ukb_birth = dplyr::if_else(state_ukb_1_birth == -1 | state_ukb_2_birth == -1, -1, 1),
                      state_ukb_enrol = dplyr::if_else(state_ukb_1_enrol == -1 | state_ukb_2_enrol == -1, -1, 1)) %>%
        # Age at first assessment
        dplyr::mutate(age_ukb_birth = pmin(age_ukb_1_birth,age_ukb_2_birth,na.rm = TRUE)) %>%
        # All population
        dplyr::right_join(pop, by = "eid") %>%
        dplyr::mutate(state_ukb = dplyr::if_else(is.na(state_ukb), 0, state_ukb),
                      state_ukb_birth = dplyr::if_else(is.na(state_ukb_birth), 0, state_ukb_birth),
                      state_ukb_enrol = dplyr::if_else(is.na(state_ukb_enrol), 0, state_ukb_enrol)) %>%
        dplyr::select("eid", "state_ukb", "state_ukb_birth","age_ukb_birth", "state_ukb_enrol")
      return(ukb_table)
    },

    "Hypertension" = {
      codes_ukb <- getCodes(outcome)$UKB_CODE
      ukb1      <- searchUKBcode(ukb, codes_ukb)

      ukb_data <- ukb1 %>%
        dplyr::right_join(pop, by = "eid") %>%
        dplyr::mutate(state_ukb = dplyr::if_else(is.na(state_ukb), 0, state_ukb),
                      state_ukb_birth = dplyr::if_else(is.na(state_ukb_birth), 0, state_ukb_birth),
                      state_ukb_enrol = dplyr::if_else(is.na(state_ukb_enrol), 0, state_ukb_enrol)) %>%
        dplyr::select("eid", "state_ukb" , "state_ukb_birth", "age_ukb_birth", "state_ukb_enrol")

      return(ukb_data)
    },

    "IS" = {
      codes_ukb <- getCodes(outcome)$UKB_CODE
      ukb <- searchUKBcode(ukb, codes_ukb)

      ukb_data <- ukb %>%
        dplyr::right_join(pop, by = "eid") %>%
        dplyr::mutate(state_ukb = dplyr::if_else(is.na(state_ukb), 0, state_ukb),
                      state_ukb_birth = dplyr::if_else(is.na(state_ukb_birth), 0, state_ukb_birth),
                      state_ukb_enrol = dplyr::if_else(is.na(state_ukb_enrol), 0, state_ukb_enrol)) %>%
        dplyr::select("eid", "state_ukb" , "state_ukb_birth", "age_ukb_birth", "state_ukb_enrol")

      return(ukb_data)
    },

    "MI" = {
      codes_ukb <- getCodes(outcome)$UKB_CODE

      ukb1 <- searchUKBcode(ukb, codes_ukb) %>%
        dplyr::rename("state_ukb_1" = "state_ukb",
                      "state_ukb_1_birth" = "state_ukb_birth",
                      "age_ukb_1_birth"   = "age_ukb_birth",
                      "state_ukb_1_enrol" = "state_ukb_enrol")

      ukb2 <- ukb %>%
        dplyr::select("eid","age_ukb_2_birth" = "f.3894.0.0") %>%
        dplyr::filter(!is.na(age_ukb_2_birth),
                      age_ukb_2_birth != -1,
                      age_ukb_2_birth != -3) %>%
        dplyr::mutate(state_ukb_2 = 1,
                      state_ukb_2_birth = 1,
                      state_ukb_2_enrol = -1) %>%
        dplyr::mutate(state_ukb_2_birth = dplyr::if_else(is.na(age_ukb_2_birth) | age_ukb_2_birth < 0,-1,state_ukb_2_birth))

      ukb_data <- ukb1 %>%
        dplyr::full_join(ukb2, by = "eid") %>%
        dplyr::mutate(state_ukb_1 = dplyr::if_else(is.na(state_ukb_1),0,state_ukb_1),
                      state_ukb_2 = dplyr::if_else(is.na(state_ukb_2),0,state_ukb_2),
                      state_ukb_1_birth = dplyr::if_else(is.na(state_ukb_1_birth),0,state_ukb_1_birth),
                      state_ukb_2_birth = dplyr::if_else(is.na(state_ukb_2_birth),0,state_ukb_2_birth),
                      state_ukb_1_enrol = dplyr::if_else(is.na(state_ukb_1_enrol),0,state_ukb_1_enrol),
                      state_ukb_2_enrol = dplyr::if_else(is.na(state_ukb_2_enrol),0,state_ukb_2_enrol)) %>%
        dplyr::mutate(state_ukb = 1,
                      state_ukb = dplyr::if_else(state_ukb_1 == -1 | state_ukb_2 == -1, -1, state_ukb),
                      state_ukb_birth = 1,
                      state_ukb_birth = dplyr::if_else(state_ukb_1_birth == -1 | state_ukb_2_birth == -1, -1, state_ukb_birth),
                      state_ukb_enrol = 1,
                      state_ukb_enrol = dplyr::if_else(state_ukb_1_enrol == -1 | state_ukb_2_enrol == -1, -1, state_ukb_enrol)) %>%
        dplyr::mutate(age_ukb_birth = pmin(age_ukb_1_birth, age_ukb_2_birth,na.rm = TRUE)) %>%
        dplyr::right_join(pop, by = "eid") %>%
        dplyr::mutate(state_ukb = dplyr::if_else(is.na(state_ukb),0,state_ukb),
                      state_ukb_birth = dplyr::if_else(is.na(state_ukb_birth), 0, state_ukb_birth),
                      state_ukb_enrol = dplyr::if_else(is.na(state_ukb_enrol), 0, state_ukb_enrol)) %>%
        dplyr::select("eid", "state_ukb", "state_ukb_birth", "age_ukb_birth", "state_ukb_enrol")

      return(ukb_data)
    },

    "T2DM" = {
      codes_ukb <- getCodes(outcome)$UKB_CODE2
      ukb1      <- searchUKBcode(ukb, codes_ukb) %>%
        dplyr::rename("state_ukb_1" = "state_ukb",
                      "state_ukb_1_birth" = "state_ukb_birth",
                      "age_ukb_1_birth"   = "age_ukb_birth",
                      "state_ukb_1_enrol" = "state_ukb_enrol")

      ukb2 <-  ukb %>%
        dplyr::select("eid","f.2443.0.0","f.2976.0.0") %>%
        dplyr::filter(!is.na(f.2443.0.0), f.2443.0.0 != -1, f.2443.0.0 != -3, f.2443.0.0 != 0,
                      !is.na(f.2976.0.0), f.2976.0.0 != -1, f.2976.0.0 != -3) %>%
        dplyr::mutate(state_ukb_2 = 1,
                      state_ukb_2_birth = 1,
                      state_ukb_2_enrol = -1) %>%
        dplyr::select("eid", "state_ukb_2", "state_ukb_2_birth", "age_ukb_2_birth" = "f.2976.0.0", "state_ukb_2_enrol") %>%
        dplyr::mutate(state_ukb_2_birth = dplyr::if_else(is.nan(age_ukb_2_birth) | age_ukb_2_birth < 0, -1, state_ukb_2_birth))

      codes_ukb <- getCodes(outcome)$UKB_CODE3
      ukb3      <- searchUKBcode1(ukb, codes_ukb) %>%
        dplyr::rename("state_ukb_3" = "state_ukb",
                      "state_ukb_3_birth" = "state_ukb_birth",
                      "age_ukb_3_birth"   = "age_ukb_birth",
                      "state_ukb_3_enrol" = "state_ukb_enrol")

      ukb4 <- ukb %>%
        dplyr::select("eid","f.30750.0.0") %>%
        dplyr::filter(!is.na(f.30750.0.0)) %>%
        dplyr::filter(f.30750.0.0 >= 48) %>%
        dplyr::mutate(state_ukb_4 = 1,
                      state_ukb_4_birth = -1,
                      age_ukb_4_birth   = 1000,
                      state_ukb_4_enrol = -1) %>%
        dplyr::select("eid", "state_ukb_4", "state_ukb_4_birth", "age_ukb_4_birth", "state_ukb_4_enrol")

      ukb_data <- ukb1 %>%
        dplyr::full_join(ukb2, by = "eid") %>%
        dplyr::full_join(ukb3, by = "eid") %>%
        dplyr::full_join(ukb4, by = "eid") %>%
        dplyr::mutate(state_ukb_1 = dplyr::if_else(is.na(state_ukb_1),0,state_ukb_1),
                      state_ukb_2 = dplyr::if_else(is.na(state_ukb_2),0,state_ukb_2),
                      state_ukb_3 = dplyr::if_else(is.na(state_ukb_3),0,state_ukb_3),
                      state_ukb_4 = dplyr::if_else(is.na(state_ukb_4),0,state_ukb_4),
                      state_ukb_1_birth = dplyr::if_else(is.na(state_ukb_1_birth),0,state_ukb_1_birth),
                      state_ukb_2_birth = dplyr::if_else(is.na(state_ukb_2_birth),0,state_ukb_2_birth),
                      state_ukb_3_birth = dplyr::if_else(is.na(state_ukb_3_birth),0,state_ukb_3_birth),
                      state_ukb_4_birth = dplyr::if_else(is.na(state_ukb_4_birth),0,state_ukb_4_birth),
                      state_ukb_1_enrol = dplyr::if_else(is.na(state_ukb_1_enrol),0,state_ukb_1_enrol),
                      state_ukb_2_enrol = dplyr::if_else(is.na(state_ukb_2_enrol),0,state_ukb_2_enrol),
                      state_ukb_3_enrol = dplyr::if_else(is.na(state_ukb_3_enrol),0,state_ukb_3_enrol),
                      state_ukb_4_enrol = dplyr::if_else(is.na(state_ukb_4_enrol),0,state_ukb_4_enrol)) %>%
        dplyr::mutate(state_ukb = dplyr::if_else(state_ukb_1 == -1 | state_ukb_2 == -1 | state_ukb_3 == -1 | state_ukb_4 == -1,-1,1),
                      state_ukb_birth = dplyr::if_else(state_ukb_1_birth == -1 | state_ukb_2_birth == -1 | state_ukb_3_birth == -1 | state_ukb_4_birth == -1, -1, 1),
                      state_ukb_enrol = dplyr::if_else(state_ukb_1_enrol == -1 | state_ukb_2_enrol == -1 | state_ukb_3_enrol == -1 | state_ukb_4_enrol == -1, -1, 1)) %>%
        dplyr::mutate(age_ukb_birth = pmin(age_ukb_1_birth, age_ukb_2_birth, na.rm = TRUE)) %>%
        dplyr::right_join(pop, by = "eid") %>%
        dplyr::mutate(state_ukb = dplyr::if_else(is.na(state_ukb),0,state_ukb),
                      state_ukb_birth = dplyr::if_else(is.na(state_ukb_birth), 0, state_ukb_birth),
                      state_ukb_enrol = dplyr::if_else(is.na(state_ukb_enrol), 0, state_ukb_enrol)) %>%
        dplyr::select("eid", "state_ukb", "state_ukb_birth", "age_ukb_birth", "state_ukb_enrol")
    }
  )
}



searchUKBcode <- function(ukb, codes_ukb, varname = "f.20002.0", varname1 = "f.20009.0"){
  ukb2 <- NULL
  for (i in 0:33){
    eval(parse(text = paste0("ukb2 <- ukb2 %>% dplyr::union_all(ukb %>%
                                                   dplyr::select(eid, f.20002.0.",i,", age_ukb = f.20009.0.",i,") %>%
                                                   dplyr::filter(f.20002.0.",i," %in% codes_ukb) %>%
                                                   dplyr::select(eid, age_ukb_birth = age_ukb) )")))
  }

  ukb2 %>%
    dplyr::group_by(eid) %>%
    dplyr::mutate(age_ukb_birth = as.numeric(age_ukb_birth)) %>%
    dplyr::summarise(age_ukb_birth = min(age_ukb_birth)) %>%
    dplyr::mutate(age_ukb_birth = round(age_ukb_birth),
                  state_ukb_birth = 1,
                  state_ukb = 1) %>%
    dplyr::mutate(state_ukb_birth = dplyr::if_else(is.na(age_ukb_birth) | age_ukb_birth < 0, -1, state_ukb_birth)) %>%
    dplyr::mutate(state_ukb_enrol = -1)
}


searchUKBcode1 <- function(x, code, varname = "f.20003.0.") {
  x %>%
    dplyr::select("eid",contains(varname)) %>%
    dplyr::filter(
      dplyr::if_any(
        dplyr::contains(varname),
        ~ . %in% code
      )
    ) %>%
    dplyr::mutate(state_ukb = 1,
                  state_ukb_birth = -1,
                  age_ukb_birth   = 1000,
                  state_ukb_enrol = -1) %>%
    dplyr::select("eid", "state_ukb", "state_ukb_birth", "age_ukb_birth", "state_ukb_enrol")
}
