PhenotypingHes <- function(outc,hes,hesD,pop){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  hes_codes <- phen %>% 
    select("HES") %>%
    filter(!is.na(HES))
  
  h <- hes %>%
    select("eid","ins_index","diag_icd10") %>%
    filter(diag_icd10 %in% hes_codes$HES) %>% # Filter patients having a code
    inner_join(hesD %>%
                 select("eid","ins_index","epistart","admidate"),
               by = c("eid","ins_index")) %>%
    mutate(epistart = as.Date(epistart, "%d/%m/%Y")) %>%
    mutate(admidate = as.Date(admidate, "%d/%m/%Y")) %>%
    mutate(event_date = pmin(epistart,admidate, na.rm = TRUE)) %>% 
    group_by(eid) %>%
    summarise(first_fracture = min(event_date, na.rm = TRUE)) %>% # Year of the first fracture
    left_join(ukb %>% select(eid,
                             year_of_birth = "X34.0.0"),
              by = "eid") %>%
    mutate(age_hes = year(first_fracture) - year_of_birth) %>% # Appending results on variable age_hes 
    mutate(state_hes = 1) %>%
    mutate(state_hes = if_else(is.na(age_hes) | age_hes < 0, -1, state_hes)) %>%
    select(eid,state_hes,age_hes) %>%
    right_join(pop, by = "eid") %>% # All hes patients
    left_join(ukb %>% select(eid,
                             year_of_birth = "X34.0.0",
                             age_of_death = "X40007.0.0"),
              by = "eid") %>%
    mutate(age_hes = if_else(is.na(age_hes) & !is.na(age_of_death), round(age_of_death),age_hes), # Death
           age_hes = if_else(is.na(age_hes), 2021 - year_of_birth, age_hes), # Alive
           state_hes = if_else(is.na(state_hes),0,state_hes)) %>% # State
    select(eid, age_hes, state_hes)
  
  return(h)
}