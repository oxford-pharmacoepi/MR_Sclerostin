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
                             year_of_birth = "X34.0.0",
                             date_assessment = "X53.0.0"),
              by = "eid") %>%
    distinct() %>%
    mutate(date_assessment = as.Date(date_assessment,"%d/%m/%Y")) %>%
    group_by(eid) %>%
    summarise(age_hes = year(first_fracture) - year(date_assessment)) %>% # Appending results on variable age_hes 
    mutate(state_hes = 1) %>%
    right_join(pop, by = "eid") %>% # All hes patients
    left_join(ukb %>% select(eid,
                             date_assessment = "X53.0.0",
                             year_of_birth = "X34.0.0",
                             age_of_death = "X40007.0.0"),
              by = "eid") %>%
    distinct() %>%
    mutate(date_assessment = as.Date(date_assessment,"%d/%m/%Y")) %>%
    mutate(age_hes = if_else(is.na(age_hes) & !is.na(age_of_death), year_of_birth + round(age_of_death) - year(date_assessment), age_hes), # Death
           age_hes = if_else(is.na(age_hes), 2021 - year(date_assessment), age_hes), # Alive
           state_hes = if_else(is.na(state_hes),0,state_hes),
           state_hes = if_else(is.na(age_hes) | age_hes < 0, -1, state_hes)) %>% # State
    select(eid, age_hes, state_hes)
  
  return(h)
}