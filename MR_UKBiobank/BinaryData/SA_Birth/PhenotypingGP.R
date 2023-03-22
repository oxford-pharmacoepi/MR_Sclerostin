PhenotypingGP <- function(outc,gp,pop,ukb){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  gp_codes <- phen %>%
    select("GP") %>%
    filter(!is.na(GP))
  
  g <- gp %>% filter(read_2 %in% gp_codes$GP | read_3 %in% gp_codes$GP) %>% 
    distinct() %>% 
    group_by(eid) %>%
    summarise(event_date = as.Date(event_dt,"%d/%m/%Y")) %>% # Age of the event
    mutate(first_fracture = min(event_date, na.rm = TRUE)) %>%
    left_join(ukb %>%
                select("eid","year_of_birth" = "X34.0.0"),
              by = "eid") %>%
    summarise(age_gp = year(first_fracture) - year_of_birth) %>%
    mutate(state_gp = 1) %>%
    mutate(state_gp = if_else(is.na(age_gp) | age_gp < 0, -1, state_gp)) %>%
    right_join(pop, by = "eid") %>%
    mutate(state_gp = if_else(is.na(state_gp),0,state_gp)) %>%
    select("eid","state_gp","age_gp")
  
  return(g)
}