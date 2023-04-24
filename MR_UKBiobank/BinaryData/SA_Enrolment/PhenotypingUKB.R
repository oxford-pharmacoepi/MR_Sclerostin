PhenotypingUKB <- function(outc,pop,ukb){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  ukb = switch(
    outc,
    "CAD" = {
      data.frame("eid" = pop,
                 "state_ukb" = 0)},
    
    "Fracture" = {
      ukb1 <- ukb %>% 
        select("eid",
               "state_ukb1" = "X3005.0.0") %>%
        filter(state_ukb1 == -1) 
      
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb2 <- searchUKBcode(ukb, codes_ukb$UKB.CODE, varname = "X20002.0.") %>%
        rename(state_ukb2 = state_ukb)
    
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        mutate(state_ukb = -1) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "Hypertension" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb_data <- searchUKBcode(ukb, codes_ukb$UKB.CODE, varname = "X20002.0.") %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "IS" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb <- searchUKBcode(ukb, codes_ukb$UKB.CODE, varname = "X20002.0.") %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "MI" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb1 <- searchUKBcode(ukb, codes_ukb$UKB.CODE, varname = "X20002.0.") %>%
        rename(state_ukb1 = state_ukb)
      
      ukb2 <- ukb %>%
        select("eid","age_ukb2" = "X3894.0.0") %>%
        filter(!is.na(age_ukb2), age_ukb2 != -1, age_ukb2 != -3) %>%
        mutate(state_ukb2 = -1) 
      
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        mutate(state_ukb = -1) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "T2DM" = {
      codes_ukb <- phen %>%
        select("UKB.CODE2") %>%
        filter(!is.na(UKB.CODE2))
      
      ukb1 <-  searchUKBcode(ukb, codes_ukb$UKB.CODE2, varname = "X20002.0.") %>%
        rename("state_ukb1" = "state_ukb") 

      ukb2 <-  ukb %>%
        select("eid","X2443.0.0","X2976.0.0") %>%
        filter(!is.na(X2443.0.0), X2443.0.0 != -1, X2443.0.0 != -3, X2443.0.0 != 0,
               !is.na(X2976.0.0), X2976.0.0 != -1, X2976.0.0 != -3) %>%
        mutate(state_ukb2 = -1) %>%
        select(eid,state_ukb2)

      codes_ukb <- phen %>%
        select("UKB.CODE3") %>%
        filter(!is.na(UKB.CODE3))
      
      ukb3 <- searchUKBcode(ukb, codes_ukb$UKB.CODE, varname = "X20003.0.") %>%
        rename(state_ukb3 = state_ukb)
      
      ukb4 <- ukb %>%
        select("eid","X30750.0.0") %>%
        filter(!is.na(X30750.0.0)) %>%
        mutate(state_ukb4 = if_else(X30750.0.0 >= 48,-1,0)) %>%
        filter(state_ukb4 == -1) %>%
        select(eid,state_ukb4)
      
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        full_join(ukb3, by = "eid") %>%
        full_join(ukb4, by = "eid") %>%
        mutate(state_ukb = -1) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
    }
    
  )
}

searchUKBcode <- function(x, code, varname = "X20003.0.") {
  x %>%
    select("eid",contains(varname)) %>%
    filter(
      if_any(
        contains(varname),
        ~ . %in% code
      )
    ) %>%
    mutate(state_ukb = -1) %>%
    select(eid,state_ukb)
}
