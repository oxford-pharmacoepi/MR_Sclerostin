
PhenotypingUKB <- function(outc,pop,ukb){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  ukb = switch(
    outc,
    "CAD" = {
      data.frame("eid" = pop,
                       "state_ukb" = 0)},
    "Fracture" = {
      ukb1 <- ukb %>% 
        select("eid","state_ukb" = "X3005.0.0") %>%
        filter(state_ukb == 1)
      
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb2 <- searchUKBcode(ukb, codes_ukb$UKB.CODE)
        
      ukb_data <- ukb1 %>%
        union(ukb2) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
    },
    
    "Hypertension" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb_data <- searchUKBcode(ukb, codes_ukb$UKB.CODE) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
      
    "IS" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb_data <- searchUKBcode(ukb, codes_ukb$UKB.CODE) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "MI" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb1 <- searchUKBcode(ukb, codes_ukb$UKB.CODE)
      
      ukb2 <- ukb %>%
        select("eid","X3894.0.0") %>%
        filter(!is.na(X3894.0.0), X3894.0.0 != -1, X3894.0.0 != -3) %>%
        mutate(state_ukb = 1) %>%
        select(eid,state_ukb)
      
      ukb_data <- ukb1 %>%
        union(ukb2) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "T2DM" = {
      codes_ukb <- phen %>%
        select("UKB.CODE2") %>%
        filter(!is.na(UKB.CODE2))
      
      ukb1 <- searchUKBcode(ukb, codes_ukb$UKB.CODE2)
      
      codes_ukb <- phen %>%
        select("UKB.CODE3") %>%
        filter(!is.na(UKB.CODE3))
      
      ukb2 <- searchUKBcode(ukb, codes_ukb$UKB.CODE3, varname = "X20003.0.")
      
      ukb3 <- ukb %>%
        select("eid","X30750.0.0") %>%
        filter(!is.na(X30750.0.0)) %>%
        mutate(state_ukb = if_else(X30750.0.0 >= 48,1,0)) %>%
        filter(state_ukb == 1) %>%
        select(eid,state_ukb)

      ukb4 <- ukb %>%
        select("eid","state_ukb" = "X2443.0.0") %>%
        filter(state_ukb == 1) %>%
        select(eid,state_ukb)
      
      ukb_data <- ukb1 %>%
        union(ukb2) %>%
        union(ukb3) %>%
        union(ukb4) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
    }
  )
  
  return(ukb)
}

searchUKBcode <- function(x, code, varname = "X20002.0.") {
  x %>%
    select("eid",contains(varname)) %>%
    filter(
      if_any(
        contains(varname),
        ~ . %in% code
      )
    ) %>%
    mutate(state_ukb = 1) %>%
    select(eid,state_ukb)
}
