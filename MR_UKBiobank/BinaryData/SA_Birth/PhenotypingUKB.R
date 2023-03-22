PhenotypingUKB <- function(outc,pop,ukb){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  ukb = switch(
    outc,
    "CAD" = {
      data.frame("eid" = pop,
                 "state_ukb" = 0,
                 "age_ukb" = 1000)
      },
    
    "Fracture" = {
      ukb1 <- ukb %>% 
        select("eid",
               "year_assesment" = "X53.0.0",
               "year_of_birth" = "X34.0.0",
               "state_ukb1" = "X3005.0.0") %>%
        filter(state_ukb1 == 1) %>%
        mutate(year_assesment = year(as.Date(year_assesment,"%d/%m/%Y")),
               age_ukb1 = year_assesment - year_of_birth,
               state_ukb1 = as.double(state_ukb1),
               state_ukb1 = if_else(is.na(age_ukb1) | age_ukb1 < 0, -1, state_ukb1)) 
      
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb2 <- NULL
      for (i in 0:33){
        eval(parse(text = paste0("ukb2 <- ukb2 %>% union_all( ukb %>% select(eid,X20002.0.",i,", age_ukb2 = X20009.0.",i,") %>% filter(X20002.0.",i," %in% codes_ukb$UKB.CODE) %>% select(eid, age_ukb2) )")))
      }
      ukb2 <- ukb2 %>%
        group_by(eid) %>%
        summarise(age_ukb2 = min(age_ukb2)) %>%
        mutate(age_ukb2 = round(age_ukb2),
               state_ukb2 = 1) %>%
        mutate(state_ukb2 = if_else(is.na(age_ukb2) | age_ukb2 < 0, -1, state_ukb2))
      
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        mutate(state_ukb1 = if_else(is.na(state_ukb1),0,state_ukb1),
               state_ukb2 = if_else(is.na(state_ukb2),0,state_ukb2)) %>%
        mutate(state_ukb = if_else(state_ukb1 == -1 | state_ukb2 == -1, -1, 1),
               age_ukb   = pmin(age_ukb1,age_ukb2,na.rm = TRUE)) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
      },
    
    "Hypertension" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb2 <- NULL
      for (i in 0:33){
        eval(parse(text = paste0("ukb2 <- ukb2 %>% union_all( ukb %>% select(eid,X20002.0.",i,", age_ukb = X20009.0.",i,") %>% filter(X20002.0.",i," %in% codes_ukb$UKB.CODE) %>% select(eid, age_ukb) )")))
      }
      ukb_data <- ukb2 %>%
        group_by(eid) %>%
        summarise(age_ukb = min(age_ukb)) %>%
        mutate(age_ukb = round(age_ukb),
               state_ukb = 1) %>%
        mutate(state_ukb = if_else(is.na(age_ukb) | age_ukb < 0, -1, state_ukb)) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
      },
    
    
      "IS" = {
        codes_ukb <- phen %>%
          select("UKB.CODE") %>%
          filter(!is.na(UKB.CODE))
        
        ukb2 <- NULL
        for (i in 0:33){
          eval(parse(text = paste0("ukb2 <- ukb2 %>% union_all( ukb %>% select(eid,X20002.0.",i,", age_ukb = X20009.0.",i,") %>% filter(X20002.0.",i," %in% codes_ukb$UKB.CODE) %>% select(eid, age_ukb) )")))
        }
        ukb_data <- ukb2 %>%
          group_by(eid) %>%
          summarise(age_ukb = min(age_ukb)) %>%
          mutate(age_ukb = round(age_ukb),
                 state_ukb = 1) %>%
          mutate(state_ukb = if_else(is.na(age_ukb) | age_ukb < 0, -1, state_ukb)) %>%
          right_join(pop, by = "eid") %>%
          mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
      },
        
    "MI" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      ukb2 <- NULL
      for (i in 0:33){
        eval(parse(text = paste0("ukb2 <- ukb2 %>% union_all( ukb %>% select(eid,X20002.0.",i,", age_ukb1 = X20009.0.",i,") %>% filter(X20002.0.",i," %in% codes_ukb$UKB.CODE) %>% select(eid, age_ukb1) )")))
      }
      ukb1 <- ukb2 %>%
        group_by(eid) %>%
        summarise(age_ukb1 = min(age_ukb1)) %>%
        mutate(age_ukb1 = round(age_ukb1),
               state_ukb1 = 1) %>%
        mutate(state_ukb1 = if_else(is.na(age_ukb1) | age_ukb1 < 0, -1, state_ukb1))
      
      ukb2 <- ukb %>%
        select("eid","age_ukb2" = "X3894.0.0") %>%
        filter(!is.na(age_ukb2), age_ukb2 != -1, age_ukb2 != -3) %>%
        mutate(state_ukb2 = 1,
               state_ukb2 = if_else(is.na(age_ukb2) | age_ukb2 < 0,-1,state_ukb2)) 
      
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        mutate(state_ukb1 = if_else(is.na(state_ukb1),0,state_ukb1),
               state_ukb2 = if_else(is.na(state_ukb2),0,state_ukb2)) %>%
        mutate(state_ukb = 1,
               state_ukb = if_else(state_ukb1 == -1 | state_ukb2 == -1, -1, state_ukb),
               age_ukb   = pmin(age_ukb1,age_ukb2,na.rm = TRUE)) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))
      },
    
    "T2DM" = {
      codes_ukb <- phen %>%
        select("UKB.CODE2") %>%
        filter(!is.na(UKB.CODE2))
      
      ukb1 <- NULL
      for (i in 0:33){
        eval(parse(text = paste0("ukb1 <- ukb1 %>% union_all(ukb1 %>% select(eid,X20002.0.",i,", age_ukb1 = X20009.0.",i,") %>% filter(X20002.0.",i," %in% codes_ukb$UKB.CODE) %>% select(eid, age_ukb1) )")))
      }
      ukb1 <- ukb1 %>%
        group_by(eid) %>%
        summarise(age_ukb1 = min(age_ukb1)) %>%
        mutate(age_ukb1 = round(age_ukb1),
               state_ukb1 = 1) %>%
        mutate(state_ukb1 = if_else(is.na(age_ukb1) | age_ukb1 < 0, -1, state_ukb1)) 
      
      ukb2 <-  ukb %>%
        select("eid","X2443.0.0","X2976.0.0") %>%
        filter(!is.na(X2443.0.0), X2443.0.0 != -1, X2443.0.0 != -3, X2443.0.0 != 0,
               !is.na(X2976.0.0), X2976.0.0 != -1, X2976.0.0 != -3) %>%
        mutate(state_ukb2 = 1) %>%
        select(eid,state_ukb2, age_ukb2 = X2976.0.0) %>%
        mutate(state_ukb2 = if_else(is.nan(age_ukb2) | age_ukb2 < 0, -1, state_ukb2))
      
      codes_ukb <- phen %>%
        select("UKB.CODE3") %>%
        filter(!is.na(UKB.CODE3))
      
      ukb3 <- searchUKBcode(UKB, codes_ukb$UKB.CODE3) %>%
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
        mutate(state_ukb1 = if_else(is.na(state_ukb1),0,state_ukb1),
               state_ukb2 = if_else(is.na(state_ukb2),0,state_ukb2),
               state_ukb3 = if_else(is.na(state_ukb3),0,state_ukb3),
               state_ukb4 = if_else(is.na(state_ukb4),0,state_ukb4)) %>%
        mutate(state_ukb = if_else(state_ukb1 == -1 | state_ukb2 == -1 | state_ukb3 == -1 | state_ukb4 == -1,-1,1),
               age_ukb   = pmin(age_ukb1,age_ukb2,na.rm = TRUE)) %>%
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
