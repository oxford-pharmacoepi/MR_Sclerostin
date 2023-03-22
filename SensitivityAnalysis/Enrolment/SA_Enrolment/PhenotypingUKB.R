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
      
      ukb2 <- ukb %>%
        select("eid",contains("X20002.0."),contains("X20009.0.")) %>%
        filter(X20002.0.0 %in% codes_ukb$UKB.CODE | X20002.0.1 %in% codes_ukb$UKB.CODE | X20002.0.2 %in% codes_ukb$UKB.CODE |
                 X20002.0.3 %in% codes_ukb$UKB.CODE | X20002.0.4 %in% codes_ukb$UKB.CODE | X20002.0.5 %in% codes_ukb$UKB.CODE |
                 X20002.0.6 %in% codes_ukb$UKB.CODE | X20002.0.7 %in% codes_ukb$UKB.CODE | X20002.0.8 %in% codes_ukb$UKB.CODE |
                 X20002.0.9 %in% codes_ukb$UKB.CODE | X20002.0.10 %in% codes_ukb$UKB.CODE | X20002.0.11 %in% codes_ukb$UKB.CODE |
                 X20002.0.12 %in% codes_ukb$UKB.CODE | X20002.0.13 %in% codes_ukb$UKB.CODE | X20002.0.14 %in% codes_ukb$UKB.CODE |
                 X20002.0.15 %in% codes_ukb$UKB.CODE | X20002.0.16 %in% codes_ukb$UKB.CODE | X20002.0.17 %in% codes_ukb$UKB.CODE |
                 X20002.0.18 %in% codes_ukb$UKB.CODE | X20002.0.19 %in% codes_ukb$UKB.CODE | X20002.0.20 %in% codes_ukb$UKB.CODE |
                 X20002.0.21 %in% codes_ukb$UKB.CODE | X20002.0.22 %in% codes_ukb$UKB.CODE | X20002.0.23 %in% codes_ukb$UKB.CODE |
                 X20002.0.24 %in% codes_ukb$UKB.CODE | X20002.0.25 %in% codes_ukb$UKB.CODE | X20002.0.26 %in% codes_ukb$UKB.CODE |
                 X20002.0.27 %in% codes_ukb$UKB.CODE | X20002.0.28 %in% codes_ukb$UKB.CODE | X20002.0.29 %in% codes_ukb$UKB.CODE |
                 X20002.0.30 %in% codes_ukb$UKB.CODE | X20002.0.31 %in% codes_ukb$UKB.CODE | X20002.0.32 %in% codes_ukb$UKB.CODE |
                 X20002.0.33 %in% codes_ukb$UKB.CODE) %>%
        mutate(state_ukb2 = -1) 
    
      ukb_data <- ukb1 %>%
        full_join(ukb2, by = "eid") %>%
        mutate(state_ukb = -1) %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "Hypertension" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb_data <- ukb %>%
        select("eid",contains("X20002.0."),contains("X20009.0.")) %>%
        filter(X20002.0.0 %in% codes_ukb$UKB.CODE | X20002.0.1 %in% codes_ukb$UKB.CODE | X20002.0.2 %in% codes_ukb$UKB.CODE |
                 X20002.0.3 %in% codes_ukb$UKB.CODE | X20002.0.4 %in% codes_ukb$UKB.CODE | X20002.0.5 %in% codes_ukb$UKB.CODE |
                 X20002.0.6 %in% codes_ukb$UKB.CODE | X20002.0.7 %in% codes_ukb$UKB.CODE | X20002.0.8 %in% codes_ukb$UKB.CODE |
                 X20002.0.9 %in% codes_ukb$UKB.CODE | X20002.0.10 %in% codes_ukb$UKB.CODE | X20002.0.11 %in% codes_ukb$UKB.CODE |
                 X20002.0.12 %in% codes_ukb$UKB.CODE | X20002.0.13 %in% codes_ukb$UKB.CODE | X20002.0.14 %in% codes_ukb$UKB.CODE |
                 X20002.0.15 %in% codes_ukb$UKB.CODE | X20002.0.16 %in% codes_ukb$UKB.CODE | X20002.0.17 %in% codes_ukb$UKB.CODE |
                 X20002.0.18 %in% codes_ukb$UKB.CODE | X20002.0.19 %in% codes_ukb$UKB.CODE | X20002.0.20 %in% codes_ukb$UKB.CODE |
                 X20002.0.21 %in% codes_ukb$UKB.CODE | X20002.0.22 %in% codes_ukb$UKB.CODE | X20002.0.23 %in% codes_ukb$UKB.CODE |
                 X20002.0.24 %in% codes_ukb$UKB.CODE | X20002.0.25 %in% codes_ukb$UKB.CODE | X20002.0.26 %in% codes_ukb$UKB.CODE |
                 X20002.0.27 %in% codes_ukb$UKB.CODE | X20002.0.28 %in% codes_ukb$UKB.CODE | X20002.0.29 %in% codes_ukb$UKB.CODE |
                 X20002.0.30 %in% codes_ukb$UKB.CODE | X20002.0.31 %in% codes_ukb$UKB.CODE | X20002.0.32 %in% codes_ukb$UKB.CODE |
                 X20002.0.33 %in% codes_ukb$UKB.CODE) %>%
        mutate(state_ukb = -1) %>%
        select("eid","state_ukb") %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    "IS" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb <- ukb %>%
        select("eid",contains("X20002.0."),contains("X20009.0.")) %>%
        filter(X20002.0.0 %in% codes_ukb$UKB.CODE | X20002.0.1 %in% codes_ukb$UKB.CODE | X20002.0.2 %in% codes_ukb$UKB.CODE |
                 X20002.0.3 %in% codes_ukb$UKB.CODE | X20002.0.4 %in% codes_ukb$UKB.CODE | X20002.0.5 %in% codes_ukb$UKB.CODE |
                 X20002.0.6 %in% codes_ukb$UKB.CODE | X20002.0.7 %in% codes_ukb$UKB.CODE | X20002.0.8 %in% codes_ukb$UKB.CODE |
                 X20002.0.9 %in% codes_ukb$UKB.CODE | X20002.0.10 %in% codes_ukb$UKB.CODE | X20002.0.11 %in% codes_ukb$UKB.CODE |
                 X20002.0.12 %in% codes_ukb$UKB.CODE | X20002.0.13 %in% codes_ukb$UKB.CODE | X20002.0.14 %in% codes_ukb$UKB.CODE |
                 X20002.0.15 %in% codes_ukb$UKB.CODE | X20002.0.16 %in% codes_ukb$UKB.CODE | X20002.0.17 %in% codes_ukb$UKB.CODE |
                 X20002.0.18 %in% codes_ukb$UKB.CODE | X20002.0.19 %in% codes_ukb$UKB.CODE | X20002.0.20 %in% codes_ukb$UKB.CODE |
                 X20002.0.21 %in% codes_ukb$UKB.CODE | X20002.0.22 %in% codes_ukb$UKB.CODE | X20002.0.23 %in% codes_ukb$UKB.CODE |
                 X20002.0.24 %in% codes_ukb$UKB.CODE | X20002.0.25 %in% codes_ukb$UKB.CODE | X20002.0.26 %in% codes_ukb$UKB.CODE |
                 X20002.0.27 %in% codes_ukb$UKB.CODE | X20002.0.28 %in% codes_ukb$UKB.CODE | X20002.0.29 %in% codes_ukb$UKB.CODE |
                 X20002.0.30 %in% codes_ukb$UKB.CODE | X20002.0.31 %in% codes_ukb$UKB.CODE | X20002.0.32 %in% codes_ukb$UKB.CODE |
                 X20002.0.33 %in% codes_ukb$UKB.CODE) %>%
        mutate_at(vars(starts_with("X20009.0.")),funs(ifelse(. == -1, 1000,.))) %>%
        mutate(state_ukb = -1) %>%
        select("eid","state_ukb") %>%
        right_join(pop, by = "eid") %>%
        mutate(state_ukb = if_else(is.na(state_ukb),0,state_ukb))},
    
    
    "MI" = {
      codes_ukb <- phen %>%
        select("UKB.CODE") %>%
        filter(!is.na(UKB.CODE))
      
      ukb1 <- ukb %>%
        select("eid",contains("X20002.0."),contains("X20009.0.")) %>%
        filter(X20002.0.0 %in% codes_ukb$UKB.CODE | X20002.0.1 %in% codes_ukb$UKB.CODE | X20002.0.2 %in% codes_ukb$UKB.CODE |
                 X20002.0.3 %in% codes_ukb$UKB.CODE | X20002.0.4 %in% codes_ukb$UKB.CODE | X20002.0.5 %in% codes_ukb$UKB.CODE |
                 X20002.0.6 %in% codes_ukb$UKB.CODE | X20002.0.7 %in% codes_ukb$UKB.CODE | X20002.0.8 %in% codes_ukb$UKB.CODE |
                 X20002.0.9 %in% codes_ukb$UKB.CODE | X20002.0.10 %in% codes_ukb$UKB.CODE | X20002.0.11 %in% codes_ukb$UKB.CODE |
                 X20002.0.12 %in% codes_ukb$UKB.CODE | X20002.0.13 %in% codes_ukb$UKB.CODE | X20002.0.14 %in% codes_ukb$UKB.CODE |
                 X20002.0.15 %in% codes_ukb$UKB.CODE | X20002.0.16 %in% codes_ukb$UKB.CODE | X20002.0.17 %in% codes_ukb$UKB.CODE |
                 X20002.0.18 %in% codes_ukb$UKB.CODE | X20002.0.19 %in% codes_ukb$UKB.CODE | X20002.0.20 %in% codes_ukb$UKB.CODE |
                 X20002.0.21 %in% codes_ukb$UKB.CODE | X20002.0.22 %in% codes_ukb$UKB.CODE | X20002.0.23 %in% codes_ukb$UKB.CODE |
                 X20002.0.24 %in% codes_ukb$UKB.CODE | X20002.0.25 %in% codes_ukb$UKB.CODE | X20002.0.26 %in% codes_ukb$UKB.CODE |
                 X20002.0.27 %in% codes_ukb$UKB.CODE | X20002.0.28 %in% codes_ukb$UKB.CODE | X20002.0.29 %in% codes_ukb$UKB.CODE |
                 X20002.0.30 %in% codes_ukb$UKB.CODE | X20002.0.31 %in% codes_ukb$UKB.CODE | X20002.0.32 %in% codes_ukb$UKB.CODE |
                 X20002.0.33 %in% codes_ukb$UKB.CODE) %>%
        mutate(state_ukb1 = -1) %>%
        select("eid","state_ukb1")
      
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
      
      ukb1 <-  ukb %>%
        select("eid",contains("X20002.0."),contains("X20009.0.")) %>%
        filter(X20002.0.0 %in% codes_ukb$UKB.CODE2 | X20002.0.1 %in% codes_ukb$UKB.CODE2 | X20002.0.2 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.3 %in% codes_ukb$UKB.CODE2 | X20002.0.4 %in% codes_ukb$UKB.CODE2 | X20002.0.5 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.6 %in% codes_ukb$UKB.CODE2 | X20002.0.7 %in% codes_ukb$UKB.CODE2 | X20002.0.8 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.9 %in% codes_ukb$UKB.CODE2 | X20002.0.10 %in% codes_ukb$UKB.CODE2 | X20002.0.11 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.12 %in% codes_ukb$UKB.CODE2 | X20002.0.13 %in% codes_ukb$UKB.CODE2 | X20002.0.14 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.15 %in% codes_ukb$UKB.CODE2 | X20002.0.16 %in% codes_ukb$UKB.CODE2 | X20002.0.17 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.18 %in% codes_ukb$UKB.CODE2 | X20002.0.19 %in% codes_ukb$UKB.CODE2 | X20002.0.20 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.21 %in% codes_ukb$UKB.CODE2 | X20002.0.22 %in% codes_ukb$UKB.CODE2 | X20002.0.23 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.24 %in% codes_ukb$UKB.CODE2 | X20002.0.25 %in% codes_ukb$UKB.CODE2 | X20002.0.26 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.27 %in% codes_ukb$UKB.CODE2 | X20002.0.28 %in% codes_ukb$UKB.CODE2 | X20002.0.29 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.30 %in% codes_ukb$UKB.CODE2 | X20002.0.31 %in% codes_ukb$UKB.CODE2 | X20002.0.32 %in% codes_ukb$UKB.CODE2 |
                 X20002.0.33 %in% codes_ukb$UKB.CODE2) %>%
        mutate(state_ukb1 = -1) %>%
        select("eid","state_ukb1") 

      ukb2 <-  ukb %>%
        select("eid","X2443.0.0","X2976.0.0") %>%
        filter(!is.na(X2443.0.0), X2443.0.0 != -1, X2443.0.0 != -3, X2443.0.0 != 0,
               !is.na(X2976.0.0), X2976.0.0 != -1, X2976.0.0 != -3) %>%
        mutate(state_ukb2 = -1) %>%
        select(eid,state_ukb2)

      codes_ukb <- phen %>%
        select("UKB.CODE3") %>%
        filter(!is.na(UKB.CODE3))
      
      ukb3 <- ukb %>%
        select("eid",contains("X20003.0.")) %>%
        filter(X20003.0.0 %in% codes_ukb$UKB.CODE3 | X20003.0.1 %in% codes_ukb$UKB.CODE3 | X20003.0.2 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.3 %in% codes_ukb$UKB.CODE3 | X20003.0.4 %in% codes_ukb$UKB.CODE3 | X20003.0.5 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.6 %in% codes_ukb$UKB.CODE3 | X20003.0.7 %in% codes_ukb$UKB.CODE3 | X20003.0.8 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.9 %in% codes_ukb$UKB.CODE3 | X20003.0.10 %in% codes_ukb$UKB.CODE3 | X20003.0.11 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.12 %in% codes_ukb$UKB.CODE3 | X20003.0.13 %in% codes_ukb$UKB.CODE3 | X20003.0.14 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.15 %in% codes_ukb$UKB.CODE3 | X20003.0.16 %in% codes_ukb$UKB.CODE3 | X20003.0.17 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.18 %in% codes_ukb$UKB.CODE3 | X20003.0.19 %in% codes_ukb$UKB.CODE3 | X20003.0.20 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.21 %in% codes_ukb$UKB.CODE3 | X20003.0.22 %in% codes_ukb$UKB.CODE3 | X20003.0.23 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.24 %in% codes_ukb$UKB.CODE3 | X20003.0.25 %in% codes_ukb$UKB.CODE3 | X20003.0.26 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.27 %in% codes_ukb$UKB.CODE3 | X20003.0.28 %in% codes_ukb$UKB.CODE3 | X20003.0.29 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.30 %in% codes_ukb$UKB.CODE3 | X20003.0.31 %in% codes_ukb$UKB.CODE3 | X20003.0.32 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.33 %in% codes_ukb$UKB.CODE3 | X20003.0.34 %in% codes_ukb$UKB.CODE3 | X20003.0.35 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.36 %in% codes_ukb$UKB.CODE3 | X20003.0.37 %in% codes_ukb$UKB.CODE3 | X20003.0.38 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.39 %in% codes_ukb$UKB.CODE3 | X20003.0.40 %in% codes_ukb$UKB.CODE3 | X20003.0.41 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.42 %in% codes_ukb$UKB.CODE3 | X20003.0.43 %in% codes_ukb$UKB.CODE3 | X20003.0.44 %in% codes_ukb$UKB.CODE3 |
                 X20003.0.45 %in% codes_ukb$UKB.CODE3 | X20003.0.46 %in% codes_ukb$UKB.CODE3 | X20003.0.47 %in% codes_ukb$UKB.CODE3) %>%
        mutate(state_ukb3 = -1) %>%
        select(eid,state_ukb3)
      
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