PhenotypingHes <- function(outc,hes){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"),
                    sheetName = outc)
  
  hes_codes <- phen %>% 
    select("HES") %>%
    filter(!is.na(HES))
  
  h <- hes %>%
    select("eid","diag_icd10") %>%
    filter(diag_icd10 %in% hes_codes$HES) %>% # Filter patients having a code
    select(eid) %>%
    distinct() %>% # No repeated patients
    mutate(state_hes = 1) %>%
    right_join(
      hes %>% 
        select("eid") %>%
        distinct(), # Add those patients without outcome
      by = "eid"
    ) %>%
    mutate(state_hes = if_else(is.na(state_hes),0,state_hes)) %>%
    select(eid,state_hes)
  
  return(h)
}

