PhenotypingGP <- function(outc,gp,pop){
  phen <- read.xlsx(here("MR_UKBiobank","BinaryData","PhenotypingR.xlsx"), sheetName = outc)
  
  gp_codes <- phen %>% select("GP") %>%
    filter(!is.na(GP))
  
  g <- gp %>% filter(read_2 %in% gp_codes$GP | read_3 %in% gp_codes$GP) %>% 
    select(eid) %>% 
    distinct() %>%
    mutate(state_gp = 1) %>%
    right_join(pop, 
               by = "eid") %>%
    mutate(state_gp = if_else(is.na(state_gp),0,state_gp))
  
  return(g)
}
