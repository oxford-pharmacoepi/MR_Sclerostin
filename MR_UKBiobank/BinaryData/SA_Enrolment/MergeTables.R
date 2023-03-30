MergeTables <- function(ukb_data,hes_data,gp_data){
  t <- hes_data %>% 
    left_join(gp_data, by = "eid") %>%
    left_join(ukb_data %>%
                select("eid","state_ukb"), by = "eid") %>%
    filter(state_hes == 1 | state_hes == 0,
           state_gp  == 1 | state_gp  == 0,
           state_ukb == 0) %>%
    mutate(
      state = if_else(state_hes == 1 | state_gp == 1,1,0),
      age   = if_else(state == 1, pmin(age_hes, age_gp,  na.rm = TRUE), age_hes)
    ) %>%
    distinct()
  
  return(t)
}