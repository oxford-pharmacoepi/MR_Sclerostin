MergeTables <- function(ukb_data,hes_data,gp_data){
  nuk <- nrow(ukb_data)
  nhs <- nrow(hes_data)
  ngp <- nrow(gp_data)
  if (nrow(ukb_data))
  t <- hes_data %>% 
    left_join(gp_data, by = "eid") %>%
    left_join(ukb_data %>%
                select("eid","state_ukb"), by = "eid") %>%
    mutate(
      state = if_else(state_hes == 1 | state_gp == 1 | state_ukb == 1,1,0)
    )
  nt <- nrow(t)
  if (nt != nuk) {stop("different number rows")}
  if (nt != nhs) {stop("different number rows")}
  if (nt != ngp) {stop("different number rows")}
  
  return(t)
}