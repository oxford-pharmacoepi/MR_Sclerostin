# ============================================================================ #
#                                     Phenotyping                              #
#                               (Logistic regression)                          #
#                            2023 - Marta Alcalde-Herraiz                      #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Phenotyping the binary outcomes using UK Biobank patient level data. The     #
# outcomes will be later used to perform a logistic regression.                #
# See:                                                                         #
# here("MR_UKBiobank","BinaryData","Phenotyping.xlsx"                          #
# for more information about the codes.                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok")))
ukb  <- as_tibble(read.delim(paste0(pathData,"UKB\\ukb669864_logistic.csv"), sep = ","))
pop <- read_delim(paste0(pathData,'UKB\\cohort.csv')) %>% select('eid')
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')
for (i in 1:length(outc)){
  if (i != 7){
    source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingHES.R"))
    hes_data <- PhenotypingHes(outc[i],hes,pop)
    
    source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingGP.R"))
    gp_data <- PhenotypingGP(outc[i],gp,pop)
  }
  
  source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingUKB.R"))
  ukb_data <- PhenotypingUKB(outc[i],pop,ukb)
  
  if (i != 7){
    source(here("MR_UKBiobank","BinaryData","Logistic","MergeTables.R"))
    t <- MergeTables(ukb_data,hes_data,gp_data)
  }else{
    t <- ukb_data
  }
  # t %>% group_by(state_gp, state_hes, state_ukb) %>% tally()
  write.csv(t,paste0(pathData,"MR_UKBiobank\\BinaryData\\Logistic\\Phenotype_",outc[i],".csv"))
}








