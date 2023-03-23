# ============================================================================ #
#                                 Phenotyping                                  #
#                       (Survival analysis - Age of birth)                     #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Phenotyping the binary outcomes using UK Biobank patient level data. The     #
# outcomes will be later used to perform a cox regression.                     #
# See:                                                                         #
# here("MR_UKBiobank","BinaryData","Phenotyping.xlsx"                          #
# for more information about the codes.                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok")))
ukb  <- as_tibble(read.delim(paste0(pathData, "UKB\\ukb669864_birth.csv"), sep = ","))

pop <- hes %>% select(eid) %>% distinct()
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')

for (i in c(1:6)){
  source(here("MR_UKBiobank","BinaryData","SA_Birth","PhenotypingHES.R"))
  hes_data <- PhenotypingHes(outc[i],hes,hesD)
  
  source(here("MR_UKBiobank","BinaryData","SA_Birth","PhenotypingGP.R"))
  gp_data <- PhenotypingGP(outc[i],gp,pop,ukb)
  
  source(here("MR_UKBiobank","BinaryData","SA_Birth","PhenotypingUKB.R"))
  ukb_data <- PhenotypingUKB(outc[i],pop,ukb)
  
  source(here("MR_UKBiobank","BinaryData","SA_Birth","MergeTables.R"))
  t <- MergeTables(ukb_data,hes_data,gp_data)
  
  # t %>% group_by(state_gp, state_hes, state_ukb) %>% tally()
  write.csv(t,here("MR_UKBiobank","BinaryData","SA_Birth",paste0("Phenotype_",outc[i],".csv")))
}
