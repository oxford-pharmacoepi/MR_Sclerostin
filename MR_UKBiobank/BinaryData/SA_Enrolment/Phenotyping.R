# ============================================================================ #
#                                     Phenotyping                              #
#                          (Survival analysis - Age of enrolment)              #
# See:                                                                         #
# here("MR_UKBiobank","BinaryData","Phenotyping.xlsx"                          #
# for more information about the codes.                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","pathData","tok","outputFolder")))
ukb  <- as_tibble(read.delim(paste0(pathData,"UKB\\ukb669864_birth.csv"), sep = ","))
pop <- read_delim(paste0(pathData,"UKB\\cohort.csv")) %>% select('eid')
outc <- c('Fracture','CAD','MI','IS','Hypertension','T2DM')
for (i in 1:6){
  source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingHES.R"))
  hes_data <- PhenotypingHes(outc[i],hes,hesD,pop)
  
  source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingGP.R"))
  gp_data <- PhenotypingGP(outc[i],gp,pop,ukb)
  
  source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingUKB.R"))
  ukb_data <- PhenotypingUKB(outc[i],pop,ukb)
  
  source(here("MR_UKBiobank","BinaryData","SA_Enrolment","MergeTables.R"))
  t <- MergeTables(ukb_data,hes_data,gp_data)
  
  write.csv(t,paste0(pathData,"MR_UKBiobank\\BinaryData\\SA_Enrolment\\Phenotype_",outc[i],".csv"))
}



