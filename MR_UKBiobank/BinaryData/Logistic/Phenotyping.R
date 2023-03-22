# ============================================================================ #
#                                     Phenotyping                             #
#                               (Logistic regression)                          #
# See:                                                                         #
# here("MR_UKBiobank","BinaryData","Phenotyping.xlsx"                          #
# for more information about the codes.                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("hes","hesD","gp","ukb")))
library(pacman)
pacman::p_load('tibble','dplyr','readxl','lubridate')

# hes  <- as_tibble(read.delim(here("Data","hesin_diag.txt"),sep = "\t", quote = ""))
# hesD <- as_tibble(read.delim(here("Data","hesin.txt"), sep = "\t", quote = ""))
# gp   <- as_tibble(read.delim(here("Data","gp_clinical.txt"),sep = "\t", quote = ""))
# ukb  <- as_tibble(read.delim(here("Data","UKB","ukb669864_logistic.csv"), sep = ","))

outc <- "MI"

pop <- hes %>% select(eid) %>% distinct()

source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingHES.R"))
hes_data <- PhenotypingHes(outc,hes)

source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingGP.R"))
gp_data <- PhenotypingGP(outc,gp,pop)

source(here("MR_UKBiobank","BinaryData","Logistic","PhenotypingUKB.R"))
ukb_data <- PhenotypingUKB(outc,pop,ukb)

source(here("MR_UKBiobank","BinaryData","Logistic","MergeTables.R"))
t <- MergeTables(ukb_data,hes_data,gp_data)

t %>% group_by(state_gp, state_hes, state_ukb) %>% tally()
write.csv(t,here("MR_UKBiobank","BinaryData","Logistic",paste0("Phenotype_",outc,".csv")))

