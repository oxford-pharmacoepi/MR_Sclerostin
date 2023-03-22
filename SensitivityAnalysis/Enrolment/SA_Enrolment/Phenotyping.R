# ============================================================================ #
#                                     Phenotyping                              #
#                          (Survival analysis - Age of enrolment)              #
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
# ukb  <- as_tibble(read.delim(here("Data","UKB","ukb669864_birth.csv"), sep = ","))

outc <- "Fracture"

pop <- unique(hes %>% select(eid))

source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingHES.R"))
hes_data <- PhenotypingHes(outc,hes,hesD)

source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingGP.R"))
gp_data <- PhenotypingGP(outc,gp,pop,ukb)

source(here("MR_UKBiobank","BinaryData","SA_Enrolment","PhenotypingUKB.R"))
ukb_data <- PhenotypingUKB(outc,pop,ukb)


# Flow chart information
# t <- hes_data %>%
#   left_join(gp_data, by = "eid") %>%
#   left_join(ukb_data %>%
#               select("eid","state_ukb"), by = "eid") %>%
#   distinct() %>%
#   inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
#                select("eid",
#                       "Caucasian" = "22006-0.0"), # Caucasian
#              by = "eid") %>%
#   filter(Caucasian == 1) %>%
#   inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
#                select(-"...1"),
#              by = "eid") %>%
#   inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") %>% # 10 PRINCIPAL COMPONENTS
#   mutate(state = if_else(state_hes == -1 | state_gp == -1 | state_ukb == -1,-1,0))
# 
# sum(t$state == -1)
# nrow(t %>% filter(state != -1))


source(here("SensitivityAnalysis","Enrolment","SA_Enrolment","MergeTables.R"))
t <- MergeTables(ukb_data,hes_data,gp_data)

write.csv(t,here("SensitivityAnalysis","Enrolment","SA_Enrolment",paste0("Phenotype_",outc,".csv")))
