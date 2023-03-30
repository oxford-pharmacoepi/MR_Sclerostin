rm(list = ls())
library(pacman)
pacman::p_load('here','dplyr','readr','trio','tidyverse','xlsx','tableone')

# Genetics ---------------------------------------------------------------------
# Reading PED files
ped <- as_tibble(read.pedfile(paste0(pathData,"SNPs\\selected_snps_c17.ped"))) %>%
  select("famid",
         "rs6503468.1" = "SNP2.1","rs6503468.2" = "SNP2.2",
         "rs9910625.1" = "SNP4.1","rs9910625.2" = "SNP4.2",
         "rs7220711.1" = "SNP5.1","rs7220711.2" = "SNP5.2",
         "rs66838809.1" = "SNP6.1","rs66838809.2" = "SNP6.2",
         "rs7213935.1" = "SNP7.1","rs7213935.2" = "SNP7.2",
         "rs80107551.1" = "SNP8.1","rs80107551.2" = "SNP8.2")

gen <- as_tibble(data.frame(ped$famid,matrix(0,length(ped$famid),6))) %>%
  rename("eid" = "ped.famid","rs6503468" = "X1","rs9910625" = "X2","rs7220711" = "X3",
         "rs66838809" = "X4","rs7213935" = "X5","rs80107551" = "X6") %>%
  mutate(rs6503468 = if_else(ped$rs6503468.1 == "T" & ped$rs6503468.2 == "T", 1,rs6503468),
         rs6503468 = if_else(ped$rs6503468.1 == "T" & ped$rs6503468.2 == "C", 1,rs6503468),
         rs6503468 = if_else(ped$rs6503468.1 == "C" & ped$rs6503468.2 == "T", 1,rs6503468),

         rs9910625 = if_else(ped$rs9910625.1 == "T" & ped$rs9910625.2 == "T", 1,rs9910625),
         rs9910625 = if_else(ped$rs9910625.1 == "T" & ped$rs9910625.2 == "C", 1,rs9910625),
         rs9910625 = if_else(ped$rs9910625.1 == "C" & ped$rs9910625.2 == "T", 1,rs9910625),

         rs7220711 = if_else(ped$rs7220711.1 == "A" & ped$rs7220711.2 == "A", 1,rs7220711),
         rs7220711 = if_else(ped$rs7220711.1 == "A" & ped$rs7220711.2 == "G", 1,rs7220711),
         rs7220711 = if_else(ped$rs7220711.1 == "G" & ped$rs7220711.2 == "A", 1,rs7220711),

         rs66838809 = if_else(ped$rs66838809.1 == "A" & ped$rs66838809.2 == "A", 1,rs66838809),
         rs66838809 = if_else(ped$rs66838809.1 == "A" & ped$rs66838809.2 == "G", 1,rs66838809),
         rs66838809 = if_else(ped$rs66838809.1 == "G" & ped$rs66838809.2 == "A", 1,rs66838809),

         rs7213935 = if_else(ped$rs7213935.1 == "A" & ped$rs7213935.2 == "A", 1,rs7213935),
         rs7213935 = if_else(ped$rs7213935.1 == "A" & ped$rs7213935.2 == "G", 1,rs7213935),
         rs7213935 = if_else(ped$rs7213935.1 == "G" & ped$rs7213935.2 == "A", 1,rs7213935),

         rs80107551 = if_else(ped$rs80107551.1 == "T" & ped$rs80107551.2 == "T", 1,rs80107551),
         rs80107551 = if_else(ped$rs80107551.1 == "T" & ped$rs80107551.2 == "C", 1,rs80107551),
         rs80107551 = if_else(ped$rs80107551.1 == "C" & ped$rs80107551.2 == "T", 1,rs80107551))

write.csv(gen,here("Tables","Table2_genetics.csv"))


# ------------------------------------------------------------------------------
ukb_tab <- read.csv(paste0(pathData,"UKB\\ukb669864_Table2.csv"))


t    <- read_delim(here("MR_UKBiobank","BinaryData","Logistic","Phenotype_T2DM.csv")) %>%
  select("eid",
         "outcome" = "state") %>%
  inner_join(read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>%
               select("eid",
                      "Caucasian" = "22006-0.0"), # Caucasian
             by = "eid") %>%
  filter(Caucasian == 1) %>%
  inner_join(read_delim(here("Tables","Table2_genetics.csv")) %>% # SNPs expression
               select(-"...1"),
             by = "eid") %>%
  inner_join(read_delim(paste0(pathData,"UKB\\PC.csv")), by = "eid") %>%# 10 PRINCIPAL COMPONENTS
  select("eid","rs6503468","rs9910625","rs7220711","rs66838809","rs7213935","rs80107551") %>%
  inner_join(ukb_tab %>%
               select(eid,
                      sex = X31.0.0,
                      age = X21003.0.0,
                      IMDE = X26410.0.0,
                      IMDS = X26427.0.0,
                      IMDW = X26426.0.0,
                      BMI = X21001.0.0,
                      ETH = X21000.0.0),
             by = "eid") %>%
  mutate(IMD = IMDE,
         IMD = if_else(is.na(IMD),IMDS,IMD),
         IMD = if_else(is.na(IMD),IMDW,IMD)) %>%
  mutate(sex = if_else(sex == 0,'Female','Male'),
         ETH = if_else(ETH %in% c(1,1001,1002,3001,4001),'White','Other ethnic group'))

snp_id <- c("rs6503468","rs9910625","rs7220711","rs66838809","rs7213935","rs80107551")

for (i in 1:length(snp_id)){
  snp_id1 <- snp_id[i]
  tab1 <- CreateTableOne(vars = c('age','sex','IMD','BMI','ETH'), data = t, strata = snp_id1, test = FALSE)
  tab3Mat <- print(tab1, smd = TRUE, quote = FALSE, noSpaces = TRUE, showAllLevels = TRUE, printToggle = FALSE)
  
  ## Save to a CSV file
  write.xlsx(tab3Mat, file = here("Tables","Table2.xlsx"), sheetName = snp_id1, append = TRUE)
}
