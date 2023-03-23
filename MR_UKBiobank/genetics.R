# ============================================================================ #
#                              GENETIC DATA                                    #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# It creates the file "genetics.csv" containing the UKBiobank data for each    #
# one of the instruments.                                                      # 
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

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
  mutate(rs6503468 = if_else(ped$rs6503468.1 == "T" & ped$rs6503468.2 == "T", 2,rs6503468),
         rs6503468 = if_else(ped$rs6503468.1 == "T" & ped$rs6503468.2 == "C", 1,rs6503468),
         rs6503468 = if_else(ped$rs6503468.1 == "C" & ped$rs6503468.2 == "T", 1,rs6503468),
         
         rs9910625 = if_else(ped$rs9910625.1 == "T" & ped$rs9910625.2 == "T", 2,rs9910625),
         rs9910625 = if_else(ped$rs9910625.1 == "T" & ped$rs9910625.2 == "C", 1,rs9910625),
         rs9910625 = if_else(ped$rs9910625.1 == "C" & ped$rs9910625.2 == "T", 1,rs9910625),
         
         rs7220711 = if_else(ped$rs7220711.1 == "A" & ped$rs7220711.2 == "A", 2,rs7220711),
         rs7220711 = if_else(ped$rs7220711.1 == "A" & ped$rs7220711.2 == "G", 1,rs7220711),
         rs7220711 = if_else(ped$rs7220711.1 == "G" & ped$rs7220711.2 == "A", 1,rs7220711),
         
         rs66838809 = if_else(ped$rs66838809.1 == "A" & ped$rs66838809.2 == "A", 2,rs66838809),
         rs66838809 = if_else(ped$rs66838809.1 == "A" & ped$rs66838809.2 == "G", 1,rs66838809),
         rs66838809 = if_else(ped$rs66838809.1 == "G" & ped$rs66838809.2 == "A", 1,rs66838809),
         
         rs7213935 = if_else(ped$rs7213935.1 == "A" & ped$rs7213935.2 == "A", 2,rs7213935),
         rs7213935 = if_else(ped$rs7213935.1 == "A" & ped$rs7213935.2 == "G", 1,rs7213935),
         rs7213935 = if_else(ped$rs7213935.1 == "G" & ped$rs7213935.2 == "A", 1,rs7213935),
         
         rs80107551 = if_else(ped$rs80107551.1 == "T" & ped$rs80107551.2 == "T", 2,rs80107551),
         rs80107551 = if_else(ped$rs80107551.1 == "T" & ped$rs80107551.2 == "C", 1,rs80107551),
         rs80107551 = if_else(ped$rs80107551.1 == "C" & ped$rs80107551.2 == "T", 1,rs80107551))

write.csv(gen,paste0(pathData,"genetics.csv"))
