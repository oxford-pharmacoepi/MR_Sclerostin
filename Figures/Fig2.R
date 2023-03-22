# UK Biobank flowchart
hes  <- as_tibble(read.delim(here("Data","hesin_diag.txt"),sep = "\t", quote = ""))
ukb  <- ukb  <- as_tibble(read.delim(here("Data","UKB","ukb669864_logistic.csv"), sep = ","))
gp   <- gp   <- as_tibble(read.delim(here("Data","gp_clinical.txt"),sep = "\t", quote = ""))

gp <- gp %>% select(eid) %>% distinct() # gp participants
hes <- hes %>% select(eid) %>% distinct() # ukb participants

gp <- gp %>% filter(eid %in% hes$eid) # gp participants in hes
ukb <- ukb %>% filter(eid %in% hes$eid) # ukb participants in hes

# Caucasian genetic ancestry
gp <- gp %>% inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
                          select("eid",
                                 "Caucasian" = "22006-0.0"), # Caucasian
                        by = "eid") %>%
  filter(Caucasian == 1)

ukb <- ukb %>% inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
                           select("eid",
                                  "Caucasian" = "22006-0.0"), # Caucasian
                         by = "eid") %>%
  filter(Caucasian == 1)

hes <- hes %>% inner_join(read_delim(here("Data","UKB","ukb669864_1.csv")) %>%
                            select("eid",
                                   "Caucasian" = "22006-0.0"), # Caucasian
                          by = "eid") %>%
  filter(Caucasian == 1)


# Genetic data
outc <- "Fracture"

gp <- gp %>% inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
                           select(-"...1"),
                         by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS

hes <- hes %>% inner_join(read_delim(here("MR_UKBiobank","genetics.csv")) %>% # SNPs expression
                          select(-"...1"),
                        by = "eid") %>%
  inner_join(read_delim(here("Data","UKB","PC.csv")), by = "eid") # 10 PRINCIPAL COMPONENTS



