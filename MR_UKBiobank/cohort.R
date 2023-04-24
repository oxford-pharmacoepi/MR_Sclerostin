# ============================================================================ #
#                       Create cohort for UK Biobank data                      #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Create the cohort for UK Biobank data                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

t <- read_delim(paste0(pathData,"hesin_diag.txt")) %>% # Participants with records in HES inpatient diagnoses dataset
  select(eid) %>% distinct() %>%
  left_join(
    read_delim(paste0(pathData,"UKB\\ukb669864_1.csv")) %>%
      select("eid",
             "caucasian" = "22006-0.0"),
    by = "eid") %>%
  left_join(
    read_delim(paste0(pathData,"UKB\\ukb669864_cov.csv")) %>%
      select("eid",
             'sex_chromosome_aneuploidy' = '22019-0.0',
             'kinship_to_other_participants' = '22021-0.0',
             'genetic_sex' = '22001-0.0'),
    by = "eid") %>%
  left_join(
    read_delim(paste0(pathData,'genetics.csv')),
    by = 'eid') %>%
  left_join(
    read_delim(paste0(pathData,"UKB\\PC.csv")), 
    by = "eid")
  
t <- t %>% 
  filter(caucasian == 1) %>%
  filter(sex == genetic_sex) %>%
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  filter(kinship_to_other_participants == 0) %>%
  filter(!is.na(PC1),!is.na(PC2),!is.na(PC3),!is.na(PC4),!is.na(PC5),
         !is.na(PC6),!is.na(PC7),!is.na(PC8),!is.na(PC9),!is.na(PC10)) %>%
  filter(!is.na(rs6503468), !is.na(rs9910625), !is.na(rs7220711), !is.na(rs66838809),
         !is.na(rs7213935), !is.na(rs80107551)) %>%
  select('eid',
         contains('rs'),
         contains('PC'),
         'sex')

gp <- gp %>% filter(event_dt != "01/01/1900",
                    event_dt != "01/01/1901",
                    event_dt != "02/02/1902",
                    event_dt != "03/03/1903",
                    event_dt != "07/07/2037")
write_csv(t,paste0(pathData,"UKB\\cohort.csv"))

  