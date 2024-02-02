# ============================================================================ #
#                                 BuildCohort                                  #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #


# Load databases
ukb <- tibble::as_tibble(read.table(paste0(pathData,'UKBiobank/ukb673925_20231210.tab'), header=TRUE, sep="\t"))
ukb$f.53.0.0 <- as.Date(ukb$f.53.0.0)
ukb$f.53.1.0 <- as.Date(ukb$f.53.1.0)
ukb$f.53.2.0 <- as.Date(ukb$f.53.2.0)
ukb$f.53.3.0 <- as.Date(ukb$f.53.3.0)

hesin      <- tibble::as_tibble(read.table(paste0(pathData,'UKBiobank/hesin.txt'), header=TRUE, sep="\t"))
hesin_diag <- tibble::as_tibble(read.table(paste0(pathData,'UKBiobank/hesin_diag.txt'), header=TRUE, sep="\t"))
gp         <- tibble::as_tibble(read.delim(paste0(pathData,'UKBiobank/gp_clinical.txt'), header=TRUE, sep = "\t", quote = ""))

# Remove invalid dates
gp <- gp %>% dplyr::filter(event_dt != "01/01/1900",
                           event_dt != "01/01/1901",
                           event_dt != "02/02/1902",
                           event_dt != "03/03/1903",
                           event_dt != "07/07/2037")


# Remove events from after 01/01/2020
gp <- gp %>%
  dplyr::mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
  dplyr::filter(event_dt < as.Date("01/01/2020", format = "%d/%m/%Y"))
hesin <- hesin %>%
  dplyr::mutate(epistart = as.Date(epistart, format = "%d/%m/%Y"),
                admidate = as.Date(admidate, format = "%d/%m/%Y")) %>%
  dplyr::filter(epistart < as.Date("01/01/2020", format = "%d/%m/%Y")) %>%
  dplyr::filter(admidate < as.Date("01/01/2020", format = "%d/%m/%Y"))

# Cohort quality control
t <- ukb %>%
  dplyr::rename("eid" = "f.eid",
                "sex" = "f.31.0.0",
                "genetic_ethnic_grouping"       = "f.22006.0.0",
                "sex_chromosome_aneuploidy"     = "f.22019.0.0",
                "kinship_to_other_participants" = "f.22021.0.0",
                "genetic_sex" = "f.22001.0.0",
                "batch" = "f.22000.0.0") %>%
  dplyr::filter(
    genetic_ethnic_grouping == 1, # European ancestry
    sex == genetic_sex, # Same registered sex and genetic sex
    is.na(sex_chromosome_aneuploidy), # No sex chromosome aneuploidy
    kinship_to_other_participants == 0 # No kinship to other participants
  )

# Genetics
# Genetic data
map <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c17.map"),col_names = FALSE) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 2)) %>%
  dplyr::group_by(X2) %>%
  dplyr::mutate(SNP = paste0(X2,"_",dplyr::row_number())) %>%
  dplyr::pull(SNP)

ped <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c17.ped"),col_names = FALSE) %>%
  dplyr::select(-"X2":-"X6")
colnames(ped) <- c("eid", map)

alleles <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c17.bim"),col_names = FALSE) %>%
  dplyr::select("SNP" = "X2", "ALLELE0" = "X5") %>%
  tidyr::pivot_wider(names_from = "SNP", values_from = "ALLELE0") %>%
  dplyr::slice(rep(1:dplyr::n(), each = nrow(ped)))
snps <- colnames(alleles)
alleles <- alleles %>%
  dplyr::rename_with(~paste0(.x,"_0")) %>%
  dplyr::mutate(eid = ped$eid) %>%
  dplyr::inner_join(ped)

for (snp in snps) {
  alleles <- alleles %>%
    dplyr::mutate(
      !!snp := dplyr::if_else(
        .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_1")]], 1, 0
      ) + dplyr::if_else(
        .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_2")]], 1, 0
      ) + dplyr::if_else(
        .data[[paste0(snp, "_1")]] == "0", NA, 0
      )
    )
}

genetic_data <- alleles %>%
  dplyr::select(eid, all_of(snps), all_of(paste0(snps,"_0"))) %>%
  dplyr::inner_join(t %>%
                      dplyr::select("eid" , "sex",
                                    "PC1" = "f.22009.0.1", "PC2" = "f.22009.0.2", "PC3" = "f.22009.0.3",
                                    "PC4" = "f.22009.0.4", "PC5" = "f.22009.0.5", "PC6" = "f.22009.0.6",
                                    "PC7" = "f.22009.0.7", "PC8" = "f.22009.0.8", "PC9" = "f.22009.0.9",
                                    "PC10" = "f.22009.0.10", "batch"),
                    by = "eid")

write.table(t, paste0(pathResults,"Cohort/ukb.txt"), row.names = FALSE)
write.table(hesin %>% dplyr::filter(eid %in% t$eid), paste0(pathResults,"Cohort/hesin.txt"), row.names = FALSE)
write.table(hesin_diag %>% dplyr::filter(eid %in% t$eid) %>% dplyr::inner_join(hesin %>% dplyr::select("eid","ins_index")), paste0(pathResults,"Cohort/hesin_diag.txt"), row.names = FALSE)
write.table(gp %>% dplyr::filter(eid %in% t$eid), paste0(pathResults,"Cohort/gp.txt"), row.names = FALSE)
write.table(genetic_data, paste0(pathResults,"Cohort/genetic_data.txt"), row.names = FALSE)


