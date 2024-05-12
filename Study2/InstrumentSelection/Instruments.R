# ============================================================================ #
#                                  Instruments                                 #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

listOfFiles <- list.files(path = paste0(pathResults,"InstrumentSelection/Instruments/"),
                          full.names = TRUE,
                          pattern = "r2")
iv <- readr::read_delim(listOfFiles) %>%
  dplyr::distinct(SNP)

readr::write_delim(iv, paste0(pathData,"Results/Study2/InstrumentSelection/instruments.txt"))
