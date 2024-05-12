outcome <- "hf"

list_of_files <- list.files(path = paste0(pathResults,"MendelianRandomisation/TemporaryFiles"),
                            recursive  = TRUE,
                            pattern    = outcome,
                            full.names = TRUE)

list_of_files_fixed  <- tibble::tibble(list_of_files) %>% dplyr::filter(grepl("Fixed",list_of_files)) %>% dplyr::filter(!grepl("PCA", list_of_files))
list_of_files_random <- tibble::tibble(list_of_files) %>% dplyr::filter(!grepl("Fixed",list_of_files))%>% dplyr::filter(!grepl("PCA", list_of_files))


readr::read_delim(list_of_files_fixed %>% dplyr::pull(list_of_files)) %>% dplyr::filter(Instruments %in% c(1,2))
