tablesList <- list()

# Manhattan plot: instrument selection


# --------------------- Instrument selection pruning --------------------------#
for(metaMethod_i in c("Fixed", "Random")){
  list_of_files <- list.files(path = paste0(pathResults,"InstrumentSelection/Instruments"),
                              recursive  = TRUE,
                              pattern    = paste0("iv_",metaMethod_i,"_gene"),
                              full.names = TRUE)
  for(i in 1:length(list_of_files)){
    file <- list_of_files[i]

    # Read file
    t <- readr::read_delim(file)

    # Extract r2 value
    r2 <- gsub(".*r2","",file)
    r2 <- gsub("_clumpw.*","",r2)

    tab <- tibble::tibble(
      "type" = "Pruning",
      "r2"   = as.numeric(r2),
      "SNP"  = t$SNP,
      "EA"   = t$effect_allele.exposure,
      "EAF"  = round(t$eaf.exposure, 2),
      "Beta" = round(t$beta.exposure, 2),
      "SE"   = round(t$se.exposure, 2),
      "PVal" = format(t$pval.exposure, scientific = TRUE, digits = 2)
    )


    if(i == 1){
      tableInstruments <- tab
    }else{
      tableInstruments <- tableInstruments %>%
        dplyr::union_all(tab)
    }
  }

  tablesList[[paste0("PruningInstruments_",metaMethod_i)]] <- tableInstruments %>%
    dplyr::arrange(r2) %>%
    dplyr::mutate(type = dplyr::if_else(r2 == 1e-6, "Single variant", type)) %>%
    dplyr::mutate(r2 = as.character(r2)) %>%
    dplyr::mutate(r2 = dplyr::if_else(r2 == "1e-06", " ", r2)) %>%
    dplyr::group_by(r2) %>%
    dplyr::mutate(type = dplyr::if_else(dplyr::row_number() >= 2, " ", type)) %>%
    dplyr::mutate(r2 = dplyr::if_else(dplyr::row_number() >= 2, " ", r2)) %>%
    flextable::flextable() %>%
    flextable::align(align = "center", part = 'all') %>%
    flextable::bold(part = "header")
}






