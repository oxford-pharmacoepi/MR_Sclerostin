# ============================================================================ #
#                             MendelianRandomisation                           #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This script performs generalised IVW when there is more than one instrument, #
# and Wald Ratio when there is only one instrument.                            #
# ============================================================================ #
source(here::here("Functions","loadOutcomeGwas.R"))
source(here::here("Functions","harmonisingData.R"))
source(here::here("Functions","ld_matrix_local.R"))
source(here::here("Functions","mendelianRandomisation.R"))
source(here::here("Functions","mendelianRandomisationPCA.R"))

# Exposure:
list_of_files <- list.files(path = paste0(pathResults,"/InstrumentSelection/Instruments/"),
                            recursive  = TRUE,
                            pattern    = "^iv_",
                            full.names = TRUE)

# Safe time not performing the analysis for all the instruments if they are the same
gwasID_i <- list(
   "bmd"  = "ebi-a-GCST006979",
   "hf"   = "ebi-a-GCST90161240",
   "cad"  = "ebi-a-GCST90132314",
   "mi"   = "ebi-a-GCST011365",
   "is"   = "ebi-a-GCST90104540",
   "hypertension" = "ukb-b-14057",
   "t2dm" = "ebi-a-GCST90132184",
   "hba1c"   = "NULL",
   "glucose" = "NULL",
   "ldl"  = "NULL",
   "hdl"  = "NULL"
)

# -------------------------------------------------------------------------
list_of_files_clumping <- tibble::tibble(list_of_files) %>% dplyr::filter(!grepl("PCA",list_of_files))
list_of_files_pca      <- tibble::tibble(list_of_files) %>% dplyr::filter(grepl("PCA", list_of_files))

resultsMR <- list()
datMR     <- list()

for(i in list_of_files_clumping$list_of_files){
  file <- i

  # Extract meta-analysis method (fixed or random)
  metaMethod <- gsub(".*iv_","",file)
  metaMethod <- gsub("_gene.*","",metaMethod)

  # Extract r2 value from the file name
  r2 <- gsub(".*r2","",file)
  r2 <- gsub("_.*","",r2)

  # Extract clumping window
  clumpw <- gsub(".*clumpw","",file)
  clumpw <- gsub(".txt*.","", clumpw)

  # Read exposure file
  exposure <- readr::read_delim(file) # Clumping method ------------------------

  for(name_i in names(gwasID_i)){
    print(name_i)
    # Clumping outcome
    outcome  <- loadOutcomeGwas(outcome = name_i, exposure = exposure)

    # Mendelian randomization using clumping
    mr_res  <- mendelianRandomisation(exposure, outcome, gwasID = gwasID_i[[name_i]], binary = (!name_i %in% c("bmd","hba1c","glucose","hdl","ldl")), name = name_i)

    varname <- paste0(name_i,"_",metaMethod,"_",r2,"_",clumpw)
    resultsMR[[varname]] <- mr_res$Results
    datMR[[varname]]     <- mr_res$Dat
  }
}


for(i in list_of_files_pca$list_of_files){
  file <- i

  # Extract meta-analysis method (fixed or random)
  metaMethod <- gsub(".*iv_","",file)
  metaMethod <- gsub("_gene.*","",metaMethod)

  # Read exposure file
  exposurePCA <- readr::read_delim(paste0(
    pathResults,"InstrumentSelection/Instruments/iv_",metaMethod)) # PCA method --

  for(name_i in names(gwasID_i)){
    # PCA outcome
    outcomePCA  <- loadOutcomeGwas(outcome = name_i, exposure = exposurePCA)

    # Mendelian randomization using PCA (>99.9% variance)
    mr_res_pca  <- mendelianRandomisationPCA(exposure = exposurePCA, outcome = outcomePCA, gwasID = gwasID_i[[name_i]], binary = (!name_i %in% c("bmd","hba1c","glucose","hdl","ldl")), name = name_i, pca_threshold = 99.9)

    varname <- paste0(name_i,"_999_",metaMethod)
    resultsMR[[varname]] <- mr_res_pca$Results
  }
}

dir.create(paste0(pathResults,"MendelianRandomisation/Results"))
for(list_name_i in 1:length(resultsMR)){
  name_i <- names(resultsMR)
  readr::write_tsv(resultsMR[[list_name_i]],paste0(pathResults,"MendelianRandomisation/Results/",name_i[list_name_i],".txt"))
}

dir.create(paste0(pathResults,"MendelianRandomisation/Data"))
name_i <- names(resultsMR)
for(list_name_i in 1:length(datMR)){
  readr::write_tsv(datMR[[list_name_i]],paste0(pathResults,"MendelianRandomisation/Data/",name_i[list_name_i],".txt"))
}

