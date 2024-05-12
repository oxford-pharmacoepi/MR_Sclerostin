# ============================================================================ #
#                             InstrumentSelection                              #
#                            Marta Alcalde-Herraiz                             #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This scripts selects the instrumental variables to be used for the MR        #
# analysis, based on LD-Clumping (see function                                 #
# "instrumentSelection_byClumping.R" and "ld_clumping.R")                      #
# ============================================================================ #
source(here::here('Functions/instrumentSelection_byClumping.R'))
source(here::here('Functions/instrumentSelection_byPCA.R'))

dat_input <- c("Fixed","Random")
clump_kb_val <- c(250e3)
clump_r2_val <- c(0.001,0.1,0.3,0.5,0.8)

for(dat_input_i in dat_input){
  for(clump_kb_val_i in clump_kb_val){
    for(clump_r2_val_i in clump_r2_val){

      iv <- instrumentSelection_byClumping(dat_input = dat_input_i,
                                           # Clump parameters
                                           clump_kb = clump_kb_val_i,
                                           clump_r2 = clump_r2_val_i,
                                           pval = 1e-6,
                                           clump_p  = 0.99)
    }
  }

  ivPCA <- instrumentSelection_byPCA(dat_input = dat_input_i, pval = 1e-5)
}


dir.create(paste0(pathResults,"InstrumentSelection/Instruments"))
listOfFiles <- list.files(path = paste0(pathResults,"InstrumentSelection/"),
                          recursive  = TRUE,
                          pattern    = "^iv",
                          full.names = TRUE)
file.copy(from = listOfFiles,
          to   = paste0(pathResults,'InstrumentSelection/Instruments/'))
file.remove(listOfFiles)

# Check if is correct:
# TwoSampleMR::clump_data(dat = TwoSampleMR::format_data(
#   readr::read_delim(paste0(pathResults,"InstrumentSelection/FixedMapped.txt")) %>%
#     dplyr::filter(pval < 1e-5) %>%
#     dplyr::filter(eaf > 0.01) %>%
#     dplyr::filter(pos >= 43753738-100e3-1,
#                   pos <= 43758791+100e3+1), type = "exposure"),
#   clump_kb = 250e3,
#   clump_r2 = 0.5,
#   pop = "EUR")





