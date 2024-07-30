# ============================================================================ #
#                                  CodeToRun                                   #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

rm(list = ls())
library('magrittr')
library("rlang")
# readr, tibble, dplyr, base, numbers, biomaRt, tidyr, here, flextable, 'ftExtra'
pathData <- "D:/Projects/MR_Sclerostin/" # to be specified !!!!!!!!!!!!!!!!!!!!!!!

# Study 1 ----------------------------------------------------------------------
dir.create(paste0(pathData,'Results'))
dir.create(paste0(pathData,'Results/Study1'))
pathResults <- paste0(pathData,'Results/Study1/')

source(here::here("Study1/MetaAnalysis/Preprocessing.R"))

gwama_directory <- paste0(pathData,"OtherSoftware/") # to be specified!!!!
metal_directory <- paste0(pathData,"OtherSoftware/")

file.copy(from = paste0(pathResults,"MetaAnalysis/Ferkingstad_NatureGenetics_GWAS.txt"),
          to   = paste0(metal_directory,"Ferkingstad_NatureGenetics_GWAS.txt"),
          overwrite = TRUE)
file.copy(from = paste0(pathResults,"MetaAnalysis/Pietzner_Science_GWAS.txt"),
          to   = paste0(metal_directory,"Pietzner_Science_GWAS.txt"),
          overwrite = TRUE)
file.copy(from = paste0(pathResults,"MetaAnalysis/Sun_Nature_GWAS.txt"),
          to   = paste0(metal_directory,"Sun_Nature_GWAS.txt"),
          overwrite = TRUE)
dir.create(here::here("Study1/MetaAnalysis/ResultsMetaAnalysis"))

source(here::here("Study1/MetaAnalysis/runGWAMA.R"))
source(here::here("Study1/MetaAnalysis/runMETAL.R"))

file.copy(from = paste0(metal_directory,"METAANALYSIS1.TBL"),
          to   = here::here("Study1/MetaAnalysis/ResultsMetaAnalysis/METAL_Fixed.TBL"))

source(here::here("Study1/MetaAnalysis/Processing.R"))
source(here::here("Study1/MetaAnalysis/Heterogeneity.R"))

source(here::here("Study1/InstrumentSelection/ObtainPositions.R"))
source(here::here("Study1/InstrumentSelection/InstrumentSelection.R"))
source(here::here("Study1/InstrumentSelection/Validation.R"))

source(here::here("Study1/MendelianRandomisation/MendelianRandomisation.R"))

# Study 2
# Select instruments to be seen
dir.create(paste0(pathData,'Results/Study2'))
dir.create(paste0(pathData,'Results/Study2/InstrumentSelection'))
source(here::here("Study2/InstrumentSelection/Instruments.R"))

pathResults <- paste0(pathData,"Results/Study2/")

# Build cases and controls
dir.create(paste0(pathData,'Results/Study2/Cohort'))
source(here::here("Study2/Cohort/BuildCohort.R"))

dir.create(paste0(pathResults,"Phenotyping"))
source(here::here("Study2/Cohort/Phenotyping.R"))

# Mendelian randomisation
dir.create(paste0(pathResults,"MendelianRandomisation"))
source(here::here("Study2/MendelianRandomisation/MendelianRandomisation.R"))

# Colocalisation
dir.create(paste0(pathData,"Results/Colocalisation"))
pathResults <- paste0(pathData,"Results/Colocalisation")
source(here::here("Colocalisation/Colocalisation.R"))

# Sensitivity analysis
dir.create(paste0(pathData,"Results/SensitivityAnalyses"))
dir.create(paste0(pathData,"Results/SensitivityAnalyses/SurvivalPackage"))
source(here::here("SurvivalPackage.R"))

# Tables
dir.create(paste0(pathData,"Results/Tables"))
pathResults <- paste0(pathData,"Results/Tables")



