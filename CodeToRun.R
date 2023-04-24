# ============================================================================ #
#                                CodeToRun                                     #
#                       2023 - Marta Alcalde-Herraiz                           #
# ============================================================================ #
rm(list=ls())
pacman::p_load('here','tibble','dplyr','readr','tidyverse','TwoSampleMR','MendelianRandomization',
               'LDlinkR','xlsx','trio','survival','coloc','ggplot2','readxl','tableone')

# Specify the following variables
pathData <- "D:\\MR_Sclerostin_Data\\" # Path to the Data directory
tok <- "93873e604a3e" # token for the LDmatrix. You can ask for one here: https://analysistools.cancer.gov/LDlink/?tab=apiaccess

source(here('Metaanalysis','preprocessing.R'))

# Instructions to use GWAMA:                                                   
#  1. First install gwama (can be found on this directory).                    
#  2. Create a .in file with the following lines:                              
#       sc.txt                                                                 
#       nn.txt                                                                 
#       ng.txt                                                                 
#     It will be read by GWAMA software to know which files it has to read.    
#  3. Open command prompt in this directory                                    
# Fixed effect method:                                                         
#  4. Write in the command line the following:                                 
# gwama -o Fixed_results/Fixed -qt
# Random effect method:                                                        
#  5. Write in the command line the following:                                 
# gwama -o Random_results/Random -qt -r                                       
# Instructions to use METAL                                                    
#  6. Click on metal.exe to open the command window.                           
#  7. Write the following instructions:                                        
#     SCHEME STDERR                                                            
#     MARKER MARKERNAME                                                        
#     ALLELELABELS EA NEA                                                      
#     PVALUELABEL pval                                                         
#     FREQLABEL EAF                                                            
#     AVERAGEFREQ ON                                                           
#     WEIGHTLABEL N                                                            
#     STDERRLABEL SE                                                           
#     EFFECTLABEL BETA                                                         
#     PROCESS sc.txt                                                           
#     PROCESS nn.txt                                                           
#     PROCESS ng.txt                                                           
#     ANALYZE                                                                  
#    The file: METAANALYSIS1.TBL.info will contain the results                 

source(here('Metaanalysis','processing.R'))
source(here('Pruning','pruning.R'))
source(here('MR_gwas','MR.R'))
source(here('MR_UKBiobank','genetics.R'))
source(here('MR_UKBiobank','cohort.R'))
source(here('MR_UKBiobank','ContinuousData','MR.R'))

hes  <- as_tibble(read.delim(paste0(pathData,"hesin_diag.txt"),sep = "\t", quote = ""))
hesD <- as_tibble(read.delim(paste0(pathData,"hesin.txt"), sep = "\t", quote = ""))
gp   <- as_tibble(read.delim(paste0(pathData,"gp_clinical.txt"),sep = "\t", quote = ""))
gp <- gp %>% filter(event_dt != "01/01/1900",
                    event_dt != "01/01/1901",
                    event_dt != "02/02/1902",
                    event_dt != "03/03/1903",
                    event_dt != "07/07/2037")

source(here('MR_UKBiobank','BinaryData','Logistic','Phenotyping.R'))
source(here('MR_UKBiobank','BinaryData','SA_Birth','Phenotyping.R'))
source(here('MR_UKBiobank','BinaryData','SA_Enrolment','Phenotyping.R'))
source(here('MR_UKBiobank','BinaryData','MR_LR.R'))
source(here('MR_UKBiobank','BinaryData','MR_SA.R'))
source(here('MR_UKBiobank','BinaryData','SA_Enrolment','MR_SA.R'))

# COLOCALIZATION ---------------------------------------------------------------
source(here('Colocalization','colocalization.R'))

# SENSITIVITY ANALYSIS ---------------------------------------------------------
# Single SNP
source(here('SensitivityAnalysis','SingleSNP','exposure_data.R'))
source(here('SensitivityAnalysis','SingleSNP','SingleSNP_gwas.R'))
source(here('SensitivityAnalysis','SingleSNP','SingleSNP_continuous.R'))
source(here('SensitivityAnalysis','SingleSNP','SingleSNP_BinaryLR.R'))
source(here('SensitivityAnalysis','SingleSNP','SingleSNP_BinarySA.R'))

# Results
source(here('Tables','UntitledR.R'))
p <- here('Figures','AllFigures.m')
system(paste0("matlab -nodisplay -r \"run('",p,"'); exit\""))

