# ============================================================================ #
#                                CodeToRun                                     #
#                       2023 - Marta Alcalde-Herraiz                           #
# ============================================================================ #
rm(list=ls())
library(pacman)
p_load('here','tibble','dplyr','readr','tidyverse','TwoSampleMR','MendelianRandomization',
       'LDlinkR','xlsx','trio','survival','coloc','ggplot2')

pathData <- "D:\\MR_Sclerostin_Data\\" # Path to the Data directory
tok <- "93873e604a3e" # token for the LDmatrix
outputFolder <- "D:\\MR_Sclerostin_Data\\" # Path to the Results directory
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
source(here('MR_UKBiobank','ContinuousData','MR.R'))

hes  <- as_tibble(read.delim(paste0(pathData,"hesin_diag.txt"),sep = "\t", quote = ""))
hesD <- as_tibble(read.delim(paste0(pathData,"hesin.txt"), sep = "\t", quote = ""))
gp   <- as_tibble(read.delim(paste0(pathData,"gp_clinical.txt"),sep = "\t", quote = ""))
ukb  <- as_tibble(read.delim(paste0(pathData,"UKB\\ukb669864_logistic.csv"), sep = ","))

source(here('MR_UKBiobank','BinaryData','Logistic','Phenotyping.R'))
source(here('MR_UKBiobank','BinaryData','SA_Birth','Phenotyping.R'))
source(here('MR_UKBiobank','BinaryData','MR_LR.R'))
source(here('MR_UKBiobank','BinaryData','MR_SA.R'))

# Survival analysis since UK Biobank enrolment
source(here('MR_UKBiobank','BinaryData','SA_Enrolment','Phenotyping.R'))
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

# Results folder
file.copy(here("Metaanalysis","Fixed_vs_random.csv"), paste0(outputFolder,"Results\\Metaanalysis\\ResultsMetaanalysis.csv"), overwrite = TRUE)
file.copy(here("Pruning","exposure_data.csv"), paste0(outputFolder,"Results\\Metaanalysis\\Instruments.csv"),overwrite = TRUE)
file.copy(here("MR_gwas","gMR_res.xlsx"),paste0(outputFolder,"Results\\MR_gwas\\MR_GWAS_results.xlsx"),overwrite = TRUE)
file.copy(here("MR_UKBiobank","ContinuousData","gMR.xlsx"),paste0(outputFolder,"Results\\MR_UKBiobank\\ContinuousData\\gMR.xlsx"),overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","Logistic","gMR_res.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\LogisticRegression\\gMR_results.xlsx"), overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","Logistic","gMR_dat.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\LogisticRegression\\gMR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","SA_Birth","gMR_res.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\CoxRegression\\gMR_results.xlsx"), overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","SA_Birth","gMR_dat.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\CoxRegression\\gMR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","SA_Enrolment","gMR_res.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\CoxRegression\\gMR_results.xlsx"), overwrite = TRUE)
file.copy(here("MR_UKBiobank","BinaryData","SA_Enrolment","gMR_dat.xlsx"), paste0(outputFolder,"Results\\MR_UKBiobank\\BinaryData\\CoxRegression\\gMR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_res_gwas.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\GWAS\\MR_results.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_dat_gwas.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\GWAS\\LogisticRegression\\MR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_res_BinaryLR.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\MR_UKBiobank\\Binary\\LogisticRegression\\MR_results.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_dat_BinaryLR.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\MR_UKBiobank\\Binary\\LogisticRegression\\MR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_res_BinarySA.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\MR_UKBiobank\\Binary\\CoxRegression\\Birth\\MR_results.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","SingleSNP","MR_dat_BinarySA.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\SingleSNP\\MR_UKBiobank\\Binary\\CoxRegression\\Birth\\MR_SNP_info.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","LeaveOneOut","Gwas","loo.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\LeaveOneOut\\GWAS\\LeaveOneOut.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","LeaveOneOut","ContinuousData","loo.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\LeaveOneOut\\MR_UKBiobank\\ContinuousData\\LeaveOneOut.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","LeaveOneOut","BinaryData_LR","loo.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\LeaveOneOut\\MR_UKBiobank\\BinaryData\\LogisticRegression\\LeaveOneOut.xlsx"), overwrite = TRUE)
file.copy(here("SensitivityAnalysis","LeaveOneOut","BinaryData_SA","loo.xlsx"), paste0(outputFolder,"Results\\SensitivityAnalysis\\LeaveOneOut\\MR_UKBiobank\\BinaryData\\CoxRegression\\LeaveOneOut.xlsx"), overwrite = TRUE)
file.copy(here("Colocalization","data.xlsx"), paste0(outputFolder,"Results\\Colocalization\\data.xlsx"), overwrite = TRUE)
file.copy(here("Colocalization","coloc.xlsx"), paste0(outputFolder,"Results\\Colocalization\\ColocalizationResults.xlsx"), overwrite = TRUE)
