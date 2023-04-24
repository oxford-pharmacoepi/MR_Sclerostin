
rm(list = setdiff(ls(),c("pathData","tok")))

unlink(here("Results","SupplementaryTables_1.xlsx"))
unlink(here("Results","SupplementaryTables_2.xlsx"))
unlink(here("Results","SupplementaryTables_3.xlsx"))
unlink(here("Results","SupplementaryTables_4.xlsx"))
unlink(here("Results","Tables.xlsx"))

# Table 1 ----------------------------------------------------------------------
rm(list = setdiff(ls(), c("pathData","tok")))
t1 <- read_delim(here('Pruning','exposure_data.csv')) %>%
  left_join(
    read_delim(here('Metaanalysis','sc.txt')) %>% 
      select('SNP' = 'MARKERNAME','Pos'), 
    by = 'SNP') %>%
  left_join(
    read_delim(here('Metaanalysis','Fixed_vs_random.csv')) %>%
      select('SNP','i2','q_statistic','q_p.value'),
    by = 'SNP'
  )
tab1 <- data.frame('SNP'  = t1$SNP,
                   'EA'   = t1$effect_allele.exposure,
                   'OA'   = t1$other_allele.exposure,
                   'EAF'  = round(t1$eaf.exposure,2),
                   'Beta' = paste0(round(t1$beta.exposure,2),
                                   ' (',
                                   round(t1$beta.exposure - 1.96*t1$se.exposure,2),
                                   ',',
                                   round(t1$beta.exposure + 1.96*t1$se.exposure,2),
                                   ')'),
                   'N'    = t1$N.exposure,
                   'Pos(GRCh37/hg19)' = paste0('17:',
                                               t1$Pos),
                   'I2'   = round(as.numeric(t1$i2),2),
                   'Q'    = round(t1$q_statistic,2),
                   'Q_P.val' = round(t1$q_p.value,2)
                   )
write.xlsx(tab1,here("Results","Tables.xlsx"),sheetName = 'Table1',append = TRUE)


# Table 2 ----------------------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok")))
outc    <- excel_sheets(here('MR_gwas','gMR_dat.xlsx'))
N       <- matrix(0,length(outc),1)
Beta    <- matrix(0,length(outc),1)
j       <- 1
for (i in outc){
  t1 <- read_xlsx(here('MR_gwas','gMR_dat.xlsx'), sheet = i) %>%
    select(outcome,samplesize.outcome) %>% distinct()
  N[j]    <- t1$samplesize.outcome
  if (i == 'eBMD'){
    t2 <- read_xlsx(here('MR_gwas','gMR_res.xlsx'), sheet = i) %>%
      select(Estimate, CI_LOW, CI_HIGH)
    Beta[j] <- paste0(round(t2$Estimate,2),
                      ' (',
                      round(t2$CI_LOW,2),
                      ', ',
                      round(t2$CI_HIGH,2),
                      ')')
  }else{
    t2 <- read_xlsx(here('MR_gwas','gMR_res.xlsx'), sheet = i) %>%
      select(OR, CI_LOW_OR, CI_HIGH_OR)
    Beta[j] <- paste0(round(t2$OR,2),
                      ' (',
                      round(t2$CI_LOW_OR,2),
                      ', ',
                      round(t2$CI_HIGH_OR,2),
                      ')')
  }
  j <- j + 1
}
tab2_continuous_data <- data.frame('Outcome' = outc, 'N' = N, 'Beta(95CI)' = Beta)
write.xlsx(tab2_continuous_data,here("Results","Tables.xlsx"),sheetName = 'Table2',append = TRUE)


# Table 3 ----------------------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok")))
outc <- excel_sheets(here('MR_UKBiobank','ContinuousData','gMR_dat.xlsx'))
tab2 <- list()
outcome <- matrix(0,length(outc),1)
N       <- matrix(0,length(outc),1)
Mean    <- matrix(0,length(outc),1)
Beta    <- matrix(0,length(outc),1)
j       <- 1
for (i in outc){
  t1 <- read_xlsx(here('MR_UKBiobank','ContinuousData','gMR_dat.xlsx'), sheet = i) %>%
    select(outcome,n.outcome, mean.outcome, std.outcome) %>% distinct()
  t2 <- read_xlsx(here('MR_UKBiobank','ContinuousData','gMR_res.xlsx'), sheet = i) %>%
    select(Estimate, CI_LOW, CI_HIGH)
  
  outcome[j] <- t1$outcome
  N[j]       <- t1$n.outcome
  Mean[j]    <- paste0(as.character(round(t1$mean.outcome, 2)),
                       ' (',as.character(round(t1$std.outcome,2)),
                       ')')
  Beta[j]    <- paste0(as.character(round(t2$Estimate,2)),
                       ' (',
                       as.character(round(t2$CI_LOW,2)),
                       ', ',
                       as.character(round(t2$CI_HIGH,2)),
                       ')')
  j <- j + 1
}


tab2_continuous_data <- data.frame('Outcome' = outcome, 'N' = N, 'Mean(SD)' = Mean, 'Beta(95CI)' = Beta)
write.xlsx(tab2_continuous_data,here("Results","Tables.xlsx"),sheetName = 'Table3',append = TRUE)

# Table 4 ----------------------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok",'outputFolder')))
outc <- excel_sheets(here('MR_UKBiobank','BinaryData','Logistic','gMR_dat.xlsx'))
tab2 <- list()

outcome  <- matrix(0,length(outc),1)
N_or     <- matrix(0,length(outc),1)
Cases_or <- matrix(0,length(outc),1)
OR95CI   <- matrix(0,length(outc),1)
N_hr     <- matrix(0,length(outc),1)
Cases_hr <- matrix(0,length(outc),1)
HR95CI   <- matrix(0,length(outc),1)
j       <- 1
for (i in outc){
  t1 <- read.xlsx(here('MR_UKBiobank','BinaryData','Logistic','gMR_dat.xlsx'), sheetName = i) %>%
    select(outcome,n.outcome, case.outcome) %>% distinct()
  t2 <- read.xlsx(here('MR_UKBiobank','BinaryData','Logistic','gMR_res.xlsx'), sheetName = i) %>%
    select(OR, CI_LOW_OR, CI_HIGH_OR)
  t3 <- read.xlsx(here('MR_UKBiobank','BinaryData','SA_Birth','gMR_dat.xlsx'), sheetName = i) %>%
    select(outcome, n.outcome, cases.outcome) %>% distinct()
  t4 <- read.xlsx(here('MR_UKBiobank','BinaryData','SA_Birth','gMR_res.xlsx'), sheetName = i) %>%
    select(OR, CI_LOW_OR, CI_HIGH_OR)
  
  outcome[j]  <- t1$outcome
  N_or[j]     <- t1$n.outcome
  Cases_or[j] <- paste0(t1$case.outcome,' (',round(t1$case.outcome/t1$n.outcome*100,0),')')
  OR95CI[j]   <- paste0(round(t2$OR,2), ' (',round(t2$CI_LOW_OR,2),', ',round(t2$CI_HIGH_OR,2),')')
  N_hr[j]     <- t3$n.outcome
  Cases_hr[j] <- paste0(t3$cases.outcome,' (',round(t3$cases.outcome/t3$n.outcome*100,0),')')
  HR95CI[j]   <- paste0(round(t4$OR,2),' (',round(t4$CI_LOW_OR,2),', ',round(t4$CI_HIGH_OR,2),')')
  j <- j + 1
}
tab4_binary_data <- data.frame('Outcome'  = outcome,
                               'N.OR'     = N_or,
                               'Cases_OR' = Cases_or,
                               'OR95CI'   = OR95CI,
                               'N.HR'     = N_hr,
                               'Cases_HR' = Cases_hr,
                               'HR95CI'   = HR95CI)
write.xlsx(tab4_binary_data,here("Results","Tables.xlsx"),sheetName = 'Table4',append = TRUE)


# Supplementary table 1 --------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok")))
stab_1 <- read_delim(here('Metaanalysis','Fixed_vs_random.csv')) %>%
  select('SNP',
         'Effect_allele' = 'EA.fixed',
         'Other_allele'  = 'OA.fixed',
         'Effect.fixed',
         'Effect.random',
         'SE.fixed'  = 'se.fixed',
         'SE.random' = 'se.random',
         'Heterogeneity.i2' = 'i2',
         'Heterogeneity.q_statistic' = 'q_statistic',
         'Heterogeneity.q_p.value'   = 'q_p.value',
         'N_samples'  = 'n_samples')

# Supplementary table 2 --------------------------------------------------------
stab_2 <- read_delim(here('Pruning','correlation.csv')) %>% rename(SNP = RS_number)

write_xlsx(list(stab_1,stab_2),here("Results","SupplementaryTables_1.xlsx"))
# Supplementary table 3 --------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok",'outputFolder')))

stab_3 <- read_delim(paste0(pathData,'UKB\\cohort.csv')) %>% 
  select('eid') %>%
  left_join(
    read_delim(paste0(pathData,'genetics.csv')),
    by = 'eid'
  ) %>%
  select(-'...1') %>%
  left_join(
    read.csv(paste0(pathData,"UKB\\ukb669864_Table2.csv")) %>%
      select(eid,
             sex  = X31.0.0,
             age  = X21003.0.0,
             IMDE = X26410.0.0,
             IMDS = X26427.0.0,
             IMDW = X26426.0.0,
             BMI  = X21001.0.0,
             ETH  = X21000.0.0),
    by = "eid") %>%
  mutate(IMD = IMDE,
         IMD = if_else(is.na(IMD),IMDS,IMD),
         IMD = if_else(is.na(IMD),IMDW,IMD)) %>%
  mutate(sex = if_else(sex == 0,'Female','Male'),
         ETH = if_else(ETH %in% c(1,1001,1002,3001,4001),'White','Other ethnic group')) %>%
  select('eid',contains('rs'),'sex','age','BMI','IMD','ETH')
names <- colnames(stab_3 %>% select(contains('rs'))) 

for (i in names){
  stab <- stab_3 %>%
    select('eid','snp' = i,'sex','age','BMI','IMD','ETH') %>%
    mutate(snp = if_else(snp == 2,
                         1,
                         snp))
  
  tab1 <- CreateTableOne(vars = c('age','sex','IMD','BMI','ETH'), data = stab, strata = 'snp', test = FALSE)
  tab1 <- print(tab1, smd = TRUE, quote = FALSE, noSpaces = TRUE, showAllLevels = TRUE, printToggle = FALSE)
  
  write.xlsx(tab1,here("Results","SupplementaryTables_2.xlsx"),sheetName = paste0('ST_3_',i),append = TRUE)
}

# Supplementary table 4 --------------------------------------------------------
rm(list = setdiff(ls(),c("pathData","tok",'outputFolder')))
outc <- excel_sheets(here('MR_UKBiobank','BinaryData','SA_Enrolment','gMR_dat.xlsx'))
outcome <- matrix(0,length(outc),1)
N       <- matrix(0,length(outc),1)
Cases   <- matrix(0,length(outc),1)
HR95CI  <- matrix(0,length(outc),1)
j       <- 1

for (i in outc){
  t1 <- read.xlsx(here('MR_UKBiobank','BinaryData','SA_Enrolment','gMR_dat.xlsx'), sheetName = i) %>%
    select('outcome','n.outcome','cases.outcome') %>% distinct()
  t2 <- read.xlsx(here('MR_UKBiobank','BinaryData','SA_Enrolment','gMR_res.xlsx'), sheetName = i) %>%
    select('OR','CI_LOW_OR','CI_HIGH_OR')
  
  outcome[j] <- t1$outcome
  N[j]       <- t1$n.outcome
  Cases[j]   <- paste0(t1$cases.outcome,
                       ' (',
                       round(t1$cases.outcome/t1$n.outcome*100,0),
                       ')')
  HR95CI[j]     <- paste0(round(t2$OR,2),
                       ' (',
                       round(t2$CI_LOW_OR,2),
                       ', ',
                       round(t2$CI_HIGH_OR,2),
                       ')')
  j <- j + 1
}

stab_4 <- data.frame('Outcome' = outcome,
                     'N'       = N,
                     'Cases'   = Cases,
                     'HR95CI'  = HR95CI)
write.xlsx(stab_4,here("Results","SupplementaryTables_3.xlsx"),sheetName = 'ST_4',append = TRUE)


# Supplementary table 5 --------------------------------------------------------
outc  <- excel_sheets(here('Colocalization','coloc.xlsx'))
nsnps <- matrix(0,length(outc),1)
H0    <- matrix(0,length(outc),1)
H1    <- matrix(0,length(outc),1)
H2    <- matrix(0,length(outc),1)
H3    <- matrix(0,length(outc),1)
H4    <- matrix(0,length(outc),1)
j     <- 1
for (i in outc){
  t <- read.xlsx(here('Colocalization','coloc.xlsx'), sheetName = i)
  nsnps[j] <- t$nsnps
  H0[j]    <- t$H0
  H1[j]    <- t$H1
  H2[j]    <- t$H2
  H3[j]    <- t$H3
  H4[j]    <- t$H4
  j <- j + 1
}
stab_5 <- data.frame('Outcome' = outc,
                     'N_snps'  = nsnps,
                     'P1'      = 0.0001*100,
                     'P2'      = 0.0001*100,
                     'P12'     = 0.00001*100,
                     'H0'      = H0*100,
                     'H1'      = H1*100,
                     'H2'      = H2*100,
                     'H3'      = H3*100,
                     'H4'      = H4*100)
write.xlsx(stab_5,here("Results","SupplementaryTables_4.xlsx"),sheetName = 'ST_5',append = TRUE)


# Supplementary table 6 --------------------------------------------------------
outc <- excel_sheets(here('SensitivityAnalysis','LeaveOneOut','Gwas','loo.xlsx'))
for (i in outc){
  stab_6 <- read.xlsx(here('SensitivityAnalysis','LeaveOneOut','Gwas','loo.xlsx'), sheetName = i) %>%
    select('SNPS',
           'MR Effect' = 'Estimate',
           'CI_95LOW'  = 'CI_LOW',
           'CI_95HIGH' = 'CI_HIGH') %>%
    mutate('Outcome' = i) %>%
    relocate(Outcome)
  
  write.xlsx(stab_6,here("Results","SupplementaryTables_4.xlsx"),sheetName = paste0('ST_6_Gwas_',i),append = TRUE)
}

outc <- excel_sheets(here('SensitivityAnalysis','LeaveOneOut','ContinuousData','loo.xlsx'))
for (i in outc){
  stab_6 <- read.xlsx(here('SensitivityAnalysis','LeaveOneOut','ContinuousData','loo.xlsx'), sheetName = i) %>%
    select('SNPS',
           'MR Effect' = 'Estimate',
           'CI_95LOW'  = 'CI_LOW',
           'CI_95HIGH' = 'CI_HIGH') %>%
    mutate('Outcome' = i) %>%
    relocate(Outcome)
  
  write.xlsx(stab_6,here("Results","SupplementaryTables_4.xlsx"),sheetName = paste0('ST_6_Continuous_',i),append = TRUE)
}

outc <- excel_sheets(here('SensitivityAnalysis','LeaveOneOut','BinaryData_LR','loo.xlsx'))
for (i in outc){
  stab_6 <- read.xlsx(here('SensitivityAnalysis','LeaveOneOut','BinaryData_LR','loo.xlsx'), sheetName = i) %>%
    select('SNPS',
           'MR Effect' = 'Estimate',
           'CI_95LOW'  = 'CI_LOW',
           'CI_95HIGH' = 'CI_HIGH') %>%
    mutate('Outcome' = i) %>%
    relocate(Outcome)
  
  write.xlsx(stab_6,here("Results","SupplementaryTables_4.xlsx"),sheetName = paste0('ST_6_Binary_',i),append = TRUE)
}

outc <- excel_sheets(here('SensitivityAnalysis','LeaveOneOut','BinaryData_SA','loo.xlsx'))
for (i in outc){
  stab_6 <- read.xlsx(here('SensitivityAnalysis','LeaveOneOut','BinaryData_SA','loo.xlsx'), sheetName = i) %>%
    select('SNPS',
           'MR Effect' = 'Estimate',
           'CI_95LOW'  = 'CI_LOW',
           'CI_95HIGH' = 'CI_HIGH') %>%
    mutate('Outcome' = i) %>%
    relocate(Outcome)
  
  write.xlsx(stab_6,here("Results","SupplementaryTables_4.xlsx"),sheetName = paste0('ST_6_Binary_SA_',i),append = TRUE)
}

# Supplementary table 7 --------------------------------------------------------
names <- c('gwas','continuous','BinaryLR','BinarySA')
for (i in names){
  outc <- excel_sheets(here('SensitivityAnalysis','SingleSNP',paste0('MR_dat_',i,'.xlsx')))
  effect  <- matrix(0,length(outc),1)
  C1      <- matrix(0,length(outc),1)
  C2      <- matrix(0,length(outc),1)
  k       <- 1
  for (j in outc){
    t <- read.xlsx(here('SensitivityAnalysis','SingleSNP',paste0('MR_res_',i,'.xlsx')), sheetName = j)
    if (i %in% c('gwas','BinaryLR','BinarySA')){
      if (j == 'eBMD'){
        effect[k]  <- t$Estimate
        C1[k]      <- t$CI_LOW
        C2[k]      <- t$CI_HIGH
      }else{
        effect[k] <- t$OR
        C1[k]     <- t$CI_LOW_OR
        C2[k]     <- t$CI_HIGH_OR
      }
    }else{
      effect[k]  <- t$Estimate
      C1[k]      <- t$CI_LOW
      C2[k]      <- t$CI_HIGH
    }
    k <- k + 1
  }
  
  stab_7 <- data.frame('Outcome' = outc,
                       'Method'  = 'Wald Ratio',
                       'N SNPs'  = 1,
                       'Effect/OR' = effect,
                       'C95_LOW'   = C1,
                       'C95_HIGH'  = C2)
  write.xlsx(stab_7,here("Results","SupplementaryTables_4.xlsx"),sheetName = paste0('ST_7_',i),append = TRUE)
}

