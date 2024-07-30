# ============================================================================ #
#                                   Tables                                     #
# ============================================================================ #

rm(list = ls())
library(magrittr)
pathData <- "D:/Projects/MR_Sclerostin/" # to be specified !!!!!!!!!!!!!!!!!!!!!!!
source(here::here("Functions/loadOutcomeGwas.R"))

gwas_outcomes <- list("bmd" = "Heel bone mineral density",
                      "hf"  = "Hip fracture",
                      "ldl" = "LDL cholesterol",
                      "hdl" = "HDL cholesterol",
                      "glucose" = "Fasting glucose",
                      "hba1c" = "HbA1c",
                      "cad" = "Coronary artery disease",
                      "mi"  = "Myocardial infarction",
                      "is"  = "Ischaemic stroke",
                      "hypertension" = "Hypertension",
                      "t2dm" = "Type 2 diabetes")

ukb_outcomes <- list("cholesterol" = "Cholesterol",
                     "ldl" = "LDL cholesterol",
                     "hdl" = "HDL cholesterol",
                     "triglycerides" = "Triglycerides",
                     "apoA"  = "Apolipoprotein-A",
                     "apoB"  = "Apolipoprotein-B",
                     "crp"   = "C-Reactive protein",
                     "lipoprotein" = "Lipoprotein (a)",
                     "glucose"     = "Glucose",
                     "hba1c" = "HbA1c",
                     "CAD"   = "Coronary artery disease",
                     "MI"    = "Myocardial infarction",
                     "IS"    = "Ischaemic stroke",
                     "Hypertension" = "Hypertension",
                     "T2DM"  = "Type 2 diabetes")

iv  <- read.table(paste0(pathData,
                         "Results/Study1/InstrumentSelection/Instruments/iv_Fixed_gene500000kb_r20.3_clumpw250000.txt"),
                  header = TRUE)

gwas_files <- list.files(path = paste0(pathData,"Results/Study1/MendelianRandomisation/Results/"),
                         recursive  = FALSE,
                         pattern    = "Fixed_0.3",
                         full.names = TRUE)

cont_files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation/"),
                         recursive  = FALSE,
                         pattern    = "Fixed_0.3_Logistic",
                         full.names = TRUE)

cat_files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation/"),
                        recursive  = FALSE,
                        pattern    = "Fixed_0.3_Logaritmic",
                        full.names = TRUE)

surv_files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation/"),
                         recursive  = FALSE,
                         pattern    = "Fixed_0.3_Birth",
                         full.names = TRUE)

# Table 1 ----------------------------------------------------------------------
t1 <- read.table(paste0(pathData,"Results/Study1/InstrumentSelection/Instruments/iv_Fixed_gene500000kb_r20.3_clumpw250000.txt"), header = TRUE) %>%
  dplyr::mutate(`Beta (95% CI)` = paste0(round(beta.exposure, digits = 2),
                                         " (",
                                         round(beta.exposure-1.96*se.exposure, digits = 2),
                                         ", ",
                                         round(beta.exposure+1.96*se.exposure, digits = 2),
                                         ")")) %>%
  dplyr::mutate(`Pos (GRCh38/hg38)` = paste0(chr.exposure,":",pos.exposure)) %>%
  dplyr::inner_join(
    readr::read_table(paste0(pathData,"Results/Study1/MetaAnalysis/Fixed.txt")) %>%
      dplyr::select("SNP", "samplesize", "i2", "Q" = "q_statistic", `Q P-Value` = "q_pvalue")
  ) %>%
  dplyr::select("SNP" = "SNP",
                "EA"  = "effect_allele.exposure",
                "OA"  = "other_allele.exposure",
                "EAF" = "eaf.exposure",
                `Beta (95% CI)`,
                "P-Value" = "pval.exposure",
                "N"   = "samplesize",
                "Pos" = `Pos (GRCh38/hg38)`,
                "i2", "Q", `Q P-Value`,
                "beta.exposure","se.exposure") %>%
  dplyr::mutate(EAF = round(EAF, digits = 2),
                i2  = round(as.numeric(i2),  digits = 2),
                Q   = round(Q,   digits = 3),
                `Q P-Value` = round(`Q P-Value`, digits = 2),
                `P-Value`   = formatC(`P-Value`, digits = 2, format = "e")) %>%
  dplyr::mutate(
    r2 = (2*EAF*(1-EAF)*(beta.exposure)^2)/(2*EAF*(1-EAF)*(beta.exposure)^2+2*EAF*(1-EAF)*N*se.exposure^2)
  ) %>%
  dplyr::mutate(
    `F statistic` = (r2*(N-2))/(1-r2)
  ) %>%
  dplyr::select(-"beta.exposure", -"se.exposure", -"r2")


# Table 2 ----------------------------------------------------------------------
# Gwas results
rm("outcome_table")
# Study1
gwas_data_files <- gwas_files %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(value = gsub("MendelianRandomisation/Results/","MendelianRandomisation/Data/", value))

for(i in gwas_data_files$value){
  # Outcome name
  outcome_name <- gsub(".*Data/","",i)
  outcome_name <- gsub("_.*","",outcome_name)

  t <- read.table(i, header = TRUE) %>%
    dplyr::mutate(Outcome = gwas_outcomes[[outcome_name]],
                  Source  = "GWAS",
                  N       = samplesize.outcome,
                  Type    = "OR",
                  Estimate = exp(beta.outcome)) %>%
    dplyr::select("Outcome", "Source", "SNP", "N", "effect_allele", "eaf.outcome",
                  "Type", "Estimate", "se.outcome", "pval.outcome")

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

files <- cont_files %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(value = gsub("MendelianRandomisation/","MendelianRandomisation/Data/", value)) %>%
  dplyr::union_all(
    cat_files %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(value = gsub("MendelianRandomisation/","MendelianRandomisation/Data/", value))
  ) %>%
  dplyr::union_all(
    surv_files %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(value = gsub("MendelianRandomisation/","MendelianRandomisation/Data/", value))
  ) %>%
  dplyr::filter(!grepl("ebmd", value) & !grepl("Fracture", value))

for(i in files$value){
  # regression
  regression <- gsub(".*3_","",i)
  regression <- gsub("_.*","",regression)

  # Outcome name
  outcome_name <- gsub(paste0(".*",regression,"_"),"",i)
  outcome_name <- gsub(".txt","",outcome_name)

  t <- read.table(i, header = TRUE) %>%
    dplyr::mutate(Outcome = ukb_outcomes[[outcome_name]],
                  Source  = "UK Biobank",
                  N       = samplesize.outcome,
                  Type    = "OR",
                  Estimate = exp(beta.outcome)) %>%
    dplyr::select("Outcome", "Source", "SNP", "N", "effect_allele", "eaf.outcome",
                  "Type", "Estimate", "se.outcome", "pval.outcome")

  if(regression == "Birth"){
    t <- t %>%
      dplyr::mutate(Type = "HR",
                    Source = "UK Biobank (survival outcome)")
  }else if (regression == "Logaritmic"){
    t <- t %>%
      dplyr::mutate(Source = "UK Biobank (categorical outcome)")
  }

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

t2 <- unlist(gwas_outcomes) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate("Type" = "GWAS") %>%
  dplyr::union_all(
    unlist(ukb_outcomes) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate("Type" = dplyr::if_else(value %in% c("Coronary artery disease", "Myocardial infarction",
                                                         "Ischaemic stroke", "Hypertension", "Type 2 diabetes"),
                                            "UK Biobank (categorical outcome)",
                                            "UK Biobank"))
  ) %>%
  dplyr::union_all(
    unlist(ukb_outcomes) %>%
      dplyr::as_tibble() %>%
      dplyr::filter(value %in% c("Coronary artery disease", "Myocardial infarction",
                                 "Ischaemic stroke", "Hypertension", "Type 2 diabetes")) %>%
      dplyr::mutate("Type" = "UK Biobank (survival outcome)")
  ) %>%
  dplyr::rename("Outcome" = "value", "Source" = "Type") %>%
  dplyr::inner_join(outcome_table, by = c("Outcome", "Source")) %>%
  dplyr::mutate(eaf.outcome = round(eaf.outcome, digits = 2),
                Estimate = round(Estimate, digits = 3),
                se.outcome = round(se.outcome, digits = 3),
                pval.outcome = formatC(pval.outcome, digits = 2, format = "e"))

write.csv(t2,paste0(pathData,"Results/Tables/t2.csv"))


# UK Biobank results
rm("outcome_table")
for(type in c("Logistic","Logaritmic","Birth")){
  if(type == "Logistic"){
    for(i in c("cholesterol","ldl","hdl","triglycerides","apoA","apoB","crp",
               "lipoprotein","glucose","hba1c")){
      t <- read.table(paste0(pathData,"Results/Study2/MendelianRandomisation/Data/Fixed_0.3_",type,"_",i,".txt"), header = TRUE) %>%
        dplyr::mutate_at(c("se.outcome", "beta.outcome"), round, digits = 3) %>%
        dplyr::mutate(eaf.outcome   = round(eaf.outcome, digits = 2)) %>%
        dplyr::mutate(pval.outcome = formatC(pval.outcome, format = "e", digits = 2 )) %>%
        dplyr::mutate(beta.outcome = beta.outcome) %>%
        dplyr::mutate(Outcome = ukb_outcomes[[i]]) %>%
        tibble::as_tibble() %>%
        dplyr::select(-"effect_allele", -"other_allele", -"id") %>%
        dplyr::rename_with(.cols = dplyr::ends_with(".outcome"), .fn = gsub, pattern = ".outcome", replacement = "") %>%
        dplyr::select(-dplyr::ends_with(".exposure")) %>%
        dplyr::select("Outcome", "SNP", "N" = "samplesize", "Effect allele" = "effect_allele",
                      "Other allele" = "other_allele", "Effect allele frequency" = "eaf",
                      "Beta/OR/HR" = "beta", "Standard error" = "se", "P value" = "pval",
        ) %>%
        dplyr::mutate(Type = type) %>%
        dplyr::relocate("Type", .after = N)

      if("outcome_table" %in% ls()){
        outcome_table <- outcome_table %>%
          dplyr::union_all(t)
      }else{
        outcome_table <- t
      }
    }
  }else{
    for(i in c("CAD","MI","IS","Hypertension","T2DM")){
      t <- read.table(paste0(pathData,"Results/Study2/MendelianRandomisation/Data/Fixed_0.3_",type,"_",i,".txt"), header = TRUE) %>%
        dplyr::mutate_at(c("eaf.outcome", "beta.outcome"), round, digits = 2) %>%
        dplyr::mutate(se.outcome   = round(se.outcome, digits = 3)) %>%
        dplyr::mutate(pval.outcome = formatC(pval.outcome, format = "e", digits = 2 )) %>%
        dplyr::mutate(beta.outcome = exp(beta.outcome)) %>%
        dplyr::mutate(Outcome = ukb_outcomes[[i]]) %>%
        tibble::as_tibble() %>%
        dplyr::select(-"effect_allele", -"other_allele", -"id") %>%
        dplyr::rename_with(.cols = dplyr::ends_with(".outcome"), .fn = gsub, pattern = ".outcome", replacement = "") %>%
        dplyr::select(-dplyr::ends_with(".exposure")) %>%
        dplyr::select("Outcome", "SNP","N" = "samplesize", "Effect allele" = "effect_allele",
                      "Other allele" = "other_allele", "Effect allele frequency" = "eaf",
                      "Beta/OR/HR" = "beta", "Standard error" = "se", "P value" = "pval",
        ) %>%
        dplyr::mutate(Type = type) %>%
        dplyr::relocate("Type", .after = N)

      if("outcome_table" %in% ls()){
        outcome_table <- outcome_table %>%
          dplyr::union_all(t)
      }else{
        outcome_table <- t
      }

    }
  }
}

t4 <- tibble::tibble("Outcome" = ukb_outcomes) %>%
  dplyr::mutate(Outcome = as.character(Outcome)) %>%
  dplyr::left_join(outcome_table, by = "Outcome")

write.csv(t4, paste0(pathData,"Results/Tables/t4.csv"))

# ============================================================================ #
#                             Supplementary tables                             #
# ============================================================================ #

# Supplementary table 4 and 5---------------------------------------------------
# Results of the Mendelian randomisation for the GWAS outcomes. -----
rm("outcome_table_cont")
rm("outcome_table_cat")
for(i in gwas_files){
  file <- i

  # Outcome name
  outcome_name <- gsub("_Fixed.*","",file)
  outcome_name <- gsub(".*Results/","",outcome_name)

  t <- read.table(i, header = TRUE) %>%
    dplyr::mutate("Outcome" = gwas_outcomes[outcome_name]) %>%
    dplyr::relocate(U95CI, .after = L95CI)

  if(outcome_name %in% c("bmd","glucose","hba1c","hdl","ldl")){
    if("outcome_table_cont" %in% ls()){
      outcome_table_cont <- outcome_table_cont %>% dplyr::union_all(t)
    }else{
      outcome_table_cont <- t
    }
  }else{
    if("outcome_table_cat" %in% ls()){
      outcome_table_cat <- outcome_table_cat %>% dplyr::union_all(t)
    }else{
      outcome_table_cat <- t
    }
  }
}

st5_continuous <- tibble::tibble("Outcome" = c("Heel bone mineral density", "LDL cholesterol",
                                               "HDL cholesterol", "Fasting glucose", "HbA1c")) %>%
  dplyr::inner_join(outcome_table_cont %>% dplyr::mutate(Outcome = as.character(Outcome)))

st5_categorical <- tibble::tibble("Outcome" = c("Hip fracture", "Coronary artery disease",
                                                "Myocardial infarction", "Ischaemic stroke",
                                                "Hypertension", "Type 2 diabetes mellitus")) %>%
  dplyr::inner_join(outcome_table_cat %>% dplyr::mutate(Outcome = as.character(Outcome)))


write.csv(st5_continuous, paste0(pathData,"Results/Tables/st5_continuous.csv"))
write.csv(st5_categorical, paste0(pathData,"Results/Tables/st5_categorical.csv"))

# Supplementary table 6 --------------------------------------------------------
# Results of the Mendelian randomisation for the UK Biobank outcomes.  -----
rm("outcome_table") #------- Continuous data -------- #
for(i in cont_files){
  # Logistic
  file <- i

  # Outcome name
  outcome_name <- gsub(".*Logistic_","",file)
  outcome_name <- gsub(".txt.*","",outcome_name)

  # Load outcome name
  t <- read.table(i, header = TRUE) %>%
    dplyr::mutate("Outcome" = ukb_outcomes[outcome_name]) %>%
    dplyr::relocate(U95CI, .after = L95CI)

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

st6 <- tibble::tibble(Outcome = c("Cholesterol", "LDL cholesterol",
                                  "HDL cholesterol",  "Apolipoprotein-A",
                                  "Apolipoprotein-B", "Lipoprotein (a)",
                                  "C-Reactive protein", "Triglycerides",
                                  "HbA1c", "Glucose", "Heel bone mineral density")) %>%
  dplyr::inner_join(
    outcome_table %>%
      dplyr::mutate(Outcome = as.character(Outcome)),
    by = "Outcome"
  )
write.csv(st6, paste0(pathData,"Results/Tables/st6_continuous.csv"))

rm("outcome_table") # ------- Categorical data -------- #
for(i in cat_files){
  file <- i

  # Outcome name
  outcome_name <- gsub(".*Logaritmic_","",file)
  outcome_name <- gsub(".txt.*","",outcome_name)

  # Load outcome name
  t <- read.table(file, header = TRUE) %>%
    dplyr::mutate("Outcome" = gsub("_.*","",Outcome)) %>%
    dplyr::mutate("Outcome" = ukb_outcomes[outcome_name]) %>%
    dplyr::relocate(U95CI, .after = L95CI)

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

st6 <- tibble::tibble(Outcome = c("Coronary artery disease","Myocardial infarction",
                                  "Ischaemic stroke","Hypertension","Type 2 diabetes")) %>%
  dplyr::inner_join(
    outcome_table %>% dplyr::mutate(Outcome = as.character(Outcome)),
    by = "Outcome"
  )
write.csv(st6, paste0(pathData,"Results/Tables/st6_categorical.csv"))


rm("outcome_table") # ------- Survival data -------- #
for(i in surv_files){
  file <- i

  # Outcome name
  outcome_name <- gsub(".*Birth_","",file)
  outcome_name <- gsub(".txt.*","",outcome_name)

  # Load outcome name
  t <- read.table(file, header = TRUE) %>%
    dplyr::mutate("Outcome" = gsub("_.*","",Outcome)) %>%
    dplyr::mutate("Outcome" = ukb_outcomes[outcome_name]) %>%
    dplyr::relocate(U95CI, .after = L95CI) %>%
    dplyr::rename("HR" = "OR")

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

st6 <- tibble::tibble(Outcome = c("Coronary artery disease","Myocardial infarction",
                                  "Ischaemic stroke","Hypertension","Type 2 diabetes")) %>%
  dplyr::inner_join(
    outcome_table %>% dplyr::mutate(Outcome = as.character(Outcome)),
    by = "Outcome"
  )
write.csv(st6, paste0(pathData,"Results/Tables/st6_survival.csv"))

# Supplementary table 7 --------------------------------------------------------
# Colocalization results within ±20kb SOST gene.
st7 <- read.table(paste0(pathData,"Results/Colocalization/colocalisation.txt"), header = TRUE) %>%
  dplyr::filter(exposure_method == "Fixed") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(H0 = H0*100) %>%
  dplyr::mutate(H1 = H1*100) %>%
  dplyr::mutate(H2 = H2*100) %>%
  dplyr::mutate(H3 = H3*100) %>%
  dplyr::mutate(H4 = H4*100) %>% dplyr::distinct()
write.csv(st7, paste0(pathData,"Results/Tables/st7.csv"))

# Supplementary table 8 --------------------------------------------------------
# Sensitivity analyses results
order <- unlist(gwas_outcomes) %>%
  dplyr::as_tibble() %>%
  dplyr::rename("Outcome" = "value") %>%
  dplyr::mutate("Source" = "GWAS") %>%
  dplyr::union_all(
    unlist(ukb_outcomes) %>%
      dplyr::as_tibble() %>%
      dplyr::rename("Outcome" = value) %>%
      dplyr::mutate("Source"  = "UK Biobank") %>%
      dplyr::mutate("Source"  = dplyr::if_else(Outcome %in% c("Coronary artery disease","Myocardial infarction",
                                                              "Ischaemic stroke","Hypertension", "Type 2 diabetes"),
                                               "UK Biobank (categorical)",
                                               Source)) %>%
      dplyr::union_all(
        unlist(ukb_outcomes) %>%
          dplyr::as_tibble() %>%
          dplyr::rename("Outcome" = value) %>%
          dplyr::mutate("Source"  = "UK Biobank") %>%
          dplyr::filter(Outcome %in% c("Coronary artery disease","Myocardial infarction",
                                       "Ischaemic stroke","Hypertension", "Type 2 diabetes")) %>%
          dplyr::mutate("Source"  = "UK Biobank (survival)")
      )
  )

# Step-wise pruning - Gwas
files <- list.files(path = paste0(pathData,"Results/Study1/MendelianRandomisation/Results"),
                    recursive  = FALSE,
                    pattern    = "Fixed",
                    full.names = TRUE)
files <- files %>%
  dplyr::as_tibble() %>%
  dplyr::filter(!grepl("Enrol", files) & !grepl("0.3", files) & !grepl("PCA", files)) %>% print(n = 200)


rm("outcome_table")
for(i in files$value){
  # r2 threshold
  r2 <- gsub(".*_Fixed_","",i)
  r2 <- gsub("_.*","",r2)

  # Outcome name
  outcome_name <- gsub(".*Results/","",i)
  outcome_name <- gsub("_.*", "",outcome_name)

  outc <- read.table(i, header = TRUE) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Sensitivity" = "Step-wise pruning",
                  "Outcome"     = as.character(gwas_outcomes[[outcome_name]])) %>%
    dplyr::mutate("Source" = "GWAS") %>%
    dplyr::mutate("r2" = r2) %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2", "Instruments", "Beta", "SE", "Pval")

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>% dplyr::union_all(outc)
  }else{
    outcome_table <- outc
  }
}

# Step-wise pruning - UK Biobank
files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation/"),
                    recursive  = FALSE,
                    pattern    = "Fixed",
                    full.names = TRUE)
files <- files %>%
  dplyr::as_tibble() %>%
  dplyr::filter(!grepl("Enrol", files) & !grepl("0.3", files) & !grepl("Fracture", files) & !grepl("ebmd", files))

for(i in files$value){
  # r2 threshold
  r2 <- gsub(".*Fixed_","",i)
  r2 <- gsub("_.*","",r2)

  r2_val <- as.numeric(r2)

  # Regression
  regression <- gsub(paste0('.*',r2,'_'), "",i)
  regression <- gsub("_.*", "",regression)

  # Outcome name
  outcome_name <- gsub(paste0('.*',regression,"_"),"",i)
  outcome_name <- gsub(".txt", "",outcome_name)

  outc <- read.table(i, header = TRUE) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Sensitivity" = "Step-wise pruning",
                  "Outcome"     = as.character(ukb_outcomes[[outcome_name]]),
                  "Type"        = regression) %>%
    dplyr::mutate("Source" = dplyr::case_when(
      Type == "Birth" ~ "UK Biobank (survival)",
      Type == "Logaritmic" ~ "UK Biobank (categorical)",
      Type == "Logistic" ~ "UK Biobank")
    ) %>%
    dplyr::mutate("r2" = r2) %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2", "Instruments", "Beta", "SE", "Pval")

  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>% dplyr::union_all(outc)
  }else{
    outcome_table <- outc
  }
}

st8_1 <- order %>%
  dplyr::inner_join(outcome_table, by = c("Outcome", "Source")) %>%
  dplyr::relocate(Sensitivity)

# Random effects method - gwas
files <- list.files(path = paste0(pathData,"Results/Study1/MendelianRandomisation/Results"),
                    recursive  = FALSE,
                    pattern    = "Random",
                    full.names = TRUE)
rm("outcome_table")
for(i in files){
  # Outcome name
  outcome_name <- gsub(".*tion/Results/","",i)
  outcome_name <- gsub("_Random.*","",outcome_name)

  # Load outcome name
  t <- read.table(i, header = TRUE) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Outcome" = as.character(gwas_outcomes[[outcome_name]]),
                  "Sensitivity" = "Random-effects method",
                  "Source" = "GWAS",
                  "r2" = "-") %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2","Instruments", "Beta","SE","Pval")


  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation"),
                    recursive  = FALSE,
                    pattern    = "Random",
                    full.names = TRUE)
files <- files %>%
  dplyr::as_tibble() %>%
  dplyr::filter(!grepl("Enrol", files) & !grepl("Fracture", files) & !grepl("ebmd", files))

for(i in files$value){
  # Regression
  regression  <- gsub(".*0.3_","",i)
  regression  <- gsub("_.*","",regression)

  # Outcome name
  outcome_name <- gsub(paste0(".*",regression,"_"),"",i)
  outcome_name <- gsub(".txt","",outcome_name)

  # Load outcome name
  t <- read.table(i, header = TRUE) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Outcome" = as.character(ukb_outcomes[[outcome_name]]),
                  "Sensitivity" = "Random-effects method",
                  "r2" = "-",
                  "Type" = regression) %>%
    dplyr::mutate("Source" = dplyr::case_when(
      Type == "Birth" ~ "UK Biobank (survival)",
      Type == "Logaritmic" ~ "UK Biobank (categorical)",
      Type == "Logistic" ~ "UK Biobank")
    ) %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2","Instruments", "Beta","SE","Pval")

    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
}

st8_2 <- order %>%
  dplyr::inner_join(outcome_table %>% dplyr::mutate(Outcome = as.character(Outcome)), by = c("Outcome", "Source")) %>%
  dplyr::relocate(Sensitivity)

# Principal components
files <- list.files(path = paste0(pathData,"Results/Study1/MendelianRandomisation/Results/"),
                    recursive  = FALSE,
                    pattern    = "PCA",
                    full.names = TRUE)
rm("outcome_table")
for(i in files){
  # Outcome name
  outcome_name <- gsub(".*tion/Results/","",i)
  outcome_name <- gsub("_.*","",outcome_name)

  # Load outcome name
  t <- readr::read_delim(i) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Outcome" = as.character(gwas_outcomes[outcome_name]),
                  "Sensitivity" = "Principal components analysis",
                  "Source" = "GWAS",
                  "r2" = "-",
                  "Instruments" = 76) %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2","Instruments", "Beta","SE","Pval")


  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}

st8_3 <- order %>%
  dplyr::filter(Source == "GWAS") %>%
  dplyr::inner_join(outcome_table, by = c("Outcome", "Source")) %>%
  dplyr::relocate(Sensitivity)

# Since UK Biobank enrolment
files <- list.files(path = paste0(pathData,"Results/Study2/MendelianRandomisation"),
                    recursive  = FALSE,
                    pattern    = "Fixed_0.3",
                    full.names = TRUE)
files <- files %>%
  dplyr::as_tibble() %>%
  dplyr::filter(grepl("Enrol", files) & !grepl("Fracture",files))

rm("outcome_table")
for(i in files$value){
  # Regression
  regression  <- gsub(".*0.3_","",i)
  regression  <- gsub("_.*","",regression)

  # Outcome name
  outcome_name <- gsub(paste0(".*",regression,"_"),"",i)
  outcome_name <- gsub(".txt","",outcome_name)

  # Load outcome name
  t <- read.table(i, header = TRUE) %>%
    dplyr::rename_with(~gsub("OR","Beta", .x), any_of("OR")) %>%
    dplyr::mutate("Outcome" = as.character(ukb_outcomes[outcome_name]),
                  "Sensitivity" = "Survival analysis since UK Biobank enrolment",
                  "r2" = "-",
                  "Source" = "UK Biobank (survival)") %>%
    dplyr::select("Sensitivity", "Outcome", "Source", "r2","Instruments", "Beta","SE","Pval")


  if("outcome_table" %in% ls()){
    outcome_table <- outcome_table %>%
      dplyr::union_all(t)
  }else{
    outcome_table <- t
  }
}
st8_4 <- order %>%
  dplyr::filter(Source == "UK Biobank (survival)") %>%
  dplyr::inner_join(outcome_table, by = c("Outcome", "Source")) %>%
  dplyr::relocate(Sensitivity)


st8 <- st8_1 %>%
  dplyr::union_all(st8_2) %>%
  dplyr::union_all(st8_3) %>%
  dplyr::union_all(st8_4) %>% print(n = 200)

st8 <- st8 %>%
  dplyr::mutate(N = dplyr::case_when(
    Outcome == "Heel bone mineral density" & Source == "GWAS" ~ "426,824",
    Outcome == "Hip fracture"              & Source == "GWAS" ~ "735,354",
    Outcome == "Coronary artery disease"   & Source == "GWAS" ~ "1,165,690",
    Outcome == "Myocardial infarction"     & Source == "GWAS" ~ "~639,000",
    Outcome == "Ischaemic stroke"          & Source == "GWAS" ~ "1,847,683",
    Outcome == "Hypertension"              & Source == "GWAS" ~ "462,933",
    Outcome == "Type 2 diabetes"           & Source == "GWAS" ~ "933,970",
    Outcome == "LDL cholesterol"           & Source == "GWAS" ~ "1,228,508",
    Outcome == "HDL cholesterol"           & Source == "GWAS" ~ "1,244,439",
    Outcome == "Fasting glucose"           & Source == "GWAS" ~ "178,455",
    Outcome == "HbA1c"                     & Source == "GWAS" ~ "132,400",
    Outcome == "Cholesterol"               & Source == "UK Biobank" ~ "263,285",
    Outcome == "LDL cholesterol"           & Source == "UK Biobank" ~ "262,794",
    Outcome == "HDL cholesterol"           & Source == "UK Biobank" ~ "241,065",
    Outcome == "Triglycerides"             & Source == "UK Biobank" ~ "263,067",
    Outcome == "Apolipoprotein-A"          & Source == "UK Biobank" ~ "239,707",
    Outcome == "Apolipoprotein-B"          & Source == "UK Biobank" ~ "261,986",
    Outcome == "C-Reactive protein"        & Source == "UK Biobank" ~ "262,728",
    Outcome == "Lipoprotein (a)"           & Source == "UK Biobank" ~ "209,333",
    Outcome == "Glucose"                   & Source == "UK Biobank" ~ "240,909",
    Outcome == "Cholesterol"               & Source == "UK Biobank" ~ "263,285",
    Outcome == "HbA1c"                     & Source == "UK Biobank" ~ "263,079",
    Outcome == "Coronary artery disease"   & Source == "UK Biobank (categorical)" ~ "276,172",
    Outcome == "Myocardial infarction"     & Source == "UK Biobank (categorical)" ~ "276,065",
    Outcome == "Ischaemic stroke"          & Source == "UK Biobank (categorical)" ~ "276,172",
    Outcome == "Hypertension"              & Source == "UK Biobank (categorical)" ~ "271,359",
    Outcome == "Type 2 diabetes"           & Source == "UK Biobank (categorical)" ~ "266,268",
    Outcome == "Coronary artery disease"   & Source == "UK Biobank (survival)" ~ "276,172",
    Outcome == "Myocardial infarction"     & Source == "UK Biobank (survival)" ~ "276,172",
    Outcome == "Ischaemic stroke"          & Source == "UK Biobank (survival)" ~ "276,172",
    Outcome == "Hypertension"              & Source == "UK Biobank (survival)" ~ "276,172",
    Outcome == "Type 2 diabetes"           & Source == "UK Biobank (survival)" ~ "276,172"
  ))

write.csv(st8, paste0(pathData,"Results/Tables/st8.csv"))

# Supplementary table 9 --------------------------------------------------------
# Study information of the cohorts used for the GWAS meta-analysis of sclerostin.
st9 <- tibble::tibble(
  "Article reference" = c("Sun et al., 2018", "Ferkingstad et al., 2021", "Pietzner et al., 2021"),
  "Outcome"           = rep(c("Sclerostin"),3),
  "N"                 = c("3,301", "35,559", "10,708"),
  "Age (SD)"          = c("43 (14)", "55 (17)", "49 (7)"),
  "% Female"          = c("49", "57","53"),
  "Study population"  = c("European (England)", "European (Iceland)", "European (England)")
)
write.csv(st9, paste0(pathData,"Results/Tables/st9.csv"))

# Supplementary table 10 -------------------------------------------------------
# Study information of the cohorts used in the published GWAS summary statistics of the outcomes.
means <- list("bmd" = "0.54 (0.12) [g/cm2]",
              "hf"  = "11,516 (2)",
              "cad" = "181,522 (16)" ,
              "mi"  = "~62,000 (10)",
              "is"  = "612,875 (33)",
              "hypertension" = "119,731 (26)",
              "t2dm"= "80,154 (9)",
              "hdl" = "51.9 [mg/dL]",
              "ldl" = "131 [mg/dL]",
              "glucose" = "-",
              "hba1c" = "5.23 [mmol/l]")
rm("outcome_table")
for(i in gwas_files){
  file <- i

  # Outcome name
  outcome_name <- gsub("_Fixed.*","",file)
  outcome_name <- gsub(".*Results/","",outcome_name)

  # Load outcome name
  outcome_gwas <- loadOutcomeGwas(outcome_name, iv)

  if(outcome_name %in% c("bmd","hdl","ldl","glucose","hba1c")){
    t <- read.table(file, header = TRUE) %>%
      dplyr::mutate(`Mean (SD)/N (%)` = means[[outcome_name]],
                    "N"       = max(outcome_gwas$samplesize.outcome %>% unique()),
                    "Outcome" = gwas_outcomes[[outcome_name]]) %>%
      dplyr::select("Outcome",
                    "N",
                    `Mean (SD)/N (%)`,
                    "N SNPs" = "Instruments")

    if("outcome_table" %in% ls()){

      outcome_table <- outcome_table %>%
        dplyr::union_all(t)

    }else{
      outcome_table <- t
    }
  }else{
    t <- read.table(file, header = TRUE) %>%
      dplyr::mutate(`Mean (SD)/N (%)` = means[[outcome_name]],
                    "N"      = max(outcome_gwas$samplesize.outcome %>% unique()),
                    "Outcome" = gwas_outcomes[[outcome_name]]) %>%
      dplyr::select("Outcome",
                    "N",
                    `Mean (SD)/N (%)`,
                    "N SNPs" = "Instruments")

    if("outcome_table" %in% ls()){
      outcome_table <- outcome_table %>%
        dplyr::union_all(t)
    }else{
      outcome_table <- t
    }
  }
}

st10 <- tibble::tibble(Outcome = c("Heel bone mineral density",
                                   "Hip fracture",
                                   "Coronary artery disease",
                                   "Myocardial infarction",
                                   "Ischaemic stroke",
                                   "Hypertension",
                                   "Type 2 diabetes mellitus",
                                   "LDL cholesterol",
                                   "HDL cholesterol",
                                   "Fasting glucose",
                                   "HbA1c")) %>%
  dplyr::inner_join(
    outcome_table,
    by = "Outcome"
  )

write.csv(st10, paste0(pathData,"Results/Tables/st10.csv"))


# Supplementary table 11 -------------------------------------------------------
# Study information of the UK Biobank cohorts.
pop <- read.table(paste0(pathData,"Results/Study2/Cohort/ukb.txt"), header = TRUE) %>%
  tibble::as_tibble()

rm("outcome_table")
outcome_table <- pop %>%
  dplyr::summarise(
    "N" = dplyr::n(),
    "Age (SD)" = paste0(round(mean(f.21003.0.0, na.rm = TRUE)), " (", round(sd(f.21003.0.0, na.rm = TRUE)), ")"),
    "Male (%)" = paste0(sum(sex, na.rm = TRUE), " (", round(sum(sex, na.rm = TRUE)/dplyr::n()*100), ")"),
    "Cases (%)" = "-",
    "Outcome"   = "Whole cohort",
    "Regression" = "-",
    "Mean followup" = "-"
  )

stats <- read.table(paste0(pathData,"Results/Study2/Phenotyping/Statistics_Continuous.txt"), header = TRUE)

for(i in c("cholesterol", "ldl", "hdl", "triglycerides","apoA","apoB", "crp", "lipoprotein","glucose","hba1c")){
  stats_i <- stats %>% dplyr::filter(outcome == i)
  t <- read.table(paste0(pathData,"Results/Study2/Phenotyping/Phenotype_Logistic_",i,".txt"), header = TRUE) %>%
    dplyr::select("eid","outcome") %>%
    dplyr::inner_join(pop) %>%
    dplyr::summarise(
      "N"        = dplyr::n(),
      "Age (SD)" = paste0(round(mean(f.21003.0.0, na.rm = TRUE)), " (", round(sd(f.21003.0.0, na.rm = TRUE)), ")"),
      "Male (%)" = paste0(sum(sex, na.rm = TRUE), " (", round(sum(sex, na.rm = TRUE)/dplyr::n()*100),  ")"),
      "Cases (%)" = paste0(round(mean(stats_i$Mean, na.rm = TRUE), digits = 2), " (", round(stats_i$SD, digits = 2), ")"),
      "Mean followup" = "-",
      "Regression" = "Logistic"
    ) %>%
    dplyr::mutate("Outcome" = ukb_outcomes[i])

  outcome_table <- outcome_table %>%
    dplyr::union_all(t  %>% dplyr::mutate(Outcome = as.character(Outcome)))
}

for(type in c("Logaritmic","Birth")){
  for(i in c("CAD","MI","IS","Hypertension","T2DM")){
    t <- read.table(paste0(pathData,"Results/Study2/Phenotyping/Phenotype_",type,"_",i,".txt"), header = TRUE) %>%
      dplyr::select("eid","state") %>%
      dplyr::inner_join(pop) %>%
      dplyr::summarise(
        "N"        = dplyr::n(),
        "Age (SD)" = paste0(round(mean(f.21003.0.0, na.rm = TRUE)), " (", round(sd(f.21003.0.0, na.rm = TRUE)), ")"),
        "Male (%)" = paste0(sum(sex, na.rm = TRUE), " (", round(sum(sex, na.rm = TRUE)/dplyr::n()*100),  ")"),
        "Cases (%)" = paste0(round(sum(state, na.rm = TRUE)), " (", round(sum(state, na.rm = TRUE)/dplyr::n()*100), ")")
      ) %>%
      dplyr::mutate("Outcome" = ukb_outcomes[i],
                    "Regression" = type)

    if(type == "Logaritmic"){
      t <- t %>% dplyr::mutate("Mean followup" = "-")
    }else{
      t <- t %>% dplyr::mutate(
        "Mean followup" = read.table(paste0(pathData,"Results/Study2/Phenotyping/Phenotype_",type,"_",i,".txt"), header = TRUE) %>%
          dplyr::select("age") %>%
          dplyr::summarise("a" = paste0(round(mean(age,na.rm = TRUE)), " (", round(sd(age, na.rm = TRUE)), ")")) %>%
          dplyr::pull(a)
      )
    }

    outcome_table <- outcome_table %>%
      dplyr::union_all(t %>% dplyr::mutate(Outcome = as.character(Outcome)))
  }
}

st11 <- outcome_table %>%
  dplyr::relocate("Outcome") %>%
  dplyr::relocate("Regression", .after = "Outcome")

write.csv(st11, paste0(pathData,"Results/Tables/st11.csv"))

# Supplementary table 1 ------
fixed <- tibble::as_tibble(read.delim(paste0(pathData,"Results/Study1/MetaAnalysis/FixedMapped.txt")))
random <- tibble::as_tibble(read.delim(paste0(pathData,"Results/Study1/MetaAnalysis/Random.txt")))

fixed |>
  dplyr::select(-c(i2,samplesize,q_statistic,q_pvalue,zscore,eaf)) |>
  dplyr::rename("fixed_beta" = "beta", "fixed_se" = "se", "fixed_pval" = "pval") |>
  dplyr::inner_join(
    random |> dplyr::select(SNP, effect_allele, other_allele, "random_beta" = "beta", "random_se" = "se", "random_pval" = "pval",
                            i2, q_statistic, q_pvalue, samplesize)
  ) |>
  dplyr::select(SNP, effect_allele, other_allele, fixed_beta, random_beta, fixed_se,
                random_se, fixed_pval, random_pval, i2, q_statistic, q_pvalue, samplesize)  |>
  write.csv(paste0(pathData,"Results/Tables/st1.csv"))

R.utils::copyDirectory(paste0(pathData,"Results/Tables"),
              "C:/Users/martaa/OneDrive - Nexus365/Marta/Projects/Sclerostin/ManucriptVersions/FirstReview/Tables",overwrite = TRUE)
