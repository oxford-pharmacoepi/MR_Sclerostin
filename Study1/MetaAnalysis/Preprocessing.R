# ============================================================================ #
#                                Preprocessing                                 #
#                            Marta Alcalde-Herraiz                             #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Data curation of the GWAS insturment-exposure                                #
# ============================================================================ #
dir.create(paste0(pathResults,'MetaAnalysis'))

recordAttrition <- function(Name, Counts, Reason){
  attrition <- tibble::tibble(Counts = .env$Counts,
                              Reason = .env$Reason)
  return(attrition)
}

qc <- function(gwas, name){
  attrition <- recordAttrition(.env$name, gwas %>% nrow(), "Number of SNPs in the GWAS")

  gwas <- gwas %>% dplyr::filter(CHR == 17)
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Number of SNPs in chromosome 17"))

  gwas <- gwas %>% dplyr::filter(EA != "D") %>% dplyr::filter(NEA != "I") %>% dplyr::filter(NEA != "!")
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with alleles equal to ´D´, ´I´, or ´!´"))

  gwas <- gwas %>% dplyr::filter(base::nchar(EA) == 1, base::nchar(NEA) == 1)
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove non biallelic SNPs"))

  gwas <- gwas %>% dplyr::filter(base::substr(MARKERNAME,1,2) == "rs") # Discard those SNPs with a non valid ID
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with a non-valid ID"))

  gwas <- gwas %>% dplyr::filter(EA != NEA)
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with the same effect allele and other allele"))

  if((gwas %>% dplyr::summarise(x = is.na(EAF)) %>% dplyr::distinct() %>% dplyr::pull()) == FALSE){
    gwas <- gwas %>% dplyr::filter(EAF != 0 | EAF != 1)
    attrition <- attrition %>%
      dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with allele frequency equal to 0 or 1"))
  }

  gwas <- gwas %>% dplyr::filter(pval >= 0 | pval <= 1)
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with P-Values < 0 or > 1"))

  gwas <- gwas %>% dplyr::filter(SE > 0 | SE != Inf)
  attrition <- attrition %>%
    dplyr::union_all(recordAttrition(.env$name, gwas %>% nrow(), "Remove SNPs with standard errors <= 0 or equal to infinity"))

  attrition <- attrition %>%
    dplyr::mutate("Removed" = c(NA, (attrition[1:nrow(attrition)-1,"Counts"] %>% dplyr::pull() - attrition[2:nrow(attrition),"Counts"] %>% dplyr::pull())))

  return(list("attrition" = attrition, "gwas" = gwas))
}


# Ferkingstad - Nature genetics gwas -------------------------------------------
# Information about the GWAS:
# Article: https://www.nature.com/articles/s41588-021-00978-w
# Download GWAS: https://www.decode.com/summarydata/
# Assembly: hg38
# Positive strand

# This gwas has the data splitted between two .txt
Ferkingstad_NatureGenetics1 <- tibble::as_tibble(readr::read_delim(paste0(pathData,'SummaryStatistics/Exposure/Ferkingstad_NatureGenetics/assocvariants.annotated.txt/assocvariants.annotated.txt')))
Ferkingstad_NatureGenetics2 <- tibble::as_tibble(readr::read_delim(paste0(pathData,'SummaryStatistics/Exposure/Ferkingstad_NatureGenetics/13101_60_SOST_SOST.txt/13101_60_SOST_SOST.txt')))

Ferkingstad_NatureGenetics <- Ferkingstad_NatureGenetics1 %>%
  dplyr::select('CHR' = 'Chrom',
         'Pos' = 'Pos', # Position of the SNP (GRCh38/hg38)
         'Name', # Name of the snp in the following form:
         'MARKERNAME' = 'rsids', # SNP Id
         "EA" = "effectAllele", # Effect allele
         "NEA" = "otherAllele", # Other allele
         "EAF" = "effectAlleleFreq") %>% # Effect allele frequency
  dplyr::inner_join(
    Ferkingstad_NatureGenetics2 %>%
      dplyr::select(
        'CHR' = 'Chrom',
        'Pos' = 'Pos', # Position of the SNP (GRCh38/hg38)
        'Name', # Name of the snp in the following form:
        'MARKERNAME' = 'rsids', # SNP Id
        "EA" = "effectAllele", # Effect allele
        "NEA" = "otherAllele", # Other allele
        "BETA" = "Beta", # Effect size
        "SE" = "SE", # Standard error of the effect size
        "pval" = "Pval", # Pvalue of the effect size
        "N" = "N", # Sample size
      ),
    by = c('CHR','Pos','Name','MARKERNAME','EA','NEA')
  ) %>%
  dplyr::mutate(dplyr::across('CHR',stringr::str_replace,'chr','')) %>%
  dplyr::mutate(CHR = as.numeric(CHR))

Ferkingstad_NatureGenetics_QC <- qc(Ferkingstad_NatureGenetics, "Ferkingstad_NatureGenetics")

readr::write_tsv(tibble::as_tibble(Ferkingstad_NatureGenetics_QC[[1]]),paste0(pathResults,'MetaAnalysis/Ferkingstad_NatureGenetics_Attrition.txt'))
readr::write_tsv(Ferkingstad_NatureGenetics_QC[[2]],paste0(pathResults,'MetaAnalysis/Ferkingstad_NatureGenetics_GWAS.txt'))

# Pietzner - Science gwas ------------------------------------------------------
# Information about the GWAS:
# Article: https://www.science.org/doi/10.1126/science.abj1541
# Download GWAS: http://www.omicscience.org/data_usage_agreement.php
# Assembly: (GRCh37 coordinates)
# Positive strand
Pietzner_Science <- tibble::as_tibble(readr::read_delim(paste0(pathData,'SummaryStatistics/Exposure/Pietzner_Science/Science_GWAS.txt'))) #hg19 coordinates
Pietzner_Science <- Pietzner_Science %>%
  dplyr::select('CHR' = 'chr',
         "Pos"  = "pos", # Position of the SNP (GRCh37/hg19)
         "MARKERNAME" = "rsid", # SNP Id
         "EA"   = "Allele1", # Effect allele
         "NEA"  = "Allele2", # Other allele
         "EAF"  = "Freq1", # Effect allele frequency
         "BETA" = "Effect", # Effect size
         "SE"   = "StdErr", # Standard error of the effect size
         "pval" = "Pvalue", # Pvalue of the effect size
         "N"    = "TotalSampleSize") %>% # Sample size
  dplyr::mutate(EA = toupper(EA), NEA = toupper(NEA))

Pietzner_Science_QC <- qc(Pietzner_Science, "Pietzner_Science")

readr::write_tsv(tibble::as_tibble(Pietzner_Science_QC[[1]]),paste0(pathResults,'MetaAnalysis/Pietzner_Science_Attrition.txt'))
readr::write_tsv(Pietzner_Science_QC[[2]],paste0(pathResults,'MetaAnalysis/Pietzner_Science_GWAS.txt'))

# Sun_Nature gwas --------------------------------------------------------------
# Information about the GWAS:
# Article: https://www.nature.com/articles/s41586-018-0175-2
# Download GWAS: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90242001-GCST90243000/GCST90242722/
# Assembly: hg38
Sun_Nature <- tibble::as_tibble(readr::read_delim(paste0(pathData,"SummaryStatistics/Exposure/Sun_Nature/GCST90242722.h.tsv/SOST.13101.60.3.tsv"))) #hg38
Sun_Nature <- Sun_Nature %>%
  dplyr::select("CHR" = "hm_chrom",
         "MARKERNAME" = "hm_rsid", # SNP Id
         "EA"   = "hm_effect_allele", # Effect allele
         "NEA"  = "hm_other_allele", # Other allele
         "BETA" = "hm_beta", # Effect size
         "SE"   = "standard_error", # Standard error of the effect size
         "pval" = "p_value", # Pvalue of the effect size
         "Pos"  = "hm_pos", # Position of the SNP (GRCh37/hg19)
         "EAF"  = "hm_effect_allele_frequency") %>%   # Effect allele frequency
  dplyr::mutate(N   = 3301)

Sun_Nature_QC <- qc(Sun_Nature)

readr::write_tsv(tibble::as_tibble(Sun_Nature_QC[[1]]),paste0(pathResults,'MetaAnalysis/Sun_Nature_Attrition.txt'))
readr::write_tsv(Sun_Nature_QC[[2]],paste0(pathResults,'MetaAnalysis/Sun_Nature_GWAS.txt'))

