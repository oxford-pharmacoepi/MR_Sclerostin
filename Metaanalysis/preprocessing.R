# ============================================================================ #
#                                PREPROCESSING                                 #
# Summary:                                                                     #
# SNPs from chromosome 17 are selected from the three different GWAS           #
#   - Nature genetics gwas                                                     #
#   - Science gwas                                                             #
#   - Nature gwas                                                              #
# It creates the files that will be used in the meta-analysis                  #
# ============================================================================ #
rm(list=ls())
library(pacman)
pacman::p_load('tibble','dplyr','tidyverse','readr','here')

# Read GWAS table
sc  <- as_tibble(read.delim(here("Data","Science","Science_GWAS.txt")))
ng1 <- as_tibble(read.delim(here("Data","NatureGenetics","NatureGenetics_GWAS1.txt")))
ng2 <- as_tibble(read.delim(here("Data","NatureGenetics","NatureGenetics_GWAS2.txt")))
nn  <- as_tibble(read.delim(here("Data","Nature","Nature_GWAS.tsv")))

# NATURE GENETICS GWAS ---------------------------------------------------------
# This gwas has the data splitted between two .txt, so first we select the
# data that is of interest for the meta-analysis.
ng1 <- ng1 %>% filter(Chrom == "chr17") %>% # Select only the SNPs in chromosome 17
  select("Pos", # Position of the SNP (GRCh38/hg38)
         "MARKERNAME" = "rsids", # SNP Id
         "Name.ng" = "Name", # Name of the snp in the following form: 
         "EA" = "effectAllele", # Effect allele 
         "NEA" = "otherAllele", # Other allele
         "EAF" = "effectAlleleFreq") # Effect allele frequency

ng2 <- ng2 %>% filter(Chrom == "chr17") %>% # Select only the SNPs in chromosome 17
  select("Pos", # Position of the SNP (GRCh38/hg38)
         "BETA" = "Beta", # Effect size
         "SE" = "SE", # Standard error of the effect size
         "pval" = "Pval", # Pvalue of the effect size
         "N" = "N", # Sample size
         "Name.ng" = "Name") # Name of the snp in the following form: 

# Save the results without the preprocessing modifications
write_tsv(ng1,here("Metaanalysis","ng1_NP.txt"))
write_tsv(ng2,here("Metaanalysis","ng2_NP.txt"))

# Preprocessing:
ng <- as_tibble(read.delim(here("Metaanalysis","ng1_NP.txt"))) %>% 
  filter(nchar(EA) == 1, nchar(NEA) == 1, # Discard those SNPs that do not have a single nucleotide variation
         MARKERNAME != ".") %>% # Discard those SNPs that do not have an ID
  inner_join(read.delim(here("Metaanalysis","ng2_NP.txt")), # Read the second file and merge it with the first one
    by = c("Name.ng","Pos")
  ) %>%
  mutate(NEA = if_else(NEA == "!",substr(Name.ng, nchar(Name.ng), nchar(Name.ng)), NEA)) %>%  # Some SNPs had this symbol "!" in the NEA column. I substitute this symbol with the last letter of the variable Name.
  filter(EA != NEA) # Delete those SNPs with the same Effect allele and other allele.
write_tsv(ng,here("Metaanalysis","ng.txt")) # Save the data.


# SCIENCE GWAS  ----------------------------------------------------------------
sc1 <- sc %>% filter(chr == 17) %>% # Select only the SNPs in chromosome 17
  select("Pos" = "pos", # Position of the SNP (GRCh37/hg19)
         "MARKERNAME" = "rsid", # SNP Id
         "Name.sc" = "MarkerName", # Name of the snp in the following form: 
         "EA" = "Allele1", # Effect allele
         "NEA" = "Allele2", # Other allele
         "EAF" = "Freq1", # Effect allele frequency
         "BETA" = "Effect", # Effect size
         "SE" = "StdErr", # Standard error of the effect size
         "pval" = "Pvalue", # Pvalue of the effect size
         "N" = "TotalSampleSize") %>% # Sample size
  mutate(EA = toupper(EA), NEA = toupper(NEA)) 

# Save the results without the preprocessing:
write_tsv(sc1,here("Metaanalysis","sc_NP.txt"))

# Preprocessing
sc <- as_tibble(read.delim(here("Metaanalysis","sc_NP.txt"))) %>% 
  filter(EA != "D", NEA != "I", # Delete those SNPs with EA = "D" andNEA = "I"
         substr(MARKERNAME,1,2) == "rs") %>% # Delete those SNPs that do not have an ID written properly.
  select(-"Name.sc")
write_tsv(sc,here("Metaanalysis","sc.txt"))


# NATURE GWAS ------------------------------------------------------------------
nn1 <-  nn %>% 
  filter(chromosome == 17) %>%
  select("Pos" = "hm_pos", # Position of the SNP (GRCh37/hg19)
         "MARKERNAME" = "hm_rsid", # SNP Id
         "Name" = "hm_variant_id",  # Name of the snp in the following form: 
         "EA" = "hm_effect_allele", # Effect allele
         "NEA" = "hm_other_allele", # Other allele
         "BETA" = "hm_beta", # Effect size
         "SE" = "standard_error", # Standard error of the effect size
         "pval" = "p_value", # Pvalue of the effect size
         "EAF" = "hm_effect_allele_frequency") # Effect allele frequency

# Save the results without the preprocessing:
write_tsv(nn1,here("Metaanalysis","nn_NP.txt"))

# Preprocessing
nn1 <- as_tibble(read.delim(here("Metaanalysis","nn_NP.txt"))) %>% 
  mutate("N" = 3301) %>% # Sample size
  filter(nchar(EA) == 1, nchar(NEA) == 1) # Discard those SNPs that do not have a single nucleotide variation
write_tsv(nn1,here("Metaanalysis","nn.txt"))

