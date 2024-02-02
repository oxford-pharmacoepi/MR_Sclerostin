# ============================================================================ #
#                                ObtainPositions                               #
#                             Marta Alcalde-Herraiz                            #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This sript obtain the positions (in build GRCh38/hg38) of the list of SNPS.  #
# It uses biomaRt package.                                                     #
# ============================================================================ #

# Choose a meta-analysed dataset (whether the random effects or the fixed one)
t <- readr::read_delim(paste0(pathResults,'MetaAnalysis/Fixed.txt'))

# Download BiocManager package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Chunks of 10,000 SNPs to calculate their position
pf <- 10000
k  <- (t %>% nrow())/pf # Number of iterations
snps <- t$SNP # List of SNPs

# This loop might be slow to run it
for(i in 1:ceiling(k)){
  # List of SNPs within the chunk
  snps1 <- snps[(1:pf)+(i-1)*pf]

  # Extract posiiton of the SNPs
  snp_locations <- tibble::as_tibble(
    biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                   filters    = "snp_filter",
                   values     = snps1,
                   mart       = biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp"))
  ) %>%
    dplyr::mutate(chr_name = as.character(chr_name)) %>%
    dplyr::filter(chr_name == "17") %>%
    dplyr::rename(SNP = refsnp_id,
                  chr = chr_name,
                  pos = chrom_start)
  readr::write_tsv(snp_locations,paste0(pathResults,"InstrumentSelection/MappingGRCh38_",i*pf,".txt"))
}

listOfFiles <- list.files(path = paste0(pathResults,"InstrumentSelection/"),
                          recursive  = TRUE,
                          pattern    = "^MappingGRCh38_",
                          full.names = TRUE)
mapping <- readr::read_delim(listOfFiles)

# See which SNPs are not mapped and try it again (it might be possible that in
# that iteration the program has failed)
snpsMissing <- mapping %>%
  dplyr::select("SNP") %>%
  dplyr::anti_join(t %>% dplyr::select(SNP))

snp_locations <- tibble::as_tibble(
  biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                 filters    = "snp_filter",
                 values     = snpsMissing,
                 mart       = biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp"))) %>%
  dplyr::mutate(chr_name = as.character(chr_name)) %>%
  dplyr::filter(chr_name == "17") %>%
  dplyr::rename(SNP = refsnp_id,
                chr = chr_name,
                pos = chrom_start)
readr::write_tsv(table_locations,paste0(pathResults,"InstrumentSelection/MappingGRCh38_Missing.txt"))

# Read all the files again
listOfFiles <- list.files(path = paste0(pathResults,"InstrumentSelection/"),
                          recursive  = TRUE,
                          pattern    = "^MappingGRCh38",
                          full.names = TRUE)
mapping <- readr::read_delim(listOfFiles) %>% dplyr::distinct()
readr::write_tsv(mapping,paste0(pathResults,"InstrumentSelection/MappingGRCh38.txt"))

dir.create(paste0(pathResults,"InstrumentSelection/Mapping"))
file.copy(from = listOfFiles,
          to   = paste0(pathResults,'InstrumentSelection/Mapping/'))
file.remove(listOfFiles)

# Build a fixed effects - random effects complete tables
fixed <- readr::read_delim(here::here('Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Fixed.out')) %>%
  dplyr::mutate(
    se = (beta_95U-beta)/1.96
  ) %>%
  dplyr::select(
    "SNP" = "rs_number",
    "effect_allele" = "reference_allele",
    "other_allele", "beta",
    "pval" = `p-value`,
    "i2",
    "se",
    "sample_size" = "n_samples"
  ) %>%
  dplyr::left_join(
    mapping,
    by = "SNP"
  ) %>%
  dplyr::left_join(
    readr::read_delim(here::here("Study1/MetaAnalysis/ResultsMetaAnalysis/METAL_Fixed.TBL")) %>%
      dplyr::select("SNP" = "MarkerName",
                    "Allele1",
                    "eaf" = "Freq1")
  ) %>%
  dplyr::mutate(
    eaf = dplyr::if_else(toupper(Allele1) == effect_allele,eaf,1-eaf)
  ) %>%
  dplyr::select(
    -Allele1
  )
readr::write_tsv(fixed,paste0(pathResults,"InstrumentSelection/FixedMapped.txt"))

random <- readr::read_delim(here::here('Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Random.out')) %>%
  dplyr::mutate(
    se = (beta_95U-beta)/1.96
  ) %>%
  dplyr::select(
    "SNP" = "rs_number",
    "effect_allele" = "reference_allele",
    "other_allele", "beta",
    "pval" = `p-value`,
    "se",
    "i2",
    "sample_size" = "n_samples"
  ) %>%
  dplyr::left_join(
    mapping,
    by = "SNP"
  ) %>%
  dplyr::left_join(
    readr::read_delim(here::here("Study1/MetaAnalysis/ResultsMetaAnalysis/METAL_Fixed.TBL")) %>%
      dplyr::select("SNP" = "MarkerName",
                    "Allele1",
                    "eaf" = "Freq1")
  ) %>%
  dplyr::mutate(
    eaf = dplyr::if_else(toupper(Allele1) == effect_allele,eaf,1-eaf)
  ) %>%
  dplyr::select(
    -Allele1
  )
readr::write_tsv(random,paste0(pathResults,"InstrumentSelection/RandomMapped.txt"))
