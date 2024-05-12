# ============================================================================ #
#                                 loadOutcomeGwas                              #
#                              Marta Alcalde-Herraiz                           #
# ---------------------------------------------------------------------------- #

# Coronary artery disease:
# - Article: https://www.nature.com/articles/s41588-022-01233-6
# - Download GWAS: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/
# - gwasID: ebi-a-GCST90132314

# Myocardial infarction:
# - Article: https://academic.oup.com/eurheartj/article/42/9/919/6126843?login=true
# - Download GWAS: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011365/harmonised/
# - gwasID: ebi-a-GCST011365

# Ischaemic stroke:
# - Article: https://www.nature.com/articles/s41586-022-05165-3
# - Download GWAS: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104540/
# - gwasID: ebi-a-GCST90104540

# Hypertension
# Article: https://data.bris.ac.uk/data/dataset/pnoat8cxo0u52p6ynfaekeigi
# Download GWAS: -ukb-b-14057

# Type 2 diabetes mellitus
# - Article: https://www.nature.com/articles/s41588-022-01058-3
# - Download GWAS:
# - gwasID: ebi-a-GCST90132184

# HDL gwas
# - Article: https://www.nature.com/articles/s41586-021-04064-3
# - Download GWAS: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/

# LDL gwas
# - Article: https://www.nature.com/articles/s41586-021-04064-3
# - Download GWAS: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/

# HbA1C gwas:
# - Article: https://pubmed.ncbi.nlm.nih.gov/34059833
# - Download GWAS: https://magicinvestigators.org/downloads/

# glucose gwas:
# - Article: https://pubmed.ncbi.nlm.nih.gov/34059833
# - Download GWAS: https://magicinvestigators.org/downloads/


loadOutcomeGwas <- function(outcome,exposure){
  outcomeData = switch(
    outcome,

    "bmd" = {
      bmd <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/bmd/30598549-GCST006979-EFO_0009270.h.tsv"))
      bmd <- bmd %>%
        dplyr::select(
          SNP = hm_rsid,
          effect_allele = hm_effect_allele,
          other_allele  = hm_other_allele,
          beta = hm_beta,
          eaf  = hm_effect_allele_frequency,
          se   = standard_error,
          pval = p_value) %>%
        dplyr::mutate(
          samplesize = 426824,
          id.outcome = 'bmd',
          outcome = "bmd"
        ) %>%
        dplyr::filter(SNP %in% exposure$SNP)
      bmd <- TwoSampleMR::format_data(dat = bmd, type = "outcome")
      return(bmd)
    },

    "hf" = {
      hf <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/hipFracture/GCST90161240_buildGRCh37.tsv"))
      hf <- hf %>%
        dplyr::select(
          SNP = variant_id,
          pval = p_value,
          effect_allele, other_allele,
          eaf = effect_allele_frequency,
          beta,
          se = standard_error) %>%
        dplyr::mutate(
          samplesize = 735354,
          id.outcome = "hf",
          outcome = "hf") %>%
        dplyr::filter(SNP %in% exposure$SNP)
      hf <- TwoSampleMR::format_data(dat = hf, type = "outcome")
      return(hf)
    },

    "cad" = {
      # CAD gwas does not contain the SNPs ID, so the following steps have to be done in order to be able to map the SNPs
      # cad <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/cad/GCST90132314_buildGRCh37.tsv"))
      # cad <- cad %>%
      #   dplyr::filter(chromosome == 17) %>%
      #   dplyr::select(
      #     chr = chromosome,
      #     pos = base_pair_location,
      #     eaf = effect_allele_frequency,
      #     beta,
      #     se = standard_error,
      #     samplesize = n,
      #     pval = p_value,
      #     effect_allele,other_allele
      #   ) %>%
      #   dplyr::mutate(
      #     effect_allele = toupper(effect_allele),
      #     other_allele  = toupper(other_allele)
      #   )
      # tab <- cad %>% dplyr::mutate('chr' = paste0("chr",chr,":",pos,"-",pos)) %>% dplyr::select(chr) %>% tibble::as_tibble()
      # # Copy and paste this file:
      # # readr::write_delim(tab,paste0(pathData,"SummaryStatistics/Outcome/cad/fileToLiftOver.txt"),col_names = FALSE)
      # # https://genome.ucsc.edu/cgi-bin/hgLiftOver

      # t1 <- readr::read_delim("C:/Users/martaa/Desktop/Projects/MR_Sclerostin_NEW/SummaryStatistics/Outcome/cad/hglft_genome_2a128_fe3990.err.txt",col_names = FALSE)
      # t2 <- readr::read_delim("C:/Users/martaa/Desktop/Projects/MR_Sclerostin_NEW/SummaryStatistics/Outcome/cad/hglft_genome_2a128_fe3990.bed", col_names = FALSE)
      #
      # cad <- tab %>%
      #   tidyr::separate(col = chr, into = c("X1","X2"), sep = ":") %>%
      #   dplyr::anti_join(t1, by = c("X1","X2")) %>%
      #   tidyr::separate(col = X2, into = c("pos","X2")) %>%
      #   dplyr::mutate("GRCh38" = t2$X2) %>%
      #   tidyr::separate(col = GRCh38, into = c("GRCh38", "X3")) %>%
      #   dplyr::mutate("GRCh38" = as.numeric(GRCh38)) %>%
      #   dplyr::mutate("pos" = as.numeric(pos)) %>%
      #   dplyr::select(-X1,-X2,-X3) %>%
      #   dplyr::left_join(cad, by = 'pos')

      # readr::write_delim(cad, paste0(pathData, "SummaryStatistics/Outcome/cad/GCST90132314_buildGRCh38.txt"))

      cad <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/cad/GCST90132314_buildGRCh38.txt")) %>%
        dplyr::right_join(
          exposure %>%
            dplyr::select(SNP, GRCh38 = 'pos.exposure')
        ) %>%
        dplyr::select(-pos) %>%
        dplyr::rename(pos = GRCh38) %>%
        dplyr::mutate(outcome = "cad_gwas", id.outcome = "cad_gwas")
      cad <- TwoSampleMR::format_data(dat = cad, type = "outcome")

      return(cad)
    },


    "mi" = {
      mi <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/mi/harmonised.qc.tsv")) %>%
        dplyr::select(
          SNP = hm_rsid,
          effect_allele = hm_effect_allele,
          other_allele  = hm_other_allele,
          beta = hm_beta,
          eaf  = hm_effect_allele_frequency,
          se   = standard_error,
          pval = p_value
        ) %>%
        dplyr::mutate(
          id.outcome = 'mi_gwas',
          outcome    = "mi_gwas",
          samplesize = 638717
        ) %>%
        dplyr::filter(SNP %in% exposure$SNP)
      mi <- TwoSampleMR::format_data(dat = mi, type = "outcome")
      return(mi)
    },

    "is" = {
      # is <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/is/Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")) %>%
      #   dplyr::filter(chromosome == 17) %>%
      #   dplyr::select(
      #     chr = chromosome,
      #     pos = base_pair_location,
      #     eaf = effect_allele_frequency,
      #     beta,
      #     se = standard_error,
      #     pval = p_value,
      #     effect_allele,other_allele
      #   ) %>%
      #   dplyr::mutate(
      #     id = 'outcome'
      #   )
      # tab <- is %>% dplyr::mutate('chr' = paste0("chr",chr,":",pos,"-",pos)) %>% dplyr::select(chr) %>% tibble::as_tibble()
      # # # Copy and paste this file:
      # readr::write_delim(tab,paste0(pathData,"SummaryStatistics/Outcome/is/fileToLiftOver.txt"),col_names = FALSE)
      # # https://genome.ucsc.edu/cgi-bin/hgLiftOver
      #
      # t1 <- readr::read_delim("C:/Users/martaa/Desktop/Projects/MR_Sclerostin_NEW/SummaryStatistics/Outcome/is/hglft_genome_2b74_fef2a0.err.txt",col_names = FALSE)
      # t2 <- readr::read_delim("C:/Users/martaa/Desktop/Projects/MR_Sclerostin_NEW/SummaryStatistics/Outcome/is/hglft_genome_2b74_fef2a0.bed", col_names = FALSE)
      #
      # is <- tab %>%
      #   tidyr::separate(col = chr, into = c("X1","X2"), sep = ":") %>%
      #   dplyr::anti_join(t1, by = c("X1","X2")) %>%
      #   tidyr::separate(col = X2, into = c("pos","X2")) %>%
      #   dplyr::mutate("GRCh38" = t2$X2) %>%
      #   tidyr::separate(col = GRCh38, into = c("GRCh38", "X3")) %>%
      #   dplyr::mutate("GRCh38" = as.numeric(GRCh38)) %>%
      #   dplyr::mutate("pos" = as.numeric(pos)) %>%
      #   dplyr::select(-X1,-X2,-X3) %>%
      #   dplyr::left_join(is, by = 'pos')
      #
      # readr::write_delim(is, paste0(pathData, "SummaryStatistics/Outcome/is/Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.txt"))
      is <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/is/Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.txt")) %>%
        dplyr::right_join(
          exposure %>%
            dplyr::select(SNP, GRCh38 = 'pos.exposure')
        ) %>%
        dplyr::select(-pos) %>%
        dplyr::rename(pos = GRCh38) %>%
        dplyr::mutate(outcome = "is.gwas", id.outcome = "is_gwas", samplesize = 1847683)
      is <- TwoSampleMR::format_data(dat = is, type = "outcome")
      return(is)
    },

    "hypertension" = {
      # Since we could not download hypertension gwas, we will extract "all the potential" instruments
      # and save the results in a folder, so we do not have to access the MR server in every iteration
      # exposure <- readr::read_delim(paste0(pathResults,"InstrumentSelection/iv_Fixed_PCA.txt"))
      # hypertension <- TwoSampleMR::extract_outcome_data(snps = exposure$SNP, outcomes = 'ukb-b-14057') %>%
      #   dplyr::select(SNP,pos,beta.outcome, se.outcome, samplesize.outcome, pval.outcome, eaf.outcome, effect_allele.outcome,
      #          other_allele.outcome) %>%
      #   dplyr::mutate(outcome = "outcome", id.outcome = "outcome")
      # readr::write_tsv(hypertension,here::here(paste0(pathData,"SummaryStatistics/Outcome/hypertension/hypertension.txt")))

      hypertension <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/hypertension/hypertension.txt")) %>%
        dplyr::filter(SNP %in% exposure$SNP) %>%
        dplyr::mutate(outcome = "hypertension_gwas", id.outcome = "hypertension_gwas", samplesize.outcome = 462933)

      return(hypertension)
    },

    "t2dm" = {
      t2dm <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/t2dm/DIAMANTE-EUR.sumstat.txt")) %>%
        dplyr::select(SNP = rsID,
                      effect_allele, other_allele,
                      eaf = effect_allele_frequency,
                      beta = `Fixed-effects_beta`,
                      se = `Fixed-effects_SE`,
                      pval = `Fixed-effects_p-value`) %>%
        dplyr::mutate(outcome = "t2dm_gwas", id.outcome = "t2dm_gwas", samplesize = 933970) %>%
        dplyr::filter(SNP %in% exposure$SNP)
      t2dm <- TwoSampleMR::format_data(t2dm, type = "outcome")
      return(t2dm)
    },

    "ldl" = {
      ldl <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/ldl/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results")) %>%
        dplyr::select("SNP" = "rsID",
                      "effect_allele" = "ALT",
                      "other_allele"  = "REF",
                      "eaf"  = "POOLED_ALT_AF",
                      "beta" = "EFFECT_SIZE",
                      "se"   = "SE",
                      "pval" = "pvalue",
                      "samplesize" = "N") %>%
        dplyr::mutate(outcome = "ldl_gwas", id.outcome = "ldl_gwas") %>%
        dplyr::filter(SNP %in% exposure$SNP)

      hdl <- TwoSampleMR::format_data(ldl, type = "outcome")
      return(hdl)
    },

    "hdl" = {
      hdl <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/hdl/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results")) %>%
        dplyr::select("SNP" = "rsID",
                      "effect_allele" = "ALT",
                      "other_allele"  = "REF",
                      "eaf"  = "POOLED_ALT_AF",
                      "beta" = "EFFECT_SIZE",
                      "se"   = "SE",
                      "pval" = "pvalue",
                      "samplesize" = "N") %>%
        dplyr::mutate(outcome = "hdl_gwas", id.outcome = "hdl_gwas") %>%
        dplyr::filter(SNP %in% exposure$SNP)

      hdl <- TwoSampleMR::format_data(hdl, type = "outcome")
      return(hdl)
    },

    "glucose" = {
      glucose <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/glucose/MAGIC1000G_FG_EUR.tsv/MAGIC1000G_FG_EUR.tsv")) %>%
        dplyr::select("SNP"  = "variant",
                      "effect_allele", "other_allele", "beta",
                      "eaf"  = "effect_allele_frequency",
                      "se"   = "standard_error",
                      "pval" = "p_value",
                      "samplesize" = "sample_size") %>%
        dplyr::mutate(outcome = "glucose_gwas", id.outcome = "glucose_gwas") %>%
        dplyr::filter(SNP %in% exposure$SNP)

      glucose <- TwoSampleMR::format_data(glucose, type = "outcome")
      return(glucose)
    },

    "hba1c" = {
      hba1c <- readr::read_delim(paste0(pathData,"SummaryStatistics/Outcome/hba1c/MAGIC1000G_HbA1c_EUR.tsv/MAGIC1000G_HbA1c_EUR.tsv")) %>%
        dplyr::select("SNP"  = "variant",
                      "effect_allele", "other_allele", "beta",
                      "eaf"  = "effect_allele_frequency",
                      "se"   = "standard_error",
                      "pval" = "p_value",
                      "samplesize" = "sample_size") %>%
        dplyr::mutate(outcome = "hba1c_gwas", id.outcome = "hba1c_gwas") %>%
        dplyr::filter(SNP %in% exposure$SNP)

      hba1c <- TwoSampleMR::format_data(hba1c, type = "outcome")
    }
  )

}
