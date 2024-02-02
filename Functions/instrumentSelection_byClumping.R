# ============================================================================ #
#                         instrumentSelection_byClumping                       #
#                             Marta Alcalde-Herraiz                            #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# Selects the SNPs to be used as instruments for the mendelian randomisation   #
# ============================================================================ #
instrumentSelection_byClumping <- function(dat_input = NULL,
                                           window    = 500e3,
                                           min_eaf   = 0.01,
                                           gene_start = 43753738,
                                           gene_end   = 43758791,
                                           chr        = 17,
                                           pval_threshold = 1e-6,
                                           # Clump parameters
                                           clump_kb = NULL,
                                           clump_r2 = NULL,
                                           clump_p  = 0.99,
                                           pop = 'EUR',
                                           bfile    = paste0(pathData,'LD_ReferencePannel/EUR'),
                                           plink_bin = paste0(pathData,'Plink/plink.exe')){
  # Mapping --------------------------------------------------------------------
  dat <- readr::read_delim(paste0(pathResults,"/InstrumentSelection/",dat_input,"Mapped.txt")) %>%
    dplyr::distinct()

  # Attrition of the SNPs ------------------------------------------------------
  a <- matrix(0,15,2)
  a[1,1] <- dat %>% nrow()
  a[1,2] <- 'All'

  # SOST region
  dat1 <- dat %>% dplyr::filter(pos >= gene_start - window - 1,
                                pos <= gene_end   + window + 1)
  a[2,1] <- dat1 %>% nrow()
  a[2,2] <- paste0('SNPs within +-',window,'kb from the end/start of the SOST gene')

  # Minimum eaf
  # dat1 <- dat1 %>% dplyr::filter(eaf >= min_eaf)
  # a[3,1] <- dat1 %>% nrow()
  # a[3,2] <- paste0('SNPs with eaf higher/equal than ',min_eaf)

  # P-value lower than the threshold
  dat1 <- dat1 %>% dplyr::filter(pval <= pval_threshold) %>%
    dplyr::mutate(id = "exposure")
  a[4,1] <- dat1 %>% nrow()
  a[4,2] <- paste0('SNPs with p-value lower/equal than ', pval_threshold)

  dat1 <- TwoSampleMR::format_data(dat = dat1, type = "exposure")

  source(here::here('Functions','ld_clumping.R'))
  iv <- ld_clumping(dat = dat1,
                    clump_kb = clump_kb,
                    clump_r2 = clump_r2,
                    clump_p  = clump_p,
                    pop      = pop,
                    bfile    = bfile,
                    plink_bin = plink_bin)



  readr::write_delim(iv,paste0(pathResults,
                               "InstrumentSelection/iv_",
                               dat_input,"_",
                               "gene",format(window, scientific = FALSE),"kb_",
                               "r2",clump_r2,"_",
                               "clumpw",format(clump_kb, scientific = FALSE),
                               ".txt"))


  readr::write_delim(as.data.frame(a),paste0(pathResults,
                               "InstrumentSelection/", dat_input, "Attrition.txt"))

  return(iv)
}




