# ============================================================================ #
#                                  ld_clumping                                 #
#                             Marta Alcalde-Herraiz                            #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# We modified the function "ld_clump" from the package "ieugwasr".             #
# This function uses PLINK clumping method, where SNPs within a specified      #
# window and in LD will be pruned, and the SNP with lowest p-value will be     #
# retained.                                                                    #
# To be able to use this function, you first need to download plink:           #
# Download plink: https://www.cog-genomics.org/plink/                          #
# And also the LD reference pannel:                                            #
# http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz                                  #
# More information about the original package and the original function can be #
# found in these links:                                                        #
# Original package: https://github.com/MRCIEU/ieugwasr                         #
# Original function: https://mrcieu.github.io/ieugwasr/articles/local_ld.html#ld-clumping-1
# ============================================================================ #

ld_clumping <- function(dat = NULL,
                        clump_kb = NULL,
                        clump_r2 = NULL,
                        clump_p = 0.99,
                        pop = "EUR",
                        bfile = NULL,
                        plink_bin = NULL)
{
  stopifnot("SNP" %in% colnames(dat))
  stopifnot("pval.exposure" %in% colnames(dat))
  dat$id <- 'exposure'

  # Make textfile
  shell <- "cmd"
  fn    <- tempfile()
  write.table(data.frame(SNP=dat[["SNP"]], P=dat[["pval.exposure"]]),
              file = fn,
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --clump ", shQuote(fn,    type = shell),
    " --clump-p1 ", clump_p,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --out ", shQuote(fn, type = shell)
  )

  system(fun2) # Executes plink via R

  res <- read.table(paste0(fn,'.clumped'), header = TRUE) # Read plink file
  unlink(fn)

  y <- dat %>% dplyr::inner_join(res %>% dplyr::select(SNP))

  if(nrow(y) > 0)
  {
    message("Removing ", nrow(dat) - nrow(y), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
  }else{
    message("No variants remaining after LD-Clumping. Something went wrong")

  }
  return(y)
}




