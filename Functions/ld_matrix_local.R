# ============================================================================ #
#                               ld_matrix_local                                #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #


# Get LD matrix using local plink binary and reference dataset
ld_matrix_local <- function(variants, bfile, plink_bin, with_alleles=TRUE)
{
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)


  fun1 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --extract ", shQuote(fn, type=shell),
    " --make-just-bim ",
    " --keep-allele-order ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun1)

  bim <- read.table(paste0(fn, ".bim"), stringsAsFactors=FALSE)

  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --extract ", shQuote(fn, type=shell),
    " --r square ",
    " --keep-allele-order ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2)
  res <- read.table(paste0(fn, ".ld"), header=FALSE) %>% as.matrix
  if(with_alleles)
  {
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep="_")
  } else {
    rownames(res) <- colnames(res) <- bim$V2
  }
  return(res)
}
