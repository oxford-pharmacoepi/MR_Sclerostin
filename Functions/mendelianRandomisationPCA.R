# ============================================================================ #
#                           mendelianRandomisationPCA                          #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #

mendelianRandomisationPCA <- function(exposure, outcome, gwasID = NULL, binary = TRUE, name = NULL, pca_threshold = 99){
  ld_matrix <- ieugwasr::ld_matrix(
    exposure$SNP,
    bfile    = paste0(pathData,'LD_ReferencePannel/EUR'),
    plink_bin = paste0(pathData,'Plink/plink.exe'))

  pattern <- ld_matrix %>%
    tibble::as_tibble() %>%
    colnames() %>%
    tibble::tibble() %>%
    tidyr::separate(col = ".", into = c("SNP", "effect_allele","other_allele"), sep = "_")

  dat_harmonised <- pattern %>%
    dplyr::select(SNP) %>%
    dplyr::left_join(
      harmonisingData(exposure = exposure,
                      outcome  = outcome,
                      gwasID   = gwasID),
      by = "SNP"
    )
  dat <- TwoSampleMR::harmonise_ld_dat(dat_harmonised,ld_matrix)

  ## Principal component analysis
  betaX <- dat$x$beta.exposure
  seX   <- dat$x$se.exposure
  betaY <- dat$x$beta.outcome
  seY   <- dat$x$se.outcome
  rho   <- dat$ld
  omega <- seY%o%seY*rho
  phi <- ((betaX/seY)%o%(betaX/seY))*rho

  # Prcomp outcomes:
  # - sdev: standard deviations of the principal components (i.e., the square roots of the eigenvalues)
  # - Rotation: matrix of variable loadings (i.e., matrix whose columns contain the eigenvectors)
  PC <- prcomp(phi, scale. = FALSE)
  K  <- which((cumsum(PC$sdev^2/sum(PC$sdev^2)*100))>pca_threshold)[1] # Number of PC

  # Transformed values:
  betaX_transformed <- as.numeric(betaX%*%PC$rotation[,1:K])
  betaY_transformed <- as.numeric(betaY%*%PC$rotation[,1:K])
  omega_transformed <- t(PC$rotation[,1:K])%*%omega%*%(PC$rotation[,1:K])

  # Inverse variance weighted method using principal component analysis
  IVW <- as.numeric((betaX_transformed%*%solve(omega_transformed)%*%betaY_transformed)*solve((betaX_transformed%*%solve(omega_transformed)%*%betaX_transformed)))
  SE  <-  as.numeric(sqrt(solve(t(betaX_transformed%*%solve(omega_transformed)%*%betaX_transformed))))
  pval <-  as.numeric(exp(-0.717*abs(IVW/SE)-0.416*(IVW/SE)^2)) # calculated based on: https://www.bmj.com/content/343/bmj.d2304


  tableResults <- list()
  tableResults$Phi <- phi

  if(binary){
    tableResults$Results <- tibble::tibble(
      "Variance explained" = cumsum(PC$sdev^2/sum(PC$sdev^2)*100)[K],
      "Number of principal components" = K,
      "Method"      = "IVW (fixed)",
      "Outcome"     = .env$name,
      "OR" = exp(-IVW),
      "U95CI" = exp(-IVW+1.96*SE),
      "L95CI" = exp(-IVW-1.96*SE),
      "SE"    =  SE,
      "Pval"  = pval)
  }else{
    tableResults$Results <- tibble::tibble(
      "Variance explained" = cumsum(PC$sdev^2/sum(PC$sdev^2)*100)[K],
      "Number of principal components" = K,
      "Method"      = "IVW (fixed)",
      "Outcome"     = .env$name,
      "Beta"  = -IVW,
      "U95CI" = -IVW+1.96*SE,
      "L95CI" = -IVW-1.96*SE,
      "SE"    =  SE,
      "Pval"  = pval)
  }
  return(tableResults)
}
