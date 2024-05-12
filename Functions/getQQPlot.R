# ============================================================================ #
#                                 getQQPlot.R                                  #
#                            Marta Alcalde-Herraiz                             #
# ============================================================================ #


getQQPlot <- function(gwas,x_lim,y_lim,
                      color = "#2166AC",
                      dot_size = 0.5,
                      reduce_dataset = 0.01
){

  # QQ plot ----------------------------------------------------------------------
  n <- nrow(gwas)
  ci <- .95

  dat <- data.frame(
    observed = sort(gwas$LOG10P),
    expected = sort(-log10(ppoints(n))),
    clower   = sort(-log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))),
    cupper   = sort(-log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))))
  )

  # Reduce dataset size
  data_top <- dat %>% dplyr::filter(observed > 1)
  data_low <- dat %>% dplyr::filter(observed <= 1) %>% dplyr::sample_frac(reduce_dataset)
  dat <- data_top %>% dplyr::full_join(data_low)

  # Calculate genetic inflation factor
  gif <- as.character(round(median(gwas$zscore^2)/qchisq(0.5,1), digits = 2))

  # Customize qqplot
  plot2 <- ggplot2::ggplot(dat, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::scale_x_continuous(limits = c(0,x_lim), expand = c(0,0), breaks = seq(0,7,1)) +
    ggplot2::scale_y_continuous(limits = c(0,y_lim), expand = c(0,0), breaks = seq(0,y_lim,1)) +
    ggplot2::geom_segment(data = . %>% dplyr::filter(expected == max(expected)),
                          ggplot2::aes(x = 0, xend = x_lim, y = 0, yend = x_lim),
                          size = 0.25, color = "grey30", lineend = "round",alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    ggplot2::geom_point(color = color,  size = dot_size) +
    ggplot2::labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
                  y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
    ggplot2::annotate("text",
             x = .5,
             y = 8,
             hjust = "right",
             label = expression(paste(lambda," = ")),
             colour="black", size = 3, family = "Calibri") +
    ggplot2::annotate("text",
                      x = .5,
                      y = 8,
                      hjust = "left",
                      label = gif,
                      colour="black", size = 3, family = "Calibri") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position="none",
      panel.border = ggplot2::element_rect(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 7),
      axis.text  = ggplot2::element_text(size = 7),
      text = ggplot2::element_text(family = "Calibri")
    )

  return(plot2)
}
