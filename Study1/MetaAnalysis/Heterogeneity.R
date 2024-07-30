# ============================================================================ #
#                                  Heterogeneity                               #
#                              Marta Alcalde-Herraiz                           #
# ============================================================================ #


## QQ plot
for(i in c("Fixed")){
  dat <- read.delim(paste0(pathResults,'MetaAnalysis/',i,'.txt')) %>%
    dplyr::mutate(LOG10P = -log10(pval))

  source(here::here("Functions/getQQPlot.R"))
  p <- getQQPlot(dat,5,10,color = "#2166AC",dot_size = .5,reduce_dataset = 0.99)
  pdf("plot.pdf")
  print(p)
  dev.off()

  # ggplot2::ggsave(paste0(pathResults,"MetaAnalysis/QQPlot_",i,".pdf"), plot = p,
  #                 width = 80, height = 60, dpi = 500, units = 'mm')
}

