# Instructions to use METAL
#  6. Click on metal.exe to open the command window.
#  7. Write the following instructions:
#     SCHEME STDERR
#     MARKER MARKERNAME
#     ALLELELABELS EA NEA
#     PVALUELABEL pval
#     FREQLABEL EAF
#     AVERAGEFREQ ON
#     WEIGHTLABEL N
#     STDERRLABEL SE
#     EFFECTLABEL BETA
#     PROCESS Ferkingstad_NatureGenetics_GWAS.txt
#     PROCESS Pietzner_Science_GWAS.txt
#     PROCESS Sun_Nature_GWAS.txt
#     ANALYZE
#    The file: METAANALYSIS1.TBL.info will contain the results

shell <- "cmd"
# Fixed effects: gwama -o ResultsMetaAnalysis/GWAMA_Fixed -qt
system(paste0(
  shQuote(paste0(metal_directory,"metal.exe"), type = shell),
  " D:/Projects/MR_Sclerostin/OtherSoftware/metal.txt"))
