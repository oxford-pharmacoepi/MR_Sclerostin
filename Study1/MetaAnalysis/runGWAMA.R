# Instructions to use GWAMA:
#  1. First install gwama in the gwama_directory
#  2. Create a "gwama.in" file with the following lines:
#       Ferkingstad_NatureGenetics_GWAS.txt
#       Pietzner_Science_GWAS.txt
#       Sun_Nature_GWAS.txt
#     It will be read by GWAMA software to know which files it has to read.
#     Save the file where gwama.exe is
#  3. Run the following file:
shell <- "cmd"
# Fixed effects: gwama -o ResultsMetaAnalysis/GWAMA_Fixed -qt
system(paste0(
  shQuote(paste0(gwama_directory,"gwama.exe"), type = shell),
  "  -i ", shQuote(paste0(gwama_directory,"gwama.in"), type = shell),
  " -qt",
  " -o ",
  shQuote(here::here("Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Fixed"), type = shell)))

# Random effects: gwama -o ResultsMetaAnalysis/GWAMA_Random -qt -r
system(paste0(
  shQuote(paste0(gwama_directory,"gwama.exe"), type = shell),
  "  -i ", shQuote(paste0(gwama_directory,"gwama.in"), type = shell),
  " -r -qt",
  " -o ",
  shQuote(here::here("Study1/MetaAnalysis/ResultsMetaAnalysis/GWAMA_Random"), type = shell)
))

