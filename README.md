# MR_Sclerostin
**Background:** Sclerostin inhibitors have been shown to improve bone mineral density (BMD) in large trials and decrease fracture risk. However, regulators like the European Medicines Agency have warned of a possible cardiovascular risk associated with their use, and restricted indication to subjects with no previous cardiovascular history. While observational studies are ongoing, there is a scarcity of information on the causal effects of sclerostin inhibition on cardiovascular risk.

**Objective:** To study the association between sclerostin and cardiovascular health, including biomarkers, risk factors, and cardiovascular events by using Mendelian randomisation (MR) methods.

For more information: [LINK TO THE ARTICLE]

# Running th analysis
Previously, you must ensure that you have a directory "${folder directory}\MR_Sclerostin_Data" with the following files:
  - GWAS_HipFracture\GCST90161240_buildGRCh37.tsv: GWAS of hip fracture by Nethander et al. (2022). Can be downloaded from: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90161001-GCST90162000/GCST90161240/
  - GWAS_IschaemicStroke\Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2: GWAS of ischaemic stroke by Mishra et al (2022). Can be downloaded from: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104540/
  - Nature\Nature_GWAS.tsv: GWAS of sclerostin by Sun et al. (2018). Can be downloaded from: http://www.phpc.cam.ac.uk/ceu/proteins/
  - NatureGenetics\NatureGenetics_GWAS1.txt: GWAS summary statistics of sclerostin by Ferkingstad et al. (2021). Can be dowloaded from: https://www.decode.com/summarydata/ Corresponds to the file "assocvariants.annotated.txt.gz".
  - NatureGenetics\NatureGenetics_GWAS2.txt: GWAS summary statistics of sclerostin by Ferkingstad et al. (2021). Can be dowloaded from: https://www.decode.com/summarydata/ Corresponds to the file "SeqId_GeneName_ProteinName.txt.gz".
  - Science\Science_GWAS.txt: GWAS summary statistics of sclerostin by Pietzner et al. (2021). Can be downloaded from: https://omicscience.org/apps/pgwas
  - 
  1. Download this entirely repository (you can download as a zip folder using Code -> Download ZIP, or you can use GitHub Desktop).
  2. Open the project MR_Sclerostin.Rproj in RStudio (when inside the project, you will see its name on the top-right of your RStudio session).
  3. Open and work though the CodeToRun.R file which should be the only file that you need to interact with. Run the lines in the file. You will notice that 

