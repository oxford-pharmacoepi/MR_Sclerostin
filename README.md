# MR_Sclerostin
**Background:** Sclerostin inhibitors have been shown to improve bone mineral density (BMD) in large trials and decrease fracture risk. However, regulators like the European Medicines Agency have warned of a possible cardiovascular risk associated with their use, and restricted indication to subjects with no previous cardiovascular history. While observational studies are ongoing, there is a scarcity of information on the causal effects of sclerostin inhibition on cardiovascular risk.

**Objective:** To study the association between sclerostin and cardiovascular health, including biomarkers, risk factors, and cardiovascular events by using Mendelian randomisation (MR) methods.

For more information: [LINK TO THE ARTICLE]

# Running the analysis
Previously, you must ensure that you have a directory "${folder directory}\MR_Sclerostin_Data" with the following files:
  - GWAS_HipFracture\GCST90161240_buildGRCh37.tsv: GWAS of hip fracture by Nethander et al. (2022). Can be downloaded from: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90161001-GCST90162000/GCST90161240/
  - GWAS_IschaemicStroke\Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2: GWAS of ischaemic stroke by Mishra et al (2022). Can be downloaded from: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104540/
  - Nature\Nature_GWAS.tsv: GWAS of sclerostin by Sun et al. (2018). Can be downloaded from: http://www.phpc.cam.ac.uk/ceu/proteins/
  - NatureGenetics\NatureGenetics_GWAS1.txt: GWAS summary statistics of sclerostin by Ferkingstad et al. (2021). Can be dowloaded from: https://www.decode.com/summarydata/ Corresponds to the file "assocvariants.annotated.txt.gz".
  - NatureGenetics\NatureGenetics_GWAS2.txt: GWAS summary statistics of sclerostin by Ferkingstad et al. (2021). Can be dowloaded from: https://www.decode.com/summarydata/ Corresponds to the file "SeqId_GeneName_ProteinName.txt.gz".
  - Science\Science_GWAS.txt: GWAS summary statistics of sclerostin by Pietzner et al. (2021). Can be downloaded from: https://omicscience.org/apps/pgwas
  - SNPs\selected_snps_c17.ped: allele patient-level data extracted from UK Biobank of all the instruments. 
  - UKB\PC.csv: The following fields extracted from UK Biobank: eid (0), sex (31-0.0), PC1 (22009-0.1), PC2 (22009-0.2), PC3 (22009-0.3), PC4 (22009-0.4), PC5 (22009-0.5), PC6 (22009-0.6), PC7 (22009-0.7), PC8 (22009-0.8), PC9 (22009-0.9), PC10 (22009-0.10)
  - UKB\ukb669864_1.csv: The following fields extracted from UK Biobank: eid, 30630-0.0, 30630-1.0,	30640-0.0, 30640-1.0, 30710-0.0, 30710-1.0, 30740-0.0, 30740-1.0, 30750-0.0, 30750-1.0, 30760-0.0, 30760-1.0, 30790-0.0, 30790-1.0, 30870-0.0, 30870-1.0, 30690-0.0, 30690-1.0, 30780-0.0, 30780-1.0, 3148-0.0, 50-0.0, 22006-0.0
  - UKB\ukb669864_birth.csv: The following fields extracted from UK Biobank: eid, 53-0.0, 34-0.0, 40007-0.0, 20002-0.0, 20002-0.1, 20002-0.2, 20002-0.3, 20002-0.4, 20002-0.5, 20002-0.6, 20002-0.7, 20002-0.8, 20002-0.9, 20002-0.10, 20002-0.11, 20002-0.12, 20002-0.13, 20002-0.14, 20002-0.15, 20002-0.16, 20002-0.17, 20002-0.18, 20002-0.19, 20002-0.20, 20002-0.21, 20002-0.22, 20002-0.23, 20002-0.24, 20002-0.25, 20002-0.26, 20002-0.27, 20002-0.28, 20002-0.29, 20002-0.30, 20002-0.31, 20002-0.32, 20002-0.33, 3005-0.0, 3894-0.0, 2443-0.0, 20003-0.0, 20003-0.1, 20003-0.2, 20003-0.3, 20003-0.4, 20003-0.5, 20003-0.6, 20003-0.7, 20003-0.8, 20003-0.9, 20003-0.10, 20003-0.11, 20003-0.12, 20003-0.13, 20003-0.14, 20003-0.15, 20003-0.16, 20003-0.17, 20003-0.18, 20003-0.19, 20003-0.20, 20003-0.21, 20003-0.22, 20003-0.23, 20003-0.24, 20003-0.25, 20003-0.26, 20003-0.27, 20003-0.28, 20003-0.29, 20003-0.30, 20003-0.31, 20003-0.32, 20003-0.33, 20003-0.34, 20003-0.35, 20003-0.36, 20003-0.37, 20003-0.38, 20003-0.39, 20003-0.40, 20003-0.41, 20003-0.42, 20003-0.43, 20003-0.44, 20003-0.45, 20003-0.46, 20003-0.47, 30750-0.0, 20009-0.0, 20009-0.1, 20009-0.2, 20009-0.3, 20009-0.4, 20009-0.5, 20009-0.6, 20009-0.7, 20009-0.8, 20009-0.9, 20009-0.10, 20009-0.11, 20009-0.12, 20009-0.13, 20009-0.14, 20009-0.15, 20009-0.16, 20009-0.17, 20009-0.18, 20009-0.19, 20009-0.20, 20009-0.21, 20009-0.22, 20009-0.23, 20009-0.24, 20009-0.25, 20009-0.26, 20009-0.27, 20009-0.28, 20009-0.29, 20009-0.30, 20009-0.31, 20009-0.32, 20009-0.33, 2976-0.0
  - UKB\ukb669864_logistic.csv: eid, 3005-0.0, 20002-0.0, 20002-0.1, 20002-0.2, 20002-0.3, 20002-0.4, 20002-0.5, 20002-0.6, 20002-0.7, 20002-0.8, 20002-0.9, 20002-0.10, 20002-0.11, 20002-0.12, 20002-0.13, 20002-0.14, 20002-0.15, 20002-0.16, 20002-0.17, 20002-0.18, 20002-0.19, 20002-0.20, 20002-0.21, 20002-0.22, 20002-0.23, 20002-0.24, 20002-0.25, 20002-0.26, 20002-0.27, 20002-0.28, 20002-0.29, 20002-0.30, 20002-0.31, 20002-0.32, 20002-0.33, 3894-0.0, 20003-0.0, 20003-0.1, 20003-0.2, 20003-0.3, 20003-0.4, 20003-0.5, 20003-0.6, 20003-0.7, 20003-0.8, 20003-0.9, 20003-0.10, 20003-0.11, 20003-0.12, 20003-0.13, 20003-0.14, 20003-0.15, 20003-0.16, 20003-0.17, 20003-0.18, 20003-0.19, 20003-0.20, 20003-0.21, 20003-0.22, 20003-0.23, 20003-0.24, 20003-0.25, 20003-0.26, 20003-0.27, 20003-0.28, 20003-0.29, 20003-0.30, 20003-0.31, 20003-0.32, 20003-0.33, 20003-0.34, 20003-0.35, 20003-0.36, 20003-0.37, 20003-0.38, 20003-0.39, 20003-0.40, 20003-0.41, 20003-0.42, 20003-0.43, 20003-0.44, 20003-0.45, 20003-0.46, 20003-0.47, 30750-0.0, 2443-0.0
  - UKB\ukb669864_Table2.csv: eid	31-0.0, 21001-0.0, 21001-1.0, 21001-2.0, 21001-3.0, 21003-0.0, 26410-0.0, 26426-0.0, 26427-0.0, 21000-0.0
  - gp_clinical.txt: GP linked data from UK Biobank.
  - hesin.txt: HES linked data from UK Biobank.
  - hesin_diag.txt: HES linked data from UK Biobank (diagnostics).

Once all these files are located in the directory: 

  1. Download this entirely repository (you can download as a zip folder using Code -> Download ZIP, or you can use GitHub Desktop).
  2. Open the project MR_Sclerostin.Rproj in RStudio (when inside the project, you will see its name on the top-right of your RStudio session).
  3. Open and work though the CodeToRun.R file which should be the only file that you need to interact with. Run the lines in the file. You will notice that you have to specify the following variables:

    -> pathData <- "...": The path to the data directory.
    -> tok <- "...": The token for the LDmatrix. Can be created here: https://analysistools.cancer.gov/LDlink/?tab=apiaccess
    -> outputFolder <- "...": The path to the folder where the results from this analysis will be saved. 
    
After running, you should have a folder named "Results" with results in your output folder.


