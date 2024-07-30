# Sclerostin levels' effect on cardiovascular biomarkers, risk factors, and cardiovascular disease: a Mendelian randomisation study
**Background:** Sclerostin inhibitors have been shown to improve bone mineral density (BMD) in large trials and decrease fracture risk. However, regulators like the European Medicines Agency have warned of a possible cardiovascular risk associated with their use, and restricted indication to subjects with no previous cardiovascular history. While observational studies are ongoing, there is a scarcity of information on the causal effects of sclerostin inhibition on cardiovascular risk.

**Objective:** To study the association between sclerostin and cardiovascular health, including biomarkers, risk factors, and cardiovascular events by using Mendelian randomisation (MR) methods.

For more information see our pre-print: https://www.researchsquare.com/article/rs-3209943/v1

# Running the analysis
Previously, you must ensure that you have a directory "{folder directory}\MR_Sclerostin" with the following directories/files:
  - LD_ReferencePannel: should contain the 1000 Genomes European reference panel (EUR.bed, EUR.bim, EUR.fam). These files will be used to calculate the LD matrix and perform pruning locally. Please refer to the following website to download the files: https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  - OtherSoftware: Install Metal and GWAMA in this directory. You can find the instructions to download the softwares in these links: https://csg.sph.umich.edu/abecasis/metal/ (METAL) and https://genomics.ut.ee/en/tools (GWAMA).
  - Plink: Install PLINK in this directory (https://www.cog-genomics.org/plink/).
  - SummaryStatistics/Exposure/Ferkingstad_NatureGenetics/assocvariants.annotated.txt/assocvariants.annotated.txt: GWAS summary statistics published by Ferkinstad et al.
  - SummaryStatistics/Exposure/Ferkingstad_NatureGenetics/13101_60_SOST_SOST.txt/13101_60_SOST_SOST.txt: GWAS summary statistics published by Ferkinstad et al.
  - SummaryStatistics/Exposure/Pietzner_Science/Science_GWAS.txt: GWAS summary statistics published by Pietzner et al.
  - SummaryStatistics/Exposure/Sun_Nature/GCST90242722.h.tsv/SOST.13101.60.3.tsv: GWAS summary statistics published by Sun et al.
  - SummaryStatistics/Outcome/bmd/30598549-GCST006979-EFO_0009270.h.tsv: GWAS summary statistics of bone mineral density.
  - SummaryStatistics/Outcome/hipFracture/GCST90161240_buildGRCh37.tsv: GWAS summary statistics of hip fracture.
  - SummaryStatistics/Outcome/cad/GCST90132314_buildGRCh37.tsv: GWAS summary statistics of coronary artery disease.
  - SummaryStatistics/Outcome/mi/harmonised.qc.tsv: GWAS summary statistics of myocardial infarction.
  - SummaryStatistics/Outcome/is/Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.txt: GWAS summary statistics of ischaemic stroke.
  - SummaryStatistics/Outcome/t2dm/DIAMANTE-EUR.sumstat.txt: GWAS summary statistics of type 2 diabetes.
  - SummaryStatistics/Outcome/ldl/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results: GWAS summary statistics of LDL cholesterol.
  - SummaryStatistics/Outcome/hdl/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.result: GWAS summary statistics of HDL cholesterol.
  - SummaryStatistics/Outcome/glucose/MAGIC1000G_FG_EUR.tsv/MAGIC1000G_FG_EUR.tsv: GWAS summary statistics of fasting glucose.
  - SummaryStatistics/Outcome/hba1c/MAGIC1000G_HbA1c_EUR.tsv/MAGIC1000G_HbA1c_EUR.tsv: GWAS summary statistics of hba1c.
  - UKBiobank/gp_clinical.txt: gp clinical linked data from UK Biobank.
  - UKBiobank/hesin.txt: hes linked data from UK Biobank.
  - UKBiobank/hesin_diag.txt: hes diagnoses linked dataset from UK Biobank.
  - UKBiobank/olink_data.txt: proteomics linked dataset from UK Biobank.
  - UKBiobank/selected_snps_c17: .bed, .bim, .fam, .ped, .bed with the alleles information from UK Biobank (only the instruments).
  - UKBiobank/ukb{application_number}.r and UKBiobank/ukb{application_number}.tab: UK Biobank dataset contaning the fields indicated in `fields_ukbiobank.txt`.

Once all these files are located in the directory: 

  1. Download this entirely repository (you can download as a zip folder using Code -> Download ZIP, or you can use GitHub Desktop).
  2. Open the project MR_Sclerostin.Rproj in RStudio (when inside the project, you will see its name on the top-right of your RStudio session).
  3. Open and work though the CodeToRun.R file which should be the only file that you need to interact with. Run the lines in the file. You will notice that you have to specify the following variables:

    -> pathData <- "...": The path to the data directory.
    -> outputFolder <- "...": The path to the folder where the results from this analysis will be saved. 
    
After running, you should have a folder named "Results" with results in your output folder.


