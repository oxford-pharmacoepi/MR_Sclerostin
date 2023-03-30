# ============================================================================ #
#                                PROCESSING                                    #
#                        2023 - Marta Alcalde-Herraiz                          #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This file reads the files created after using GWAS and METAL and converts    #
# them into .csv files.                                                        #
# ============================================================================ #
rm(list = setdiff(ls(),c("pathData","tok")))

# Read GWAMA file --------------------------------------------------------------
# Read the output files
t1 <- as_tibble(read.delim(here("Metaanalysis","Fixed_results","Fixed.out"), sep = "\t"))
t2 <- as_tibble(read.delim(here("Metaanalysis","Random_results","Random.out"), sep = "\t"))
# Add SE column
t1 <- t1 %>%
  mutate(se = (beta-beta_95L)/1.96)
t2 <- t2 %>%
  mutate(se = (beta-beta_95L)/1.96)

# Save the results
write.csv(t1,here("Metaanalysis","Fixed_results","Fixed.csv"))
write.csv(t2,here("Metaanalysis","Random_results","Random.csv"))

# Read METAL file --------------------------------------------------------------
t3 <- as_tibble(read.delim(here("Metaanalysis","METAANALYSIS1.TBL"), sep = "\t"))

write.csv(t3,here("Metaanalysis","Fixed_results","Fixed_metal.csv"))

# Table to compare the results of both analysis
tab <- t1 %>%
  select("SNP" = "rs_number",
         "EA.fixed" = "reference_allele",
         "OA.fixed" = "other_allele",
         "Effect.fixed" = "beta",
         "pval.fixed" = "p.value",
         "se.fixed" = "se") %>%
  inner_join(t2 %>% 
               select("SNP" = "rs_number",
                      "EA.random" = "reference_allele",
                      "OA.random" = "other_allele",
                      "Effect.random" = "beta",
                      "pval.random" = "p.value",
                      "i2","q_statistic","q_p.value","n_samples",
                      "se.random" = "se"),
             by = "SNP") %>% 
  inner_join(t3 %>%
               select("SNP" = "MarkerName",
                      "EA.fixed.metal" = "Allele1",
                      "OA.fixed.metal" = "Allele2",
                      "Effect.fixed.metal" = "Effect",
                      "pval.fixed.metal" = "P.value",
                      "se.fixed.metal" = "StdErr"),
             by = "SNP")

# Save the results
write.csv(tab,here("Metaanalysis","Fixed_vs_random.csv"))

