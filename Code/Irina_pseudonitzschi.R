## Script for analyzing Psuedonitzchia dataset
## Written by Zach Quinlan for Irina Koester 
## Created 11/08/2019
## Modified into .R rather than .RMD 11/12/2019
## Merged with changes from Irina on 11/12/2019

# LOADING -- Libraries ----------------------------------------------------
## Load in all of the data frames and libraries
#Data manipulations
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(randomForest)

# PCoA and visualizations
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)


# LOADING -- Functions ----------------------------------------------------
flag_background <- function(data, 
                            min_val = 0.5, 
                            blank_columns = match(names(select(data, contains("blank", ignore.case = TRUE))), names(data)), 
                            sample_columns = match(names(select(data,-c(contains("blank", ignore.case = TRUE),"feature_number"))), 
                                                   names(data))) 
{
  require("tidyverse")
  data$max_blanks <- apply(data[blank_columns], 1, max)
  data$mean_samples <- apply(data[sample_columns], 1, mean, na.rm = TRUE)
  
  no_background <- data%>%
    mutate(background_features = case_when(mean_samples*min_val > max_blanks ~ "real",
                                           TRUE ~ "background"))%>%
    dplyr::select(-c(max_blanks, mean_samples))
}

flag_transient <- function(data, 
                           sample_columns = match(names(select(data, -contains("blank", ignore.case = TRUE))), names(data)), 
                           noise_level = 2E5, 
                           replication_number = 3) 
{
  require("tidyverse")
  no_transient <- data%>%
    add_column(samples_over_noise = rowSums(.[sample_columns] > noise_level), .before = 2)%>%
    mutate(transient_features = case_when(samples_over_noise >= replication_number ~ "real",
                                          TRUE ~ "transient"))%>%
    dplyr::select(-samples_over_noise)
}


map <- purrr::map
select <- dplyr::select

tidy <- broom::tidy
rename <- dplyr::rename

# LOADING -- Dataframes ---------------------------------------------------
quant_raw <- read_csv("./Raw/quant_all.csv")%>%
  select(-c(2:3))%>%
  rename(feature_number = 1)

metadata_quant <- read_tsv("./Raw/metadata_table.tsv")%>%
  mutate(filename = gsub("mzXML", "mzML", filename))

cat_df <- read_csv("./Raw/Pn_Ex2_MASTERx_Canopus_categories_probability.csv")

chl <- read_xlsx("Raw/Pn_Ex2_Chlorophyll.xlsx")

feature_info <- read_csv("./Raw/Pn_Ex2_MASTERx_elements.csv")%>%
  rename(feature_number = 1)
  select(feature_number, everything())

otu_df <- read_tsv("./Raw/Irina_2018_16s_exp_GEL.swarm.tax")

otu_samples <- read_csv("./Raw/Pn_16S_identifiers.csv")%>%
  rename("sample_name" = "SampleID")%>%
  rename("sample_code" = "OTU_name")


# CLEANING -- Removing Blanks ---------------------------------------------
culture_blanks <- (metadata_quant%>%
                     filter(SampleType == "blank_culturemedia"))$filename%>%
  as.vector()

culture_samples <- (metadata_quant%>%
                     filter(SampleType == "culture_multiplespecies"))$filename%>%
  as.vector()

quant_blanks_env <- quant_raw%>%
  flag_background(blank_columns =  match(names(select(., Blank_Fieldtrip.mzML)), names(.)))%>%
  filter(background_features == "real")%>%
  select(-background_features)

quant_culture_blanks_removed <- quant_blanks_env%>%
  select(c(feature_number, culture_blanks, culture_samples))%>%
  flag_background(blank_columns = match(names(select(., culture_blanks)), names(.)))%>%
  filter(background_features == "real")%>%
  select(-background_features)
  
quant_df <- quant_blanks_env%>%
  select(-c(culture_blanks, culture_samples))%>%
  full_join(quant_culture_blanks_removed, by = "feature_number")%>%
  flag_transient()%>%
  filter(transient_features == "real")%>%
  select(-transient_features)

# CLEANING -- Stats dataframes --------------------------------------------------------------------------------
## Cleaning all of the data
quant_stats <- quant_df%>%
  gather(sample_code, xic, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>%
  unite(sample_code, c("Experiment", "Organism", "biological_replicates"), sep = "_", remove = FALSE)%>%
  left_join(chl, by = "sample_code")%>%
  select(-sample_code)%>%
  unite(sample_code, c("Experiment", "Organism", 
                      "biological_replicates", "DOM_fil", 
                      "technical_replicates"), sep = "_")%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic),
         chl_norm = ra/chl,
         asin = asin(sqrt(chl_norm)))%>%
  ungroup()%>%
  rename("feature_number" = "SampleID")%>%
  select(-c(ra, xic,  chl_norm, chl))%>%
  # group_by(feature_number)%>%
  # filter(sum(.$asin) != 0)%>%
  separate(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")
  # ungroup()

cat_stats <- cat_df%>%
  select(1, 25:ncol(.))%>%
  gather(cat, prob, 2:ncol(.))%>%
  group_by(cat)%>%
  filter(max(prob) > 0.5)%>%
  ungroup()%>%
  mutate(asin = asin(sqrt(prob)))%>%
  select(-prob)%>%
  spread(cat, asin)%>%
  left_join(cat_df[1:24], by = "FeatureID")

# CLEANING -- OTU table -------------------------------------------------------------------
otu_clean <- otu_df%>%
  select(-c(1:2))%>%
  rownames_to_column("otu_number")%>%
  gather(sample_code, reads, 3:ncol(.))%>%
  left_join(otu_samples, ., by = "sample_code")%>%
  separate(sample_name, c("Experiment", "Organism", "biological_replicate"), sep = "_")%>%
  unite(sample_name, c("Organism", "biological_replicate"), sep = "_")%>%
  separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  mutate(Class = case_when(Class %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                           TRUE ~ as.character(Class)),
         Order = case_when(Class == "unclassified" ~ "",
                           TRUE ~ as.character(Order)),
         Family = case_when(Class == "unclassified" ~ "",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Class == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Class == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                           TRUE ~ as.character(Order)),
         Family = case_when(Order == "unclassified" ~ "",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Order == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Order == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Family == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Family == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Genus == "unclassified" ~ "",
                         OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                         TRUE ~ as.character(OTU)))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "otu_number"), sep = ";")%>%
  select(-c(OTU, Experiment, sample_code))

# PRE-STATS -- OTU TABLE --------------------------------------------------
## Making the stats dataframes for OTU, family and classes
otu_stats <- otu_clean%>%
  group_by(sample_name)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()

family_stats <- otu_clean%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "otu_number"), sep = ";")%>%
  select(-c("Genus", "otu_number"))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";")%>%
  group_by(sample_name,Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample_name)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()

# STATS -- SET SEED -------------------------------------------------------
set.seed(295034) # Setting the seed before we do any stats


# STATS ANOVA -- OTUs One-way -----------------------------------------------------
otu_aov <- otu_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  mutate(Organism = as.factor(Organism))%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)


# STATS ANOVA -- Quant Two-way -------------------------------------------------------------------
aov_pvalues <- quant_stats%>%
  group_by(feature_number)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism*DOM_fil, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

aov_organism_sigs <- (aov_pvalues%>%
               filter(!term == "DOM_fil"))$feature_number%>%
  as.factor()%>%
  unique()%>%
  as.vector()

aov_DOM_fil_sigs <- (aov_pvalues%>%
                        filter(!term == "Organism"))$feature_number%>%
  as.factor()%>%
  unique()%>%
  as.vector()

aov_all_sigs <- (aov_pvalues)$feature_number%>%
  as.factor()%>%
  unique()%>%
  as.vector()

# STATS RANDOM FOREST -- QUANT Organism ----------------------------------------------
quant_org_rf_prep <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% aov_organism_sigs)%>%
  mutate(asin = as.numeric(asin))%>%
  spread(feature_number, asin)%>%
  select(-c(Experiment, biological_replicates, DOM_fil, technical_replicates))%>%
  mutate(Organism = as.factor(Organism))

names(quant_org_rf_prep) <- make.names(names(quant_org_rf_prep))
  
quant_rf_org <- randomForest(Organism ~ ., quant_org_rf_prep, 
               importance = TRUE, proximity = TRUE, nPerm = 10,
               ntree = 50000, na.action = na.exclude)

rf_quant_org_mda <- quant_rf_org$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(MeanDecreaseAccuracy >= (top_n(., 30, MeanDecreaseAccuracy)%>%
                                                               arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30]~ "important",
                                             TRUE ~ "not important"),
         multiseries_important = case_when(multiseries >= (top_n(., 30, multiseries)%>%
                                                             arrange(-multiseries))$multiseries[30]~ "important",
                                           TRUE ~ "not important"),
         delicatissima_important = case_when(delicatissima >= (top_n(., 30, delicatissima)%>%
                                                                 arrange(-delicatissima))$delicatissima[30]~ "important",
                                             TRUE ~ "not important"),
         galaxiae_important = case_when(galaxiae >= (top_n(., 30, galaxiae)%>%
                                                       arrange(-galaxiae))$galaxiae[30]~ "important",
                                        TRUE ~ "not important"),
         hasleana_important = case_when(hasleana >= (top_n(., 30, hasleana)%>%
                                                       arrange(-hasleana))$hasleana[30]~ "important",
                                        TRUE ~ "not important"),
         subpacifica_important = case_when(subpacifica >= (top_n(., 30, subpacifica)%>%
                                                             arrange(-subpacifica))$subpacifica[30]~ "important",
                                           TRUE ~ "not important"))

write_csv(rf_quant_org_mda, "Analyzed/RF_quant_organism.csv")


# STATS RANDOM FOREST -- QUANT unfilfil  -----------------------------------------
quant_dom_rf_prep <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% aov_DOM_fil_sigs)%>%
  mutate(asin = as.numeric(asin))%>%
  spread(feature_number, asin)%>%
  select(-c(Experiment, biological_replicates, Organism, technical_replicates))%>%
  mutate(DOM_fil = as.factor(DOM_fil))

names(quant_dom_rf_prep) <- make.names(names(quant_dom_rf_prep))

quant_rf_dom <- randomForest(DOM_fil ~ ., quant_dom_rf_prep, 
                             importance = TRUE, proximity = TRUE, nPerm = 10,
                             ntree = 50000, na.action = na.exclude)

rf_quant_dom_mda <- quant_rf_dom$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(MeanDecreaseAccuracy >= (top_n(., 30, MeanDecreaseAccuracy)%>%
                                                                        arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30] ~ "important",
                                             TRUE ~ "not important"))

write_csv(rf_quant_dom_mda, "Analyzed/RF_quant_dom.csv")



# STATS RANDOM FOREST -- OTUs ---------------------------------------------
sig_otu <- otu_aov$Taxonomy%>%
  as.vector()

otu_rf_df <- otu_stats%>%
  filter(Taxonomy %in% sig_otu)%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)%>%
  separate(sample_name, c("Organism", "biological_replicate"), sep = "_")%>%
  select(-biological_replicate)%>%
  mutate(Organism = as.factor(Organism))

names(otu_rf_df) <- make.names(names(otu_rf_df))
  
otu_rf <-   randomForest(Organism ~ ., otu_rf_df, 
                         importance = TRUE, proximity = TRUE, nPerm = 10,
                         ntree = 50000, na.action = na.exclude)

otu_rf_mda <- otu_rf$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(MeanDecreaseAccuracy >= (top_n(., 30, MeanDecreaseAccuracy)%>%
                                                                        arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30]~ "important",
                                             TRUE ~ "not important"),
         multiseries_important = case_when(multiseries >= (top_n(., 30, multiseries)%>%
                                                             arrange(-multiseries))$multiseries[30]~ "important",
                                           TRUE ~ "not important"),
         delicatissima_important = case_when(delicatissima >= (top_n(., 30, delicatissima)%>%
                                                                 arrange(-delicatissima))$delicatissima[30]~ "important",
                                             TRUE ~ "not important"),
         galaxiae_important = case_when(galaxiae >= (top_n(., 30, galaxiae)%>%
                                                       arrange(-galaxiae))$galaxiae[30]~ "important",
                                        TRUE ~ "not important"),
         hasleana_important = case_when(hasleana >= (top_n(., 30, hasleana)%>%
                                                       arrange(-hasleana))$hasleana[30]~ "important",
                                        TRUE ~ "not important"),
         subpacifica_important = case_when(subpacifica >= (top_n(., 30, subpacifica)%>%
                                                             arrange(-subpacifica))$subpacifica[30]~ "important",
                                           TRUE ~ "not important"))

important_org_otu <- (otu_rf_mda%>%
                        mutate(feature = gsub("\\.", ";", feature))%>%
                        top_n(30, MeanDecreaseAccuracy))$feature%>%
  unique()%>%
  as.vector()

write_csv(otu_rf_mda, "Analyzed/Otu_rf_mda.csv")


# POST-STATS -- MINI Quant Table organism ---------------------------------------------------
important_quant_org <- (rf_quant_org_mda%>%
    mutate(feature = gsub("X", "", feature))%>%
    filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

mini_quant_org <-quant_stats%>%
  filter(feature_number %in% important_quant_org)%>%
  spread(feature_number, asin)

write_csv(mini_quant_org, "Analyzed/mini_quant_org.csv")  

# POST-STATS -- Mini Quant Table unfil ------------------------------------------
important_quant_dom <- (rf_quant_dom_mda%>%
                          mutate(feature = gsub("X", "", feature))%>%
                          filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

mini_quant_dom <- quant_stats%>%
  filter(feature_number %in% important_quant_dom)%>%
  spread(feature_number, asin)

write_csv(mini_quant_dom, "Analyzed/mini_quant_dom.csv")  


# POST-STATS -- Mini RF Table OTUs -------------------------------------
important_quant_otu <- (otu_rf_mda%>%
                          mutate(feature = gsub("X", "", feature))%>%
                          filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

mini_quant_otu <- quant_stats%>%
  filter(feature_number %in% important_quant_otu)%>%
  spread(feature_number, asin)

write_csv(mini_quant_otu, "Analyzed/mini_quant_otu.csv")  


# PRE-MATRIX QUANT AND CAT -- Organism ---------------------------------------------
rf_sd <- (rf_quant_org_mda%>%
            filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy))%>%
            mutate(feature = gsub("X", "", feature)))$feature%>%
  as.vector()

cat_clean_org <- cat_stats%>%
  filter(FeatureID %in% rf_sd)%>%
  column_to_rownames("FeatureID")%>%
  data.matrix(rownames.force = NA)

canopus_available_features_org <- rownames(cat_clean_org)%>% as.vector()

quant_binary_org <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% canopus_available_features_org)%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature)%>%
  mutate(binary = case_when(asin > 0.01*max(asin) ~ 1,
                            TRUE ~ 0))%>%
  ungroup()%>%
  select(-asin)%>%
  spread(feature_number, binary)%>%
  column_to_rownames("feature")%>%
  data.matrix(rownames.force = NA)


# PRE-MATRIX QUANT AND CAT -- DOM_Fil ---------------------------------------------
cat_clean_dom <- cat_stats%>%
  filter(FeatureID %in% aov_DOM_fil_sigs)%>%
  column_to_rownames("FeatureID")%>%
  data.matrix(rownames.force = NA)

canopus_available_features_dom <- rownames(cat_clean_dom)%>% as.vector()

quant_binary_dom <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% canopus_available_features_dom)%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature)%>%
  mutate(binary = case_when(asin > 0.01*max(asin) ~ 1,
                            TRUE ~ 0))%>%
  ungroup()%>%
  select(-asin)%>%
  spread(feature_number, binary)%>%
  column_to_rownames("feature")%>%
  data.matrix(rownames.force = NA)



# MATRIX MULTIPLICATION --  Organism--------------------------------------------
matrix_multiplied_org <- quant_binary_org%*%cat_clean_org%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample_code")%>%
  gather(category, mult, 2:ncol(.))%>%
  filter(category != "DBE-O")%>%
  mutate(log10 = log10(mult + 1))%>%
  select(-mult)%>%
  spread(category, log10)

multi_matrix_tidy_org <- matrix_multiplied_org%>%
  gather(category, mult, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

# MATRIX MULTIPLICATION -- DOM_Fil--------------------------------------------
matrix_multiplied_dom <- quant_binary_dom%*%cat_clean_dom%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample_code")%>%
  gather(category, mult, 2:ncol(.))%>%
  filter(category != "DBE-O")%>%
  mutate(log10 = log10(mult + 1))%>%
  select(-mult)%>%
  spread(category, log10)

multi_matrix_tidy_dom <- matrix_multiplied_dom%>%
  gather(category, mult, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")


# STATS ANOVA -- org matrix -----------------------------------------------
aov_matrix_org <- multi_matrix_tidy_org%>%
  group_by(category)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(mult ~ Organism, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

matrix_aov_org_sig <- aov_matrix_org$category%>%
  as.vector()

write_csv(aov_matrix_org, "Analyzed/anova_pvals_matrix_org.csv")

# STATS ANOVA -- dom matrix -----------------------------------------------
aov_matrix_dom <- multi_matrix_tidy_dom%>%
  group_by(category)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(mult ~ DOM_fil, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

matrix_aov_dom_sig <- aov_matrix_dom$category%>%
  as.vector()

write_csv(aov_matrix_org, "Analyzed/anova_pvals_matrix_dom.csv")

# STATS RANDOM FOREST -- Matrix  Organism-------------------------------------------
multi_matrix_random_forest_df <- multi_matrix_tidy_org%>%
  filter(category %in% matrix_aov_org_sig)%>%
  spread(category, mult)%>%
  select(c(Organism, 7:ncol(.)))%>%
  mutate(Organism = as.factor(Organism))


names(multi_matrix_random_forest_df) <- make.names(names(multi_matrix_random_forest_df))

rf_matrix <- randomForest(Organism ~ ., multi_matrix_random_forest_df, 
                          importance = TRUE, proximity = TRUE, nPerm = 10,
                          ntree = 50000, na.action = na.exclude)

top30_org <- (rf_matrix$importance%>% 
                  as.data.frame()%>%
                  rownames_to_column("feature")%>%
                  top_n(10, MeanDecreaseAccuracy))$feature%>%
  as.vector()


rf_matrix_mda_org <- rf_matrix$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(feature %like any% top30_org ~ "important",
                                             TRUE ~ "not important"),
         multiseries_important = case_when(multiseries >= (top_n(., 10, multiseries)%>%
                                                             arrange(-multiseries))$multiseries[10]~ "important",
                                           TRUE ~ "not important"),
         delicatissima_important = case_when(delicatissima >= (top_n(., 10, delicatissima)%>%
                                                             arrange(-delicatissima))$delicatissima[10]~ "important",
                                           TRUE ~ "not important"),
         galaxiae_important = case_when(galaxiae >= (top_n(., 10, galaxiae)%>%
                                                             arrange(-galaxiae))$galaxiae[10]~ "important",
                                           TRUE ~ "not important"),
         hasleana_important = case_when(hasleana >= (top_n(., 10, hasleana)%>%
                                                             arrange(-hasleana))$hasleana[10]~ "important",
                                           TRUE ~ "not important"),
         subpacifica_important = case_when(subpacifica >= (top_n(., 10, subpacifica)%>%
                                                          arrange(-subpacifica))$subpacifica[10]~ "important",
                                           TRUE ~ "not important"))

write_csv(rf_matrix_mda_org,"./Analyzed/RF_matrix_organism_mda.05.csv")

ggplot(rf_matrix_mda_org, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity")

# STATS RANDOM FOREST -- Matrix  Fil_Unfil -------------------------------------------
multi_matrix_random_forest_UnfilFil_df <- multi_matrix_tidy_dom%>%
  filter(category %in% matrix_aov_dom_sig)%>%
  spread(category,mult)%>%
  select(c(DOM_fil, 7:ncol(.)))%>%
  mutate(DOM_fil = as.factor(DOM_fil))

names(multi_matrix_random_forest_UnfilFil_df) <- make.names(names(multi_matrix_random_forest_UnfilFil_df))

rf_matrix_UnfilFil <- randomForest(DOM_fil ~ ., multi_matrix_random_forest_UnfilFil_df,
                          importance = TRUE, proximity = TRUE,
                          ntree = 50000, na.action=na.exclude)

top30_unfil <- (rf_matrix_UnfilFil$importance%>%
                  as.data.frame()%>%
                  rownames_to_column("feature")%>%
                  top_n(15, MeanDecreaseAccuracy))$feature%>%
  as.vector()

rf_matrix_UnfilFil_mda <- rf_matrix_UnfilFil$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(feature %like any% top30_unfil ~ "important",
                                             TRUE ~ "not important"),
         dom_mda_important = case_when(DOM >= (top_n(., 15, DOM)%>%
                                                   arrange(-DOM))$DOM[15]~ "important",
                                       TRUE ~ "not important"),
         filt_mda_important = case_when(Unfil >= (top_n(., 15, Unfil)%>%
                                                 arrange(-Unfil))$Unfil[15]~ "important",
                                       TRUE ~ "not important"))

write_csv(rf_matrix_UnfilFil_mda,"./Analyzed/RF_matrix_UnfilFil_mda.05.csv")

ggplot(rf_matrix_UnfilFil_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity")



# STATS PERMANOVA - org and unfilfil  ---------------------------------------------------
matrix_permanova_org <- matrix_multiplied_org%>%
  gather(category, mult, 2:ncol(.))%>%
  mutate(mult = mult +1)%>%
  spread(category, mult)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

permanova_org <- matrix_permanova_org%>%
  # column_to_rownames("sample_code")%>%
  adonis(.[7:ncol(.)] ~ Organism*DOM_fil, ., perm = 1000, method = "bray", na.rm = TRUE) 

permanova_org

matrix_permanova_dom <- matrix_multiplied_dom%>%
  gather(category, mult, 2:ncol(.))%>%
  mutate(mult = mult +1)%>%
  spread(category, mult)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

permanova_dom <- matrix_permanova_dom%>%
  # column_to_rownames("sample_code")%>%
  adonis(.[7:ncol(.)] ~ Organism*DOM_fil, ., perm = 1000, method = "bray", na.rm = TRUE) 

permanova_dom


# POST-STATS -- mini-matrix organism -----------------------------------------------
important_org_compounds <- (rf_matrix_mda_org%>%
                              mutate(feature = gsub("X", "", feature))%>%
                              top_n(30, MeanDecreaseAccuracy))$feature%>%
  unique()%>%
  as.vector()

mini_matrix_org <- matrix_multiplied_org%>%
  gather(feature, val, 2:ncol(.))%>%
  mutate(feature = gsub("[[:space:]]", ".", feature))%>%
  mutate(feature = gsub("-", ".", feature))%>%
  filter(feature %in% important_org_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_org, "Analyzed/mini_matrix_important_org.csv")

# STATS - T-TEST Important features org ---------------------------------------
feature_info_test <- feature_info%>%
  gather(variable, response, 2:ncol(.))%>%
  mutate(importance = case_when(feature_number %in% important_quant_org ~ "important",
                                TRUE ~ "not"),
         importance = as.factor(importance))%>%
  group_by(variable)%>%
  nest()%>%
  mutate(data = map(data, ~ t.test(.x$response ~ .x$importance, alternative = "greater")),
                    p_value = map(data, ~ .x["p.value"][[1]]))%>%
  select(-data)%>%
  ungroup()%>%
  mutate(p_value = as.numeric(p_value),
         FDR = p.adjust(p_value, method = "BH"))

write_csv(feature_info_test, "Analyzed/Ttest_elements_org.csv")

# STATS - T-TEST Important features DOM ---------------------------------------
feature_info_test_dom <- feature_info%>%
  gather(variable, response, 2:ncol(.))%>%
  mutate(importance = case_when(feature_number %in% important_quant_dom ~ "important",
                                TRUE ~ "not"),
         importance = as.factor(importance))%>%
  group_by(variable)%>%
  nest()%>%
  mutate(data = map(data, ~ t.test(.x$response ~ .x$importance, alternative = "greater")),
         p_value = map(data, ~ .x["p.value"][[1]]))%>%
  select(-data)%>%
  ungroup()%>%
  mutate(p_value = as.numeric(p_value),
         FDR = p.adjust(p_value, method = "BH"))

write_csv(feature_info_test_dom, "Analyzed/Ttest_elements_dom.csv")


# STATS Correlation analysis ----------------------------------------------
## Correlation analysis
## Doing this between OTU and multiplied matrix
correlation_matrix <- matrix_multiplied_org%>%
  gather(compound, val, 2:ncol(.))%>%
  mutate(compound = gsub("[[:space:]]", ".", compound))%>%
  mutate(compound = gsub("-", ".", compound))%>%
  filter(compound %in% important_org_compounds)%>%
  spread(compound, val)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil",
                          "technical_replicates"), sep = "_")%>%
  select(-c(Experiment, technical_replicates))%>%
  group_by(Organism, biological_replicates, DOM_fil)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

correlation_microbe <- otu_stats%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)

correlation_table <- correlation_matrix%>%
  unite(sample_name, c("Organism", "biological_replicates"), sep = "_")%>%
  group_by(DOM_fil)%>%
  nest()%>%
  mutate(data = map(data, ~ left_join(.x, correlation_microbe, by = "sample_name")%>%
                      gather(microbe, microbe_asin, contains(";"))%>%
                      gather(category, category_asin, 2:31)))


correlation_pvals <- correlation_table%>%
  unnest(data)%>%
  ungroup()%>%
  group_by(DOM_fil, microbe, category)%>%
  filter(sum(category_asin) > 0)%>%
  nest()%>%
  mutate(corr = map(data, ~ cor.test(.x$category_asin, .x$microbe_asin, method = "pearson")%>%
                      broom::tidy()))%>%
  dplyr::select(-data)%>%
  unnest(corr)%>%
  ungroup()%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))

write_csv(correlation_pvals, "Analyzed/correlation_analysis.csv")


# POST-STATS -- mini-matrix unfilfil -----------------------------------------------
important_unfil_compounds <- (rf_matrix_UnfilFil_mda%>%
                              mutate(feature = gsub("X", "", feature))%>%
                                top_n(30, MeanDecreaseAccuracy))$feature%>%
  unique()%>%
  as.vector()

mini_matrix_dom <- matrix_multiplied_dom%>%
  gather(feature, val, 2:ncol(.))%>%
  mutate(feature = gsub("[[:space:]]", ".", feature))%>%
  filter(feature %in% important_unfil_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_dom, "Analyzed/mini_matrix_important_dom.csv")


# POST-STATS -- mini-matrix both org unfilfil important features ----------
mini_matrix_all <- mini_matrix_org%>%
  left_join(mini_matrix_dom, by = "sample_code")

write_csv(mini_matrix_all, "Analyzed/mini_matrix_important_all.csv")

# POST STATS -- matrix for HC ---------------------------------------------
otu_hc <- otu_stats%>%
  filter(Taxonomy %in% important_org_otu)%>%
  group_by(Taxonomy)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  select(-c(reads, ra, asin))%>%
  spread(Taxonomy, zscore)%>%
  rename("sample_code" = "sample_name")

hc_matrix <- mini_matrix_org%>%
  gather(category, asin, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism",
                          "biological_replicate", "DOM_fil",
                          "technical_replicate"), sep = "_")%>%
  group_by(category)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(category, Organism, biological_replicate, DOM_fil)%>%
  select(-technical_replicate)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  unite(sample_code, c("Organism", "biological_replicate"), sep = "_")%>%
  spread(category, zscore)%>%
  left_join(otu_hc, by = "sample_code")
  
  
write_csv(hc_matrix, "Analyzed/hc_matrix.csv")

# VISUALIZATION -- hc for quant features ----------------------------------
quant_hc <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% canopus_available_features_org)%>%
  unite(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  select(-asin)%>%
  spread(feature_number, zscore)

write_csv(quant_hc, "Analyzed/hc_features.csv")

# VISUALIZATION -- PCoA org and unfilfil -------------------------------------------------
##Quant all
pcoa_quant <- quant_stats%>%
  unite(sample, c("Experiment", "Organism", "biological_replicates", "DOM_fil", "technical_replicates"), sep = "_")%>%
  spread(feature_number, asin)%>%
  column_to_rownames("sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pcoa_quant$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## otu
pcoa_otu <- otu_stats%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)%>%
  column_to_rownames("sample_name")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pcoa_otu$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## Organism Matrix
pcoa_org <- matrix_multiplied_org%>%
  gather(cat, val, 2:ncol(.))%>%
  mutate(val = val+ min(val) +1)%>%
  spread(cat, val)%>%
  column_to_rownames("sample_code")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()
  

pcoa_org$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## Unfil vs Fil matrix
pcoa_dom <- matrix_multiplied_dom%>%
  gather(cat, val, 2:ncol(.))%>%
  mutate(val = val+ min(val) +1)%>%
  spread(cat, val)%>%
  column_to_rownames("sample_code")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()


pcoa_dom$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

#Settings for Pcoas
pcoa_settings <- function(x) {
  ggplot(x, aes(Axis.1, Axis.2, color = Organism, shape = DOM_fil)) +
  geom_point(stat = "identity") +
    scale_color_manual(values = c("#75d648", "#ae2da9", "#2d67c7", "#f27304", "#64d6f7")) +
    scale_shape_manual(values=c(1, 16))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(color="black"))
}


## Plotting PCoAs
pdf("Plots/PCoA_all.pdf", width = 7, height = 5)  
pcoa_quant$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>% 
  pcoa_settings() +
  ylab(str_c("Axis 2", " (", round(pcoa_quant$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_quant$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ggtitle("quant all")

pcoa_otu$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Organism", "biological_replicates"), sep = "_")%>% 
  ggplot(aes(Axis.1, Axis.2, color = Organism)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_color_manual(values = c("#75d648", "#ae2da9", "#2d67c7", "#f27304", "#64d6f7")) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(color="black"))+
  ylab(str_c("Axis 2", " (", round(pcoa_otu$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_otu$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ggtitle("otus")

pcoa_org$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>%  
  pcoa_settings() +
  ylab(str_c("Axis 2", " (", round(pcoa_org$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_org$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ggtitle("Organism 0.05")

pcoa_dom$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>%  
  pcoa_settings() +
  ylab(str_c("Axis 2", " (", round(pcoa_dom$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_dom$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ggtitle("UnfilFil 0.05")
dev.off()


# VISUALIZATION -- Mean Decrease Accuracy scatterplots -------------------------
mda_theme <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                   axis.title.x = element_text(size=14, face="bold"),
                   axis.title.y = element_text(size=14, face="bold"),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   panel.grid.major = element_blank(), # get rid of major grid
                   panel.grid.minor = element_blank(), # get rid of minor grid
                   legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                   legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                   axis.line = element_line(color="black"))

pdf("./Plots/Mean_Decrease_Accuracy.pdf", height = 5, width = 7)
ggplot(rf_quant_org_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("Organism QUANT Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (mean(rf_quant_org_mda$MeanDecreaseAccuracy + sd(rf_quant_org_mda$MeanDecreaseAccuracy))),
             col = "red") +
  mda_theme

ggplot(rf_quant_dom_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("DOM Quant Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(rf_quant_dom_mda, 30, MeanDecreaseAccuracy)%>%
                             arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme

ggplot(otu_rf_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("OTU Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(otu_rf_mda, 30, MeanDecreaseAccuracy)%>%
                             arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme

ggplot(rf_matrix_mda_org, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("Organism Matrix Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(rf_matrix_mda_org, 30, MeanDecreaseAccuracy)%>%
               arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme

ggplot(rf_matrix_UnfilFil_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("FilUnfil Matrix Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(rf_matrix_UnfilFil_mda, 30, MeanDecreaseAccuracy)%>%
               arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme
dev.off()





