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


rename_sample_codes_ms <- function(x) {
  new <- x%>%
    mutate(sample_code_ms = gsub(".mzML", "", sample_code_ms),
           sample_code_ms = gsub(".mzXML", "", sample_code_ms))%>%
    left_join(sample_rename%>%
                select(1:2), by = "sample_code_ms")%>%
    select(-sample_code_ms)
}

map <- purrr::map
select <- dplyr::select

tidy <- broom::tidy
rename <- dplyr::rename



# LOADING -- Dataframes ---------------------------------------------------
sample_rename <- read_csv("./Raw/Rename_MS_SampleIDs.csv")%>%
  rename(sample_code_ms = ID_MS,
         sample_code = ID_new)

quant_raw <- read_csv("./Raw/quant_all.csv")%>%
  select(-c(2:3))%>%
  rename(feature_number = 1)%>%
  gather(sample_code_ms, xic, 2:ncol(.))%>%
  rename_sample_codes_ms()%>%
  spread(sample_code, xic)


metadata_quant <- read_tsv("./Raw/Pn_metadata_table.tsv")%>%
  mutate(filename = gsub("mzXML", "mzML", filename))%>%
  rename(sample_code_ms = filename)%>%
  rename_sample_codes_ms()

lib_id <- read_tsv("Raw/GNPS_LibIds.tsv")%>%
  rename(feature_number = '#Scan#')%>%
  mutate(feature_number = as.character(feature_number))

cat_df <- read_csv("./Raw/Input_Canopus.csv")

chl <- read_xlsx("Raw/Pn_Ex2_Chlorophyll.xlsx")

feature_info <- read_csv("./Raw/ZODIAC_ElementalComposition_input.csv")%>%
  rename(feature_number = 1)%>%
  select(feature_number, everything())

sample_rename_16S <- read_csv("./Raw/Rename_16S_SampleIDs.csv")%>%
  rename(sample_code_16S_old = ID_16S,
         sample_code_16S = ID_16S_new)

otu_taxonomy <- read_csv("./Raw/Pn_16S_no-plastids_taxonomy_fixed.csv")%>%
  rename("#OTU ID" = "Feature ID")

otu_df <- read_csv("./Raw/no-plastids-dada2-table.csv")%>%
  gather(sample_code_16S_old, reads, 2:ncol(.))%>%
  right_join(sample_rename_16S, by = "sample_code_16S_old")%>%
  select(-sample_code_16S_old)%>%
  spread(sample_code_16S, reads)%>%
  left_join(otu_taxonomy,by = "#OTU ID")%>%
  mutate(Taxon = gsub('D_[0-9]__', '', Taxon))%>%
  select(-c("#OTU ID", "Confidence"))%>%
  select("Taxon",everything())
  


# CLEANING -- Removing Blanks ---------------------------------------------
field_blanks <- (metadata_quant%>%
                      filter(SampleType == "blank_extraction"))$sample_code%>%
  as.vector()

culture_blanks <- (metadata_quant%>%
                     filter(SampleType == "blank_culturemedia",
                            sample_code != "Exp1_MediaBlank_A_DOM_100"))$sample_code%>%
  as.vector()

culture_samples <- (metadata_quant%>%
                     filter(ATTRIBUTE_Experiment == "Exp2_Culture"))$sample_code%>%
  as.vector()

quant_blanks_env <- quant_raw%>%
  flag_background(blank_columns =  match(names(select(., field_blanks)), names(.)))%>%
  filter(background_features == "real")%>%
  select(-background_features)

quant_culture_blanks_removed <- quant_blanks_env%>%
  select(c(feature_number, culture_blanks, culture_samples))%>%
  flag_background(blank_columns = match(names(select(., culture_blanks)), names(.)))%>%
  filter(background_features == "real")%>%
  select(-background_features)
  
quant_df <- quant_blanks_env%>%
  select(-c(culture_blanks, culture_samples))%>%
  right_join(quant_culture_blanks_removed, by = "feature_number")

# CLEANING -- Stats dataframes --------------------------------------------------------------------------------
## Cleaning all of the data
quant_stats <- quant_df%>%
  gather(sample_code, xic, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_", remove = FALSE)%>%
  filter(Experiment == "Exp2")%>%
  unite(sample_code, c("Organism", "biological_replicates"), sep = "_", remove = FALSE)%>%
  left_join(chl%>%
              mutate(sample_code = gsub("Pn_", "Pn-", sample_code)), by = "sample_code")%>%
  select(-sample_code)%>%
  filter(Organism != 'MediaBlank')%>%
  unite(sample_code, c("Experiment", "Organism", 
                      "biological_replicates", "DOM_fil", 
                      "technical_replicates"), sep = "_")%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic),
         chl_norm = ra/chl,
         asin = asin(sqrt(chl_norm)))%>%
  ungroup()%>%
  select(-c(ra, xic,  chl_norm, chl))%>%
  # group_by(feature_number)%>%
  # filter(sum(.$asin) != 0)%>%
  separate(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")
  # ungroup()

cat_stats <- cat_df%>%
  filter(ZodiacScore > 0.98)%>%
  select(1, 24:ncol(.))%>%
  gather(cat, prob, 2:ncol(.))%>%
  group_by(cat)%>%
  filter(max(prob) > 0.5)%>%
  ungroup()%>%
  mutate(asin = asin(sqrt(prob)))%>%
  select(-prob)%>%
  spread(cat, asin)

# CLEANING -- OTU table -------------------------------------------------------------------
otu_clean <- otu_df%>%
  rownames_to_column("otu_number")%>%
  gather(sample_code_16S, reads, 3:ncol(.))%>%
  unite(Taxonomy, c("otu_number", "Taxon"), sep = ";")


# PRE-STATS -- OTU TABLE --------------------------------------------------
## Making the stats dataframes for OTU, family and classes
otu_stats <- otu_clean%>%
  separate(sample_code_16S, c("Experiment", "Organism", "biological_replicate", "technical_replicates"), sep = "_")%>%
  filter(Experiment == "Exp2")%>%
  unite(sample_code_16S, c("Experiment","Organism", "biological_replicate", "technical_replicates"), sep = "_")%>%
  group_by(Taxonomy)%>%
  filter(sum(reads) != 0)%>%
  filter(!Taxonomy %like any% c('%Eukaryota%', '%Cyanobacteria%'))%>%
  ungroup()%>%
  group_by(sample_code_16S)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))
  
  

family_stats <- otu_stats%>%
  separate(Taxonomy, c("otu_number", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")%>%
  select(-c("Genus", "otu_number"))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";")%>%
  mutate(Taxonomy = gsub(';NA', '', Taxonomy))%>%
  group_by(sample_code_16S,Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample_code_16S)%>%
  mutate(ra = reads/sum(reads, na.rm = TRUE),
  asin = asin(sqrt(ra)))%>%
  group_by(Taxonomy)%>%
  mutate(mean = mean(ra, na.rm = TRUE),
         sd = sd(ra, na.rm = TRUE))%>%
  #select(-c(reads, ra))%>%
  unique()


class_stats <- otu_stats%>%
  separate(Taxonomy, c("otu_number", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")%>%
  select(-c("Order", "Family", "Genus", "otu_number"))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class"), sep = ";")%>%
  mutate(Taxonomy = gsub(';NA', '', Taxonomy))%>%
  group_by(sample_code_16S,Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample_code_16S)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))%>%
  group_by(Taxonomy)%>%
  mutate(mean = mean(ra, na.rm = TRUE),
         sd = sd(ra, na.rm = TRUE))%>%
  #select(-c(reads, ra))%>%
  unique()


###write file for OTU table (all in Exp 2)

otu_tableSI <- otu_stats%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  filter(Experiment == "Exp2")%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()%>%
  mutate(Organism = as.factor(Organism))%>%
  group_by(Experiment, Taxonomy)%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates",
                       "technical_replicates"), sep = "_")%>%
  select(-c(ra, asin))%>%
  spread(sample_name, reads)

all_otu <- otu_tableSI$Taxonomy%>%
  as.vector()

write_csv(otu_tableSI, "Analyzed/OTU_Reads_all.csv")

# STATS -- SET SEED -------------------------------------------------------
set.seed(295034) # Setting the seed before we do any stats

# STATS PERMANOVA - ASV  ---------------------------------------------------
# STATS PERMANOVA - Quant org and unfilfil  ---------------------------------------------------
quant_permanova_org <- quant_stats%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                  "technical_replicates"), sep = "_")%>%
  spread(feature_number, asin)

  
  gather(category, mult, 2:ncol(.))%>%
  mutate(mult = mult +1)%>%
  spread(category, mult)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

##still needs to be finished....


# STATS ANOVA -- OTUs One-way -----------------------------------------------------
otu_aov <- otu_stats%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                          "biological_replicates", "technical_replicates"), sep = "_")%>%
  filter(Experiment == "Exp2")%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()%>%
  mutate(Organism = as.factor(Organism))%>%
  group_by(Experiment, Taxonomy)%>%
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
  rownames_to_column("feature")

important_quant_org <- (rf_quant_org_mda%>%
                          mutate(feature = gsub("X", "", feature))%>%
                          filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

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
  rownames_to_column("feature")

write_csv(rf_quant_dom_mda, "Analyzed/RF_quant_dom.csv")

important_quant_dom <- (rf_quant_dom_mda%>%
                          mutate(feature = gsub("X", "", feature))%>%
                          filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

# STATS RANDOM FOREST -- OTUs ---------------------------------------------
sig_otu <- otu_aov$Taxonomy%>%
  as.vector()

otu_rf_df <- otu_stats%>%
  filter(Taxonomy %in% sig_otu)%>%
  select(-c(reads, ra))%>%
  mutate(Taxonomy = gsub(";", "_", Taxonomy),
         Taxonomy = gsub(" ", "SPACE", Taxonomy),
         Taxonomy = gsub("-", "LINE", Taxonomy))%>%
  spread(Taxonomy, asin)%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  filter(Experiment == "Exp2")%>%
  select(-c(Experiment, biological_replicates, technical_replicates))%>%
  mutate(Organism = as.factor(Organism))

names(otu_rf_df) <- make.names(names(otu_rf_df))
  
otu_rf <- randomForest(Organism ~ ., otu_rf_df, 
                         importance = TRUE, proximity = TRUE, nPerm = 10,
                         ntree = 50000, na.action = na.exclude)

otu_rf_mda <- otu_rf$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")
  

important_org_otu <- (otu_rf_mda%>%
                        mutate(feature = gsub("_", ";", feature),
                               feature = gsub("SPACE", " ", feature),
                               feature = gsub("LINE", "-", feature))%>%
                        filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  unique()%>%
  as.vector()

write_csv(otu_rf_mda, "Analyzed/Otu_rf_mda.csv")


# POST-STATS -- MINI Quant Table organism ---------------------------------------------------
mini_quant_org <-quant_stats%>%
  filter(feature_number %in% important_quant_org)%>%
  spread(feature_number, asin)

write_csv(mini_quant_org, "Analyzed/mini_quant_org.csv")  


# POST-STATS -- Mini Quant Table unfil ------------------------------------------
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

# STATS - T-TEST Important features elements org ---------------------------------------
feature_info_test <- feature_info%>%
  filter(ZodiacScore > 0.98)%>%
  inner_join(quant_culture_blanks_removed, by = "feature_number")%>%
  select(1:20)%>%
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

# STATS - T-TEST Important features elements DOM ---------------------------------------
feature_info_test_dom <- feature_info%>%
  filter(ZodiacScore > 0.98)%>%
  inner_join(quant_culture_blanks_removed, by = "feature_number")%>%
  select(1:20)%>%
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

# STATS - T-TEST Important features org canopus ---------------------------------------
feature_info_test_canopus <- cat_stats%>%
  rename('feature_number' = 'FeatureID')%>%
  inner_join(quant_culture_blanks_removed, by = "feature_number")%>%
  select(1:532)%>%
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

important_canopus_org <- (feature_info_test_canopus%>%
                            filter(FDR <= 0.05))$variable%>%
  as.vector()

write_csv(feature_info_test_canopus, "Analyzed/Ttest_canopus_org.csv")

# STATS - T-TEST Important features dom canopus ---------------------------------------
feature_info_test_DOM_canopus <- cat_stats%>%
  rename('feature_number' = 'FeatureID')%>%
  inner_join(quant_culture_blanks_removed, by = "feature_number")%>%
  select(1:532)%>%
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

important_canopus_dom <- (feature_info_test_DOM_canopus%>%
                            filter(FDR <= 0.05))$variable%>%
  as.vector()

write_csv(feature_info_test_DOM_canopus, "Analyzed/Ttest_canopus_dom.csv")


# PRE-MATRIX QUANT AND CAT -- Organism ---------------------------------------------
cat_clean_org <- cat_stats%>%
  filter(FeatureID %in% important_quant_org)%>%
  column_to_rownames("FeatureID")%>%
  select(important_canopus_org)%>%
  data.matrix(rownames.force = NA)


canopus_available_features_org <- rownames(cat_clean_org)%>% as.vector()

quant_binary_org <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% canopus_available_features_org)%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature_number)%>%
  mutate(binary = case_when(asin > mean(asin) ~ 1,
                            TRUE ~ 0))%>%
  ungroup()%>%
  select(-asin)%>%
  spread(feature_number, binary)%>%
  column_to_rownames("feature")%>%
  data.matrix(rownames.force = NA)

write.csv(quant_binary_org,"./Analyzed/quant_binary_org.csv")


# PRE-MATRIX QUANT AND CAT -- DOM_Fil ---------------------------------------------
cat_clean_dom <- cat_stats%>%
  filter(FeatureID %in% important_quant_dom)%>%
  column_to_rownames("FeatureID")%>%
  select(important_canopus_dom)%>%
  data.matrix(rownames.force = NA)

canopus_available_features_dom <- rownames(cat_clean_dom)%>% as.vector()

quant_binary_dom <- quant_stats%>%  ## Okay so here we are first making the data "tidy"
  filter(feature_number %in% canopus_available_features_dom)%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature)%>%
  mutate(binary = case_when(asin > mean(asin) ~ 1,
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
  mutate(log10 = log10(mult + 1))%>%
  select(-mult)%>%
  spread(category, log10)

write_csv(matrix_multiplied_org, "Analyzed/matrix_multiplied_org.csv")

multi_matrix_tidy_org <- matrix_multiplied_org%>%
  gather(category, mult, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

# MATRIX MULTIPLICATION -- DOM_Fil--------------------------------------------
matrix_multiplied_dom <- quant_binary_dom%*%cat_clean_dom%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample_code")%>%
  gather(category, mult, 2:ncol(.))%>%
  mutate(log10 = log10(mult + 1))%>%
  select(-mult)%>%
  spread(category, log10)

multi_matrix_tidy_dom <- matrix_multiplied_dom%>%
  gather(category, mult, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

# STATS PERMANOVA - Matrix org and unfilfil  ---------------------------------------------------
matrix_permanova_org <- matrix_multiplied_org%>%
  gather(category, mult, 2:ncol(.))%>%
  mutate(mult = mult +1)%>%
  spread(category, mult)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

permanova_matrix_org <- matrix_permanova_org%>%
  # column_to_rownames("sample_code")%>%
  adonis(.[7:ncol(.)] ~ Organism*DOM_fil, ., perm = 1000, method = "bray", na.rm = TRUE) 

permanova_matrix_org

matrix_permanova_dom <- matrix_multiplied_dom%>%
  gather(category, mult, 2:ncol(.))%>%
  mutate(mult = mult +1)%>%
  spread(category, mult)%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")

permanova_matrix_dom <- matrix_permanova_dom%>%
  # column_to_rownames("sample_code")%>%
  adonis(.[7:ncol(.)] ~ Organism*DOM_fil, ., perm = 1000, method = "bray", na.rm = TRUE) 

permanova_matrix_dom


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
                  top_n(30, MeanDecreaseAccuracy))$feature%>%
  as.vector()


rf_matrix_mda_org <- rf_matrix$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(feature %like any% top30_org ~ "important",
                                             TRUE ~ "not important"),
         multiseries_important = case_when(`Pn-multiseries` >= (top_n(., 10, `Pn-multiseries`)%>%
                                                             arrange(-`Pn-multiseries`))$`Pn-multiseries`[10]~ "important",
                                           TRUE ~ "not important"),
         delicatissima_important = case_when(`Pn-delicatissima` >= (top_n(., 10, `Pn-delicatissima`)%>%
                                                             arrange(-`Pn-delicatissima`))$`Pn-delicatissima`[10]~ "important",
                                           TRUE ~ "not important"),
         galaxiae_important = case_when(`Pn-galaxiae` >= (top_n(., 10, `Pn-galaxiae`)%>%
                                                             arrange(-`Pn-galaxiae`))$`Pn-galaxiae`[10]~ "important",
                                           TRUE ~ "not important"),
         hasleana_important = case_when(`Pn-hasleana` >= (top_n(., 10, `Pn-hasleana`)%>%
                                                             arrange(-`Pn-hasleana`))$`Pn-hasleana`[10]~ "important",
                                           TRUE ~ "not important"),
         subpacifica_important = case_when(`Pn-subpacifica` >= (top_n(., 10, `Pn-subpacifica`)%>%
                                                          arrange(-`Pn-subpacifica`))$`Pn-subpacifica`[10]~ "important",
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
                  top_n(30, MeanDecreaseAccuracy))$feature%>%
  as.vector()


rf_matrix_UnfilFil_mda <- rf_matrix_UnfilFil$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(feature %like any% top30_unfil ~ "important",
                                             TRUE ~ "not important"),
         dom_mda_important = case_when(DOM >= (top_n(., 30, DOM)%>%
                                                   arrange(-DOM))$DOM[30]~ "important",
                                       TRUE ~ "not important"),
         filt_mda_important = case_when(Unfil >= (top_n(., 30, Unfil)%>%
                                                 arrange(-Unfil))$Unfil[30]~ "important",
                                       TRUE ~ "not important"))

write_csv(rf_matrix_UnfilFil_mda,"./Analyzed/RF_matrix_UnfilFil_mda.05.csv")

ggplot(rf_matrix_UnfilFil_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity")




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
  mutate(feature = gsub(",", ".", feature))%>%
  filter(feature %in% important_org_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_org, "Analyzed/mini_matrix_important_org.csv")


# STATS Correlation analysis ----------------------------------------------
## Correlation analysis
## Doing this between OTU all and multiplied matrix
correlation_matrix_classes <- matrix_multiplied_org%>%
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

correlation_microbe_all <- otu_stats%>%
  filter(sample_code_16S %like% "%Exp2%")%>%
  rename('sample_name' = 'sample_code_16S')%>%
  filter(Taxonomy %like any% all_otu)%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)

correlation_table_classes <- correlation_matrix_classes%>%
  filter(DOM_fil != 'Unfil')%>%
  mutate(Experiment = 'Exp2')%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates"), sep = "_")%>%
  left_join(correlation_microbe_all%>%
              mutate(sample_name = gsub('.{2}$', "", sample_name)), by = "sample_name")%>%
  gather(microbe, microbe_asin, contains(";"))%>%
  gather(category, category_asin, 3:32)


correlation_pvals_classes <- correlation_table_classes%>%
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

write_csv(correlation_pvals_classes, "Analyzed/correlation_analysis_otuAll_classes.csv")

## Doing this between OTU 29 sig and multiplied matrix
correlation_matrix_classes <- matrix_multiplied_org%>%
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

correlation_microbe_Top29 <- otu_stats%>%
  filter(sample_code_16S %like% "%Exp2%")%>%
  rename('sample_name' = 'sample_code_16S')%>%
  filter(Taxonomy %like any% sig_otu)%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)

correlation_table_Top29_classes <- correlation_matrix_classes%>%
  filter(DOM_fil != 'Unfil')%>%
  mutate(Experiment = 'Exp2')%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates"), sep = "_")%>%
  left_join(correlation_microbe_Top29%>%
              mutate(sample_name = gsub('.{2}$', "", sample_name)), by = "sample_name")%>%
  gather(microbe, microbe_asin, contains(";"))%>%
  gather(category, category_asin, 3:32)


correlation_pvals_Top29_classes <- correlation_table_Top29_classes%>%
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

write_csv(correlation_pvals_Top29_classes, "Analyzed/correlation_analysis_otu29_classes.csv")

## Doing this between OTU families and multiplied matrix
correlation_microbe_family <- family_stats%>%
  filter(sample_code_16S %like% "%Exp2%")%>%
  filter(sum(asin) != 0)%>%
  rename('sample_name' = 'sample_code_16S')%>%
  select(-c(reads, ra, mean, sd))%>%
  spread(Taxonomy, asin)

correlation_table_family_classes <- correlation_matrix_classes%>%
  filter(DOM_fil != 'Unfil')%>%
  mutate(Experiment = 'Exp2')%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates"), sep = "_")%>%
  left_join(correlation_microbe_family%>%
              mutate(sample_name = gsub('.{2}$', "", sample_name)), by = "sample_name")%>%
  gather(microbe, microbe_asin, contains(";"))%>%
  gather(category, category_asin, 3:32)


correlation_pvals_family_classes <- correlation_table_family_classes%>%
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

write_csv(correlation_pvals_family_classes, "Analyzed/correlation_analysis_family_classes.csv")

## Doing this between OTU classes and multiplied matrix

correlation_microbe_class <- class_stats%>%
  filter(sample_code_16S %like% "%Exp2%")%>%
  filter(sum(asin) != 0)%>%
  rename('sample_name' = 'sample_code_16S')%>%
  select(-c(reads, ra, mean, sd))%>%
  spread(Taxonomy, asin)

correlation_table_class_classes <- correlation_matrix_classes%>%
  filter(DOM_fil != 'Unfil')%>%
  mutate(Experiment = 'Exp2')%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates"), sep = "_")%>%
  left_join(correlation_microbe_class%>%
              mutate(sample_name = gsub('.{2}$', "", sample_name)), by = "sample_name")%>%
  gather(microbe, microbe_asin, contains(";"))%>%
  gather(category, category_asin, 3:32)


correlation_pvals_class_classes <- correlation_table_class_classes%>%
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

write_csv(correlation_pvals_class_classes, "Analyzed/correlation_analysis_class_classes.csv")

## Correlation analysis
## Doing this between OTU and LibIDs
correlation_quant_LibID <- quant_stats%>%
  unite(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                   "technical_replicates"), sep = "_")%>%
  group_by(feature)%>%
  ungroup()%>%
  spread(feature_number, asin)%>%
  separate(feature, c("Experiment", "Organism", "biological_replicates", "DOM_fil",
                          "technical_replicates"), sep = "_")%>%
  select(-c(Experiment, technical_replicates))%>%
  group_by(Organism, biological_replicates, DOM_fil)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(Organism, biological_replicates, DOM_fil, "1795", "1855", "3988", "4471", "4582", "864", "5432", "5563", "6693", "3667", "4075", "7249", "4580", "8081")


correlation_microbe_all <- otu_stats%>%
  filter(sample_code_16S %like% "%Exp2%")%>%
  rename('sample_name' = 'sample_code_16S')%>%
  filter(Taxonomy %like any% all_otu)%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)

correlation_table_LibIDs <- correlation_quant_LibID%>%
  filter(DOM_fil != 'Unfil')%>%
  mutate(Experiment = 'Exp2')%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates"), sep = "_")%>%
  left_join(correlation_microbe_all%>%
              mutate(sample_name = gsub('.{2}$', "", sample_name)), by = "sample_name")%>%
  gather(microbe, microbe_asin, contains(";"))%>%
  gather(category, category_asin, 3:16)


correlation_pvals_LibIDs <- correlation_table_LibIDs%>%
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

write_csv(correlation_pvals_LibIDs, "Analyzed/correlation_analysis_LibIds.csv")


# POST-STATS -- mini-matrix unfilfil -----------------------------------------------

important_unfil_compounds <- (rf_matrix_UnfilFil_mda%>%
                              mutate(feature = gsub("X", "", feature))%>%
                              top_n(30, MeanDecreaseAccuracy))$feature%>%
  unique()%>%
  as.vector()

mini_matrix_dom <- matrix_multiplied_dom%>%
  gather(feature, val, 2:ncol(.))%>%
  mutate(feature = gsub("[[:space:]]", ".", feature))%>%
  mutate(feature = gsub("-", ".", feature))%>%
  mutate(feature = gsub(",", ".", feature))%>%
  filter(feature %in% important_unfil_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_dom, "Analyzed/mini_matrix_important_dom.csv")

# POST STATS -- matrix for HC ---------------------------------------------
otu_hc <- otu_stats%>%
  filter(Taxonomy %in% sig_otu)%>%
  group_by(Taxonomy)%>%
  filter(sample_code_16S %like% "Exp2%")%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  select(-c(reads, ra, asin))%>%
  ungroup()%>%
  spread(Taxonomy, zscore)
  # rename("sample_code" = "sample_code_16S")

compound_org_hc <- mini_matrix_org%>%
  gather(category, asin, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism",
                          "biological_replicate", "DOM_fil",
                          "technical_replicate"), sep = "_")%>%
  group_by(category)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(category, Organism, biological_replicate, DOM_fil)%>%
  select(-technical_replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  #filter(DOM_fil == "DOM")%>%
  #select(-c(asin, DOM_fil))%>%
  select(-c(asin))%>%
  #unite(sample, c("Organism", "biological_replicate"), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicate", "DOM_fil"), sep = "_")%>%
  spread(category, zscore)
  # left_join(otu_hc, by = "sample_code")

compound_org_hc_DOM <- mini_matrix_org%>%
  gather(category, asin, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism",
                          "biological_replicate", "DOM_fil",
                          "technical_replicate"), sep = "_")%>%
  group_by(category)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(category, Organism, biological_replicate, DOM_fil)%>%
  select(-technical_replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  #unite(sample, c("Organism", "biological_replicate"), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicate"), sep = "_")%>%
  spread(category, zscore)
# left_join(otu_hc, by = "sample_code")

compound_dom_hc <- mini_matrix_dom%>%
  gather(category, asin, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism",
                          "biological_replicate", "DOM_fil",
                          "technical_replicate"), sep = "_")%>%
  group_by(category)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(category, Organism, biological_replicate, DOM_fil)%>%
  select(-technical_replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(asin))%>%
  unite(sample, c("Organism", "biological_replicate", "DOM_fil"), sep = "_")%>%
  spread(category, zscore)
# left_join(otu_hc, by = "sample_code")
  
feature_org_hc <- mini_quant_org%>%
  gather(feature_number, asin, 6:ncol(.))%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(feature_number, Organism, biological_replicates, DOM_fil)%>%
  select(-technical_replicates)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  left_join(lib_id%>%
              select(feature_number, Compound_Name), by = "feature_number")%>%
  unite(feature, c(feature_number, Compound_Name), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicates"), sep = "_")%>%
  spread(feature, zscore)

feature_dom_hc <- mini_quant_dom%>%
  gather(feature_number, asin, 6:ncol(.))%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(feature_number, Organism, biological_replicates, DOM_fil)%>%
  select(-technical_replicates)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(-c(asin))%>%
  unite(sample, c("Organism", "biological_replicates", "DOM_fil"), sep = "_")%>%
  left_join(lib_id%>%
              select(feature_number, Compound_Name), by = "feature_number")%>%
  unite(feature, c(feature_number, Compound_Name), sep = "_")%>%
  spread(feature, zscore)
  
write_csv(compound_org_hc, "Analyzed/compound_org_hc.csv")
write_csv(compound_org_hc_DOM, "Analyzed/compound_org_hc_DOM.csv")
write_csv(compound_dom_hc, "Analyzed/compound_dom_hc.csv")
write_csv(feature_org_hc, "Analyzed/feature_org_hc.csv")
write_csv(feature_dom_hc, "Analyzed/feature_dom_hc.csv")


# VISUALIZATION -- Cytoscape ----------------------------------------------
cyto_base <- quant_df%>%
  gather(sample_code, xic, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_", remove = FALSE)%>%
  unite(sample_code, c("Organism", "biological_replicates"), sep = "_", remove = FALSE)%>%
  left_join(chl%>%
              mutate(sample_code = gsub("Pn_", "Pn-", sample_code)), by = "sample_code")%>%
  select(-sample_code)%>%
  filter(!Organism %like% '%PPLBlank%',
         !Organism %like% '%MediaBlank%')%>%
  unite(sample_code, c("Experiment", "Organism", 
                          "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_", remove = FALSE)%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic),
         chl_norm = ra/chl,
         asin = asin(sqrt(chl_norm)))

cyto_exp2_dom <- cyto_base%>%
  group_by(DOM_fil, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'ra'))%>%
  filter(Experiment == "Exp2")%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  spread(DOM_fil, chl_norm)%>%
  select(-4)

cyto_exp2_org <- cyto_base%>%
  group_by(Organism, DOM_fil, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'ra'))%>%
  filter(Experiment == "Exp2")%>%
  unite(sample, c(Organism, DOM_fil), sep = "_")%>%
  spread(sample, chl_norm)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

cyto_exp2_org2 <- cyto_base%>%
  group_by(Organism, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'ra'))%>%
  filter(Experiment == "Exp2")%>%
  unite(sample, c(Organism), sep = "_")%>%
  spread(sample, chl_norm)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

cyto_exp1_org2 <- cyto_base%>%
  group_by(Organism, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'ra'))%>%
  filter(Experiment == "Exp1")%>%
  unite(sample, c(Experiment, Organism), sep = "_")%>%
  spread(sample, chl_norm)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

cyto_piers <- cyto_base%>%
  group_by(Organism, biological_replicates, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'chl_norm'))%>%
  filter(Experiment == "Piers")%>%
  unite(sample, c(Organism, biological_replicates), sep = "_")%>%
  spread(sample, ra)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

cyto_CCE <- cyto_base%>%
  group_by(Organism, biological_replicates, feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'chl_norm'))%>%
  filter(Experiment == "CCE-P1706")%>%
  unite(sample, c(Organism, biological_replicates), sep = "_")%>%
  spread(sample, ra)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

cyto_piers_mean <- cyto_base%>%
  group_by(feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'chl_norm'))%>%
  filter(Experiment == "Piers")%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  rename(mean_piers = 2)

cyto_CCE_mean <- cyto_base%>%
  group_by(feature_number)%>%
  select(-c('chl', 'xic', 'asin', 'chl_norm'))%>%
  filter(Experiment == "CCE-P1706")%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  rename(mean_CCE = 2)

cyto_xic_sum <- cyto_base%>%
  group_by(feature_number)%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm'))%>%
  filter(Experiment == "Exp2")%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  rename(sum_xic_exp2 = 2)

cyto_xic_sum_env <- cyto_base%>%
  group_by(feature_number)%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm'))%>%
  filter(Experiment != "Exp1")%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  rename(sum_xic_exp2env = 2)

cyto_full <- cyto_exp2_org%>%
  full_join(cyto_exp2_org2, by = "feature_number")%>%
  full_join(cyto_exp2_dom, by = "feature_number")%>%
  full_join(cyto_exp1_org2, by = "feature_number")%>%
  full_join(cyto_piers, by = "feature_number")%>%
  full_join(cyto_piers_mean, by = "feature_number")%>%
  full_join(cyto_CCE, by = "feature_number")%>%
  full_join(cyto_CCE_mean, by = "feature_number")%>%
  full_join(cyto_xic_sum, by = "feature_number")%>%
  full_join(cyto_xic_sum_env, by = "feature_number")

write_csv(cyto_full, "./Analyzed/cyto_node_table.csv")


# VISUALIZATION -- Venn Diagramfeatures-------------------------------------------
CCE_compare <- cyto_base%>%
  filter(Experiment == "CCE-P1706")%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm',"Experiment", "Organism", 
            "biological_replicates", "DOM_fil", 
            "technical_replicates"))%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic))%>%
  ungroup()%>%
  group_by(feature_number)%>%
  filter(max(ra) >= 0.001)%>%
  ungroup()%>%
  select(-ra)%>%
  filter(feature_number %in% important_quant_org)%>%
  group_by(feature_number)%>%
  filter(sum(xic) != 0)%>%
  ungroup()%>%
  spread(sample_code, xic)

length(important_quant_org)
length(CCE_compare$feature_number)


CCE_features <- cyto_base%>%
  filter(Experiment == "CCE-P1706")%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm',"Experiment", "Organism", 
            "biological_replicates", "DOM_fil", 
            "technical_replicates"))%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic))%>%
  ungroup()%>%
  group_by(feature_number)%>%
  filter(max(ra) >= 0.001)%>%
  ungroup()%>%
  select(-ra)%>%
  group_by(feature_number)%>%
  filter(sum(xic) != 0)%>%
  ungroup()%>%
  spread(sample_code, xic)



Piers_compare <- cyto_base%>%
  filter(Experiment == "Piers")%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm',"Experiment", "Organism", 
            "biological_replicates", "DOM_fil", 
            "technical_replicates"))%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic))%>%
  ungroup()%>%
  group_by(feature_number)%>%
  filter(max(ra) >= 0.001)%>%
  ungroup()%>%
  select(-ra)%>%
  filter(feature_number %in% important_quant_org)%>%
  group_by(feature_number)%>%
  filter(sum(xic) != 0)%>%
  ungroup()%>%
  spread(sample_code, xic)

length(important_quant_org)
length(Piers_compare$feature_number)

Piers_features <- cyto_base%>%
  filter(Experiment == "Piers")%>%
  select(-c('chl', 'ra', 'asin', 'chl_norm',"Experiment", "Organism", 
            "biological_replicates", "DOM_fil", 
            "technical_replicates"))%>%
  group_by(sample_code)%>%
  mutate(ra = xic/sum(xic))%>%
  ungroup()%>%
  group_by(feature_number)%>%
  filter(max(ra) >= 0.001)%>%
  ungroup()%>%
  select(-ra)%>%
  group_by(feature_number)%>%
  filter(sum(xic) != 0)%>%
  ungroup()%>%
  spread(sample_code, xic)

Env_compare <- inner_join(Piers_features, CCE_features, by = 'feature_number')

length(Env_compare$feature_number)


# VISUALIZATION -- Venn Diagram ASV  --------------------------------------
otu_vis <- otu_stats%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()

asv_env <- otu_vis%>%
  filter(Experiment == "CCE-P1706")%>%
  unite(sample_code, c("Experiment", "Organism", 
                       "biological_replicates", "technical_replicates"), sep = "_")%>%
  group_by(Taxonomy)%>%
  filter(max(ra) >= 0.0001)%>%
  ungroup()%>%
  select(-c(reads, asin))%>%
  spread(sample_code, ra)

asv_exp2 <- otu_vis%>%
  filter(Experiment == "Exp2")%>%
  unite(sample_code, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  group_by(Taxonomy)%>%
  filter(max(ra) >= 0.0001)%>%
  ungroup()%>%
  select(-c(reads, asin))%>%
  spread(sample_code, ra)

length(asv_exp2$Taxonomy)
length(asv_env$Taxonomy)

asv_overlap <- inner_join(asv_env, asv_exp2, by = 'Taxonomy')

length(asv_overlap$Taxonomy)

# VISUALIZATIONS -- Class/Family stacked bar chart--------------------------------------
#Defining the colors for the Class/Family plot
#First 3 are Alpha's = orangy red
#1 Bacterioplankton = green
#9 Gammaproteobacteria = purple:blue
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
brewer.pal.info

colors_taxonomy <- c("#f27304", "#F68F1D", "#FAAC35",
                     "#2d67c7",
                     "#7C28A8", "#ae2da9",
                     "#75D648","#68C946","#5BBC43","#4EAF41","#40A13E","#33943C","#268739")
##greybrown757761, purplemix 952BA9
otu_vis%>%  
  filter(Experiment == "Exp2")%>%
  group_by(Taxonomy)%>%
  filter(sum(reads) > 0)%>%
  ungroup()%>%
  group_by(Organism, Taxonomy)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  separate(Taxonomy, c("otu_number", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")%>%
  unite(tax, c('Class', 'Family'), sep = '; ', remove = FALSE)%>%
  mutate(Organism = gsub('Pn-', '', Organism))%>%
  ggplot(aes(Organism, ra, fill = tax)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = colors_taxonomy) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_text(size=14, face="bold", color = 'black'),
        axis.title.y = element_text(size=14, face="bold", color = 'black'),
        axis.text.x = element_text(face="bold", angle = 60, hjust = 1, size=14, color = 'black'),
        axis.text.y = element_text(face="bold", size=14, color = 'black'),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(color="black"))


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

## otu Exp1 and 2 del, sub, mul

pcoa_otu_bothExp <- otu_stats%>%
  filter(sample_code_16S != "Exp1_Pn-subpacifica_A_I", 
         sample_code_16S != "Exp1_Pn-subpacifica_A_II",
         sample_code_16S != "Exp1_Pn-delicatissima_B_I",
         sample_code_16S != "Exp1_Pn-delicatissima_B_II",)%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  filter(Experiment != "CCE-P1706", 
         Experiment != "Piers")%>%
  filter(Organism != "Pn-hasleana", 
         Organism != "Pn-galaxiae")%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates",
                              "technical_replicates"), sep = "_")%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)%>%
  column_to_rownames("sample_name")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pcoa_otu_bothExp$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## otu Exp2

pcoa_otu_Exp2 <- otu_stats%>%
  separate(sample_code_16S, c("Experiment", "Organism", 
                              "biological_replicates", "technical_replicates"), sep = "_")%>%
  filter(Experiment == "Exp2")%>%
  group_by(Taxonomy)%>%
  filter(sum(asin) != 0)%>%
  ungroup()%>%
  unite(sample_name, c("Experiment", "Organism", "biological_replicates",
                       "technical_replicates"), sep = "_")%>%
  select(-c(reads, ra))%>%
  spread(Taxonomy, asin)%>%
  column_to_rownames("sample_name")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pcoa_otu_Exp2$values[1:10,]%>%
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
  geom_point(stat = "identity", size = 4) +
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

pcoa_otu_bothExp$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "technical_replicates"), sep = "_")%>%
  ggplot(aes(Axis.1, Axis.2, color = Organism, shape = Experiment)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#75d648", "#f27304", "#64d6f7")) +
  scale_shape_manual(values=c(17,16))+
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
  ylab(str_c("Axis 2", " (", round(pcoa_otu_bothExp$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_otu_bothExp$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ggtitle("otus")

pcoa_otu_Exp2$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "technical_replicates"), sep = "_")%>%
  ggplot(aes(Axis.1, Axis.2, color = Organism, shape = Experiment)) +
  geom_point(stat = "identity", aes(size = 0.1)) +
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
  ylab(str_c("Axis 2", " (", round(pcoa_otu_Exp2$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_otu_Exp2$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
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
  geom_hline(yintercept = (mean(rf_quant_org_mda$MeanDecreaseAccuracy + sd(rf_quant_org_mda$MeanDecreaseAccuracy))),
             col = "red") +
  mda_theme

ggplot(otu_rf_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("OTU Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
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





