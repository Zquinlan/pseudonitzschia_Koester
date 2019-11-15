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

map <- purrr::map
select <- dplyr::select

tidy <- broom::tidy
rename <- dplyr::rename

# LOADING -- Dataframes ---------------------------------------------------
quant_df <- read_csv("./Raw/Pn_Ex2_MASTERx_quant_raw.csv")
cat_df <- read_csv("./Raw/Pn_Ex2_MASTERx_Canopus_categories_probability.csv")
chl <- read_xlsx("Raw/Pn_Ex2_Chlorophyll.xlsx")

otu_df <- read_tsv("./Raw/Irina_2018_16s_exp_GEL.swarm.tax")
otu_samples <- read_csv("./Raw/Pn_16S_identifiers.csv")%>%
  rename("sample_name" = "SampleID")%>%
  rename("sample_code" = "OTU_name")

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
  gather(sample_code, reads, 2:ncol(.))%>%
  left_join(otu_samples, ., by = "sample_code")%>%
  separate(sample_name, c("Experiment", "Organism", "biological_replicate"), sep = "_")%>%
  unite(sample_name, c("Organism", "biological_replicate"), sep = "_")%>%
  select(-c(sample_code, Experiment))%>%
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
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  group_by(sample_name, Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")

# PRE-STATS -- OTU TABLE --------------------------------------------------
## Making the stats dataframes for OTU, family and classes
otu_stats <- otu_clean%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  group_by(sample_name)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))


family_stats <- otu_clean%>%
  select(-c(Genus, OTU))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";")%>%
  group_by(sample_name, Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample_name)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))

class_stats <- otu_clean%>%
  select(-c(Order, Family, Genus, OTU))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class"), sep = ";")%>%
  group_by(sample_name, Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample_name)%>%
  mutate(ra = reads/sum(reads),
         asin = asin(sqrt(ra)))  


# STATS -- SET SEED -------------------------------------------------------
set.seed(295034) # Setting the seed before we do any stats


# STATS ANOVA -- OTUs One-way -----------------------------------------------------
otu_aov <- otu_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

family_aov <- family_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

class_aov <- class_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)


# STATS TUKEY -- Family -----------------------------------------------------
family_tukey <- family_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      TukeyHSD(p.adjust.methods = "BH")%>%
                      tidy()))%>%
  unnest(data)%>%
  filter(adj.p.value < 0.05)

class_tukey <- class_stats%>%
  separate(sample_name, c("Organism", "Replicate"), sep = "_")%>%
  group_by(Taxonomy)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Organism, .x)%>%
                      TukeyHSD(p.adjust.methods = "BH")%>%
                      tidy()))%>%
  unnest(data)%>%
  filter(adj.p.value < 0.05)

write_csv(family_tukey, "Analyzed/family_tukey.csv")
write_csv(class_tukey, "Analyzed/class_tukey.csv")

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

# STATS RANDOM FOREST -QUANT NOT BINARY TEST ----------------------------------------------
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

top30_org_quant <- (quant_rf_org$importance%>% 
                as.data.frame()%>%
                rownames_to_column("feature")%>%
                top_n(30, MeanDecreaseAccuracy))$feature%>%
  as.vector()


rf_quant_org <- quant_rf_org$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(mean_decrease_important = case_when(feature %like any% top30_org_quant ~ "important",
                                             TRUE ~ "not important"))

write_csv(rf_quant_org, "Analyzed/RF_quant_organism.csv")

# POST-STATS -- MINI QUANT TABLE TEST ---------------------------------------------------
important_quant <- (rf_quant_org%>%
    mutate(feature = gsub("X", "", feature))%>%
    filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy) + sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

mini_quant <-quant_stats%>%
  filter(feature_number %in% important_quant)%>%
  spread(feature_number, asin)

write_csv(mini_quant, "Analyzed/mini_quant.csv")  

# PRE-MATRIX QUANT AND CAT -- Organism ---------------------------------------------
cat_clean_org <- cat_stats%>%
  filter(FeatureID %in% aov_organism_sigs)%>%
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

# STATS RANDOM FOREST -- Matrix  Organism-------------------------------------------
multi_matrix_random_forest_df <- multi_matrix_tidy_org%>%
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
                                             TRUE ~ "not important"))

write_csv(rf_matrix_mda_org,"./Analyzed/RF_matrix_organism_mda.05.csv")

ggplot(rf_matrix_mda_org, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity")

# STATS RANDOM FOREST -- Matrix  Fil_Unfil -------------------------------------------
multi_matrix_random_forest_UnfilFil_df <- multi_matrix_tidy_dom%>%
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
                                             TRUE ~ "not important"))

write_csv(rf_matrix_UnfilFil_mda,"./Analyzed/RF_matrix_UnfilFil_mda.05.csv")

ggplot(rf_matrix_UnfilFil_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity")


# STATS Correlation analysis ----------------------------------------------
## Correlation analysis
## Doing this between OTU and multiplied matrix
correlation_matrix <- matrix_multiplied_org%>%
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
                      gather(category, category_asin, 2:447)))


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
  mutate(FDR = p.adjust(p.value, method = "BH"))

write_csv(correlation_pvals, "Analyzed/correlation_analysis.csv")


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
  as.vector()

mini_matrix_org <- matrix_multiplied_org%>%
  gather(feature, val, 2:ncol(.))%>%
  mutate(feature = gsub("[[:space:]]", ".", feature))%>%
  filter(feature %in% important_org_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_org, "Analyzed/mini_matrix_important_org.csv")

# POST-STATS -- mini-matrix unfilfil -----------------------------------------------
important_unfil_compounds <- (rf_matrix_UnfilFil_mda%>%
                              mutate(feature = gsub("X", "", feature))%>%
                                top_n(30, MeanDecreaseAccuracy))$feature%>%
  as.vector()

mini_matrix_dom <- matrix_multiplied_dom%>%
  gather(feature, val, 2:ncol(.))%>%
  mutate(feature = gsub("[[:space:]]", ".", feature))%>%
  filter(feature %in% important_unfil_compounds)%>%
  spread(feature, val)

write_csv(mini_matrix_dom, "Analyzed/mini_matrix_important_dom.csv")



# POST STATS -- matrix for HC ---------------------------------------------
otu_hc <- otu_stats%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  select(-c(reads, ra, asin))%>%
  spread(Taxonomy, zscore)%>%
  rename("sample_code" = "sample_name")

hc_matrix <- mini_matrix_org%>%
  gather(category, asin, 2:ncol(.))%>%
  separate(sample_code, c("Experiment", "Organism",
                          "biological_replicate", "DOM_fil",
                          "technical_replicate"), sep = "_")%>%
  group_by(category, Organism, biological_replicate, DOM_fil)%>%
  select(-technical_replicate)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  unite(sample_code, c("Organism", "biological_replicate"), sep = "_")%>%
  spread(category, zscore)%>%
  left_join(otu_hc, by = "sample_code")
  
  
write_csv(hc_matrix, "Analyzed/hc_matrix.csv")


# VISUALIZATION -- PCoA org and unfilfil -------------------------------------------------
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

## Plotting PCoAs
pdf("Plots/PCoA_org_unfilfil_matrix_05.pdf", width = 7, height = 5)  
pcoa_org$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>%  
  ggplot(aes(Axis.1, Axis.2, color = Organism, shape = DOM_fil)) +
  geom_point(stat = "identity") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) +
  ylab(str_c("Axis 2", " (", round(pcoa_org$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_org$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
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
  ggtitle("Organism 0.05")

pcoa_dom$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample_code")%>%
  separate(sample_code, c("Experiment", "Organism", "biological_replicates", "DOM_fil", 
                          "technical_replicates"), sep = "_")%>%  
  ggplot(aes(Axis.1, Axis.2, color = Organism, shape = DOM_fil)) +
  geom_point(stat = "identity") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) +
  ylab(str_c("Axis 2", " (", round(pcoa_dom$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  xlab(str_c("Axis 1", " (", round(pcoa_dom$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
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
ggplot(rf_matrix_mda_org, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("Organism Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(rf_matrix_mda_org, 30, MeanDecreaseAccuracy)%>%
               arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme

ggplot(rf_matrix_UnfilFil_mda, aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("FilUnfil Mean Decrease Accuracy pval = 0.05") +
  xlab("Features (decreasing mda)") +
  ylab("Mean Decrease Accuracy") +
  geom_hline(yintercept = (top_n(rf_matrix_UnfilFil_mda, 30, MeanDecreaseAccuracy)%>%
               arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  mda_theme
dev.off()

