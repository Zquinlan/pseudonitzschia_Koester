analog <- read_tsv("Raw/Pn_Analogs0.8.tsv")%>%
  rename(feature_number = '#Scan#')%>%
  mutate(feature_number = as.character(feature_number))

lib_id <- read_tsv("Raw/Pn_LibIds0.7.tsv")%>%
  rename(feature_number = '#Scan#')%>%
  mutate(feature_number = as.character(feature_number))

feature_org_LibId <- mini_quant_org%>%
  gather(feature_number, asin, 6:ncol(.))%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(feature_number, Organism, biological_replicates, DOM_fil)%>%
  select(-technical_replicates)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  #filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  inner_join(lib_id%>%
              select(feature_number, Compound_Name), by = "feature_number")%>%
  unite(feature, c(feature_number, Compound_Name), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicates"), sep = "_")

feature_org_analog <- mini_quant_org%>%
  gather(feature_number, asin, 6:ncol(.))%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(feature_number, Organism, biological_replicates, DOM_fil)%>%
  select(-technical_replicates)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  #filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  inner_join(analog%>%
               select(feature_number, Compound_Name), by = "feature_number")%>%
  unite(feature, c(feature_number, Compound_Name), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicates"), sep = "_")
  
feature_org_combined <- mini_quant_org%>%
  gather(feature_number, asin, 6:ncol(.))%>%
  group_by(feature_number)%>%
  mutate(zscore = (asin - mean(asin))/sd(asin))%>%
  ungroup()%>%
  group_by(feature_number, Organism, biological_replicates, DOM_fil)%>%
  select(-technical_replicates)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  #filter(DOM_fil == "DOM")%>%
  select(-c(asin, DOM_fil))%>%
  inner_join(combined_classyfire%>%
               select(feature_number,"Parent Level 1"), by = "feature_number")%>%
  unite(feature, c(feature_number, "Parent Level 1"), sep = "_")%>%
  unite(sample, c("Organism", "biological_replicates"), sep = "_")
