#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)
library(biclustermd)
library(webchem)
library(classyfireR)
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(gtable)


#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename


##read in files from GNPS: LibIds and Analogs

LibIDs <- read_csv("./Raw/Pn_LibIds0.7.csv")%>%
  gather(col, val, 2:ncol(.))%>%
  mutate(var = "library")%>%
  unite(col, c(col, var), sep = "_")%>%
  spread(col, val)

Analogs <- read_csv("./Raw/Pn_Analog0.8.csv")%>%
  gather(col, val, 2:ncol(.))%>%
  mutate(var = "analog")%>%
  unite(col, c(col, var), sep = "_")%>%
  spread(col, val)

##combine Inchi
###if LibId keep, if not use analog
LibIDs_combined <- LibIDs%>%
  full_join(Analogs, by = "#Scan#")%>%
  replace_na(list(INCHI_library = "N/A"))%>%
  replace_na(list(INCHI_analog = "N/A"))%>%
  mutate(inchi_combined = case_when(INCHI_library != "N/A" ~ INCHI_library,
                                    INCHI_library == "N/A" ~ INCHI_analog))

#if it doesn't start with "InChi=" add "InChi=" to fix the few wrong Inchi which would fail in Classyfire
LibIDs_combined_inchifixed <- LibIDs_combined%>%
  mutate(inchi_combined = case_when(!inchi_combined %like% "InChI=%" ~ paste("InChI=", inchi_combined, sep = ""),
                           TRUE ~ as.character(inchi_combined)))

#convert inchi of combined table to inchikey for Classyfire
Combined_classyfire <- LibIDs_combined_inchifixed%>%
  select(`#Scan#`, inchi_combined)%>%
  filter(inchi_combined != "InChI=N/A",inchi_combined != "InChI=n/a")%>%
  group_by(`#Scan#`)%>%
  filter(row_number(inchi_combined) == 1 )%>%
  mutate(inchi_key = cs_inchi_inchikey(inchi_combined))%>%
  ungroup()

## write file to use for Classyfire at https://cfb.fiehnlab.ucdavis.edu/

write_csv(Combined_classyfire,"./Raw/Pn_GNPS_combined_Inchikeys.csv")

## read file exported from Classyfire
Classyfire_download <- read_csv("./Raw/Pn_classyfire_download.csv")%>%
  rename(inchi_key = InChIKey)

#merge results from Classyfire back to feature numbers
### CAUTION: used bin not joined command, means it's merge by row number, check if correct
Classyfire <- bind_cols(Combined_classyfire, Classyfire_download)

write_csv(Classyfire,"./Raw/Pn_combined_Classyfire.csv")


