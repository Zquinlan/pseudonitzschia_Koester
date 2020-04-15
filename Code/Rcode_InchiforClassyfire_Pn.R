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

LibIDs <- read_tsv("./Raw/Pn_LibIds0.7.tsv")%>%
  gather(col, val, 2:ncol(.))%>%
  mutate(var = "library")%>%
  unite(col, c(col, var), sep = "_")%>%
  spread(col, val)

Analogs <- read_tsv("./Raw/Pn_Analogs0.8.tsv")%>%
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

##replace VZFRNCSOCOPNDB-AOKDLOFSSA-N
##with InChI=1S/C15H21NO6/c1-8(4-3-5-9(2)14(19)20)11-7-16-13(15(21)22)10(11)6-12(17)18/h3-5,9-11,13,16H,6-7H2,1-2H3,(H,17,18)(H,19,20)(H,21,22)/b5-3+,8-4-/t9-,10+,11-,13+/m1/s1
#mutate(inchi_combined = gsub("VZFRNCSOCOPNDB-AOKDLOFSSA-N", "InChI=1S/C15H21NO6/c1-8(4-3-5-9(2)14(19)20)11-7-16-13(15(21)22)10(11)6-12(17)18/h3-5,9-11,13,16H,6-7H2,1-2H3,(H,17,18)(H,19,20)(H,21,22)/b5-3+,8-4-/t9-,10+,11-,13+/m1/s1", inchi_combined))


#if it doesn't start with "InChi=" add "InChi=" (except N/A)

#convert inchi of combined table to inchikey for Classyfire
LibIds_inchikeys <- LibIDs_combined%>%
  select(`#Scan#`, inchi_combined)%>%
  filter(inchi_combined != "N/A")%>%
  group_by(`#Scan#`)%>%
  filter(row_number(inchi_combined) == 1 )%>%
  mutate(inchi_key = cs_inchi_inchikey(inchi_combined))%>%
  ungroup()

## write file to use for Classyfire at https://cfb.fiehnlab.ucdavis.edu/

write_csv(LibIds_inchikeys,"Pn_GNPS_combined_Inchikeys.csv")





##example
inchi <-  "InChI=1xxS/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-
2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1"
inchi <-  "1xxS/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-
2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1"
# convert InChI to CSID
cs_inchi_inchikey(inchi)


write_csv(LibIDs_combined,"./Analyzed/Pn_GNPS_Inchikeys_combined.csv")

###for just one table it would be:

failednaming <- read_csv("./Analyzed/Inchi_failednaming.csv")