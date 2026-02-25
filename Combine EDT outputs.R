##############################
# EDT Outputs Organization
#
# This script was created to take multiple excel files of EDT outputs and
# combine so all species and scenarios are in one file which can be the single 
# source of inputs into the LCM
#
# 2/12/26
# B Hamilton
##############################
setwd("C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs")

<<<<<<< HEAD

library(readxl)
library(dplyr)
library(stringr)
library(purrr)
=======
https://kleinschmidtgroup.sharepoint.com/:f:/r/sites/projects/Chehalis/Shared%20Documents/Internal%20Documents/26-27%20Biennium/KA_LCM/Data?csf=1&web=1&e=rSGvEu

C:\Users\Jason.Romine\OneDrive - Kleinschmidt Associates\Chehalis_LCM_devo\EDT_results

#===  Functions  

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


libs <- c("readxl", "dplyr", "stringr", "purrr")
ipak(libs)
>>>>>>> c9f8a3644c118351a7287ed0a3391ac6156bb48b

#For KS EDT
data_dir <- "C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs"

<<<<<<< HEAD
#For SEPA EDT
data_dir <- "C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/SEPA LCM/Chinook/data"
=======
data_dir <- "C:/Users/Jason.Romine/OneDrive - Kleinschmidt Associates/Chehalis_LCM_devo/EDT_results/Fall Chinook"
>>>>>>> c9f8a3644c118351a7287ed0a3391ac6156bb48b


#Path to files
files <- list.files(
  path = data_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)
files

#Make function to read and combine files
read_one <- function(f) {
  
  # --- pick correct sheet ---
  sheets <- excel_sheets(f)
  data_sheet <- sheets[!tolower(sheets) %in% 
                         c("properties","metadata","info")][1]
  
  df <- read_excel(f, sheet = 'Values')
  
  # --- parse filename ---
  fname <- tools::file_path_sans_ext(basename(f))
  parts <- stringr::str_split(fname, " MEDT ", simplify = TRUE)
  
  species  <- parts[1]
  scenario <- parts[2]
  
  # --- add + reorder columns ---
  df %>%
    dplyr::mutate(
      Species = Species,               # these are uppercase in SEPA, lowercase in KS
      Scenario = Scenario              #
    ) %>%
    dplyr::select(Species, Scenario, everything())
}

combined_df <- map_dfr(files, read_one)



#Save
write.csv(
  combined_df,
  file = file.path(data_dir, "Master SEPA EDT Outputs.csv"),
  row.names = FALSE
)
