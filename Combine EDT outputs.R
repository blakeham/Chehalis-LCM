##############################
# EDT Outputs Organization
#
# This script was created to take multiple excel files of EDT outputs and
# combine so all species and scenarios are in one file which can be the single 
# source of inputs into the LCM
#
# 2/12/26
#
##############################
setwd("C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs")

install.packages(c("readxl", "dplyr", "stringr", "purrr"))
library(readxl)
library(dplyr)
library(stringr)
library(purrr)

data_dir <- "C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs"

#Path to files
files <- list.files(
  path = data_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

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
      Species = species,
      Scenario = scenario
    ) %>%
    dplyr::select(Species, Scenario, everything())
}

combined_df <- map_dfr(files, read_one)

#Save
write.csv(
  combined_df,
  file = file.path(data_dir, "Master EDT Outputs.csv"),
  row.names = FALSE
)
