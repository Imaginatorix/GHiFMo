library(dplyr)
library(readr)

important_properties = c("smiles")

molecules <- read_delim("data/molecules.csv") %>% 
  select(important_properties)

cleaned <- distinct(molecules, smiles)
  # na.omit(cleaned)

write.table(cleaned, "data/cleaned.csv", row.names=FALSE, quote=FALSE)
