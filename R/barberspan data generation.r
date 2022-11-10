# Saving data for the barberspan location

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

rm(list=ls())

source("R/cwac_functions.r")

# load information regarding bird species
robertscsv <- read_excel('data/Roberts DB WRJD.xls')
migrantStatus <- read_excel("data/barberspanIDs.xlsx")

# ----------------------------------------------------------------------------
# BARBERSPAN data generation

# getting species list
species <- get_local_species('north%20west', '26352535')

# getting counts of species
counts <- species_counts('north%20west', '26352535', species$id)
counts <- as.data.frame(counts)

# save counts
save(counts, file = "data/barberspan_counts.RData")


# joining dataset with roberts DB using scientific name
robertsDB <- left_join(species, robertscsv, by = 'Scientific name')
robertsDB <- robertsDB %>% select('id', `species`, `Common name`, `Scientific name`, `Nest site`, `Foraging substratum`)
robertsDB <- left_join(robertsDB, migrantStatus, by = 'id')


# fill in missing nest site data
robertsDB[which(robertsDB$`Common name` == 'Eurasian Curlew'),]$`Nest site` <- 'Grass'
robertsDB[which(robertsDB$`Common name` == 'African Darter'),]$`Nest site` <- 'Diverse'
robertsDB[which(robertsDB$`Common name` == 'Bar-tailed Godwit'),]$`Nest site` <- 'Shrub'
robertsDB[which(robertsDB$`Common name` == 'Black-tailed Godwit'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Common Greenshank'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Common Black-headed Gull'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Montagus Harrier'),]$`Nest site` <- 'Diverse'
robertsDB[which(robertsDB$`Common name` == 'Black-headed Heron'),]$`Nest site` <- 'Tree'
robertsDB[which(robertsDB$`Common name` == 'African Sacred Ibis'),]$`Nest site` <- 'Tree'
robertsDB[which(robertsDB$`Common name` == 'Osprey'),]$`Nest site` <- 'Diverse'
robertsDB[which(robertsDB$`Common name` == 'Pink-backed Pelican'),]$`Nest site` <- 'Tree'
robertsDB[which(robertsDB$`Common name` == 'Caspian Plover'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Common Ringed Plover'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Greater Sand Plover'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Grey Plover'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Black-winged Pratincole'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Ruff'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Sanderling'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Buff-breasted Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Common Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Curlew Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Green Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Marsh Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Terek Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Wood Sandpiper'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Little Stint'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Caspian Tern'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'White-winged Tern'),]$`Nest site` <- 'Diverse'
robertsDB[which(robertsDB$`Common name` == 'Ruddy Turnstone'),]$`Nest site` <- 'Ground'
robertsDB[which(robertsDB$`Common name` == 'Grey Wagtail'),]$`Nest site` <- 'Diverse'
robertsDB[which(robertsDB$`Common name` == 'Yellow Wagtail'),]$`Nest site` <- 'Grass'
robertsDB[which(robertsDB$`Common name` == 'Cape Cormorant'),]$`Nest site` <- 'Diverse'
robertsDB$`Nest site`[is.na(robertsDB$`Nest site`)] <- ''

robertsDB[which(robertsDB$`Scientific name` == 'Phalacrocorax carbo'),]$`Common name` <- 'Great Cormorant'
robertsDB[which(robertsDB$`Scientific name` == 'Anas hybrid'),]$`Common name` <- 'Anas hybrid'
robertsDB[which(robertsDB$`Scientific name` == 'Alopochen aegyptiacus'),]$`Common name` <- 'Egyptian Goose'

robertsDB <- robertsDB[-c(16,96,101),]

# classifying species according to how many different types of foriging substratum they possess

robertsDB <- robertsDB %>%
  mutate(numForage = (str_count(`Foraging substratum`, '[0-9]')))

save(robertsDB, file = "data/robertsDB.RData")




