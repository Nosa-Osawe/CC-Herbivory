# Load libraries ----
rm(list = ls())

library(tidyverse)
library(jsonlite)
library(daymetr)
library(nlme)
library(performance)
library(corrplot)
library(nlme)
library(lme4)
 
# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))



# this was maually selected after taking a look at plants with broad spatial coverage using the Tableau software.
generaGoodSpatial = c("Acer", "Aesculus", "Amelanchier", "Betula", "Carpinus", "Carya", 
            "Celtis", "Cercis", "Cornus", "Crataegus", "Diospyros", "Elaeagnus", 
            "Fagus", "Fraxinus", "Hamamelis", "Ilex", "Juglans", "Kalmia", 
            "Lindera", "Liquidambar", "Liriodendron", "Lonicera", "Magnolia", 
            "Malus", "Morella", "Morus", "Nyssa", "Platanus", "Prunus", 
            "Quercus", "Rhododendron", "Sassafras", "Tilia", "Ulmus", 
            "Vaccinium", "Viburnum")



min.Genus.nSurv = 6
julianWindow = 120:230
min.nJulianWeekYearSite = 7
min.nSurvWeekYearSite = 20
minLatRange = 6

goodGenusSurv = fullDataset %>% 
  filter( julianweek %in% julianWindow,
    !HerbivoryScore %in% c(-128, -1),
    plantGenus %in% generaGoodSpatial
  ) %>% 
  group_by(Name, Latitude, Longitude, Year, ObservationMethod, julianweek, plantGenus) %>% 
  summarise(nSurv.genus = n_distinct(ID)) %>% 
  filter(nSurv.genus >= min.Genus.nSurv) %>% 
  as.data.frame()


goodGenusYear = goodGenusSurv %>% 
  group_by(Name, Latitude, Longitude, ObservationMethod, Year, plantGenus) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite) %>% 
  as.data.frame()


goodGenus = goodGenusYear %>% 
  group_by(plantGenus) %>% 
  summarise(latMin = min(Latitude),
            latMax = max(Latitude),
            latRange = round((latMax - latMin), 2)) %>% 
  filter(latRange >= minLatRange) # 1/2 A dosen good plant!



goodGenusYear %>% 
  inner_join(goodGenus, by = "plantGenus") %>% 
  group_by(plantGenus) %>% 
  summarise(nSite = n_distinct(Name),
            minLat = min(Latitude),
            maxLat = max(Latitude)) 

   # Seems like just a few plants can be used to study single genus herb-Lat patterns
