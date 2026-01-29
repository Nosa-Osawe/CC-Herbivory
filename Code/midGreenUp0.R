
# Load libraries ----
rm(list = ls())

library(tidyverse)
library(jsonlite)
library(terra)
library(sf)
library(dplyr)
library(purrr)


# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))



# Define the directory where your files are stored
data_dir = "Data/greenup"

#  Get list of all Greenup_0 files
file_list = list.files(data_dir, pattern = "Greenup_0.*\\.tif$", full.names = TRUE)




sites = fullDataset %>% 
  filter(Longitude >= -96) %>% 
  select(Name, Longitude, Latitude) %>% 
  group_by(Name, Longitude, Latitude) %>% 
  summarise(n()) %>% 
  select(Name, Longitude, Latitude ) %>%  as.data.frame()



#  Prepare your base sites (do this once outside the loop)
pts = st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)
# Define the target CRS string
target_crs = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +datum=WGS84"
pts_aea = st_transform(pts, target_crs)
buffers_vect = vect(st_buffer(pts_aea, dist = 1000)) # Convert to terra 'SpatVector' once (1km2)

# 4. The Loop (using map or lapply)
results_list = lapply(file_list, function(f) {
  
  # Extract year from filename (assuming format ...doyYYYY...)
  # This regex looks for 'doy' followed by 4 digits
  year_str = regmatches(f, regexpr("doy\\d{4}", f))
  year_val = gsub("doy", "", year_str)
  
  # Load and project
  r = rast(f)
  r_aea = project(r, target_crs)
  
  # Extract mean value for buffers
  # Extract returns a dataframe; we want the 2nd column (the values)
  ext_values = extract(r_aea, buffers_vect, fun = mean, na.rm = TRUE)[, 2]
  
  # Create a dataframe for this specific year
  # We convert to DOY relative to the specific year
  origin_date = as.Date(paste0(year_val, "-01-01"))
  
  df_year = data.frame(
    ext_values
  )
  
  # Rename column to include the year (e.g., Greenup_2018)
  colnames(df_year) = paste0("Greenup_", year_val)
  
  return(df_year)
})

# Combine results back to the original sites dataframe
all_greenup_data = do.call(cbind, results_list)
sites_greenup = cbind(sites, all_greenup_data)



#  Calculate DOY
# We create the Date object first, then subtract Jan 1st of that year
greenup_doy = sites_greenup %>%
  pivot_longer(
    cols = starts_with("Greenup_"), 
    names_to = "Year", 
    values_to = "Raw_Value",
    names_prefix = "Greenup_"
  ) %>%
  mutate(
    Year = as.numeric(Year),
    DOY = as.numeric(
      as.Date(Raw_Value, origin = "1970-01-01") - 
        as.Date(paste0(Year, "-01-01")) + 1)) 



# I have to find a way to deal with the NANs.
# we may approach this problem by estimations from sites which we already have data for.


greenupLm = lm(DOY    ~ Year + Latitude + Longitude , data = greenup_doy ,
               na.action = na.exclude)
summary(greenupLm)

greenup_doy$DOY_pred =predict(greenupLm, newdata = greenup_doy)


greenup_doy$DOY_filled = ifelse(
  is.na(greenup_doy$DOY),
  greenup_doy$DOY_pred,
  greenup_doy$DOY
)



summary(greenup_doy$DOY_filled[is.na(greenup_doy$DOY)])


greenup_doy$DOY_source = ifelse(
  is.na(greenup_doy$DOY),
  "predicted",
  "observed"
)

greenup_doy %>% 
  write.csv(file = "Data/greenup_doy.csv", row.names = FALSE)

##############################################################################################

greenup_doy = read.csv("Data/greenup_doy.csv")

greenup_doy %>% 
  filter(is.na(DOY)) %>% 
  write.csv(file = "Data/NANs_greenup_doy.csv", row.names = FALSE)
