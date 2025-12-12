library(tidyverse)
library(daymetr)



Herb1 = fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1)
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
  summarise(
    Herbcount = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  mutate(
    Herb_category = case_when(
      Herb == 0 ~ "Herb_0",
      Herb == 1 ~ "Herb_1",
      Herb == 2 ~ "Herb_2",
      Herb == 3 ~ "Herb_3",
      Herb == 4 ~ "Herb_4",
    )
  ) %>% 
  pivot_wider(
    names_from = Herb_category,
    values_from = Herbcount
  ) %>% 
  left_join(
    fullDataset %>%   
      filter(
        !HerbivoryScore %in% c(-128, -1)
      ) %>% 
      group_by(Name, Year, julianweek) %>% 
      summarise(nSurv = n_distinct(ID)), by = c("Name", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  mutate(Herb_0 = Herb_0/ nSurv,
         Herb_1 = Herb_1/ nSurv,
         Herb_2 = Herb_2/ nSurv,
         Herb_3 = Herb_3/ nSurv,
         Herb_4 = Herb_4/ nSurv) %>% 
  rename(H0.prop = Herb_0,
         H1.prop = Herb_1,
         H2.prop = Herb_2,
         H3.prop = Herb_3,
         H4.prop = Herb_4)


herbUsingMean= fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
  ) %>% 
  group_by(Name, Latitude, Year, julianweek, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Latitude, Year, julianweek, Herb) %>% 
  summarise(
    Herbcount = n(),
    nSurv = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Herb,
    values_from = Herbcount
  ) %>% 
  mutate(
    H0 = `0` * 0,
    H1 = `1` * 3,
    H2 = `2` * 7,
    H3 = `3` * 17,
    H4 = `4` * 62
  ) %>% 
  mutate(
    totalHerb = rowSums(across(H0:H4), na.rm = TRUE)
  ) %>% 
  select(Name, Latitude, Year, julianweek, nSurv, totalHerb) %>% 
  as.data.frame() %>% 
  group_by(Name, Latitude, Year, julianweek) %>% 
  summarise(nSurv = sum(nSurv),
            totalHerb = sum(totalHerb)) %>% 
  mutate(Herb_standardized = totalHerb/nSurv) 

herbUsingMean2= herbUsingMean%>% 
  left_join(
    herbUsingMean %>% 
      group_by(Name) %>% 
      summarise(nYear = n_distinct(Year)),
    by = c("Name")) %>%
  mutate(Name = reorder(Name, Latitude))

fullHerb = Herb1 %>% 
  left_join(herbUsingMean2 %>% select(-nSurv),
                    by = c("Name", "Year", "julianweek")) %>% 
  rename(nSurvHerb = nSurv)



prop_fullDataset = fullDataset %>%
  filter( # potentially filter for (1) julian window and (2) sites
         WetLeaves == 0) %>% 
  group_by(Name,ObservationMethod, Year, julianweek, ID) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
            grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
            fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
            daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0)) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(caterpillar_prop = mean(caterpillar),
            spider_prop = mean(spider),
            beetle_prop = mean(beetle),
            truebug_prop = mean(truebug),
            hopper_prop  = mean(hopper),
            ant_prop = mean(ant),
            grasshopper_prop = mean(grasshopper),
            fly_prop = mean(fly),
            daddylonglegs_prop = mean(daddylonglegs),
            nSurv = n_distinct(ID))  



dens_fullDataset = fullDataset %>%
  filter( # potentially filter for (1) julian window and (2) sites
    WetLeaves == 0) %>% 
  group_by(Name, ObservationMethod, Latitude, Longitude, Year, julianweek) %>%
  summarize(nSurv = n_distinct(ID),
            caterpillar_density = sum(Group == 'caterpillar', na.rm = TRUE)/nSurv,
            spider_density = sum(Group == 'spider', na.rm = TRUE)/nSurv,
            beetle_density = sum(Group == 'beetle', na.rm = TRUE)/nSurv,
            truebug_density = sum(Group == 'truebugs', na.rm = TRUE)/nSurv,
            hopper_density = sum(Group == 'leafhopper', na.rm = TRUE)/nSurv,
            ant_density = sum(Group == 'ant', na.rm = TRUE)/nSurv,
            grasshopper_density = sum(Group == "grasshopper", na.rm = TRUE)/nSurv,
            fly_density = sum(Group == "fly", na.rm = TRUE)/nSurv,
            daddylonglegs_density = sum(Group == "daddylonglegs", na.rm = TRUE)/nSurv) 







prop.dens.arthropod = prop_fullDataset %>% 
  left_join(dens_fullDataset %>% select(-nSurv), 
            by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  rename(nSurvArthropod = nSurv) %>% 
  left_join( 
  dens_fullDataset %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYearArthropod = n_distinct(Year)), 
  by = c("Name", "ObservationMethod"))

sum(is.na(prop.dens.arthropod)) # all good
 


Herb.Arthropod = fullHerb %>% select(-Latitude) %>% 
  left_join(prop.dens.arthropod, 
            by = c("Name", "Year", "julianweek"),
            relationship = "many-to-many") %>% # due to multiple observation methods in one site
  data.frame() %>% rename(nYearHerb = nYear)



# Calculate herbivory anomaly-----

# Caterpillar density anomaly-----



Herb.Arthropod %>% 
  group_by(Name, ObservationMethod, Year, nYearArthropod) %>%
  summarise(
    maxOcc = max(caterpillar_density, na.rm = TRUE),
    CentroidOccurrence = mean(caterpillar_density, na.rm = TRUE),
    Centroid = sum(julianweek * caterpillar_density, na.rm = TRUE) / sum(caterpillar_density, na.rm = TRUE),
    .groups = "drop"
  ) %>% view()

Herb.Arthropod$nYearHerb 

# Caterpillar occurrence anomaly-----






#  Temperature data retrieval----


# -- First, retrieve the good sites based on the those that are good for the actual phenometrics
# 
 

# Call the good sites GoodSiteYearLatLon, which  should have Name, Latitude, Longitude, and Year columns.


tmp_file = tempfile(fileext = ".csv")

GoodSiteYearLatLon %>%
  rename(
    site = Name,
    lat = Latitude,
    lon = Longitude
  ) %>%
  write.csv(tmp_file, row.names = FALSE)

# pass temp CSV to function
TempSiteData = download_daymet_batch(
  file_location = tmp_file,
  start = min(GoodSiteYear$Year),
  end = 2024, # this is the most recent available in daymetr
  internal = TRUE
)

# remove temporary file 
unlink(tmp_file)

# The TempSiteData is a list of list. check TempSiteData[[1]]$
# so this function would extract the information and add the site information


TempSiteData_clean <- lapply(TempSiteData, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)
})

AllTempData = bind_rows(TempSiteData_clean)



## -- Combine centroids with Temperature----
TempAnomalPhenoCentroid = AllTempData %>% 
  filter(yday %in% TempDayWindow) %>% 
  group_by(site, year) %>% 
  summarise(meanTmin = mean(tmin..deg.c.),
            meanTmax = mean(tmax..deg.c.),
            meanPreci = mean(prcp..mm.day.)) %>% 
  left_join(
    AllTempData %>% 
      filter(yday %in% TempDayWindow) %>% 
      group_by(site) %>% 
      summarise(AllmeanTmin = mean(tmin..deg.c.),
                AllmeanTmax = mean(tmax..deg.c.),
                AllmeanPreci = mean(prcp..mm.day.)),
    by = c("site")) %>% 
  mutate(AnomalTmin = meanTmin - AllmeanTmin,
         AnomalTmax = meanTmax - AllmeanTmax,
         AnomalPreci = meanPreci - AllmeanPreci)%>% 
  right_join(anomalPhenoCentroid %>%  filter(Year != "2025"), # because daymetr has no 2025 yet.
             by = c("site" = "Name", "year" = "Year")) %>% 
  left_join(GoodSiteYearLatLon, by = c("site" = "Name")) %>% 
  mutate(Group = factor(Group, levels =c("Caterpillar",
                                         "Ant",
                                         "Spider")))  


