# Load libraries ----

library(tidyverse)
library(jsonlite)
library(daymetr)



# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))



# proportion of herbivory per herb class----


# convert ordinal leaf damage scores into a weighted estimate of tissue loss, 
# aggregate them per site × week × year, and standardize by survey effort to obtain a continuous, 
# effort-corrected herbivory intensity suitable for phenological modeling

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



# Herb score using the mean of the scale ----

# here herbUsingMean is just using the mean/median of the herbivory scale to standardise the data

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
  mutate(totalHerbS = totalHerb/nSurv) # totalHerbS is the standardized total heribory (to account for survey effort)

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



# Arthropod data: prop of surv per week----
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


# Arthropod density: density per week----
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
 


# Herbivory + Arthropod data
Herb.Arthropod = fullHerb %>% select(-Latitude) %>% 
  left_join(prop.dens.arthropod, 
            by = c("Name", "Year", "julianweek"),
            relationship = "many-to-many") %>% # due to multiple observation methods in one site
  data.frame() %>% rename(nYearHerb = nYear) %>% filter(!is.na(ObservationMethod))




# Calculate herbivory anomaly-----





julianWindow = 120:230
min.nJulianWeekYearSite = 5
min.nSurvWeekYearSite = 10
min.nYear = 3



# Two different data cleaning pipeline leads to sightly different outcome:

# 1. Define a good JulianWeek
# 2. Define good nSurv
# 3. Define good survey years

goodJulianWeek = Herb.Arthropod %>% 
  filter(julianweek %in% julianWindow) %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite)

 

goodnSurv = fullDataset %>% 
  right_join(goodJulianWeek, by = c("Name", "ObservationMethod", "Year")) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(nSurv = n_distinct(ID)) %>% 
  filter(nSurv >= min.nSurvWeekYearSite)



goodnYear = Herb.Arthropod %>% 
  right_join(goodnSurv, by = c("Name","ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear = n_distinct(Year)) %>% 
  filter(nYear >= min.nYear) %>% 
  filter(!is.na(ObservationMethod))


goodHerbData = Herb.Arthropod %>% 
  inner_join(goodJulianWeek, by = c("Name", "ObservationMethod", "Year")) %>% 
  inner_join(goodnSurv, by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  inner_join(goodnYear, by = c("Name", "ObservationMethod"))

dim(goodHerbData)
#------------------------------------------------------------------------------------------

# A. Define what counts as herbivory (0-4). Already done in preparation of Herb.Arthropod dataframe
# 1. Define a good nSurv
# 2. Define a good JulianWindow
# 3. Define a good nJulianWeek
# 4. Define a good nYearHerb



goodnSurv2 = fullDataset %>% 
  right_join(Herb.Arthropod  %>% filter(julianweek %in% julianWindow),
             by = c("Name",  "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name,  ObservationMethod, Year, julianweek) %>% 
  summarise(nSurv = n_distinct(ID)) %>% 
  filter(nSurv >= min.nSurvWeekYearSite) %>% data.frame()


goodJulianWeek2 = Herb.Arthropod %>% 
  right_join(goodnSurv2, by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite)



goodnYear2 = Herb.Arthropod %>% 
  right_join(goodJulianWeek2, by = c("Name", "ObservationMethod", "Year")) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYearHerb = n_distinct(Year)) %>% 
  filter(nYearHerb >= min.nYear)



goodHerbData2 = Herb.Arthropod %>% 
  inner_join(goodnYear2, by  = c("Name", "ObservationMethod")) %>% # keep good years
  inner_join(goodJulianWeek2, by = c("Name", "ObservationMethod", "Year")) %>% # keep good weeks
  inner_join(goodnSurv2, by = c("Name", "ObservationMethod", "Year", "julianweek")) # Keep good surveys

dim(goodHerbData2)



anti_join(goodHerbData, goodHerbData2, 
          by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  select(Name, ObservationMethod, Year, julianweek, nYearHerb)



  
#######################################################################################################

# Fit linear regression models for each julian week  as predictor for total herbivory

goodHerbData2_collapsed = goodHerbData2 %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod)) %>% 
  select(siteObserv, Year, julianweek, totalHerbS) %>%
  group_by(siteObserv, Year, julianweek) %>%
  summarise(totalHerbS = mean(totalHerbS), .groups = "drop") # use mean here due to repeated record

modelGoodHerb = goodHerbData2_collapsed %>% 
  right_join(
    goodHerbData2_collapsed %>% 
      group_by(siteObserv, Year) %>% 
      summarise(nJulianweek = n_distinct(julianweek)) %>% 
      filter(nJulianweek >= min.nJulianWeekYearSite), # This may not be necessary anymore
    by = c("siteObserv", "Year")
  ) %>% 
  as.data.frame() %>% 
  mutate(
    julianweek_C = scale(julianweek),
    julianweek_N = (julianweek - min(julianweek)) / (max(julianweek) - min(julianweek))
  )



#####################################################################################


fit_herb_model =function(df) {
  
  mod <- lm(totalHerbS ~ julianweek_N, data = df)
  sm  <- summary(mod)
  
  coef_ci <- confint(mod)
  
  tibble(
    intercept = coef(mod)[1],
    effect    = round(coef(mod)[2], 2),
    r2        = round(sm$r.squared, 2),
    intercept_lwr = coef_ci[1, 1],
    intercept_upr = coef_ci[1, 2],
    effect_lwr    = coef_ci[2, 1],
    effect_upr    = coef_ci[2, 2]
  )
}


herbModelOutput = modelGoodHerb %>%
  group_by(siteObserv, Year) %>%
  group_modify(~ fit_herb_model(.x)) %>% # passes each group's time series into the model
  ungroup()


herbModelOutput.Lat = herbModelOutput %>% 
  separate(siteObserv, into = c("Name", "ObservationMethod"), sep = "_") %>% 
  left_join(fullDataset %>% select(Name, Latitude, Longitude) %>% 
              group_by(Name) %>% 
              summarise(Latitude = mean(Latitude),
                        Longitude = mean(Longitude)),
            by = c("Name"))

#########################################################################################



# Herbivory + caterpillar Anomaly ----

centroidHerbCat=  goodHerbData2 %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod)) %>% 
  group_by(siteObserv, Year) %>% 
  summarise(maxHerb = max(totalHerbS, na.rm = TRUE),
            centroidweek = sum(julianweek * totalHerbS, na.rm = TRUE)/sum(totalHerbS, na.rm = TRUE),
            maxH0.prop = max(H0.prop, na.rm = TRUE),
            centroidH0 = sum(julianweek * H0.prop, na.rm = TRUE)/sum(H0.prop, na.rm = TRUE),
            maxH1.prop = max(H1.prop, na.rm = TRUE),
            centroidH1 = sum(julianweek * H1.prop, na.rm = TRUE)/sum(H1.prop, na.rm = TRUE),
            maxH2.prop = max(H2.prop, na.rm = TRUE),
            centroidH2 = sum(julianweek * H2.prop, na.rm = TRUE)/sum(H2.prop, na.rm = TRUE),
            maxH3.prop = max(H3.prop, na.rm = TRUE),
            centroidH3 = sum(julianweek * H3.prop, na.rm = TRUE)/sum(H3.prop, na.rm = TRUE),
            maxH4.prop = max(H4.prop, na.rm = TRUE),
            centroidH4 = sum(julianweek * H4.prop, na.rm = TRUE)/sum(H4.prop, na.rm = TRUE),
            maxcaterpillar_prop = max(caterpillar_prop, na.rm = TRUE),
            centroidcaterpillar_prop = sum(julianweek * caterpillar_prop, na.rm = TRUE)/sum(caterpillar_prop, na.rm = TRUE),
            maxcaterpillar_density = max(caterpillar_density, na.rm = TRUE),
            centroidcaterpillar_density = sum(julianweek * caterpillar_density, na.rm = TRUE)/sum(caterpillar_density, na.rm = TRUE)) %>% 
  as.data.frame()


HerbCatAnomaly = centroidHerbCat %>% 
  left_join(
    centroidHerbCat %>% 
      group_by(siteObserv) %>% 
      summarise(
        mean.maxHerb = mean(maxHerb, na.rm = TRUE),
        mean.centroidweek = mean(centroidweek, na.rm = TRUE),
        mean.maxH0.prop = mean(maxH0.prop, na.rm = TRUE),
        mean.centroidH0 = mean(centroidH0, na.rm = TRUE),
        mean.maxH1.prop = mean(maxH1.prop, na.rm = TRUE),
        mean.centroidH1 = mean(centroidH1, na.rm = TRUE),
        mean.maxH2.prop = mean(maxH2.prop, na.rm = TRUE),
        mean.centroidH2 = mean(centroidH2, na.rm = TRUE),
        mean.maxH3.prop = mean(maxH3.prop, na.rm = TRUE),
        mean.centroidH3 = mean(centroidH3, na.rm = TRUE),
        mean.maxH4.prop = mean(maxH4.prop, na.rm = TRUE),
        mean.centroidH4 = mean(centroidH4, na.rm = TRUE),
        mean.maxcaterpillar_prop = mean(maxcaterpillar_prop, na.rm = TRUE),
        mean.centroidcaterpillar_prop = mean(centroidcaterpillar_prop, na.rm = TRUE),
        mean.maxcaterpillar_density = mean(maxcaterpillar_density, na.rm = TRUE),
        mean.centroidcaterpillar_density = mean(centroidcaterpillar_density, na.rm = TRUE)
      ),
    by = c("siteObserv")
  ) %>% 
  mutate(
    maxHerbAnomaly = maxHerb - mean.maxHerb,
    centroidweekAnomaly = centroidweek - mean.centroidweek,
    maxH0.propAnomaly = maxH0.prop - mean.maxH0.prop,
    centroidH0Anomaly = centroidH0 - mean.centroidH0,
    maxH1.propAnomaly = maxH1.prop - mean.maxH1.prop,
    centroidH1Anomaly = centroidH1 - mean.centroidH1,
    maxH2.propAnomaly = maxH2.prop - mean.maxH2.prop,
    centroidH2Anomaly = centroidH2 - mean.centroidH2,
    maxH3.propAnomaly = maxH3.prop - mean.maxH3.prop,
    centroidH3Anomaly = centroidH3 - mean.centroidH3,
    maxH4.propAnomaly = maxH4.prop - mean.maxH4.prop,
    centroidH4Anomaly = centroidH4 - mean.centroidH4,
    maxcaterpillar_propAnomaly = maxcaterpillar_prop - mean.maxcaterpillar_prop,
    centroidcaterpillar_propAnomaly = centroidcaterpillar_prop - mean.centroidcaterpillar_prop,
    maxcaterpillar_densityAnomaly = maxcaterpillar_density - mean.maxcaterpillar_density,
    centroidcaterpillar_densityAnomaly = centroidcaterpillar_density - mean.centroidcaterpillar_density
  ) %>% as.data.frame()


HerbCatAnomaly_herbModel = left_join(HerbCatAnomaly,
                                     herbModelOutput, 
                                     by = c("siteObserv", "Year")) %>% data.frame()




