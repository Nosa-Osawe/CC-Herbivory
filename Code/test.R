
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
library(ggpmisc)


# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))






greenup_doy = read.csv("Data/greenup_doy.csv")

# Determine time steps away from mid greenup DOY ----
JN1.steps = 60
JN2.steps = JN1.steps + 70


# INTERESTING CC SITES.... (but feel free to check others!)

# [1] "Acadia NP - Alder"                    "Acadia NP - Sundew"                   "BHRC - Hemlock Draw"                 
# [4] "BHRC Pan Hollow Uplands"              "Beare Swamp in Rouge Park"            "Belle Isle State Park, Lancaster, VA"
# [7] "Birmingham Zoo AL"                    "Brookbanks Park"                      "Highlands Biological Station"        
# [10] "Middle Mill"                          "NC Botanical Garden"                  "Prairie Ridge Ecostation"            
# [13] "RVCC"                                 "Riverbend Park"                       "Stage Nature Center"                 
# [16] "UNC Chapel Hill Campus"               "Virginia Zoo"                         "Walker Nature Center"    




# INTERSTING CC PLANG GENUS .... (but feel free to check others!)

# "Acer"  "Quercus"   "Prunus"    "Fagus"   "Carya"

source("Code/herbivoryFunctions.R")

unique(fullDataset$Name)
sort(unique(fullDataset$plantGenus))



h0 = 0
h1 = 3
h2 = 7
h3 = 17
h4 = 62


min.nJulianWeekYearSite = 6   # minimum number of julianWeeks for each site-year
min.nSurvWeekYearSite = 10  # minimum number of site-year-week surveys
min.nYear = 3        # minimum number of survey site-year
TempDayWindow = 120: 230   # Temperature window

# get_site is a custom function


site = get_site( "All", fullDataset)
genusPlant = get_genera(plantGenus = "Fagus", fullDataset = fullDataset) # for now, this works well for one plantGenus ONLY.

fullDataset2 = fullDataset %>% 
  filter(
    plantGenus %in% genusPlant,
    Name %in% site
  )



# Simple mid-greenup EDA
greenup_doy %>% 
  group_by(Name) %>% 
  summarise(minDOY = min(DOY_filled),
            maxDOY = max(DOY_filled)) %>% 
  as.data.frame()  


# There is no 2025 Mid-greenUp data yet!
# We estimate 2025 as same as 2024 midgreenup DOY (Its our best guess!)


greenup_doy_2025 =  greenup_doy %>%
  filter(Year == 2024) %>%          # take only 2024 rows
  mutate(Year = 2025)               # relabel as 2025

greenup_doy_extended = bind_rows(
  greenup_doy,
  greenup_doy_2025
)


fullDatasetGreen = fullDataset2 %>% 
  left_join(greenup_doy_extended, by = c("Name", "Longitude", "Latitude", "Year")) %>% 
  filter(
    julianday >= DOY_filled + JN1.steps,  
    julianday <= DOY_filled + JN2.steps  
  )






Herb1 = fullDatasetGreen %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1)
  ) %>% 
  group_by(Name, plantGenus, Year, julianweek,  ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, plantGenus, Year, julianweek, Herb) %>% 
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
    values_from = Herbcount,
    values_fill = 0 
  ) %>% 
  left_join(
    fullDatasetGreen %>%   
      filter(
        !HerbivoryScore %in% c(-128, -1)
      ) %>% 
      group_by(Name, plantGenus, Year, julianweek) %>% 
      summarise(nSurv = n_distinct(ID)), by = c("Name", "plantGenus", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  mutate(
    across(starts_with("Herb_"), ~ .x / nSurv)
  ) %>% 
  rename_with(
    ~ sub("Herb_", "H.prop", .x),
    starts_with("Herb_")
  )


# Herb score using the mean of the scale ----

# here herbUsingMean is just using the mean/median of the herbivory scale to standardize the data

herb_weights = c(
  "0" = h0,
  "1" = h1,
  "2" = h2,
  "3" = h3,
  "4" = h4
)

herbUsingMean = fullDatasetGreen %>% 
  filter(!HerbivoryScore %in% c(-128, -1)) %>% 
  group_by(Name, plantGenus, Latitude, Year, julianweek, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, plantGenus, Latitude, Year, julianweek, Herb) %>% 
  summarise(
    Herbcount = n(),
    nSurv = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Herb,
    values_from = Herbcount,
    values_fill = 0
  ) %>% 
  mutate(
    totalHerb = rowSums(
      across(
        intersect(names(herb_weights), names(.)),
        ~ .x * herb_weights[cur_column()]
      ),
      na.rm = TRUE
    )
  ) %>% 
  select(Name, plantGenus, Latitude, Year, julianweek, nSurv, totalHerb) %>% 
  group_by(Name, plantGenus, Latitude, Year, julianweek) %>% 
  summarise(
    nSurv = sum(nSurv),
    totalHerb = sum(totalHerb),
    .groups = "drop"
  ) %>% 
  mutate(
    totalHerbS = totalHerb / nSurv
  )
# totalHerbS is the standardized total heribory (to account for survey effort)




herbUsingMean2= herbUsingMean%>% 
  left_join(
    herbUsingMean %>% 
      group_by(Name, plantGenus) %>% 
      summarise(nYear = n_distinct(Year)),
    by = c("Name", "plantGenus")) %>%
  mutate(Name = reorder(Name, Latitude))

fullHerb = Herb1 %>% 
  left_join(herbUsingMean2 %>% select(-nSurv),
            by = c("Name", "plantGenus", "Year", "julianweek")) %>% 
  rename(nSurvHerb = nSurv)



# Arthropod data: prop of surv per week----
prop_fullDataset2 = fullDatasetGreen %>%
  filter( # potentially filter for (1) julian window and (2) sites
    WetLeaves == 0) %>% 
  group_by(Name, plantGenus, ObservationMethod, Year, julianweek, ID) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
            grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
            fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
            daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0)) %>% 
  group_by(Name, plantGenus,  ObservationMethod, Year, julianweek) %>% 
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

dens_fullDataset2 = fullDatasetGreen %>%
  filter(WetLeaves == 0) %>% 
  group_by(Name, Latitude, Longitude, plantGenus, ObservationMethod, Year, julianweek) %>%
  summarise(
    nSurv = n_distinct(ID),
    
    caterpillar_density =
      sum(Quantity[Group == "caterpillar"], na.rm = TRUE) / nSurv,
    
    spider_density =
      sum(Quantity[Group == "spider"], na.rm = TRUE) / nSurv,
    
    beetle_density =
      sum(Quantity[Group == "beetle"], na.rm = TRUE) / nSurv,
    
    truebug_density =
      sum(Quantity[Group == "truebugs"], na.rm = TRUE) / nSurv,
    
    hopper_density =
      sum(Quantity[Group == "leafhopper"], na.rm = TRUE) / nSurv,
    
    ant_density =
      sum(Quantity[Group == "ant"], na.rm = TRUE) / nSurv,
    
    grasshopper_density =
      sum(Quantity[Group == "grasshopper"], na.rm = TRUE) / nSurv,
    
    fly_density =
      sum(Quantity[Group == "fly"], na.rm = TRUE) / nSurv,
    
    daddylonglegs_density =
      sum(Quantity[Group == "daddylonglegs"], na.rm = TRUE) / nSurv,
    
    .groups = "drop"
  ) %>% as.data.frame()








# Proportion and density data for all interesting arthropod 

prop.dens.arthropod = prop_fullDataset2 %>% 
  left_join(dens_fullDataset2 %>% select(-nSurv), 
            by = c("Name",  "plantGenus", "ObservationMethod", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  rename(nSurvArthropod = nSurv) %>% 
  left_join( 
    dens_fullDataset2 %>% 
      group_by(Name, plantGenus, ObservationMethod, ) %>% 
      summarise(nYearArthropod = n_distinct(Year)), 
    by = c("Name", "plantGenus", "ObservationMethod"))

sum(is.na(prop.dens.arthropod)) # all good

prop_fullDataset2[!complete.cases(prop_fullDataset2), ] %>% as.data.frame() # NA's in plantGenus

# Herbivory + Arthropod data
Herb.Arthropod = fullHerb %>% select(-Latitude) %>% 
  left_join(prop.dens.arthropod, 
            by = c("Name", "plantGenus", "Year", "julianweek"),
            relationship = "many-to-many") %>% # due to multiple observation methods in one site
  data.frame() %>% rename(nYearHerb = nYear) %>% filter(!is.na(ObservationMethod))


# Data wrangling pipeline:----

# A. Define what counts as herbivory (0-4). Already done in preparation of Herb.Arthropod dataframe
# 1. Define a good nSurv
# 2. Define a good JulianWindow
# 3. Define a good nJulianWeek
# 4. Define a good nYearHerb



goodnSurv2 = fullDatasetGreen %>% 
  right_join(Herb.Arthropod,
             by = c("Name", "plantGenus", "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name, plantGenus, ObservationMethod, Year, julianweek) %>% 
  summarise(nSurv = n_distinct(ID)) %>% 
  filter(nSurv >= min.nSurvWeekYearSite) %>% data.frame()


goodJulianWeek2 = Herb.Arthropod %>% 
  right_join(goodnSurv2, by = c("Name", "plantGenus", "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name, plantGenus, ObservationMethod, Year) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite)



goodnYear2 = Herb.Arthropod %>% 
  right_join(goodJulianWeek2, by = c("Name", "plantGenus", "ObservationMethod", "Year")) %>% 
  group_by(Name, plantGenus, ObservationMethod) %>% 
  summarise(nYearHerb = n_distinct(Year)) %>% 
  filter(nYearHerb >= min.nYear)



goodHerbData2 = Herb.Arthropod %>% 
  inner_join(goodnYear2, by  = c("Name", "plantGenus", "ObservationMethod" )) %>% # keep good years
  inner_join(goodJulianWeek2, by = c("Name", "plantGenus", "ObservationMethod", "Year")) %>% # keep good weeks
  inner_join(goodnSurv2, by = c("Name", "plantGenus", "ObservationMethod", "Year", "julianweek")) # Keep good surveys


# anti_join(goodHerbData, goodHerbData2, 
#           by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
#   select(Name, ObservationMethod, Year, julianweek, nYearHerb)
# 



#######################################################################################################

#  Julian week as a predictor of totalherbivory ----
# Fit linear regression models for each julian week  as predictor for total herbivory




goodHerbData2_collapsed = goodHerbData2 %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod)) %>% 
  group_by(siteObserv, Name, plantGenus, ObservationMethod, Year, julianweek) %>%
  summarise(totalHerbS = mean(totalHerbS), .groups = "drop") # use mean here due to repeated record


modelGoodHerb = goodHerbData2_collapsed %>% 
  right_join(
    goodHerbData2_collapsed %>% 
      group_by(siteObserv, plantGenus, Year) %>% 
      summarise(nJulianweek = n_distinct(julianweek)) %>% 
      filter(nJulianweek >= min.nJulianWeekYearSite), # This may not be necessary anymore
    by = c("siteObserv", "Year")
  ) %>% 
  as.data.frame() %>% 
  mutate(
    julianweek_C = scale(julianweek),
    julianweek_N = (julianweek - min(julianweek)) / (max(julianweek) - min(julianweek))
  )


# function that makes juilianweek predictor of total herbivory at a site-year----
# This may not work well if we choose more than one plantGenus per site.

fit_herb_model =function(df) {
  
  mod <- lm(totalHerbS ~ julianweek_N, data = df)
  sm  <- summary(mod)
  
  coef_ci <- confint(mod)
  
  tibble(
    intercept = coef(mod)[1],
    effect    = round(coef(mod)[2], 2),
    r2        = round(sm$r.squared, 2),
  )
}


herbModelOutput = modelGoodHerb %>%
  group_by(siteObserv, Year) %>%
  group_modify(~ fit_herb_model(.x)) %>% # passes each group's time series into the model
  ungroup()

# Add latitude information to the 'herbivory ~ julianweek' model
herbModelOutput.Lat = herbModelOutput %>% 
  separate(siteObserv, into = c("Name", "ObservationMethod"), sep = "_") %>% 
  left_join(fullDatasetGreen %>% select(Name, Latitude, Longitude) %>% 
              group_by(Name) %>% 
              summarise(Latitude = mean(Latitude),
                        Longitude = mean(Longitude)),
            by = c("Name"))

#########################################################################################



# Herbivory + caterpillar Anomaly ----


# Herbivory + caterpillar Anomaly ----

centroidHerbCat=  goodHerbData2 %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod)) %>% 
  group_by(siteObserv, Year) %>% 
  summarise(maxHerb = max(totalHerbS, na.rm = TRUE),
            centroidweek = sum(julianweek * totalHerbS, na.rm = TRUE)/sum(totalHerbS, na.rm = TRUE),
            maxH0.prop = max(H.prop0, na.rm = TRUE),
            centroidH0 = sum(julianweek * H.prop0, na.rm = TRUE)/sum(H.prop0, na.rm = TRUE),
            maxH1.prop = max(H.prop1, na.rm = TRUE),
            centroidH1 = sum(julianweek * H.prop1, na.rm = TRUE)/sum(H.prop1, na.rm = TRUE),
            maxH2.prop = max(H.prop2, na.rm = TRUE),
            centroidH2 = sum(julianweek * H.prop2, na.rm = TRUE)/sum(H.prop2, na.rm = TRUE),
            maxH3.prop = max(H.prop3, na.rm = TRUE),
            centroidH3 = sum(julianweek * H.prop3, na.rm = TRUE)/sum(H.prop3, na.rm = TRUE),
            maxH4.prop = max(H.prop4, na.rm = TRUE),
            centroidH4 = sum(julianweek * H.prop4, na.rm = TRUE)/sum(H.prop4, na.rm = TRUE),
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




###################################################################################################


## Correlate important variables----
corrplot(cor(HerbCatAnomaly_herbModel %>%
               select(maxHerb,
                      # maxH0.prop, maxH2.prop, maxH3.prop, maxH4.prop, 
                      effect, intercept, r2, centroidweek,
                      centroidcaterpillar_density, centroidcaterpillar_prop, maxcaterpillar_density, maxcaterpillar_prop,
                      centroidweekAnomaly, centroidcaterpillar_propAnomaly)%>%
               mutate(across(everything(), ~ifelse(is.infinite(.), NA, .))),
             use = "pairwise.complete.obs"),  
         method = "color", 
         type = "upper",          
         addCoef.col = "black",    
         tl.col = "black",         
         tl.srt = 45,              
         diag = FALSE)  

## Maximum total herb ~ max caterpillar prop occurrence ----

ggplot(data = HerbCatAnomaly_herbModel,
       aes(x= maxcaterpillar_prop, y = maxHerb)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )





ggplot(data = HerbCatAnomaly_herbModel,
       aes(x= centroidcaterpillar_prop, y = centroidweek)) + geom_point(size = 3, color = "darkgrey")+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid herbivory",
                        subtitle = "Slope Est: 0.49 in Random-intercept model")


summary(lme(
  centroidweek ~ centroidcaterpillar_prop,
  random = ~ 1 | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
))


centroid= lme(
  centroidweek ~ centroidcaterpillar_prop,
  random = ~ 1 | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
)
r2(centroid)



ggplot(data = HerbCatAnomaly_herbModel,
       aes(x= centroidcaterpillar_density, y = centroidweek)) + 
  geom_point(size = 3, color = "darkgrey")+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid herbivory",
                        subtitle = "Slope Est: 0.42 in Random-intercept model")


summary(lme(
  centroidweek ~ centroidcaterpillar_density,
  random = ~ 1 | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
))

r2(lme(
  centroidweek ~ centroidcaterpillar_density,
  random = ~ 1 | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
))



ggplot(data = HerbCatAnomaly_herbModel,
       aes(y= centroidweekAnomaly, x = centroidcaterpillar_propAnomaly)) + geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid herbivory anomaly", x= "Centroid anomaly for caterpillar proportion of occurence")


summary(
  lme( # This is a random slope-only model (no random intercept!)
    centroidweekAnomaly ~ centroidcaterpillar_propAnomaly,
    random = ~ 0 + centroidcaterpillar_propAnomaly | siteObserv,
    data = HerbCatAnomaly_herbModel,
    method = "REML",
    na.action = na.omit
  )
)

r2(lme( # This is a random slope-only model (no random intercept!)
  centroidweekAnomaly ~ centroidcaterpillar_propAnomaly,
  random = ~ 0 + centroidcaterpillar_propAnomaly | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
))

ggplot(data = HerbCatAnomaly_herbModel,
       aes(y= centroidweekAnomaly, x = centroidcaterpillar_densityAnomaly)) + geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid herbivory anomaly", 
                        x= "Centroid anomaly for caterpillar density")



summary(
  lme( # This is a random slope-only model (no random intercept!)
    centroidweekAnomaly ~ centroidcaterpillar_densityAnomaly,
    random = ~ 0 + centroidcaterpillar_densityAnomaly | siteObserv,
    data = HerbCatAnomaly_herbModel,
    method = "REML",
    na.action = na.omit
  )
)

r2(lme( # This is a random slope-only model (no random intercept!)
  centroidweekAnomaly ~ centroidcaterpillar_densityAnomaly,
  random = ~ 0 + centroidcaterpillar_densityAnomaly | siteObserv,
  data = HerbCatAnomaly_herbModel,
  method = "REML",
  na.action = na.omit
))






HerbCatAnomaly_herbModel.Lat = fullDatasetGreen %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(Latitude = mean(Latitude, na.rm = TRUE),
            Longitude = mean(Longitude, na.rm = TRUE)) %>% 
  right_join( HerbCatAnomaly_herbModel %>% extract(
    siteObserv,
    into = c("Name", "ObservationMethod"),
    regex = "(.*)_(.*)",
    remove = FALSE
  ), 
  by = c("Name" = "Name", "ObservationMethod"))






#  Herb_Rate ~ Temperature ----

# -- First, retrieve the good sites based on the those that are good for the actual phenometrics


tmp_file = tempfile(fileext = ".csv")

AnomalySites = HerbCatAnomaly_herbModel.Lat[, c("Name","Year", "Latitude", "Longitude")] %>% 
  filter(Longitude >= -100) %>% 
  group_by(Name, Latitude, Longitude) %>% 
  summarise(n = n()) %>% 
  select(-n) %>% as.data.frame()   

AnomalyDaymetr = AnomalySites %>% 
  rename(
    site = Name,
    lat = Latitude,
    lon = Longitude
  ) %>%
  write.csv(tmp_file, row.names = FALSE)


# pass temp CSV to function
TempAnomaly = download_daymet_batch(
  file_location = tmp_file,
  start = min(HerbCatAnomaly_herbModel.Lat$Year),
  end = 2024, # this is the most recent available in daymetr
  internal = TRUE
)

# remove temporary file 
unlink(tmp_file)


TempAnomalyData_clean <- lapply(TempAnomaly, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)
})

AllTempAnomaly= bind_rows(TempAnomalyData_clean)






HerbTempAnomaly = AllTempAnomaly %>% 
  filter(yday %in% TempDayWindow) %>% 
  group_by(site, year) %>% 
  summarise(meanTmin = mean(tmin..deg.c.), # average site level yearly minimum temperature
            meanTmax = mean(tmax..deg.c.),
            meanPreci = mean(prcp..mm.day.)) %>% 
  left_join(
    AllTempAnomaly %>% 
      filter(yday %in% TempDayWindow) %>% 
      group_by(site) %>% 
      summarise(AllmeanTmin = mean(tmin..deg.c.), # Across all years, what is the sites average minimum temperature (?)
                AllmeanTmax = mean(tmax..deg.c.),
                AllmeanPreci = mean(prcp..mm.day.)),
    by = c("site")) %>% 
  mutate(AnomalTmin = meanTmin - AllmeanTmin,
         AnomalTmax = meanTmax - AllmeanTmax,
         AnomalPreci = meanPreci - AllmeanPreci)%>% 
  inner_join(HerbCatAnomaly_herbModel.Lat %>%  filter(Year != "2025"), # because daymetr has no 2025 yet.
             by = c("site" = "Name", "year" = "Year")) %>% data.frame()


## Correlate important variables----

corrplot(cor(HerbTempAnomaly %>%
               select(maxHerb, effect, intercept, r2, centroidweek,
                      centroidcaterpillar_density, centroidcaterpillar_prop, maxcaterpillar_prop,
                      centroidweekAnomaly, centroidcaterpillar_propAnomaly, centroidcaterpillar_densityAnomaly,
                      AnomalTmin,  AnomalTmax, AnomalPreci, Latitude)%>%
               mutate(across(everything(), ~ifelse(is.infinite(.), NA, .))),
             use = "pairwise.complete.obs"),  
         method = "color", 
         type = "upper",          
         addCoef.col = "black",    
         tl.col = "black",         
         tl.srt = 80,              
         diag = FALSE)  







ggplot(data = HerbTempAnomaly,
       aes(y= centroidweekAnomaly, x = AnomalTmax)) + geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid herbivory anomaly", x= " Max.Temperature Anomaly")




ggplot(data = HerbTempAnomaly,
       aes(y= centroidcaterpillar_propAnomaly, x = AnomalTmax)) + geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(
      after_stat(eq.label),
      after_stat(rr.label),
      after_stat(p.value.label),
      sep = "~~~"
    )),
    parse = TRUE
  ) +theme_bw() + labs( y = "Centroid caterpillar proportion anomaly", 
                        x= " Max.Temperature Anomaly")


summary( # Theoretically, the intercept is supposed to be approximmately zero for all site 
  lme( # This is a random slope-only model (no random intercept!)
    centroidcaterpillar_propAnomaly ~ AnomalTmax,
    random = ~ 0 + AnomalTmax | siteObserv,
    data = HerbTempAnomaly,
    method = "REML",
    na.action = na.omit
  )
)

AnomalTmax.CentroidProp =  lme( # This is a random slope-only model (no random intercept!)
  centroidcaterpillar_propAnomaly ~ AnomalTmax,
  random = ~ 0 + AnomalTmax | siteObserv,
  data = HerbTempAnomaly,
  method = "REML",
  na.action = na.omit
)
r2(AnomalTmax.CentroidProp) # same R2 as a fixed effect model.


####################################################################################################################



library(ggh4x)

fullDataset %>% 
  filter(Year >= 2017) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear = n_distinct(Year), .groups = "drop") %>% 
  arrange(desc(nYear)) %>% 
  slice_head(n = 100) %>% 
  inner_join(
    fullDataset %>% 
      filter(Year >= 2017) %>% 
      group_by(Name, ObservationMethod) %>% 
      summarise(nSurv = n_distinct(ID), .groups = "drop") %>% 
      arrange(desc(nSurv)) %>% 
      slice_head(n = 100),
    
    by = c("Name", "ObservationMethod")
  ) %>% 
  inner_join(
    fullDataset %>% 
      filter(Year >= 2017) %>% 
      filter(!is.na(julianday)) %>% 
      mutate(
        julianMonth = case_when(
          julianday <= 30  ~ "Jan",
          julianday <= 60  ~ "Feb",
          julianday <= 90  ~ "Mar",
          julianday <= 120 ~ "Apr",
          julianday <= 150 ~ "May",
          julianday <= 180 ~ "Jun",
          julianday <= 210 ~ "Jul",
          julianday <= 240 ~ "Aug",
          julianday <= 270 ~ "Sep",
          julianday <= 300 ~ "Oct",
          julianday <= 330 ~ "Nov",
          julianday <= 366 ~ "Dec"
        )
      ),
    by = c("Name", "ObservationMethod")
  ) %>% 
  group_by(Name, Year, julianMonth) %>% 
  summarise(nSurv = n_distinct(ID), .groups = "drop") %>% 
  mutate(
    Name = factor(Name, levels = rev(sort(unique(Name)))),
    julianMonth = factor(julianMonth, levels = month.abb)
  ) %>% 
  ggplot(aes(x = julianMonth,
             y = Name,
             fill = nSurv)) +
  geom_tile(color = "white") +
  facet_grid(. ~ Year,
            # scales = "free_x", 
             space = "free_x") +
  scale_fill_viridis_c(name = "Survey intensity") +
  labs(x = " ", y = " ") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.spacing.x = unit(0, "lines")
  )


