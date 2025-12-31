
rm(list = ls())

library(DescTools)
library(tidyverse)
library(jsonlite)
library(ineq)

library(lme4)
library(interactions)




# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))



julianWindow = 120:230 # tweak this
min.nJulianWeekYearSite = 5
min.nSurvWeekYearSite = 20
min.nYear = 1
minLat = 30
maxLat = 60  # Maybe not useful

HerbivoryGini = fullDataset %>% 
  filter( julianweek %in% julianWindow,
          Latitude > minLat,
    !HerbivoryScore %in% c(-128, -1)
  ) %>% 
  group_by(Name, Latitude, ObservationMethod, Year, julianweek, ID) %>% 
  summarise(HerbivoryScore = mean(HerbivoryScore)) %>% 
  mutate( 
    HerbActual = as.integer(case_when(
      HerbivoryScore == 0 ~ 0.001, # I use this instead of 0, because I am not sure gini should work well with 0.
      HerbivoryScore == 1 ~ 3,
      HerbivoryScore == 2 ~ 7,
      HerbivoryScore == 3 ~ 17,
      HerbivoryScore == 4 ~ 62))) %>% 
  left_join(
  fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1)
  ) %>% 
  group_by(Name, Latitude, ObservationMethod, Year, julianweek) %>% 
  summarise(
    nHerbSurv = n_distinct(ID),
    .groups = "drop"
  ), by = c("Name", "Latitude", "ObservationMethod", "Year", "julianweek")) %>% 
  filter(nHerbSurv >= min.nSurvWeekYearSite) %>% 
  group_by(Name, Latitude, ObservationMethod, Year, julianweek) %>% 
  summarise(nHerbSurv = mean(nHerbSurv),
            Gini.index = round(DescTools::Gini(HerbActual), 4),
            Lac = round (ineq::Lasym(HerbActual), 4)) %>% # Make adjustment here...
  mutate(TruncGini.index = ifelse(Gini.index < 0.01, 0.01 , Gini.index)) %>% 
  mutate(TruncGini.index = ifelse(TruncGini.index > 0.99, 0.99, TruncGini.index)) %>% 
  mutate(ilogitGini = qlogis(TruncGini.index)) %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod ))
 
summary(HerbivoryGini) # very few NAs produed in the Gini.index Calculation (13 of them, for now)

# Check: it seems that all herbscore is 0. So it may be safe to award a gini of 0 (absolute unism in herbivory level)
fullDataset %>% 
  filter( julianweek %in% julianWindow,
          !HerbivoryScore %in% c(-128, -1)
  ) %>% inner_join(
HerbivoryGini %>% filter(is.na(Gini.index)) %>% select(Name, Latitude, ObservationMethod, Year, julianweek, nHerbSurv),
by = c("Name", "Latitude", "ObservationMethod", "Year", "julianweek")) %>% 
  select(Name, ObservationMethod, Year, julianweek, HerbivoryScore, ID) %>%data.frame()


HerbivoryGini %>% filter(is.na(Lac))

HerbivoryGiniAdj = HerbivoryGini %>% 
  mutate(Gini.index = ifelse(is.na(Gini.index), 0, Gini.index),
         Lac = ifelse(is.na(Lac), 1, Lac)) %>% # change to 1 for same cause as the gini index NAs
  mutate(TruncGini.index = ifelse(Gini.index < 0.01, 0.01 , Gini.index)) %>% # can't ilogit a 0
  mutate(TruncGini.index = ifelse(TruncGini.index > 0.99, 0.99, TruncGini.index)) %>% 
  mutate(ilogitGini = qlogis(TruncGini.index)) %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod ))
  
summary(HerbivoryGiniAdj) # look good!

 
fitHerbGini =  lme(
  fixed = ilogitGini ~ Latitude,
  random = ~ 1 | siteObserv/Year,
  data = HerbivoryGiniAdj,
  method = "REML"
)

summary(fitHerbGini)
r2(fitHerbGini)

beta_fitHerbGini <- fixef(fitHerbGini)

ggplot(HerbivoryGiniAdj, aes(x = Latitude, y = ilogitGini, colour  = julianweek )) +
  geom_point(alpha = 0.3, size = 3) +
  geom_abline(
    intercept = beta_fitHerbGini["(Intercept)"],
    slope = beta_fitHerbGini["Latitude"],
    linewidth = 1.2,
    color = "red"
  ) +
  theme_classic() +
  labs(
    x = "Latitude",
    y = "ilogit Gini",
    title = "Across latituide, a few herbivivory score dominate the total herbivory",
    subtitle = "Marginal R2: 0.044"
  )






fitHerbLac_lmer <- lmer(
  Lac ~ ilogitGini * Latitude   + (1 | siteObserv/Year),
  data = HerbivoryGiniAdj,
  REML = TRUE
)

summary(fitHerbLac_lmer)
r2(fitHerbLac_lmer)

interact_plot(
  model = fitHerbLac_lmer,
  modx  = Latitude,
  pred  = ilogitGini,
  plot.points = FALSE,
  interval = FALSE,
  data = HerbivoryGiniAdj
)  




######################################################################################################33


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




# Inclusiuon criteria -----



julianWindow = 120:230 # tweak this
min.nJulianWeekYearSite = 5
min.nSurvWeekYearSite = 20
min.nYear = 1
TempDayWindow = 120: 230
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



goodHerbData2.gini =  goodHerbData2 %>%
  mutate(H_totalProp = rowSums(select(., H1.prop:H4.prop), na.rm = TRUE))

# health check
goodHerbData2.gini %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(hSum = sum(H_totalProp)) %>% 
  filter(hSum >1.001 | hSum < 0.999) # good!





herbGini = goodHerbData2.gini %>% 
  dplyr::select(Name, Latitude, ObservationMethod, Year, julianweek, Herb, H_totalProp) %>% 
  mutate(
    Herb_level = case_when(
      Herb == '0' ~ "H0",
      Herb == '1' ~ "H3",
      Herb == '2' ~ "H7",
      Herb == '3' ~ "H17",
      Herb == '4' ~ "H62"
    )
  ) %>% 
  pivot_wider(
    names_from = Herb_level,
    values_from = H_totalProp
  )%>% 
  mutate(
    across(
      .cols = c(H0, H3, H7, H17, H62),
      .fns  = ~ tidyr::replace_na(.x, 0)
    )
  ) %>%
  dplyr::select(-Herb) %>% 
  pivot_longer(
    cols = -c(Name, Latitude, ObservationMethod, Year, julianweek),
    names_to = "HerbClass",
    values_to = "Prop"
  ) %>% 
  group_by(Name, Latitude, ObservationMethod, Year, julianweek, HerbClass) %>% 
  summarise(Prop = sum(Prop)) %>% as.data.frame() %>% 
  mutate(HerbClassNum = as.integer(sub("H", "", HerbClass))) %>% 
  mutate(HerbClassNum = HerbClassNum / max(HerbClassNum)) %>%  # just to re-scale it, so 1 is the highest.
  group_by(Name, Latitude, ObservationMethod, Year, julianweek) %>% 
  summarise(
    Gini.index = DescTools::Gini(Prop),
    LAC = Lasym(Prop),
    Gini.index2 = DescTools::Gini(
      x = HerbClassNum,
      weights = round(100 * Prop) # will produce NaNs if only one weight is present and others are 0. 
    ), # This is same as a gini index of 1, so, to solve the NaN problem,
    .groups = "drop" #  I will replace all NaNs by 0.99 (so it can be ilogit transformed)
  ) %>% 
  as.data.frame() %>% 
  mutate(Gini.index2 = ifelse(is.nan(Gini.index2), 0.99, Gini.index2)) %>% 
  mutate(truncGini = ifelse((Gini.index > 0.99), 0.99, Gini.index)) %>% 
  mutate(ilogitGini.index = qlogis(truncGini),
         ilogitGini.index2 = qlogis(Gini.index2)) %>% # we compute an ilogit from 1
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod ))

summary(herbGini)



fit_ilogitGini1 =  lme(
  fixed = ilogitGini.index ~ Latitude,
  random = ~ 1 | siteObserv/Year,
  data = herbGini,
  method = "REML"
)

summary(fit_ilogitGini1)
r2(fit_ilogitGini1)

beta_fit_ilogitGini1 <- fixef(fit_ilogitGini1)

ggplot(herbGini, aes(x = Latitude, y = ilogitGini.index)) +
  geom_point(alpha = 0.3, size = 3) +
  geom_abline(
    intercept = beta_fit_ilogitGini1["(Intercept)"],
    slope = beta_fit_ilogitGini1["Latitude"],
    linewidth = 1.2,
    color = "red"
  ) +
  theme_classic() +
  labs(
    x = "Latitude",
    y = "ilogit Gini",
    title = "Variation in Herbivory (level) increases with latitude",
    subtitle = "Marginal R2: 0.044"
  )




# Julian week as predictor of gini variation in herbiivory


fit_giniWeek =   lme(
  fixed = ilogitGini.index ~ julianweek + Latitude,
  random = ~ 1 | siteObserv,
  data = herbGini2,
  method = "REML"
)

summary(fit_giniWeek)
r2(fit_giniWeek)

