 
rm(list = ls())
# Load libraries ----

library(tidyverse)
library(jsonlite)
library(daymetr)
library(ggpmisc)
library(nlme)
library(performance)
library(DescTools)
library(ineq)
library(jtools)
library(corrplot)

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



### Herbivory + Arthropod data----
Herb.Arthropod = fullHerb %>% select(-Latitude) %>% 
  left_join(prop.dens.arthropod, 
            by = c("Name", "Year", "julianweek"),
            relationship = "many-to-many") %>% # due to multiple observation methods in one site
  data.frame() %>% rename(nYearHerb = nYear) %>% filter(!is.na(ObservationMethod))




# Calculate herbivory anomaly-----

### Inclusion criteria ----

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

 
 




#######################################################################################################

# Fit linear regression models for each julian week  as predictor for total herbivory

# first make site and observation method one site.
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
      filter(nJulianweek >= min.nJulianWeekYearSite), # Earlier filtering (might) have solved this problem. See that dim (modelGoodHerb) keeps the row length the same
    by = c("siteObserv", "Year")
  ) %>% 
  as.data.frame() %>% 
  mutate(
    julianweek_C = scale(julianweek), # Z-score transformation
    julianweek_N = (julianweek - min(julianweek)) / (max(julianweek) - min(julianweek)) # Min- max normalization
  )

# The idea behind mutating with the Min-Max normalized variable is that I can calculate total herbivory per unit increase over the 
# growing season. Note that the unit here is growing season (not per julian week)
# The z transformed variables is not of real use for now. But I think of its use if I am to switch to carrying out the analysis 
# using the brms package (for Bayesian statistics)


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
  separate(siteObserv, into = c("Name", "ObservationMethod"),
           sep = "_", remove = FALSE ) %>% 
  left_join(fullDataset %>% select(Name, Latitude, Longitude) %>% 
              group_by(Name) %>% 
              summarise(Latitude = mean(Latitude),
                        Longitude = mean(Longitude)),
            by = c("Name"))

 

summary(
  lm(effect ~ Latitude + intercept,
     data = herbModelOutput.Lat,
     weights = r2+0.0001)
)

 
 

ggplot(herbModelOutput.Lat, aes(x = Latitude, y = effect, weight = (100*(r2+ 0.001)))) +
  geom_point(size = 3, aes(colour = r2), alpha = 0.9, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  ) +
  theme_bw() +
  labs(
    y = "Rate of herbivory per growing season"
  )


# This one here does not make sense, so I comment it out
# fit_lme_wt <- lme(
#   fixed = effect ~ Latitude ,
#   random = ~ 1 | siteObserv,
#   weights = varFixed(~ (1-r2 + 0.001)), # Weigh by % of variance explained.
#   data = herbModelOutput.Lat,
#   method = "REML"
# )
# summary(fit_lme_wt)
# r2(fit_lme_wt)
# 





herbModelOutput.Lat



#  Herb_Rate ~ Temperature ----

# -- First, retrieve the good sites based on the those that are good for the actual phenometrics


tmp_file = tempfile(fileext = ".csv")

HerbSites = herbModelOutput.Lat[, c("Name","Year", "Latitude", "Longitude")] %>% 
  filter(Longitude >= -100) %>% 
  group_by(Name, Latitude, Longitude) %>% 
  summarise(n = n()) %>% 
  select(-n) %>% as.data.frame()   

HerbSitesDaymetr = HerbSites %>% 
  rename(
    site = Name,
    lat = Latitude,
    lon = Longitude
  ) %>%
  write.csv(tmp_file, row.names = FALSE)


# pass temp CSV to function
TempSiteData = download_daymet_batch(
  file_location = tmp_file,
  start = min(herbModelOutput.Lat$Year),
  end = 2024, # this is the most recent available in daymetr
  internal = TRUE
)

# remove temporary file 
unlink(tmp_file)


TempSiteData_clean <- lapply(TempSiteData, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)
})

AllTempData = bind_rows(TempSiteData_clean)



HerbTemp = AllTempData %>% 
  filter(yday %in% TempDayWindow) %>% 
  group_by(site, year) %>% 
  summarise(meanTmin = mean(tmin..deg.c.), # average site level yearly minimum temperature
            meanTmax = mean(tmax..deg.c.),
            meanPreci = mean(prcp..mm.day.)) %>% 
  left_join(
    AllTempData %>% 
      filter(yday %in% TempDayWindow) %>% 
      group_by(site) %>% 
      summarise(AllmeanTmin = mean(tmin..deg.c.), # Across all years, what is the sites average minimum temperature (?)
                AllmeanTmax = mean(tmax..deg.c.),
                AllmeanPreci = mean(prcp..mm.day.)),
    by = c("site")) %>% 
  mutate(AnomalTmin = meanTmin - AllmeanTmin,
         AnomalTmax = meanTmax - AllmeanTmax,
         AnomalPreci = meanPreci - AllmeanPreci)%>% 
  inner_join(herbModelOutput.Lat %>%  filter(Year != "2025"), # because daymetr has no 2025 yet.
             by = c("site" = "Name", "year" = "Year")) %>% data.frame()
 

HerbTemp_W = HerbTemp %>% # This is the dataframe I would be working with to find the best predictor for "effect"
  select(
    siteObserv,  Longitude, year, 
    r2, intercept, effect, 
    Latitude, meanTmin, meanTmax, meanPreci,
    AnomalTmin, AnomalTmax, AnomalPreci 
    )

summary(HerbTemp_W)

  # Remove siteObserv and keep only numeric variables
  corrplot(cor(HerbTemp_W %>%
    select(-siteObserv) %>%        
    select(where(is.numeric)),     
    use = "pairwise.complete.obs"),  
  method = "color", 
  type = "upper",          
  addCoef.col = "black",    
  tl.col = "black",         
  tl.srt = 45,              
  diag = FALSE)    
 


fit_temp<- lme(
  fixed = effect ~ meanTmin ,
  random = ~ 1 | siteObserv,
  weights = ~ (1-r2+ 0.001),  
  data = HerbTemp_W,
  method = "REML"
)
summary(fit_temp)
r2(fit_temp)


 
## Anomaly on rate of herbivory increase ----


HerbRateAnomaly = herbModelOutput.Lat %>% 
  group_by(siteObserv) %>% 
  summarise(meanEffect = mean(effect, na.rm = TRUE),
            meanIntercept = mean(intercept, na.rm = TRUE),
            meanR2 = mean(r2, na.rm = TRUE)) %>% data.frame() %>% 
  right_join(herbModelOutput.Lat, by = "siteObserv") %>% 
  mutate(effectAnomal = effect - meanEffect,
         interceptAnomal = intercept - meanIntercept,
         r2Anomal = r2 - meanR2) %>% as.data.frame() %>% 
  select(-c(intercept, Name, ObservationMethod, intercept, effect, r2, 
            intercept_lwr, intercept_upr, effect_lwr, effect_upr, Latitude, Longitude)) %>% 
  inner_join(
  HerbTemp %>% select(- c( intercept_lwr, intercept_upr, effect_lwr, effect_upr)), 
  by = c("siteObserv", "Year" = "year")) %>% as.data.frame()
 

corrplot(cor(HerbRateAnomaly %>%
               select(where(is.numeric)),     
             use = "pairwise.complete.obs"),  
         method = "color", 
         type = "upper",          
         addCoef.col = "black",    
         tl.col = "black",         
         tl.srt = 45,              
         diag = FALSE)  




ggplot(HerbRateAnomaly, aes(x = AnomalTmin, y = effectAnomal))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )


ggplot(HerbRateAnomaly, aes(x = AnomalTmax, y = effectAnomal))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )



ggplot(HerbRateAnomaly, aes(x = AnomalTmin, y = r2Anomal))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )

ggplot(HerbRateAnomaly, aes(x = AnomalTmax, y = r2Anomal))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )


ggplot(HerbRateAnomaly, aes(x = AnomalPreci, y = effectAnomal))+
  geom_point()+
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
  )















