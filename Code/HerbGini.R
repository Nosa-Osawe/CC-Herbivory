
rm(list = ls())

library(DescTools)
library(tidyverse)
library(jsonlite)
library(ineq)

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

ggplot(HerbivoryGiniAdj, aes(x = Latitude, y = ilogitGini)) +
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




 
library(lme4)

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
