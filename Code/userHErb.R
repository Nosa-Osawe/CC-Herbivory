
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



### Inclusion criteria ----

julianWindow = 120:230 # tweak this
min.nJulianWeekYearSite = 5
min.nSurvWeekYearSite = 20
min.nYear = 1
TempDayWindow = 120: 230
min.User.nSurv = 10

# proportion of herbivory per herb class----


# convert ordinal leaf damage scores into a weighted estimate of tissue loss, 
# aggregate them per site × week × year, and standardize by survey effort to obtain a continuous, 
# effort-corrected herbivory intensity suitable for phenological modeling

Herb1User = fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1)
  ) %>% 
  group_by(Name, Year, julianweek, UserFKOfObserver, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Year, julianweek, UserFKOfObserver, Herb) %>% 
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
      group_by(Name, Year, julianweek, UserFKOfObserver) %>% 
      summarise(nSurv = n_distinct(ID)), 
    by = c("Name", "Year", "julianweek", "UserFKOfObserver")) %>% 
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
         H4.prop = Herb_4) %>% 
  filter(nSurv >= min.User.nSurv)



herbUsingMeanUser= fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
  ) %>% 
  group_by(Name, Latitude, Year, julianweek, UserFKOfObserver, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Latitude, Year, julianweek, UserFKOfObserver, Herb) %>% 
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
  select(Name, Latitude, Year, julianweek, UserFKOfObserver, nSurv, totalHerb) %>% 
  as.data.frame() %>% 
  group_by(Name, Latitude, Year, julianweek, UserFKOfObserver) %>% 
  summarise(nSurv = sum(nSurv),
            totalHerb = sum(totalHerb)) %>% 
  mutate(totalHerbS = totalHerb/nSurv) %>% # totalHerbS is the standardized total heribory (to account for survey effort)
  filter(nSurv >= min.User.nSurv)

herbUsingMean2User= herbUsingMeanUser%>% 
  left_join(
    herbUsingMeanUser %>% 
      group_by(Name) %>% 
      summarise(nYear = n_distinct(Year)),
    by = c("Name")) %>%
  mutate(Name = reorder(Name, Latitude))

fullHerbUser = Herb1User %>% 
  left_join(herbUsingMean2User %>% select(-nSurv),
            by = c("Name", "Year", "julianweek", "UserFKOfObserver")) %>% 
  rename(nSurvHerb = nSurv)




 
 

#------------------------------------------------------------------------------------------

# A. Define what counts as herbivory (0-4). Already done in preparation of Herb.Arthropod dataframe
# 1. Define a good nSurv
# 2. Define a good JulianWindow
# 3. Define a good nJulianWeek
# 4. Define a good nYearHerb


min.nJulianWeekYearSite

goodUserSurv = fullHerbUser %>% 
  filter(nSurvHerb >= min.User.nSurv)

goodUserJulianWeek = goodUserSurv %>% 
  filter(julianweek %in% julianWindow) %>% 
  group_by(Name, Year, UserFKOfObserver) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite)

goodUserYear = goodUserJulianWeek %>% 
  left_join(fullHerbUser, 
            by = c("Name", "Year", "UserFKOfObserver")) %>% 
  filter(nYear >= min.nYear)  # this should be based on a recalculated number of nYears.


goodUserYear =   goodUserSurv %>% 
  inner_join(goodUserJulianWeek,  by = c("Name", "Year", "UserFKOfObserver")) %>% 
  group_by(Name,  julianweek, UserFKOfObserver) %>% 
  summarise(nYear = n_distinct(Year)) %>% filter(nYear >= min.nYear)
  


userGoodHerbData  = goodUserSurv %>% select(-nYear) %>% 
  inner_join(goodUserJulianWeek, by = c("Name", "Year", "UserFKOfObserver")) %>% 
  inner_join(goodUserYear, by = c("Name", "julianweek", "UserFKOfObserver"))

dim(userGoodHerbData)




 





# first make site and observation method one site.
userGoodHerbData2_collapsed = userGoodHerbData %>% 
  select(Name, Year, julianweek, UserFKOfObserver, totalHerbS) %>%
  group_by(Name, Year, julianweek, UserFKOfObserver) %>%
  summarise(totalHerbS = mean(totalHerbS), .groups = "drop") # use mean here due to repeated record


modelGoodHerbUser = userGoodHerbData2_collapsed %>% 
  mutate(
    julianweek_N = (julianweek - min(julianweek)) / (max(julianweek) - min(julianweek)),
    UserFKOfObserver = as.factor(UserFKOfObserver)
  )


modelGoodHerbUser_filtered = modelGoodHerbUser %>%
  right_join(
    modelGoodHerbUser %>% 
      group_by(Name, Year, julianweek) %>% 
      summarise(
        nUser = n_distinct(UserFKOfObserver),
        .groups = "drop"
      ) %>% 
      filter(nUser >= 2),
    by = c("Name", "Year", "julianweek")
  )

modelGoodHerbUser_filtered = modelGoodHerbUser %>%
  group_by(Name, Year) %>%
  filter(n_distinct(UserFKOfObserver) >= 2) %>%
  ungroup()


userFit_herb_model =function(df) {
  
  mod <- lm(totalHerbS ~ julianweek_N + UserFKOfObserver, data = df)
  sm  <- summary(mod)
  
  coef_ci <- confint(mod)
  
  tibble(
    intercept = coef(mod)[1],
    effect    = round(coef(mod)[2], 2),
    r2        = round(sm$r.squared, 2)
  )
}


userHerbModelOutput = modelGoodHerbUser_filtered %>%
  group_by(Name, Year) %>%
  group_modify(~ userFit_herb_model(.x)) %>% # passes each group's time series into the model
  ungroup()
 
 
userHerbModelOutput %>%
  filter(complete.cases(.)) %>% as.data.frame()

userHerbModelOutput.Lat = userHerbModelOutput %>%
  filter(complete.cases(.)) %>% 
  left_join(fullDataset %>% select(Name, Latitude, Longitude) %>% 
              group_by(Name) %>% 
              summarise(Latitude = mean(Latitude),
                        Longitude = mean(Longitude)),
            by = c("Name"))

summary(
  lm(effect ~ Latitude,
     data = userHerbModelOutput.Lat,
     weights = r2+0.0001)
)



ggplot(userHerbModelOutput.Lat, aes(x = Latitude, y = effect, weight = r2)) +
  geom_point(size = 3, aes(colour = r2), alpha = 0.9, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  ) +
  theme_bw() +
  labs( title = "Adjusted for Observer level variance",
    subtitle = "Mixed model showd R2= .64;marginal r2= 0.45; model : fixed = effect ~ Latitude + intercept,
+   random = ~ 1 | siteObserv, Est:1.56",
    y = "Rate of herbivory per growing season"
  )


userfit_lme_wt <- lme(
  fixed = effect ~ Latitude,
  random = ~ 1 | Name,
  weights = ~ (r2+ 0.001),  
  data = userHerbModelOutput.Lat,
  method = "REML"
)
summary(userfit_lme_wt)
r2(userfit_lme_wt)















