
# Script is useful to check the herbivory, arthropod density and occurrence at a community (site level)
# see the plantGenusHerbArthopod.R script for plantGenus - site level analysis/EDA. 
# Do not forget to source the herbivoryFunction.R
# This script, as well as others in this CC-Herbivory should run independent of the others. 
# So, it is recommended that you do not forget to clear the environment every time you want to run a new script, as object names are shared and may mean different things for different scripts. "rm(list = ls())" Does this well.


# Clear R's brain
rm(list = ls())


# Load libraries ----
library(tidyverse)
library(jsonlite)

# CC data ----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))


# INTERESTING CC SITES.... (but feel free to check others!)
# [1] "Acadia NP - Alder"                    "Acadia NP - Sundew"                   "BHRC - Hemlock Draw"                 
# [4] "BHRC Pan Hollow Uplands"              "Beare Swamp in Rouge Park"            "Belle Isle State Park, Lancaster, VA"
# [7] "Birmingham Zoo AL"                    "Brookbanks Park"                      "Highlands Biological Station"        
# [10] "Middle Mill"                          "NC Botanical Garden"                  "Prairie Ridge Ecostation"            
# [13] "RVCC"                                 "Riverbend Park"                       "Stage Nature Center"                 
# [16] "UNC Chapel Hill Campus"               "Virginia Zoo"                         "Walker Nature Center"    

source("Code/herbivoryFunctions.R")

unique(fullDataset$Name) # CC cites!



h0 = 0
h1 = 3
h2 = 7
h3 = 17
h4 = 62

julianWindow = 120:300       # julainWindow (what i use typically is : 120:230)
min.nJulianWeekYearSite = 6   # minimum number of julianWeeks for each site-year
min.nSurvWeekYearSite = 20   # minimum number of site-year-week surveys
min.nYear = 3        # minimum number of survey site-year
TempDayWindow = 120: 230   # Temperature window

# get_site is a custom function
site = get_site( "Acadia NP - Alder"  , fullDataset) # 'All' selects all



fullDataset2 = fullDataset %>% 
  filter(Name %in% site)



Herb1 = fullDataset2 %>% 
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
    values_from = Herbcount,
    values_fill = 0 
  ) %>% 
  left_join(
    fullDataset2 %>%   
      filter(
        !HerbivoryScore %in% c(-128, -1)
      ) %>% 
      group_by(Name, Year, julianweek) %>% 
      summarise(nSurv = n_distinct(ID)), by = c("Name", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  mutate(
    across(starts_with("Herb_"), ~ .x / nSurv)
  ) %>% 
  rename_with(
    ~ sub("Herb_", "H.prop", .x),
    starts_with("Herb_")
  )


# Herb score using the mean of the scale ----

# here herbUsingMean is just using the mean/median of the herbivory scale to standardise the data

herb_weights = c(
  "0" = h0,
  "1" = h1,
  "2" = h2,
  "3" = h3,
  "4" = h4
)

herbUsingMean = fullDataset2 %>% 
  filter(!HerbivoryScore %in% c(-128, -1)) %>% 
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
  select(Name, Latitude, Year, julianweek, nSurv, totalHerb) %>% 
  group_by(Name, Latitude, Year, julianweek) %>% 
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
      group_by(Name) %>% 
      summarise(nYear = n_distinct(Year)),
    by = c("Name")) %>%
  mutate(Name = reorder(Name, Latitude))

fullHerb = Herb1 %>% 
  left_join(herbUsingMean2 %>% select(-nSurv),
            by = c("Name", "Year", "julianweek")) %>% 
  rename(nSurvHerb = nSurv)



# Arthropod data: prop of surv per week----
prop_fullDataset2 = fullDataset2 %>%
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

dens_fullDataset2 = fullDataset2 %>%
  filter(WetLeaves == 0) %>% 
  group_by(Name, ObservationMethod, Latitude, Longitude, Year, julianweek) %>%
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
            by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  as.data.frame() %>% 
  rename(nSurvArthropod = nSurv) %>% 
  left_join( 
    dens_fullDataset2 %>% 
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




Herb.Arthropod = Herb.Arthropod %>% 
  filter(julianweek %in% julianWindow) %>% 
  select(Name, ObservationMethod, Year, julianweek,
         totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(
    totalHerbS = mean(totalHerbS, na.rm = TRUE),
    caterpillar_prop = mean(caterpillar_prop, na.rm = TRUE) * 100,
    caterpillar_density = mean(caterpillar_density, na.rm = TRUE),
    .groups = "drop"
  ) %>% left_join(
    Herb.Arthropod %>% 
      filter(julianweek %in% julianWindow) %>% 
      select(Name, ObservationMethod, Year, julianweek,
             totalHerbS, caterpillar_prop, caterpillar_density) %>% 
      group_by(Name, ObservationMethod, Year, julianweek) %>% 
      summarise(
        totalHerbS = mean(totalHerbS, na.rm = TRUE),
        caterpillar_prop = mean(caterpillar_prop, na.rm = TRUE) * 100,
        caterpillar_density = mean(caterpillar_density, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      group_by(Name, ObservationMethod, Year) %>% 
      summarise(maxtotalHerbS = max(totalHerbS ),
                maxcaterpillar_prop = max(caterpillar_prop),
                maxcaterpillar_density = max(caterpillar_density)), by = c("Name", "ObservationMethod", "Year")
  ) %>% 
  mutate(caterpillar_prop.min.max = caterpillar_prop/ maxcaterpillar_prop,
         totalHerbS.min.max = totalHerbS / maxtotalHerbS,
         caterpillar_density.min.max = caterpillar_density/ maxcaterpillar_density) %>% 
  as.data.frame()


# I found that this filtering sometimes eliminates some weeks of survery where nSurv is below threshold. So that the highest number
# max(X/)Xmax is sometimes != 1.  For now, I think its fine that way. We can keep in mind that this is only for viz purpose.
# Data wrangling pipeline:----



# A. Define what counts as herbivory (0-4). Already done in preparation of Herb.Arthropod dataframe
# 1. Define a good nSurv
# 2. Define a good JulianWindow
# 3. Define a good nJulianWeek
# 4. Define a good nYearHerb



goodnSurv2 = fullDataset2 %>% 
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


# anti_join(goodHerbData, goodHerbData2, 
#           by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
#   select(Name, ObservationMethod, Year, julianweek, nYearHerb)
# 



#######################################################################################################

#  Julian week as a predictor of totalherbivory ----
# Fit linear regression models for each julian week  as predictor for total herbivory



goodHerbData2_collapsed = goodHerbData2 %>% 
  mutate(siteObserv = paste0(Name, sep = "_", ObservationMethod)) %>% 
  group_by(siteObserv, Name, ObservationMethod, Year, julianweek) %>%
  summarise(totalHerbS = mean(totalHerbS),
            caterpillar_prop  = mean(caterpillar_prop),
            caterpillar_density = mean(caterpillar_density),
            maxtotalHerbS = mean(maxtotalHerbS),
            maxcaterpillar_prop = mean(maxcaterpillar_prop),
            maxcaterpillar_density = mean(maxcaterpillar_density),
            caterpillar_prop.min.max = mean(caterpillar_prop.min.max),
            totalHerbS.min.max = mean(totalHerbS.min.max),
            caterpillar_density.min.max = mean(caterpillar_density.min.max),
            nYearHerb  = mean(nYearHerb),
            nJulianWeek = mean(nJulianWeek),
            nSurv = mean(nSurv),
            .groups = "drop") # use mean here due to repeated record




goodHerbData2_collapsed %>% 
  pivot_longer(
    cols = c(totalHerbS.min.max, caterpillar_prop.min.max, caterpillar_density.min.max ),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = julianweek, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    ObservationMethod ~ Year,
    scales = "free_x"
  ) +
  scale_color_manual(
    values = c( "orange", "#d95f02", "darkgreen"
    )
  ) +
  labs(
    x = "Julian week",
    y = "Min-max normalization of the herbivory and proportion of occurence (%)",
    color = "",
    title =  "",              
    subtitle = site
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

