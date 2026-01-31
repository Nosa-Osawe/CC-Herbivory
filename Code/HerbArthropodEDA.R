
rm(list = ls())


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



# This bunnch of filters would be used 


# Load functions

source("Code/herbivoryFunctions.R")

unique(fullDataset$Name)
unique(fullDataset$plantGenus)


# Herbivory scoring ----

# Change as you find logically convincing

h0 = 0
h1 = 3
h2 = 7
h3 = 17
h4 = 62

julianWindow = 120:230 # tweak this
min.nJulianWeekYearSite = 5
min.nSurvWeekYearSite = 5
min.nYear = 1
site = get_site("All", fullDataset) # 'All' selects all
plantGenus = get_genera("All", fullDataset) 


fullDataset2 = fullDataset %>% 
  filter(plantGenus %in% plantGenus,
          Name %in% site)


# proportion of herbivory per herb class----


# convert ordinal leaf damage scores into a weighted estimate of tissue loss, 
# aggregate them per site × week × year, and standardize by survey effort to obtain a continuous, 
# effort-corrected herbivory intensity suitable for phenological modeling

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
    values_from = Herbcount
  ) %>% 
  left_join(
    fullDataset2 %>% 
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

herbUsingMean= fullDataset2 %>% 
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
    H0 = `0` * h0,
    H1 = `1` * h1,
    H2 = `2` * h2,
    H3 = `3` * h3,
    H4 = `4` * h4  # this is subject to change.
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



Herb.ArthropodMinMax = Herb.Arthropod %>% 
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


# A. Define what counts as herbivory (0-4). Already done in preparation of Herb.Arthropod dataframe
# 1. Define a good nSurv
# 2. Define a good JulianWindow
# 3. Define a good nJulianWeek
# 4. Define a good nYearHerb



goodnSurv2MinMax = fullDataset2 %>% 
  right_join(Herb.ArthropodMinMax  %>% filter(julianweek %in% julianWindow),
             by = c("Name",  "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name,  ObservationMethod, Year, julianweek) %>% 
  summarise(nSurv = n_distinct(ID)) %>% 
  filter(nSurv >= min.nSurvWeekYearSite) %>% data.frame()


goodJulianWeek2MinMax = Herb.ArthropodMinMax %>% 
  right_join(goodnSurv2MinMax, by = c("Name", "ObservationMethod", "Year", "julianweek")) %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite)



goodnYear2MinMax = Herb.ArthropodMinMax %>% 
  right_join(goodJulianWeek2MinMax, by = c("Name", "ObservationMethod", "Year")) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYearHerb = n_distinct(Year)) %>% 
  filter(nYearHerb >= min.nYear)



goodHerbData2MinMax = Herb.ArthropodMinMax %>% 
  inner_join(goodnYear2MinMax, by  = c("Name", "ObservationMethod")) %>% # keep good years
  inner_join(goodJulianWeek2MinMax, by = c("Name", "ObservationMethod", "Year")) %>% # keep good weeks
  inner_join(goodnSurv2MinMax, by = c("Name", "ObservationMethod", "Year", "julianweek")) # Keep good surveys


# added siteObserv just in case analysis need to be done at the site-observ level.
goodHerbData2_collapsed_MinMax = goodHerbData2MinMax %>% 
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




goodHerbData2_collapsed_MinMax %>% 
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
    title = plantGenus,              
    subtitle = site
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )




