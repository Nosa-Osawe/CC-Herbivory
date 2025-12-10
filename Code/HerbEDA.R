library(tidyverse)
library(jsonlite)
library(daymetr)

### 1. Read in CC! data files
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))

#  Herbivory EDA ----



fullDataset %>% 
  filter(!HerbivoryScore %in% c(-128, -1)) %>% 
  group_by(Name, HerbivoryScore) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% data.frame()


fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden"),
    Year == 2024
  ) %>% 
  group_by(Name, julianweek, HerbivoryScore) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(y = count, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_grid(Name ~ HerbivoryScore)+
  theme_bw()



fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden"),
  ) %>% 
  group_by(Name, julianweek, HerbivoryScore) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(y = count, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_grid(Name ~ HerbivoryScore)+
  theme_bw()

fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden|Eno River State Park"),
    Year >= 2021
  ) %>% 
  group_by(Name, Year, julianweek, HerbivoryScore) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(x = julianweek, y = count, color = as.factor(Year), group = Year)) +
  geom_point(alpha = 0.6, size  = 1) +
  geom_line(size = 0.5) + 
  facet_grid(Name ~ HerbivoryScore) +
  theme_bw() +
  labs(color = "Year")



fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge"),
    Year >= 2024
  ) %>% 
  group_by(Name, Code, julianweek) %>% 
  summarise(nSurv = n_distinct(ID),
            herb = mean(HerbivoryScore)) %>%
  ggplot(aes(y = herb, x = julianweek)) +
  geom_point()+
  facet_wrap( ~ Code)








fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, regex("prairie", ignore_case = TRUE )),
    Year == 2017
  ) %>% 
  group_by(PlantSpecies, julianweek) %>% 
  summarise(HerbivoryScore = mean(HerbivoryScore)) %>% 
  ggplot(aes(y = HerbivoryScore, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_wrap(~PlantSpecies)



fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, 
               "Prairie Ridge|NC Botanical Garden|UNC Chapel Hill Campus|Eno River State Park|Acadia NP - Sundew|Acadia NP - Alder"),
    Year >= 2021
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb  = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
  summarise(Herbcount = n(),
            nSurv = n_distinct(ID),
            .groups = "drop") %>% 
  ggplot(aes(x = julianweek, y = Herbcount, color = as.factor(Year), group = Year)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", 
                                "#009E73", "#F0E442", "darkgrey"))+
  facet_grid(Name ~ as.factor(Herb), scales = "free_y") +
  theme_bw() +
  labs(color = "Year")





herb= fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden|UNC Chapel Hill Campus|Eno River State Park|Acadia NP - Sundew|Acadia NP - Alder"),
    Year >= 2021
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
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
    H1 = `1` * 1,
    H2 = `2` * 2,
    H3 = `3` * 3,
    H4 = `4` * 4
  ) %>% 
  mutate(
    totalHerb = rowSums(across(H0:H4), na.rm = TRUE)
  ) %>% 
  select(Name, Year, julianweek, nSurv, totalHerb) %>% 
  as.data.frame() %>% 
  group_by(Name, Year, julianweek) %>% 
  summarise(nSurv = sum(nSurv),
            totalHerb = sum(totalHerb)) %>% 
  mutate(Herb_standardized = totalHerb/nSurv)




cor(herb$nSurv, herb$totalHerb)

ggplot(data = herb, 
       aes(x = julianweek, y = totalHerb, color = as.factor(Year))) +
  geom_point(aes(size = nSurv), alpha = 0.7) +
  scale_color_manual(values =c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", size = "nSurv")


ggplot(data = herb %>% filter(nSurv >=6), 
       aes(x = julianweek, y = Herb_standardized, color = as.factor(Year))) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values =c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(color = "Year")



ggplot(data = herb %>% filter(nSurv >= 6), 
       aes(x = julianweek, 
           y = Herb_standardized, 
           color = as.factor(Year))) +
  geom_point(alpha = 0.7) +
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   # controls wiggliness
    se = TRUE,
    alpha = 0.15                 # transparency of SE ribbon
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", fill = "Year")







ggplot(data = herb , 
       aes(x = julianweek, y = totalHerb, color = as.factor(Year))) +
  stat_smooth(aes(fill = as.factor(Year)), se = FALSE, method = "loess", alpha = 0.1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(~ Name, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", fill = "Year")




ggplot(data =  herb %>% filter(nSurv >=6),
       aes(x = julianweek, y = Herb_standardized, color = as.factor(Year))) +
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,   
    alpha = 0.1
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 3, scales = "free_y") +  # <- 3 columns, rows determined automatically
  theme_bw() +
  labs(color = "Year", fill = "Year")

fullDataset %>% 
  filter(Name %in% c("Acadia NP - Alder", "Acadia NP - Sundew", "NC Botanical Garden", "Prairie Ridge Ecostation", "UNC Chapel Hill Campus"),
         Year>2020) %>% 
  group_by(Name, PlantSpecies) %>% 
  summarise(count = n_distinct(ID)) %>% data.frame()








ggplot(data= herb, aes(x = julianweek, y = count, color = as.factor(Year), group = Year)) +
  geom_point(alpha = 0.6, size  = 1) +
  geom_line(size = 0.5) + 
  facet_grid(Name ~ HerbivoryScore) +
  theme_bw() +
  labs(color = "Year")




Herbivory= fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb = mean(HerbivoryScore), .groups = "drop") %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
  summarise(
    Herbcount = n_distinct(ID),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Herb,
    values_from = Herbcount
  ) %>% 
  mutate(
    H0 = `0` * 0,
    H1 = `1` * 1,
    H2 = `2` * 2,
    H3 = `3` * 3,
    H4 = `4` * 4
  ) %>% 
  mutate(
    totalHerb = rowSums(across(H0:H4), na.rm = TRUE)
  ) %>% 
  select(Name, Year, julianweek, Herbcount, totalHerb) %>% 
  as.data.frame() %>% 
  group_by(Name, Year, julianweek) %>% 
  summarise(nSurv = sum(nSurv),
            totalHerb = sum(totalHerb)) %>% 
  mutate(Herb_standardized = totalHerb/nSurv) %>% 
  rename(TotalnSurv = nSurv)



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
         Herb_4 = Herb_4/ nSurv)



Herb1 %>% 
  filter(Name %in% c(
    "Acadia NP - Alder", 
    "Acadia NP - Sundew",
    "NC Botanical Garden", 
    "Prairie Ridge Ecostation",
    "UNC Chapel Hill Campus"
  ),
  Year>2020,
  nSurv>5) %>% 
  pivot_longer(cols = -c(Name, Year, julianweek, Herb, nSurv),
               names_to = "HerbCat", 
               values_to = "HerbStand") %>% 
  ggplot(aes(x = julianweek, y = HerbStand, color = as.factor(Year), group = Year)) +
  geom_point()+
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   # controls wiggliness
    se = TRUE,
    alpha = 0.15                 # transparency of SE ribbon
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(Name ~ as.factor(HerbCat), scales = "free_y") +
  theme_bw() +
  guides(fill = "none")+
  labs(color = "Year",
       x= "Day", y  = "Proportion of herbivory Category")



Herb1 %>%  
  filter(
    Name %in% (
      Herb1 %>% 
        group_by(Name) %>% 
        summarise(meanH1 = mean(Herb_1, na.rm = TRUE)) %>% 
        arrange(desc(meanH1)) %>% 
        slice_head(n = 8) %>% 
        pull(Name)
    ),
    Year>2020,
    nSurv>5) %>% 
  pivot_longer(cols = -c(Name, Year, julianweek, Herb, nSurv),
               names_to = "HerbCat", 
               values_to = "HerbStand") %>% 
  ggplot(aes(x = julianweek, y = HerbStand, color = as.factor(Year), group = Year)) +
  geom_point()+
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),    
    se = TRUE,
    alpha = 0.15                  
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(Name ~ as.factor(HerbCat), scales = "free_y") +
  theme_bw() +
  guides(fill = "none")+
  labs(color = "Year",
       x= "Day", y  = "Proportion of herbivory Category")

Herb1 %>%
  group_by(Name) %>% 
  summarise(nSurv = sum(nSurv)) %>% arrange(desc(nSurv))




Herb1 %>%  
  filter(
    Name %in% (
      Herb1 %>%
        group_by(Name) %>% 
        summarise(nSurv = sum(nSurv)) %>% arrange(desc(nSurv)) %>% 
        slice_head(n = 12) %>% 
        pull(Name)
    ),
    Year>2020,
    nSurv>5) %>% 
  filter(!Name %in% c("Acadia NP - Sundew", "Acadia NP - Alder", "Stage Nature Center")) %>% 
  pivot_longer(cols = -c(Name, Year, julianweek, Herb, nSurv),
               names_to = "HerbCat", 
               values_to = "HerbStand") %>% 
  ggplot(aes(x = julianweek, y = HerbStand, color = as.factor(Year), group = Year)) +
  geom_point()+
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   
    se = TRUE,
    alpha = 0.15                  
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(Name ~ as.factor(HerbCat), scales = "free_y") +
  theme_bw() +
  guides(fill = "none")+
  labs(color = "Year",
       x= "Day", y  = "Proportion of herbivory Category")








Herb1 %>%  
  filter(
    Name %in% (
      Herb1 %>%
        group_by(Name) %>% 
        summarise(nSurv = sum(nSurv)) %>% arrange(desc(nSurv)) %>% 
        slice_head(n = 6) %>% 
        pull(Name)
    ),
    Year>2020,
    nSurv>5,
    julianweek>= 140 & julianweek <= 200) %>% 
  filter(!Name %in% c("Acadia NP - Sundew", "Acadia NP - Alder", "Stage Nature Center")) %>% 
  pivot_longer(cols = -c(Name, Year, julianweek, Herb, nSurv),
               names_to = "HerbCat", 
               values_to = "HerbStand") %>% 
  mutate(Year = as.factor(Year)) %>% 
  ggplot(aes(y = HerbStand, x = Year)) +
  geom_boxplot(aes(fill = Year))+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(Name ~ as.factor(HerbCat), scales = "free_y") +
  theme_bw() +
  guides(fill = "none")+
  labs( x= "Year", y  = "Proportion of herbivory Category")




fullDataset %>% 
  group_by(PlantSpecies) %>% 
  summarise(nSurv = n_distinct(ID)) %>% 
  arrange(desc(nSurv)) %>% 
  print(n= 10)




RedMaple = fullDataset %>% 
  filter(
    PlantSpecies=="Red maple",
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
         Herb_4 = Herb_4/ nSurv)



RedMaple %>%  
  filter(
    Name %in% (
      RedMaple %>%
        group_by(Name) %>% 
        summarise(nSurv = sum(nSurv)) %>% arrange(desc(nSurv)) %>% 
        slice_head(n = 6) %>% 
        pull(Name)
    ),
    Year>2020,
    nSurv>5,
    julianweek>= 140 & julianweek <= 200) %>% 
  pivot_longer(cols = -c(Name, Year, julianweek, Herb, nSurv),
               names_to = "HerbCat", 
               values_to = "HerbStand") %>% 
  ggplot(aes(x = julianweek, y = HerbStand, color = as.factor(Year), group = Year)) +
  geom_point()+
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   
    se = TRUE,
    alpha = 0.15                  
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_grid(Name ~ as.factor(HerbCat), scales = "free_y") +
  theme_bw() +
  guides(fill = "none")+
  labs(color = "Year",
       x= "Day", y  = "Proportion of herbivory Category",
       subtitle = "Red maple")


# Total herbivory ----

# Herbivory scores are as follows: 0 - none, 1 - trace (<5%), 2 - light (5-10%), 3 – moderate (10-25%), 4 - heavy (>25%).

# using their mean
# none = 0, Trace = 3, light = 7, Moderate = 17, Heavy = 62


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
  mutate(Herb_standardized = totalHerb/nSurv) %>% 
  left_join(
    herbUsingMean %>% 
      group_by(Name) %>% 
      summarise(nYear = n_distinct(Year)),
    by = c("Name")) %>%
  mutate(Name = reorder(Name, Latitude))




herbUsingMean %>% 
  filter(Name %in% c(
    "Acadia NP - Alder", 
    "Acadia NP - Sundew",
    "NC Botanical Garden", 
    "Prairie Ridge Ecostation",
    "UNC Chapel Hill Campus",
    "Brookbanks Park"), 
    Year>2020, nSurv >= 6) %>% 
  ggplot(aes(x = julianweek, 
             y = Herb_standardized, 
             color = as.factor(Year))) +
  geom_point(alpha = 0.7) +
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   # controls wiggliness
    se = TRUE,
    alpha = 0.15                 # transparency of SE ribbon
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", fill = "Year")




herbUsingMean %>% 
  filter( 
    nYear >= 4,
    Year>2020, 
    nSurv >= 6,
  ) %>% 
  ggplot(aes(x = julianweek, 
             y = Herb_standardized, 
             color = as.factor(Year))) +
  geom_point(alpha = 0.7) +
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),   # controls wiggliness
    se = TRUE,
    alpha = 0.15                 # transparency of SE ribbon
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "red", "black", "purple")) +
  facet_wrap(~ Name, ncol = 4, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", fill = "Year")


herbUsingMean <- herbUsingMean %>%
  mutate(Name = reorder(Name, latitude))


herbUsingMean %>% 
  filter( 
    nYear >= 4,
    Year >= 2020, 
    nSurv >= 10
  ) %>% 
  ggplot(aes(x = julianweek, 
             y = Herb_standardized, 
             color = as.factor(Year))) +
  geom_point(alpha = 0.7) +
  stat_smooth(
    aes(fill = as.factor(Year)),
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    alpha = 0.15
  ) +
  scale_color_manual(values = c( "purple", "#E69F00", "#56B4E9", "#009E73", "red", "black")) +
  scale_fill_manual(values = c("purple", "#E69F00", "#56B4E9", "#009E73", "red", "black")) +
  facet_wrap(~ Name, ncol = 5, scales = "free_y") +
  theme_bw() +
  labs(color = "Year", fill = "Year")

