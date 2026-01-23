##############
# greenup ~ GDD
##############

# specify run date --------------------------------------------------------

run_date <- '2024-07-17'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)

# read in greenup phenometrics and gdd data -----------------------------------

gdd <- data.table::fread('data/gdd_calcs.txt')
gdd = data.table::fread('gdd_calcs.txt')

greenup <- readRDS("data/pheno-data-2020-08-25.rds") %>% 
  filter(!is.na(gr_mn)) %>% 
  distinct(cell,year, .keep_all = TRUE)

arrive <- greenup %>% 
  filter(!is.na(arr_GAM_mean))

na <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"), returnclass = "sf")

#ggplot() +
#  geom_sf(na, mapping = aes(), fill = NA) +
#  geom_point(greenup, mapping = aes(x = cell_lng, y = cell_lat)) +
#  theme_bw() +
#  coord_sf(xlim = c(-100, -55),
#           ylim = c(25,60))
#
#ggplot() +
#  geom_sf(na, mapping = aes(), fill = NA) +
#  geom_point(arrive, mapping = aes(x = cell_lng, y = cell_lat)) +
#  theme_bw() +
#  coord_sf(xlim = c(-100, -55),
#           ylim = c(25,60))

# reduce data areas east of 95 degrees W and 
# reducedata to areas North of 24 degrees N (Youngflesh et al., 2021) & S of 50 degrees N
greenup <- greenup %>% 
  filter(cell_lat > 24 & cell_lat < 50) %>% 
  filter(cell_lng > -95)
  

# process data -------------------------------------------------------------

#join fledge and GDD
# select limited number of data columns
gu_df <- ungroup(greenup) %>%
  dplyr::left_join(gdd) %>% 
  dplyr::select(cell, year, cell_lat, spring.gdd, gr_mn)

#ggplot(gu_df, aes(x = spring.gdd)) +
#  geom_histogram() +
#  facet_wrap(cell~species)

# center gdd for each cell
gu_df <- gu_df %>% 
  filter(year >= 2002 & year <= 2017) %>%  # restrict to same years
  group_by(cell) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()

# center lat across all cells
gu_df <- ungroup(gu_df) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) # center latitude across cell


# only include cells with at least three years of data
yearsGrouping <- gu_df %>% 
  group_by(cell) %>% 
  summarise(nYears = length(unique(year)))
CellToKeep <- filter(yearsGrouping, nYears >= 3)
gu_df <- gu_df %>% 
  filter(cell %in% CellToKeep$cell)

gu_df <- gu_df %>% 
  mutate(cellFactor = cell)

gu_df$cellFactor <- factor(gu_df$cellFactor)
gu_df$cellFactor <- droplevels(gu_df$cellFactor)
gu_df$cellFactor <- as.integer(gu_df$cellFactor)

#order by cell to make sure lat, cell match
gu_df <- gu_df %>%
  arrange(cell, year)

# data for model ----------------------------------------------------------

DATA <- list(y = as.vector(gu_df$gr_mn),
             gdd = as.vector(gu_df$spring_gdd_center),
             Cell = as.vector(gu_df$cellFactor),
             lat = as.vector(distinct(gu_df, 
                                      cell, lat_center))$lat_center[,1],
             N = length(gu_df$gr_mn),
             M = length(unique(gu_df$cellFactor)))


# fit model ---------------------------------------------------------------

options(mc.cores = parallel::detectCores())
#set_cmdstan_path('c:/Users/Mike/Documents/.cmdstan/cmdstan-2.33.1/')
# DELTA <- 0.92
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 4000

#compile model - this will tell you if the model syntax is valid before running the model
mod <- cmdstanr::cmdstan_model('scripts/stan/greenup.stan',
                               force_recompile = TRUE)

#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA
# step_size = STEP_SIZE


#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('greenup-results-', run_date),
                  dir = paste0('results'),
                  # mkdir = paste0('bird-ml-phylo-vint-measured-', run_date),
                  mkdir = paste0('greenup-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('greenup-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('greenup-data-', run_date),
                  cp_file = c(paste0('scripts/stan/greenup.stan'), 
                              paste0('scripts/greenup.R')),
                  cp_file_names = c(paste0('greenup-', run_date, '.stan'),
                                    paste0('greenup-', run_date, '.R')))

MCMCvis::MCMCplot(fit, 
                  params = 'beta')

MCMCvis::MCMCplot(fit, 
                  params = 'gamma')

MCMCvis::MCMCplot(fit, 
                  params = 'theta')

