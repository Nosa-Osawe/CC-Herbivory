half_area_point <- function(beta0, beta1, xmax) {
  
  # Zero-crossing (start of positive herbivory)
  x0 <- -beta0 / beta1
  
  # If curve never becomes positive
  if (x0 >= xmax) return(NA)
  
  # Total positive area
  A_total <- (beta0 * xmax + 0.5 * beta1 * xmax^2) -
    (beta0 * x0   + 0.5 * beta1 * x0^2)
  
  # Solve for x* where cumulative area = 0.5 * total area
  # Quadratic: (beta1/2)x^2 + beta0 x + C = 0
  C <- -(beta0 * x0 + 0.5 * beta1 * x0^2 + 0.5 * A_total)
  
  roots <- polyroot(c(C, beta0, beta1 / 2))
  roots <- Re(roots[abs(Im(roots)) < 1e-8])
  
  # Keep the biologically meaningful root
  x_star <- roots[roots >= x0 & roots <= xmax]
  
  if (length(x_star) == 0) return(NA)
  
  return(x_star[1])
}


#############################################################################################



library(dplyr)
library(purrr)

half_area_time <- function(df,
                           x = "julianweek_N",
                           y = "totalHerbS",
                           n_grid = 500) {
  
  # Fit linear model
  fit <- lm(reformulate(x, y), data = df)
  
  # Prediction grid
  x_seq <- seq(min(df[[x]], na.rm = TRUE),
               max(df[[x]], na.rm = TRUE),
               length.out = n_grid)
  
  pred <- predict(fit, newdata = setNames(data.frame(x_seq), x))
  
  # Truncate at zero
  pred_pos <- pmax(pred, 0)
  
  # Numerical integration (trapezoidal rule)
  dx <- diff(x_seq)
  area_segments <- dx * (head(pred_pos, -1) + tail(pred_pos, -1)) / 2
  
  total_area <- sum(area_segments)
  
  if (total_area == 0) {
    return(NA_real_)  # no positive herbivory
  }
  
  cum_area <- cumsum(area_segments)
  
  # Find x where cumulative area reaches 50%
  x_mid <- x_seq[which(cum_area >= 0.5 * total_area)[1] + 1]
  
  return(x_mid)
}



herb_midpoints <- modelGoodHerb %>%
  group_by(siteObserv, Year) %>%
  summarise(
    half_area_julianweek_N =
      half_area_time(cur_data()),
    .groups = "drop"
  )

herb_midpoints %>% data.frame()



############################################################################################################################################
julianWindow = 120:230
# EDA herbivory VS cat occurence and density ----
Herb.Arthropod %>% 
  filter(Name == "Acadia NP - Alder") %>% 
  filter(julianweek %in% julianWindow) %>% 
  filter(nSurvHerb >= min.nSurvWeekYearSite) %>% 
  select(Name, Year, julianweek, totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, Year, julianweek) %>% 
  summarise(totalHerbS = mean(totalHerbS),
            caterpillar_prop = mean(caterpillar_prop) * 100,
            caterpillar_density = mean(caterpillar_density))%>%
  pivot_longer(
    cols = c(totalHerbS, caterpillar_prop),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = julianweek, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Year, scales = "free_x") +
  scale_color_manual(
    values = c("totalHerbS" = "#1b9e77",
               "caterpillar_prop" = "#d95f02"),
    # labels = c("Total herbivory (%)",
    #            "Caterpillar proportion")
  ) +
  labs(
    x = "Julian week",
    y = "Percentage (%)",
    color = "",
    title = "Acadia NP - Alder",
    subtitle = "Seasonal dynamics of herbivory and caterpillar occurrence"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )


Herb.Arthropod %>% 
  filter(Name == "Acadia NP - Sundew") %>% 
  filter(julianweek %in% julianWindow) %>% 
  filter(nSurvHerb >= min.nSurvWeekYearSite) %>% 
  select(Name, Year, julianweek, totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, Year, julianweek) %>% 
  summarise(totalHerbS = mean(totalHerbS),
            caterpillar_prop = mean(caterpillar_prop) * 100,
            caterpillar_density = mean(caterpillar_density))%>%
  pivot_longer(
    cols = c(totalHerbS, caterpillar_prop),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = julianweek, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Year, scales = "free_x") +
  scale_color_manual(
    values = c("totalHerbS" = "#1b9e77",
               "caterpillar_prop" = "#d95f02"),
    # labels = c("Total herbivory (%)",
    #            "Caterpillar proportion")
  ) +
  labs(
    x = "Julian week",
    y = "Percentage (%)",
    color = "",
    title = "Acadia NP - Alder",
    subtitle = "Seasonal dynamics of herbivory and caterpillar occurrence"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )


Herb.Arthropod %>% 
  filter(Name == "Prairie Ridge Ecostation") %>% 
  filter(nSurvHerb >= min.nSurvWeekYearSite) %>% 
  filter(julianweek %in% julianWindow) %>% 
  select(Name, ObservationMethod, Year, julianweek, totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, Year, julianweek) %>% 
  summarise(totalHerbS = mean(totalHerbS),
            caterpillar_prop = mean(caterpillar_prop) * 100,
            caterpillar_density = mean(caterpillar_density))%>%
  pivot_longer(
    cols = c(totalHerbS, caterpillar_prop),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = julianweek, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Year, scales = "free_x") +
  scale_color_manual(
    values = c("totalHerbS" = "#1b9e77",
               "caterpillar_prop" = "#d95f02")) +
  labs(
    x = "Julian week",
    y = "Percentage (%)",
    color = "",
    title = "Acadia NP - Alder",
    subtitle = "Seasonal dynamics of herbivory and caterpillar occurrence"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )



Herb.Arthropod %>% 
  filter(Name == "Prairie Ridge Ecostation") %>% 
  filter(julianweek %in% julianWindow) %>% 
  filter(nSurvHerb >= min.nSurvWeekYearSite) %>% 
  filter(!Year %in% c(2014, 2015)) %>% 
  select(Name, ObservationMethod, Year, julianweek,
         totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(
    totalHerbS = mean(totalHerbS, na.rm = TRUE),
    caterpillar_prop = mean(caterpillar_prop, na.rm = TRUE) * 100,
    caterpillar_density = mean(caterpillar_density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(totalHerbS, caterpillar_prop),
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
    values = c(
      "totalHerbS" = "#1b9e77",
      "caterpillar_prop" = "#d95f02"
    )
  ) +
  labs(
    x = "Julian week",
    y = "Percentage (%)",
    color = "",
    title = "Prairie Ridge Ecostation",
    subtitle = "Seasonal dynamics of herbivory and caterpillar occurrence"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )





Herb.Arthropod %>% 
  filter(Name == "NC Botanical Garden") %>% 
  filter(julianweek %in% julianWindow) %>% 
  filter(!Year %in% c(2014, 2015, 2020)) %>% 
  select(Name, ObservationMethod, Year, julianweek,
         totalHerbS, caterpillar_prop, caterpillar_density) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>% 
  summarise(
    totalHerbS = mean(totalHerbS, na.rm = TRUE),
    caterpillar_prop = mean(caterpillar_prop, na.rm = TRUE) * 100,
    caterpillar_density = mean(caterpillar_density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(totalHerbS, caterpillar_prop),
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
    values = c(
      "totalHerbS" = "#1b9e77",
      "caterpillar_prop" = "#d95f02"
    )
  ) +
  labs(
    x = "Julian week",
    y = "Percentage (%)",
    color = "",
    title = "NC Botanical Garden",
    subtitle = "Seasonal dynamics of herbivory and caterpillar occurrence"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
