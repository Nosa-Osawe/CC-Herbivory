
################################################## Bayesian modelling ##############################################

# Run the provius code in Herb_arthropod analysis first else this would not work.



modelGoodHerb %>% 
  filter(siteObserv == "Prairie Ridge Ecostation_Visual",
         Year == 2024) %>% 
  ggplot(aes( x = julianweek_N, y = totalHerbS))+
  geom_point(size =10)


prairie2024_prior = brm(
  totalHerbS ~ julianweek_N,
  data = modelGoodHerb %>% 
    filter(siteObserv == "Prairie Ridge Ecostation_Visual",
           Year == 2024),
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "Intercept", lb = 0),
    prior(normal(0, 10), class = "b")),
  sample_prior = "only",
  iter = 1000, chains = 2, warmup = 500,
  cores = 2, backend = "cmdstanr"
)

pp_check(prairie2024_prior, ndraws = 100) +
  ggtitle("Prior predictive check: Prairie Ridge 2024")+
  coord_cartesian(xlim = c(-10, 10)) +
  theme_minimal()

yrep <- posterior_predict(prairie2024_prior, draws = 100)
range(yrep)


conditional_effects(
  prairie2024_prior,
  effects = "julianweek_N",
  spaghetti = TRUE,
  ndraws = 50
) 





praire2024 = brm(
  totalHerbS ~ julianweek_C,
  data = modelGoodHerb %>% 
    filter(siteObserv == "Prairie Ridge Ecostation_Visual",
           Year == 2024),
  family = gaussian(),
  prior = c(
    prior(normal(10, 10), class = "Intercept", lb = 0),
    prior(normal(0,3), class = "b")
  ),
  iter = 1000, warmup = 500, chains = 2, 
  cores = 2, backend = "cmdstanr"
)

summary(praire2024)
plot(praire2024)

plot(conditional_effects(praire2024), 
     points=TRUE, 
     point_args=c(alpha=0.5, size =10),
     ask=FALSE)

########################################################################################


 

 

library(brms)

fit_herb_model_brms <- function(df) {
  
  fit <- brm(
    totalHerbS ~ julianweek_N,
    data = df,
    family = gaussian(),
    prior = c(
      prior(normal(1, 1), class = "Intercept", lb = 0),
      prior(normal(0, 10), class = "b")
    ),
    iter = 1000, warmup = 500,
    chains = 2, cores = 2,
    backend = "cmdstanr",
    silent = TRUE, refresh = 0
  )
  
  post <- posterior_summary(fit, probs = c(0.025, 0.975))
  
  r2 <- bayes_R2(fit)
  
  tibble(
    intercept     = post["Intercept", "Estimate"],
    effect        = post["b_julianweek_N", "Estimate"],
    intercept_lwr = post["Intercept", "Q2.5"],
    intercept_upr = post["Intercept", "Q97.5"],
    effect_lwr    = post["b_julianweek_N", "Q2.5"],
    effect_upr    = post["b_julianweek_N", "Q97.5"],
    r2            = mean(r2)   # posterior mean R²
  )
}

coef_table <- modelGoodHerb %>%  # replace with desired dataframe
  filter(siteObserv == "Acadia NP - Alder_Visual") %>%
  group_by(siteObserv, Year) %>%
  group_modify(~ fit_herb_model_brms(.x)) %>%
  ungroup()



 
