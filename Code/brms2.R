  library(brms)
  library(bayesplot)
  library(performance)
  library(tidyverse)
  library(ecostats)
  try(dev.off())
  
  data("globalPlants")
  plot(globalPlants$lat, log(globalPlants$height))
  
   
  fit1 = brm(
    log(height) ~ lat,
    data = globalPlants,
    prior = prior(normal(-0.05, 0.02), class = "b", coef = "lat"),
    backend = "cmdstanr",
    chains = 4,
    cores = 4
  )
  
  
  summary(fit1)
  plot(fit1)
  
  plot(conditional_effects(fit1), points = T)
  
  
  # Comparing the posterior (blue) to the prior (red)
  
  mcmc_dens(fit1, pars = c("b_lat")) + # posterior
    geom_function(                    # Prior
      fun  = dnorm,
      args = list(mean = -0.05, sd = 0.02),
      colour = "red",
      linewidth = 1.5
    ) +
    xlim(-0.1, 0)
  
  
  
  
  post = as_draws_matrix(fit1)
  str(post)
  
  mcmc_hist(fit1, pars = "Intercept")
  
  pairs(
    fit1,
    variable = c("b_Intercept", "b_lat", "sigma"),
    off_diag_fun = "hex"
  )
  
  pairs( # This is what happens internally: Brms uses the mean centered intercept.
    fit1,
    variable = c("Intercept", "b_lat", "sigma"),
    off_diag_fun = "hex"
  )
  
  
  
  head(post[, 1:3])
  plot(conditional_effects(fit1), points = TRUE)
  plot(conditional_effects(fit1, prob = 0.80), points = T)
  plot(conditional_effects(fit1, prob = 0.99), points = T)
  plot(conditional_effects(fit1, spaghetti = TRUE, ndraws = 500), points = TRUE) # draws from 500 samples. (we have 4000!)
  
  
  
  
  plot(conditional_effects(fit1, method = "posterior_predict"), points = TRUE) # prediction interval 
                                                                  # that uses botht the deterministic and stochastic process
  
  fitted(fit1) %>% round(3)
  predict(fit1) %>% round(3)
  
  fitted(fit1, newdata = data.frame(lat = c(60))) %>% round(3)
  
  hypothesis(fit1, "lat<0") # the entire probability mass is negative.
  
  mu40 = posterior_epred(fit1, newdata = data.frame(lat = 40))
  mu60 = posterior_epred(fit1, newdata = data.frame(lat = 60))
  
  hist(mu40)
  hist(mu60)
  
  sum(mu40>mu60)/ length(mu40) # all difference are positive
  
  
  
  # - Posterior predictive checks ----
  
  bayes_R2(fit1)
  pp_check(fit1, ndraws = 100)
  
  # classical observed vs predicted plot
  
  pp_check(fit1, type = "scatter_avg")
  
  
  
  check_model(fit1)
  
  
  # - Prior predictive checks ----
  
  get_prior(fit1)
  
  fit1.prior=  brm(
    log(height) ~ lat,
    data = globalPlants,
    prior = prior(normal(-0.05, 0.02), class = "b", coef = "lat"),
    sample_prior = "only",
    backend = "cmdstanr",
    chains = 4,
    cores = 4
  )
  
  summary(fit1.prior, prior = TRUE)
  plot(conditional_effects(fit1.prior), points = TRUE)
  plot(conditional_effects(fit1.prior, spaghetti = TRUE, ndraws = 200), points = TRUE)
  pp_check(fit1.prior, ndraws = 100)
  
  
  fit2.prior=  brm(
    log(height) ~ lat,
    data = globalPlants,
    prior = c(prior(normal(-0.05, 0.02), class = "b", coef = "lat"),
              prior(student_t(3, 1.1, 1), class = Intercept)),
    sample_prior = "only",
    backend = "cmdstanr",
    chains = 4,
    cores = 4
  )
  
  summary(fit2.prior, priors = TRUE)
  
  plot(conditional_effects(fit2.prior), points = TRUE)
  plot(conditional_effects(fit2.prior, spaghetti = TRUE, ndraws = 100), points = TRUE)


stancode(fit1)
    


###############################################################################################################----

library(emmeans)


globalPlants$z.lat = scale(globalPlants$lat)
globalPlants$z.rain = scale(globalPlants$rain)

ggplot(globalPlants, aes(z.lat, log(height))) + geom_point(alpha = 0.5)
ggplot(globalPlants, aes(z.rain, log(height))) + geom_point(alpha = 0.5)

fit.ln.add = brm(log(height) ~ z.lat + z.rain,
                 prior = prior(normal(-1, 1), class = b, coef = z.lat) +
                   prior(normal(+1, 1), class = b, coef = z.rain),
                 backend = "cmdstanr",
                 data = globalPlants)

summary(fit.ln.add, prior = TRUE)
plot(fit.ln.add)



plot(conditional_effects(fit.ln.add),
     points = TRUE,
     point_args = c(alpha = .5),
     ask = FALSE)

# plot by variable

plot(conditional_effects(fit.ln.add, effects = "z.lat"),
     points = TRUE,
     point_args = c(alpha = .5),
     ask = FALSE)

# although the model does not contain an interaction
# prediction for 3 levels of rain (mean -1sd, mean, mean + 1 sd)


plot(conditional_effects(fit.ln.add, effects = "z.rain:z.lat", prob  = 0))
plot(conditional_effects(fit.ln.add, effect = "z.lat:z.rain", prob = .8))

# model evaluation through posterior predictive checks
pp_check(fit.ln.add, ndraws = 100)
check_model(fit.ln.add)


update.packages()
install.packages("devtools")    # install packages from github
install.packages("brms")        # our main software package
install.packages("ggplot2")     # plotting
install.packages("bayesplot")   # additional plotting tools
install.packages("sfsmisc")     # mathematical integration through data points
install.packages("performance") # model evaluation
install.packages("arm")         # model evaluation
install.packages("GGally")      # pairs plots
install.packages("emmeans")     # post-hoc analysis
install.packages("ecostats")    # some datasets
devtools::install_github("jfieberg/Data4Ecologists", force = TRUE) # more datasets

library(GGally)
library(Data4Ecologists)







