#Cleaned up code

library(statsecol)
library(zoo)
library(jagsUI)
library(MCMCvis)
library(coda)
library(dplyr)
library(formattable)
library(ggplot2)
data("wildebeest")

non_na_indices = which(!is.na(wildebeest$Nhat))
wildebeest <- na.locf(wildebeest, fromLast=T)


wildebeest$growth <- c(wildebeest$Nhat[-1]/wildebeest$Nhat[-nrow(wildebeest)], 0)




sink("project2.txt")
cat("
model{

  # Priors and constraints
  N[1] ~ dnorm(N1, 0.01)

  b0 ~ dnorm(0, 0.0001)  # Prior for the Intercept
  b1 ~ dnorm(0, 0.0001)  # Prior for the slope

  sig.n ~ dunif(0, 1)     # Standard deviation for the population
  tau.n <- pow(sig.n, -2) # Precision for the state process

  # Likelihood - State process
  for(t in 1:(nyrs-1)) {
    log(r[t]) <- log.r[t]       # Link function for the parameter
    log.r[t] <- b0 + b1*rain[t] # Linear model for the logarithm of the growth rate
    
    # We use T(0, ) to make sure that the values sampled 
    # are positively truncated. This means that if a negative
    # value is sampled, a new value will be sampled until positive.
    # We want the population counts to be positive all the time
    # but the support of the normal distribution is (-inf, inf)
    # so we need this truncation to conform with the boundaries
    # of a real population.
    N[t+1] ~ dnorm(r[t] * N[t] - catch[t], tau.n)T(0,)
  }

  # Likelihood - Observation process
  for (t in 1:nyrs) {
     y[t] ~ dnorm(N[t], tau.y[t])

  }
}
", fill = TRUE)
sink()


set.seed(5)


#Years we want to project forward
projection <- 5

# Passing in parameters from dataset
model_parameters <- list(
  # Vector of rain quantity observed each year
  # plus extra values for projection
  rain = c(wildebeest$rain,rep(mean(wildebeest$rain), projection)),
  # Catch - the amount harvested (or illegal poaching)
  # that happened on the wildebeest population including
  # values for the projection.
  catch = c(wildebeest$Catch,rep(0.04, projection)), 
  # Vector of population observations
  # We need to consider a logarithm transformation 
  # as the model is modelling the logarithm of the counts
  # in order to avoid negative values
  y = c(wildebeest$Nhat, rep(NA,projection)),
  # Number of years to include in the simulation
  # Each sample will contain nyrs of data
  nyrs = nrow(wildebeest) + projection,
  # Initial population estimate
  # We use the logarithm of the initial population
  # when the surveillance of the population started
  N1 = wildebeest$Nhat[1],
  # precision for the observation process
  tau.y = c(wildebeest$sehat^(-2), sample(wildebeest$sehat^(-2), projection)),
  # indices for which the Nhat column in the wildebeest dataset has value
  non_na_indices = non_na_indices
)


# Initial values
# For b0, b1 and sig.n, we will generate randomly
# values to make sure that multiple chains are not
# starting from the same initial configuration.
initial_values <- function() {
  list(
    b0 = runif(1, -3, 3),
    b1 = runif(1, -3, 3),
    sig.n = runif(1, 0, 1),
    N = c(wildebeest$Nhat, rep(NA, projection))
  )
}

# Parameters that we want to monitor from the model
# b0, b1, sig.n and sig.y will be monitored to see if the
# chains converge to the posterior distribution while N.est
# y.est will be used to draw conclusions. 
monitored_parameters <- c("b0", "b1", "sig.n", "N")

# Number of chains that we want to run
nc <- 3
# Number of iterations for burn-in. These
# iteration will be discarded before considering
# summaries or chain convergence diagnostics.
nb <- 100000
# Total number of iterations to run the chain.
# We started with smaller samples and increased
# the values until we got to some satisfactory
# effective sample sizes for the parameters of
# interest.
ni <- 100000 + nb

# Number of thinning samples. This is set to 10
# because a thinning of 1 shows a very high level
# of autocorrelation for all 4 parameters of interest.
# With a thinning of 10, only sig.y has autocorrelation.
nt <- 10

# Start the MCMC algorithm with the parameters provided
out <- jags(
  data = model_parameters,
  inits = initial_values,
  parameters.to.save = monitored_parameters,
  model.file = "project2.txt",
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.thin = nt
)

# Draw trace plots and density plots for the monitored parameters
MCMCtrace(
  out,                                 # the fitted model
  params = monitored_parameters[1:4],  # out parameters of interest
  iter = ni,                           # plot all iterations
  pdf = FALSE,                         # don't write to a PDF
  ind = FALSE                          # chain specific densities
)




MCMC_Param_Summary = MCMCsummary(out,
                                 params = monitored_parameters[1:3] #out parameters of interest
)



#Putting outputs into a data frame
projection_plot_df <- data.frame(
  Year = c(wildebeest$year,1990:1994),
  Mean = out$mean$N,
  Lower = out$q2.5$N,
  Upper = out$q97.5$N,
  Obs = c(wildebeest$Nhat,rep(NA,projection))
)

ggplot(data = projection_plot_df) + 
  geom_ribbon(aes(x=Year, y=Mean, ymin=Lower, ymax=Upper),
              fill="grey", alpha = 0.25) +
  geom_line(aes(x=Year, y=Mean), size=1.2, color="red") + 
  geom_line(aes(x=Year, y=Obs), color = "green") +
  geom_point(aes(x=Year, y=Obs), size=1.2) +
  theme_bw()+
  labs(title = 'Population abundance and a five year projection (1990-1994)',
       y = 'Mean Wildebeest Abundance (Millions)',
       x = 'Observation Year')+
  theme(plot.title = element_text(hjust = 0.5))


#Putting outputs into a data frame
projection_plot_df <- data.frame(
  Year = c(wildebeest$year,1990:1994),
  Mean = out$mean$N,
  Lower = out$q2.5$N,
  Upper = out$q97.5$N,
  Obs = c(wildebeest$Nhat,rep(NA,projection))
)

projection_plot_df = projection_plot_df[2:35, ]



#Formatting the table for Our Estimates
estimates_data_frame = data.frame(
  Year = c(wildebeest$year,1990:1994),
  Lower = round(out$q2.5$N,2),
  `Mean_Abundance` = round(out$mean$N,2),
  Upper = round(out$q97.5$N,2),
  stringsAsFactors = FALSE
)

#Setting Some Colors to Call for the table
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"
light_blue = '#ADD8E6'
light_green = '#DAF7A6'
light_yellow = '#FFFFE0'

formatted_CI = formattable(
  estimates_data_frame[31:35,],
  caption = "Credible Interval Table",
  align =c("l","c","c","c","c"),
  list(Year = formatter(
    "span", style = ~ style(color = "black",font.weight = "bold")),
    Mean_Abundance = color_tile(light_green, light_green),
    `Lower`= color_tile(light_blue, light_blue),
    `Upper`= color_tile(light_yellow, light_yellow)
  )
)
formatted_CI

MCMC_Param_Summary = data.frame(round(MCMC_Param_Summary,2))
colnames(MCMC_Param_Summary) = c("Mean", "SD", "2.5%", "50%", "97.5%", "Rhat", "n.eff")


formatted_params = formattable(
  MCMC_Param_Summary,
  caption = "Parameter Summaries",
  align =c("c","c","c","c","c","c","c"),
  list(mean = formatter(
    "span", style = ~ style(color = "black",font.weight = "bold")),
    Mean = color_tile(light_green, light_green),
    SD = color_tile(light_green, light_green),
    `2.5%`= color_tile(light_blue, light_blue),
    `50%`= color_tile(light_yellow, light_yellow),
    `97.5%` = color_tile(light_green, light_green),
    `Rhat`= color_tile(light_blue, light_blue),
    `n.eff`= color_tile(light_yellow, light_yellow)
  )
)
formatted_params


growth_rain <- data.frame(
  Growth = wildebeest$growth[-nrow(wildebeest)],
  Rain = wildebeest$rain[-nrow(wildebeest)]
)

library(ggthemes)

ggplot(data = growth_rain, aes(x = Growth, y = Rain))+
  geom_point(size = 2, color = 'red')+
  labs(title = "Rain vs Growth",
       y = "Rain",
       x = "Growth") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15))

growth_rain$Sqrtgrowth <- sqrt(growth_rain$Growth)
linear_model <- lm(Sqrtgrowth ~ Rain, data=growth_rain)
linear_model_df = data.frame(round(anova(linear_model), 2))
colnames(linear_model_df) = c("Df", "SSE", "MSE", "F-Val", "P-Val")

formatted_lm = formattable(
  linear_model_df,
  caption = "ANOVA Summary",
  align =c("c","c","c","c","c"),
  list(df = formatter(
    "span", style = ~ style(color = "black",font.weight = "bold")),
    Df = color_tile(light_green, light_green),
    SSE = color_tile(light_green, light_green),
    `MSE`= color_tile(light_blue, light_blue),
    `F-Val`= color_tile(light_yellow, light_yellow),
    `P-Val` = color_tile(light_green, light_green)
  )
)

formatted_lm

fit <- glm(formula = Growth ~ Rain,  family = Gamma(link="log"), data=growth_rain)
glm_df = data.frame(round(anova(fit, test="F"), 2))
glm_df
colnames(glm_df) = c("Df", "Deviance", "Resid-Df","Resid-Deviance", "F-Val", "P-Val")
rownames(glm_df) = c("Intercept", "Rain")



formatted_glm = formattable(
  glm_df,
  caption = "ANOVA Summary",
  align =c("c","c","c","c","c"),
  list(df = formatter(
    "span", style = ~ style(color = "black",font.weight = "bold")),
    Df = color_tile(light_green, light_green),
    Deviance = color_tile(light_green, light_green),
    `Resid-Df`= color_tile(light_blue, light_blue),
    `F-Val`= color_tile(light_yellow, light_yellow),
    `P-Val` = color_tile(light_green, light_green),
    `Resid-Deviance` = color_tile(light_green, light_green)
  )
)
formatted_glm









