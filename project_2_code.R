#Mock Up code for project but using jupyter for erik
library(statsecol)
library(zoo)
library(jagsUI)
library(MCMCvis)
data("wildebeest")
wildebeest
wildebeest <- na.locf(wildebeest)
wildebeest

sink("project2.txt")
cat("
model{

  # Priors and constraints
  N[1] ~ dnorm(N1, 0.01)

  b0 ~ dnorm(0, 0.0001)  # Prior for the Intercept
  b1 ~ dnorm(0, 0.0001)  # Prior for the slope

  sig.n ~ dunif(0, 1)     # Standard deviation for the population
  tau.n <- pow(sig.n, -2) # Precision for the state process

  sig.y ~ dunif(0, 1)     # Standard deviation for the observations process
  tau.y <- pow(sig.y, -2) # Precision for the observation process

  # Likelihood - State process
  for(t in 1:(nyrs-1)) {
    log(r[t]) <- log.r[t]       # Link function for the parameter
    log.r[t] <- b0 + b1*rain[t] # Linear model for the logarithm of the growth rate

    N[t+1] ~ dnorm(r[t] * N[t] - catch[t], tau.n)
  }

  # Likelihood - Observation process
  for (t in 1:nyrs) {
     y[t] ~ dnorm(N[t], tau.y)

  }
}
", fill = TRUE)
sink()

projection <- 6

# Passing in parameters from dataset
model_parameters <- list(
  # Vector of rain quantity observed each year
  rain = wildebeest$rain,
  # Catch - the amount harvested (or illegal poaching)
  # that happened on the wildebeest population
  catch = wildebeest$Catch, 
  # Vector of population observations
  # We need to consider a logarithm transformation 
  # as the model is modelling the logarithm of the counts
  # in order to avoid negative values
  y = wildebeest$Nhat,
  # Number of years to include in the simulation
  # Each sample will contain nyrs of data
  nyrs = nrow(wildebeest),
  # Initial population estimate
  # We use the logarithm of the initial population
  # when the surveillance of the population started
  N1 = wildebeest$Nhat[1]
)


# Initial values
# For the initial values, we will generate randomly
# values to make sure that multiple chains are not
# starting from the same initial configuration.
initial_values <- function() {
  list(
    b0 = runif(1, -2, 2),
    b1 = runif(1, -2, 2),
    sig.n = runif(1, 0, 1),
    sig.y = runif(1, 0, 1)
  )
}

# Parameters that we want to monitor from the model
# b0, b1, sig.n and sig.y will be monitored to see if the
# chains converge to the posterior distribution while N.est
# y.est will be used to draw conclusions. 
monitored_parameters <- c("b0", "b1",  "sig.n", "sig.y", "N")

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

# Number of thinning samples. This is set to 1
# so no thinning occurs, meaning that we consider
# the posterior samples to be sufficiently independent
# to the point that no thinning is necessary.
nt <- 1

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
