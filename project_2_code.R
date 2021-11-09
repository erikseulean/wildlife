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
  log.n1 ~ dnorm(N1, 0.01) # Initial population size
  log.N[1] <- log.n1

  b0 ~ dnorm(0, 0.00001)  # Prior for the Intercept
  b1 ~ dnorm(0, 0.00001)  # Prior for the slope

  sig.n ~ dunif(0, 1)     # Standard deviation for the population
  tau.n <- pow(sig.n, -2) # Precision for the state process

  sig.y ~ dunif(0, 1)     # Standard deviation for the observations process
  tau.y <- pow(sig.y, -2) # Precision for the observation process

  # Likelihood - State process
  for(t in 1:(nyrs-1)) {
    log.N[t+1] ~ dnorm(r[t]+ log.N[t], tau.n)
    log(r[t]) <- log.r[t]
    log.r[t] <- b0 + b1*rain[t]
  }

  # Likelihood - Observation process
  for (t in 1:nyrs) {
     log.y[t] ~ dnorm(log.N[t], tau.y)

  }

  # Derive population and observation sizes on real scale
  for (t in 1:nyrs) {
    N.est[t] <- exp(log.N[t])
    y.est[t] <- exp(log.y[t])
  }
}
", fill = TRUE)
sink()

projection <- 6

# Passing in parameters from dataset
wildebeest_model_parameters <- list(
  rain = wildebeest$rain, # vector of rain quantity observed each year
  y = log(wildebeest$Nhat), # vector of population observations
  nyrs = nrow(wildebeest), # number of years included in the data
  N1 = log(wildebeest$Nhat[1]) # initial population estimate
)


# Initial values
wildebeest_inits <- function() {
  list(b0 = runif(1, -2, 2),
       b1 = runif(1, -2, 2),
       sig.n = runif(1, 0, 1),
       sig.y = runif(1, 0, 1)
  )
}

parms <- c("b0", "b1", "N.est", "y.est")

nc <- 3
nb <- 100000
ni <- 100000 + nb
nt <- 1

out <- jags(data = wildebeest_model_parameters,
            inits = wildebeest_inits,
            parameters.to.save = parms,
            model.file = "project2.txt",
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb,
            n.thin = nt)
MCMCtrace(out,                      #the fitted model
          params = parms[1:2],      #out parameters of interest
          iter = ni,                #plot all iterations
          pdf = FALSE,              #DON'T write to a PDF
          ind = FALSE)              #chain specific densities
