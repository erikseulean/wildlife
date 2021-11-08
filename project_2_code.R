#Mock Up code for project but using jupyter for erik
library(statsecol)
library(zoo)
library(jagsUI)
library(MCMCvis)
data("wildebeest")
wildebeest
wildebeest <- na.locf(wildebeest)
wildebeest

#Prior for beta 0
# b0 <- rnorm(0, 0.00001)
# b1 <- rnorm(0, 0.00001)
# sigma.n <- runif(0, 1)
# sigma.y <- runif(0, 1)

sink("project2.txt")
cat("
model{

  # Priors and constraints
  log.n1 ~ dnorm(0.92333, 0.01) # Initial population size
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

#passing in parameters from dataset
wildebeest_dta <- list(
  rain = c(wildebeest$rain, rep(NA, projection)),
  y = c(log(wildebeest$Nhat), rep(NA, projection)),
  nyrs = nrow(wildebeest) + projection
)


# Initial values
wildebeest_inits <- function() {
  list(b0 = runif(1, -2, 2),
       b1 = runif(1, -2, 2),
       sig.y = 1,
       sig.n = 1)
}

parms <- c("b0", "b1")

nc <- 3
nb <- 1000
ni <- 5000 + nb
nt <- 1

out <- jags(data = wildebeest_dta,
            inits = wildebeest_inits,
            parameters.to.save = parms,
            model.file = "project2.txt",
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb,
            n.thin = nt)
out
par(mar = (c(0, 0, 10, 100)))
MCMCtrace(out,                 #the fitted model
          params = parms,      #out parameters of interest
          iter = ni,           #plot all iterations
          pdf = FALSE,         #DON'T write to a PDF
          ind = FALSE)         #chain specific densities
