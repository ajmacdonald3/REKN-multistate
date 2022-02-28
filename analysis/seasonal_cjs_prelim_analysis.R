################################################################################
# RED KNOT SEASONAL SURVIVAL - CJS MODEL
#
# using m-array format
#
# single state seasonal survival before moving to multistate
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)
library(cowplot)
library(viridis)
library(R2ucare)

set.seed(42)
theme_set(theme_bw())

# GOF testing
rekn <- read_inp("./processed-data/rekn-cjs-enchist.inp")

rekn_hist <- rekn$encounter_histories
rekn_freq <- rekn$sample_size

test3sr(rekn_hist, rekn_freq) # transience
test3sm(rekn_hist, rekn_freq) # resighting rates
test2ct(rekn_hist, rekn_freq) # equal resightability
test2cl(rekn_hist, rekn_freq) # equal resightability (before and after)
overall_CJS(rekn_hist, rekn_freq)

# create separate m-arrays for first capture and subsequent recaptures
enchist_cjs <- readRDS("./processed-data/rekn-cjs-enchist.rds")
CH <- enchist_cjs[,2:39]

cap <- apply(CH, 1, sum)
ind <- which(cap >= 2)
CH.R <- CH[ind,] # first period CH recaptured at least once
CH.N <- CH[-ind,] # first period CH never recaptured
# remove first capture
first <- numeric()
for (i in 1:dim(CH.R)[1]){
  first[i] <- min(which(CH.R[i,]==1))
}
CH.R1 <- CH.R
for (i in 1:dim(CH.R1)[1]){
  CH.R1[i,first[i]] <- 0
}
# create m-array of those recaptured at least once
CH.marray <- marray(CH.R1)
# create CH matrix for first period, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.R1)[1]){
  second[i] <- min(which(CH.R1[i,]==1))
}
CH.R2 <- matrix(0, nrow = dim(CH.R)[1], ncol = dim(CH.R)[2])
for (i in 1:dim(CH.R)[1]){
  CH.R2[i, first[i]] <- 1
  CH.R2[i, second[i]] <- 1
}
# create m-array for these
CH.R.marray <- marray(CH.R2)
# the last column ought to show the number of transients not recaptured
# again and should all be zeros, since all of them are released as "residents"
CH.R.marray[,dim(CH)[2]] <- 0
# create the m-array for transients never recaptured and addit to the previous m-array
CH.N.marray <- marray(CH.N)
CH.T.marray <- CH.R.marray + CH.N.marray

# Specify model in BUGS language

sink("cjs-seasonal.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
logit(phi.trans[t]) <- mu.t[season[t]] + epsilon.t[t]
epsilon.t[t] ~ dnorm(0, tau.t[season[t]])T(-15,15) # Range restriction
phi.t[t] <- 1/(1+exp(-mu.t[season[t]]-epsilon.t[t]))

logit(phi.app[t]) <- mu.phi[season[t]] + epsilon.phi[t]
epsilon.phi[t] ~ dnorm(0, tau.phi[season[t]])T(-15,15) # Range restriction
phi[t] <- 1/(1+exp(-mu.phi[season[t]]-epsilon.phi[t]))

logit(psight[t]) <- mu.p[season[t]] + epsilon.p[t]
epsilon.p[t] ~ dnorm(0, tau.p[season[t]])T(-15,15) # Range restriction
p[t] <- 1/(1+exp(-mu.p[season[t]]-epsilon.p[t]))

}

for (s in 1:4){

mu.t[s] <- log(mean.phi.t[s] / (1-mean.phi.t[s]))
mean.phi.t[s] ~ dunif(0, 1) # Prior for mean transience
sigma.t[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.t[s] <- pow(sigma.t[s], -2)
sigma.t2[s] <- pow(sigma.t[s], 2)

mu.phi[s] <- log(mean.phi[s] / (1-mean.phi[s]))
mean.phi[s] ~ dunif(0, 1) # Prior for mean survival
sigma.phi[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.phi[s] <- pow(sigma.phi[s], -2)
sigma2.phi[s] <- pow(sigma.phi[s], 2)

mu.p[s] <- log(mean.p[s] / (1-mean.p[s]))
mean.p[s] ~ dunif(0, 1) # Prior for mean survival
sigma.p[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.p[s] <- pow(sigma.p[s], -2)
sigma2.p[s] <- pow(sigma.p[s], 2)

}

# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
marr.t[t,1:n.occasions] ~ dmulti(pr.t[t,], rel.t[t])
marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
}

# Define the cell probabilities of the m-arrays
# Main diagonal
for (t in 1:(n.occasions-1)){
q[t] <- 1-p[t] # Probability of non-recapture
pr.t[t,t] <- phi.trans[t]*p[t]
pr[t,t] <- phi[t]*p[t]

# Above main diagonal
for (j in (t+1):(n.occasions-1)){
pr.t[t,j] <- phi.trans[t]*prod(phi[(t+1):j])*prod(q[t:(j-1)])*p[j]
pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
} # j

# Below main diagonal
for (j in 1:(t-1)){
pr.t[t,j] <- 0
pr[t,j] <- 0
} # j
} # t

# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
pr.t[t,n.occasions] <- 1-sum(pr.t[t,1:(n.occasions-1)])
pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
} # t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      expmarr.t[t,j] <- rel.t[t]*pr.t[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.org.t[t,j] <- pow((pow(marr.t[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   marr.t.new[t,1:n.occasions] ~ dmulti(pr.t[t, ], rel.t[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.t.new[t,j] <- pow((pow(marr.t.new[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
fit <- sum(E.org[,])
fit.t <- sum(E.org.t[,])
fit.new <- sum(E.new[,])
fit.t.new <- sum(E.t.new[,])
}
",fill = TRUE)
sink()

# sink("cjs-seasonal.jags")
# cat("
# model {
# 
# # Priors and constraints
# for (t in 1:(n.occasions-1)){
#    logit(phi[t]) <- mu.phi[season[t]] + eps.phi[t]
#    eps.phi[t] ~ dnorm(0, tau.phi[season[t]])T(-10,10)
#    logit(p[t]) <- mu.p[season[t]] + eps.p[t] 
#    eps.p[t] ~ dnorm(0, tau.p[season[t]])T(-10,10)
#    }
# 
# for(s in 1:4){
#    mean.phi[s] ~ dunif(0, 1)             # Prior for mean survival
#    mu.phi[s] <- logit(mean.phi[s])       # Logit transformation
#    sigma.phi[s] ~ dunif(0, 10)               # Prior for standard deviation
#    tau.phi[s] <- pow(sigma.phi[s], -2)
#    sigma2.phi[s] <- pow(sigma.phi[s], 2)
#    
#    mean.p[s] ~ dunif(0, 1)             # Prior for mean resighting
#    mu.p[s] <- logit(mean.p[s])         # Logit transformation
#    sigma.p[s] ~ dunif(0, 10)               # Prior for standard deviation
#    tau.p[s] <- pow(sigma.p[s], -2)
#    sigma2.p[s] <- pow(sigma.p[s], 2)
# }
# 
# # Temporal variance on real scale
# # sigma2.real <- sigma2 * pow(mean.phi, 2) * pow((1-mean.phi), 2) 
# 
# 
# # Define the multinomial likelihood
# for (t in 1:(n.occasions-1)){
#    marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
#    }
# 
# # Define the cell probabilities of the m-array:
# # Main diagonal
# for (t in 1:(n.occasions-1)){
#    q[t] <- 1-p[t]
#    pr[t,t] <- phi[t]*p[t]	
#    # Above main diagonal
#    for (j in (t+1):(n.occasions-1)){
#       pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
#       } #j	
#    # Below main diagonal
#    for (j in 1:(t-1)){
#       pr[t,j]<-0
#       } #j
#    } #t
# 
# # Last column: probability of non-recapture
# for (t in 1:(n.occasions-1)){
#    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
#    } # t
# 
# # Assess model fit using Freeman-Tukey statistic
# # Compute fit statistics for observed data
# for (t in 1:(n.occasions-1)){
#    for (j in 1:n.occasions){
#       expmarr[t,j] <- rel[t]*pr[t,j]
#       E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
#       }
#    }
# # Generate replicate data and compute fit stats from them
# for (t in 1:(n.occasions-1)){
#    marr.new[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
#    for (j in 1:n.occasions){
#       E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
#       }
#    }
# fit <- sum(E.org[,])
# fit.new <- sum(E.new[,])
# 
# }
# ",fill = TRUE)
# sink()

# Bundle data
jags.data <- list(marr.t = CH.T.marray, marr = CH.marray, n.occasions = dim(CH)[2], rel.t = rowSums(CH.T.marray),
                  rel = rowSums(CH.marray), season = rep(c(1,2,3,4), length.out = dim(CH.marray)[2]-1))

# Initial values
inits <- function(){list(mean.phi = runif(4, 0, 1), sigma.phi = runif(4, 0, 5),
                         mean.phi.t = runif(4, 0, 1), sigma.phi.t = runif(4, 0, 5),
                         mean.p = runif(4, 0, 1), sigma.p = runif(4, 0, 5))}  

# Parameters monitored
parameters <- c("mean.phi.t", "mean.phi", "mean.p", "phi.t", "phi", "p",
                "sigma.2t", "sigma2.phi", "sigma2.p", "fit.t", "fit.t.new", "fit", "fit.new")

# MCMC settings
ni <- 500000
nt <- 30
nb <- 200000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.seasonal <- jags(jags.data, inits, parameters, "cjs-seasonal.jags", n.chains = nc, n.thin = nt,
                     n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.seasonal, digits = 3)

sims.list <- cjs.seasonal$sims.list
sims.list <- as.data.frame(sims.list)

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

ppcheck

mean(sims.list$fit.t.new > sims.list$fit.t)

ppcheck.t <- ggplot(sims.list, aes(x = fit.t, y = fit.t.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

ppcheck.t
