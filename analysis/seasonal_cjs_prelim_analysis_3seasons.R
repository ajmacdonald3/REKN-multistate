################################################################################
# PREPARE CJS ENCOUNTER HISTORIES
#
# Season 1 = Apr-Jun
# Season 2 = Jul-Sep
# Season 3 = Oct-Mar
#
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)
library(cowplot)
library(viridis)
library(R2ucare)

################################################################################

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

################################################################################

# create separate m-arrays for first capture and subsequent recaptures

# load data
rekn_enchist <- readRDS("./processed-data/rekn-cjs-enchist.rds")

# convert BirdID to row names so a numeric matrix can be generated
rownames(rekn_enchist) <- NULL

enchist <- rekn_enchist %>% 
  column_to_rownames(var = "BirdID")

enchist <- sapply(enchist, as.numeric)

# function to create a m-array based on capture-histories (CH)
m_array <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# Create separate m-arrays for first capture and subsequent recaptures
CH <- enchist

cap <- apply(CH, 1, sum)
ind <- which(cap >= 2)
CH.R <- CH[ind,] # First period CH recaptured at least once
CH.N <- CH[-ind,] # First period CH never recaptured
# Remove first capture
first <- numeric()
for (i in 1:dim(CH.R)[1]){
  first[i] <- min(which(CH.R[i,]==1))
}
CH.R1 <- CH.R
for (i in 1:dim(CH.R)[1]){
  CH.R1[i,first[i]] <- 0
}
# Create m-array of those recaptured at least once
CH.marray <- m_array(CH.R1)
# Create CH matrix for first period, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.R1)[1]){
  second[i] <- min(which(CH.R1[i,]==1))
}
CH.R2 <- matrix(0, nrow = dim(CH.R)[1], ncol = dim(CH.R)[2])
for (i in 1:dim(CH.R)[1]){
  CH.R2[i,first[i]] <- 1
  CH.R2[i,second[i]] <- 1
}
# Create m-array for these
CH.R.marray <- m_array(CH.R2)
# The last column ought to show the number of transients not recaptured
# again and should all be zeros, since all of them are released as "residents"
CH.R.marray[,dim(CH)[2]] <- 0
# Create the m-array for transients never recaptured and add it to the
# previous m-array
CH.N.marray <- m_array(CH.N)
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

for (s in 1:3){

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

# Derived parameters
phi.nb <- sqrt(mean.phi[3])

for (t in 1:(n.occasions-1)){
  phi.nb.t[t] <- sqrt(phi[t])
}

}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.t = CH.T.marray, marr = CH.marray, n.occasions = dim(CH)[2], rel.t = rowSums(CH.T.marray),
                  rel = rowSums(CH.marray), season = rep(c(1,2,3), length.out = dim(CH.marray)[2]-1))

# Initial values
inits <- function(){list(mean.phi = runif(3, 0, 1), sigma.phi = runif(3, 0, 5),
                         mean.phi.t = runif(3, 0, 1), sigma.phi.t = runif(3, 0, 5),
                         mean.p = runif(3, 0, 1), sigma.p = runif(3, 0, 5))}  

# Parameters monitored
parameters <- c("mean.phi.t", "mean.phi", "mean.p", "phi.t", "phi", "p", "phi.nb", "phi.nb.t",
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

png(filename = "figures/seasonal-cjs-prelim/ppcheck.png", width = 8, height = 8,
    units = "in", res = 600)

print(ppcheck)

dev.off()

mean(sims.list$fit.t.new > sims.list$fit.t)

ppcheck.t <- ggplot(sims.list, aes(x = fit.t, y = fit.t.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/ppcheck.t.png", width = 8, height = 8,
    units = "in", res = 600)

print(ppcheck.t)

dev.off()

saveRDS(cjs.seasonal$summary, "./analysis-output/seasonal_cjs_prelim_analysis_3seasons_summary2.rds")
write.csv(cjs.seasonal$summary, "./analysis-output/seasonal_cjs_prelim_analysis_3seasons_summary2.csv")

saveRDS(cjs.seasonal$sims.list, "./analysis-output/seasonal_cjs_prelim_analysis_3seasons_simslist2.rds")

################################################################################

plot.data <- read_csv("./analysis-output/seasonal_cjs_prelim_analysis_3seasons_plotdata2.csv")
plot.means <- read_csv("./analysis-output/seasonal_cjs_prelim_analysis_3seasons_plotmeans2.csv")

# set custom theme for all plots
theme_cust <- function() {
  theme_classic() %+replace%
    theme(axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=12, colour = "black"),
          axis.title.y = element_text(size=14, angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          strip.text.x = element_text(size=14, face = "bold"),
          legend.text = element_text(size=12),
          #legend.key.height = unit(1, "mm"),
          plot.title = element_text(size = 12, hjust = 0, vjust = 1.5),
          #panel.border = element_rect(size =0.5, fill = "transparent"),
          plot.margin = margin(10, 10, 10, 15))
}

dodge <- position_dodge(.2)

# apparent seasonal survival
phi_plot <- ggplot() +
  # geom_rect(data = plot.means,
  #           aes(xmin=-Inf, xmax=Inf, ymin=lcl, ymax=ucl, fill = season), alpha=0.5) +
  geom_hline(data = plot.means,
             aes(yintercept=mean, colour = season), size = 1) +
  geom_errorbar(data = plot.data, aes(x=as.factor(year), ymin=lcl, ymax=ucl, group = season),
                width=0, size=0.5, colour="black", linetype=1, position = dodge) +
  # geom_line(data = plot.data, aes(x=as.factor(year), y=mean, group = var),
  #           linetype="dashed", size=0.5) +
  geom_point(data = plot.data, aes(x=as.factor(year), y=mean, fill = season, group = season),
             size=3, shape = 21, position = dodge) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  #scale_x_continuous(labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%d %b")) +
  #ylim(0, 1) +
  ylab("seasonal survival") +
  theme_cust() +
  theme(legend.position = c(0.85,0.15),
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.x = element_blank())

png(filename = paste0("figures/seasonal-cjs-prelim-phi-plot2.png"),
    width=6, height=4, units="in", res=600)

plot(phi_plot)

dev.off()

phi_plot2 <- ggplot() +
  geom_rect(data = plot.means,
            aes(xmin=-Inf, xmax=Inf, ymin=lcl, ymax=ucl), fill = "grey90", alpha=0.5) +
  geom_hline(data = plot.means,
             aes(yintercept=mean), colour = "black", size = 0.5, linetype = "dashed") +
  geom_errorbar(data = plot.data, aes(x=as.factor(year), ymin=lcl, ymax=ucl, group = season),
                width=0, size=0.5, colour="black", linetype=1) +
  # geom_line(data = plot.data, aes(x=as.factor(year), y=mean, group = var),
  #           linetype="dashed", size=0.5) +
  geom_point(data = plot.data, aes(x=as.factor(year), y=mean, fill = season, group = season),
             size=2, shape = 21) +
  facet_wrap(. ~ season, nrow = 3) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  #scale_x_continuous(labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%d %b")) +
  #ylim(0, 1) +
  xlab("year") +
  ylab("seasonal survival") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

png(filename = paste0("figures/seasonal-cjs-prelim/phi.point2.png"),
    width=6, height=4, units="in", res=600)

plot(phi_plot2)

dev.off()

mean.phi.plot2 <- ggplot() +
  geom_errorbar(data = plot.means, aes(x=season, ymin=lcl, ymax=ucl, group = season),
                width=0, size=0.5, colour="black", linetype=1) +
  geom_point(data = plot.means, aes(x=season, y=mean, fill = season, group = season),
             size=6, shape = 21) +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/mean.phi.point2.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.phi.plot2)

dev.off()

# plot posterior distributions
sims.list <- readRDS("./analysis-output/seasonal_cjs_prelim_analysis_3seasons_simslist.rds")
sims.list <- as.data.frame(sims.list)

# format output
mean.phi.mod <- sims.list %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.phi.1", "prebreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.2", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.3", "nonbreeding"))

mean.t.mod <- sims.list %>% 
  select(mean.phi.t.1, mean.phi.t.2, mean.phi.t.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.phi.t.1", "prebreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.t.2", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.t.3", "nonbreeding"))

mean.p.mod <- sims.list %>% 
  select(mean.p.1, mean.p.2, mean.p.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.p.1", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.p.2", "nonbreeding")) %>%
  mutate(season = str_replace(season, "mean.p.3", "prebreeding"))

phi.mod <- sims.list %>% 
  select(phi.1, phi.2, phi.3, phi.4, phi.5, phi.6, phi.7, phi.8, phi.9, phi.10,
         phi.11, phi.12, phi.13, phi.14, phi.15, phi.16, phi.17, phi.18, phi.19, phi.20,
         phi.21, phi.22, phi.23, phi.24, phi.25, phi.26, phi.27, phi.28) %>% 
  pivot_longer(cols = 1:28, names_to = "parameter", values_to = "estimate") %>% 
  arrange(match(parameter,
                c("phi.1", "phi.2", "phi.3", "phi.4", "phi.5", "phi.6", "phi.7", "phi.8", "phi.9", "phi.10",
                  "phi.11", "phi.12", "phi.13", "phi.14", "phi.15", "phi.16", "phi.17", "phi.18", "phi.19", "phi.20",
                  "phi.21", "phi.22", "phi.23", "phi.24", "phi.25", "phi.26", "phi.27", "phi.28")),
          desc(estimate)) %>% 
  mutate(season = rep(c("prebreeding", "postbreeding", "nonbreeding"), each = 30000, length.out = 840000)) %>%   
  mutate(year = rep(2009:2018, each = 90000, length.out = 840000))

# make plots
mean.phi.plot <- ggplot() +
  geom_violin(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  #geom_point(mean.phi.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/mean.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.phi.plot)

dev.off()

mean.t.plot <- ggplot() +
  geom_violin(mean.t.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.t.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.t.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  #geom_point(mean.t.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Survival/transience probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/mean.t.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.t.plot)

dev.off()

mean.p.plot <- ggplot() +
  geom_violin(mean.p.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.p.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.p.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  #geom_point(mean.p.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/mean.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.p.plot)

dev.off()

phi.plot <- ggplot() +
  geom_violin(phi.mod, mapping = aes(x = as.factor(year), y = estimate, group = year, fill = season), alpha = 0.6) +
  geom_boxplot(phi.mod, mapping = aes(x = as.factor(year), y = estimate, group = year),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phi.mod, mapping = aes(x = as.factor(year), y = estimate, group = year),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  #geom_point(phi.sim, mapping = aes(x = year, y = value), size = 2, colour = "black") +
  facet_wrap(. ~ season, nrow = 3) +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Year") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-prelim/phi.png", width = 8, height = 5,
    units = "in", res = 600)

print(phi.plot)

dev.off()

