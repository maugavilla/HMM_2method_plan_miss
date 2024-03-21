

library(mHMMbayes)
library(depmixS4)
library(simsem)
library(portableParallelSeeds)
library(fpc)



# simulating data for 10 subjects with each 100 categorical observations
n_t <- 20 ## time points
n <- 10 ## subjects
m <- 3 ## states
n_dep <- 2 ## dependent variables
q_emiss <- c(4,4) ## number of categories per dep var

## transition matrix
gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.7, 0.1,
                  0.2, 0.2, 0.6), ncol = m, byrow = TRUE)

## emission probabilities for each n_dep. row=m, col=q_emiss
emiss_distr <- list(matrix(c(0.5, 0.5, 0.0, 0.0,
                             0.1, 0.1, 0.8, 0.0,
                             0.0, 0.0, 0.1, 0.9), nrow = m, 
                           ncol = q_emiss, byrow = TRUE),
                    matrix(c(0.6, 0.4, 0.0, 0.0,
                             0.2, 0.0, 0.8, 0.0,
                             0.0, 0.1, 0.1, 0.8), nrow = m, 
                           ncol = q_emiss, byrow = TRUE))

data1 <- sim_mHMM(n_t = n_t, 
                  n = n, 
                  gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                  gamma = gamma, 
                  emiss_distr = emiss_distr, 
                  var_gamma = 0, # variance pf transition matrix 
                  var_emiss = c(0,0), # variance between emission matrices n_dep length
                  return_ind_par = T)
head(data1$obs)
head(data1$states)
data1$subject_gamma


dat <- data.frame(data1$obs)
colnames(dat) <- c("id","v1","v2")
head(dat)

### missing data
#?simsem::impose

dat_mis <- imposeMissing(dat, pmMCAR = .2, ignoreCols = "id")
head(dat_mis)
summary(dat_mis)


###
### run with depmix
###
mod <- depmix(list(v1~1,v2~1),data=dat,
              nstates=3,
              family=list(multinomial("identity"),
                          multinomial("identity")),
              ntimes= as.numeric(table(dat$id)) )
# print the model, formulae and parameter values
mod

# fit the model by calling fit
#fm <- fit(mod)
#summary(fm)

fm2 <- multistart(mod, nstart=10, initIters=50000)
summary(fm2)

sts <- viterbi(fm2)


# run with missing

mod_m <- depmix(list(v1~1,v2~1),data=dat_mis,
              nstates=3,
              family=list(multinomial("identity"),
                          multinomial("identity")),
              ntimes= as.numeric(table(dat$id)) )
# print the model, formulae and parameter values
mod_m

# fit the model by calling fit
#fm_m <- fit(mod_m)
#summary(fm_m)

fm2_m <- multistart(mod_m, nstart=10, initIters=50000)
summary(fm2_m)

sts_m <- viterbi(fm2_m)

dface <- dist(dat[,-1])
cluster.stats(d=dface,
              clustering = sts_m$state, # data1$states[,"state"]
              alt.clustering =  data1$states[,"state"])[c("entropy","corrected.rand","vi")]

table(data1$states[,"state"], 
      sts_m$state)


### impose missing by logical matrix

head(dat)
dim(dat)
summary(dat)

## create logical matrix, where TRUE equals missing
miss_mat <- matrix(0, ncol = 2, nrow = nrow(dat))
miss_mat[,1] <- imposeMissing(cbind(miss_mat[,1]), pmMCAR = .1)
miss_mat[,2] <- imposeMissing(cbind(miss_mat[,2]), pmMCAR = .7)
miss_mat <- data.frame(is.na(miss_mat))

dat_mis2 <- imposeMissing(dat, ignoreCols = "id", 
                          logical = miss_mat )
head(dat_mis2)
summary(dat_mis2)



mod_m2 <- depmix(list(v1~1+cov1,v2~1),data=dat_mis2,
                nstates=3,
                family=list(multinomial("identity"),
                            multinomial("identity")),
                ntimes= as.numeric(table(dat$id)) )
# print the model, formulae and parameter values
mod_m2

# fit the model by calling fit
#fm_m2 <- fit(mod_m2)
#summary(fm_m2)

fm2_m2 <- multistart(mod_m2, nstart=10, initIters=50000)
summary(fm2_m2)

sts_m2 <- viterbi(fm2_m2)

dface <- dist(dat[,-1])
cluster.stats(d=dface,
              clustering = sts_m2$state, # data1$states[,"state"]
              alt.clustering =  data1$states[,"state"])[c("entropy","corrected.rand","vi")]



## no missing
## pmcar: increase at LFS indicator
## structure missing on LFS indicator: longitudinal, fully excluded, mixed

## are predictions "equal" to complete indicator: 
## compare to real state and complete indicator
## out of sample pred: use subsample for analysis, predict left out

## add predictors to emiis probs for V1? as measure of bias

## classirication accuracy measures? like
ss <- subset(sts_m2, state == 1)
apply(cbind(ss[,2]), 2, function(x){c(mean(x), sd(x))})

