

library(mHMMbayes)
library(depmixS4)
library(simsem)
library(portableParallelSeeds)
library(fpc)
## cluster.stats(d=distance,
##clustering = clust_ward_OM,
##alt.clustering =  sim_class)$corrected.rand,


# simulating data for 10 subjects with each 100 categorical observations
n_t <- 100 ## time points
n <- 10 ## subjects
m <- 3 ## states
n_dep <- 2 ## dependent variables
q_emiss <- c(4,4) ## number of categories

## transition matrix
gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.7, 0.1,
                  0.2, 0.2, 0.6), ncol = m, byrow = TRUE)

## emission probabilities for each n_dep
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
?simsem::impose

dat_mis <- imposeMissing(dat, pmMCAR = .2, ignoreCols = "id")
head(dat_mis)
summary(dat_mis)


###
###
###

out_sim <- mHMM(s_data = dat,
                data_distr = 'categorical',
                gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                start_val = c(list(gamma), emiss_distr),
                #emiss_hyp_prior = manual_prior_emiss,
                mcmc = list(J = 1000, burn_in = 200))

out_sim
summary(out_sim)

obtain_gamma(out_sim, level = 'group')
obtain_gamma(out_sim, level = 'subject')

obtain_emiss(out_sim, level = 'group')
obtain_emiss(out_sim, level = 'subject')


states1 <- vit_mHMM(s_data = dat, object = out_sim)
head(states1)
summary(states1)

## run with missing
## doesnt run with missing!!

out_sim_m <- mHMM(s_data = dat_mis,
                data_distr = 'categorical',
                gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                start_val = c(list(gamma), emiss_distr),
                #emiss_hyp_prior = manual_prior_emiss,
                mcmc = list(J = 1000, burn_in = 200))

out_sim_m
summary(out_sim_m)

obtain_gamma(out_sim_m, level = 'group')
obtain_gamma(out_sim_m, level = 'subject')

obtain_emiss(out_sim_m, level = 'group')
obtain_emiss(out_sim_m, level = 'subject')


states1_m <- vit_mHMM(s_data = dat_mis, object = out_sim_m)
head(states1_m)
summary(states1_m)

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
fm <- fit(mod)
summary(fm)

fm2 <- multistart(mod, nstart=10, initIters=50000)
summary(fm2)


# run with missing

mod_m <- depmix(list(v1~1,v2~1),data=dat_mis,
              nstates=3,
              family=list(multinomial("identity"),
                          multinomial("identity")),
              ntimes= as.numeric(table(dat$id)) )
# print the model, formulae and parameter values
mod_m

# fit the model by calling fit
fm_m <- fit(mod_m)
summary(fm_m)

fm2_m <- multistart(mod_m, nstart=10, initIters=50000)
summary(fm2_m)

sts_m <- viterbi(fm2_m)

dface <- dist(dat)
cluster.stats(d=dface,
              clustering = data1$states[,"state"],
              alt.clustering =  sts_m$state)$corrected.rand

table(data1$states[,"state"], 
      sts_m$state)
