## MCMC fitness
## By country, by genotype 
library(questionr)
library(treeio)
library(ape)
library(rstan)
library(gridExtra)
library(grid)
library(gtable)
library(doParallel)  
library(loo)
rstan_options(auto_write = TRUE)

############# R(t) model
model.MCMC <- stan_model(file = 'Estimate_relative_fitness/Model_relative_fitness_genotypes_vaccine_switch_250421.stan')

############# data for MCMC ######################################################
data.MCMC = readRDS('Estimate_relative_fitness/Data_model_allcountries_refgeno5.rds')

############# parameters vaccination #############################################
data.MCMC$R_every_pre_vacc = 12
data.MCMC$number_R_pre_vacc = 1; 

data.MCMC$R_every_post_vacc = 12
data.MCMC$number_R_post_vacc = 1;

## Vaccine introduction (vector of length countries)
data.MCMC$yearF0 = data.MCMC$WCV_introduction
  
## ACV introduction (vector of length countries)
data.MCMC$yearIntroduction = data.MCMC$booster_introduction

##################################################################################
## Run MCMC 
##################################################################################
name_file = 'Estimate_relative_fitness/Output_WCV_booster_ref3'

no_cores = 3
registerDoParallel(cores=no_cores)  
cl = makeCluster(no_cores) 

seed <- floor(runif(1, min = 1, max = 1E6))
# seed = 124
print(paste0('seed = ', seed))

f0_init = function(nb_countries, nb_geno){
  res = matrix(0, ncol = nb_geno, nrow = nb_countries)
  res[,1] = rnorm(nb_countries, mean = 1-nb_geno*0.08, sd = 0.01)
  for(i in 2:nb_geno){
    res[,i] = rnorm(nb_countries, mean = 0.08, sd = 0.01)
  }
}

foreach(i = 1:3)  %dopar% {
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC, 
                        show_messages = TRUE, 
                        chains = 1, cores = 1,iter= 2000, chain_id = i,
                        control = list(adapt_delta = 0.97, max_treedepth = 13), 
                        init = list(list(alpha = rnorm(1, 0.9, 0.02), 
                                         f0 = f0_init(data.MCMC$nb_countries,6),
                                         overdispersion_inv = rnorm(1, 0.75, 0.1), 
                                         fitness_genotypes_pre_vacc = matrix(rnorm(5*data.MCMC$number_R_pre_vacc,0,0.1), ncol = data.MCMC$number_R_pre_vacc, nrow = 5),
                                         fitness_genotypes_post_vacc = matrix(rnorm(5*data.MCMC$number_R_post_vacc,0,0.1), ncol = data.MCMC$number_R_post_vacc, nrow = 5))))  # iter =10000 et plus de chaines
  
  fit = list(fit=fit_delay,
             data= data.MCMC)
  

  log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik_1), cores = 2)
  loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
  waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
  print(loo_1)
  print(waic_1)
  

  Chains=rstan::extract(fit$fit)
  saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
  saveRDS(data.MCMC, file = paste0(name_file, '_data_', i, '.rds'))
  
  
  m = monitor(fit$fit, print = F)
  fit$monitor = m
  saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
}

print('Reading 1')
fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
print('Reading 2')
fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
print('Reading 3')
fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))

print('Fit')
fit = NULL
fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))

fit$data = fit1$data

print('Chains')
Chains = rstan::extract(fit$fit)

print('Writing fit')
saveRDS(fit, paste0(name_file, '_fit_all.rds'))

print('Writing chains')
saveRDS(Chains, paste0(name_file, '_chains_all.rds'))

################################################################################
## COmpute WAIC and LOO of the model
################################################################################
library(loo)
log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
print(waic_1)
################################################################################