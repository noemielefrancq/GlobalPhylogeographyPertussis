#######################################################
## Figure 2: Proportion of pairs with MRCA<2y
##   as a function of location population size
#######################################################
### Author: Noemie Lefrancq
### Last modification: 28/01/2022
#######################################################


## Load necessary packages
library(ape)
library(stringr)
library(lubridate)
library(REdaS)
library(questionr)



## Load tree
tree = read.nexus('Data/Beast_tree_MCC_26012022.tree')



## Create the global master dataframes
case.data_all = read.csv(file = 'Data/Metadata_analyses_26012022.csv', sep = ';')
colnames(case.data_all)[2] = 'Isolates'
tree_tmp = tree
tree_tmp = drop.tip(tree_tmp, tip = 'Tohama-I')
name.isolates = sapply(tree_tmp$tip.label, function(x)strsplit(x,split="_")[[1]][[1]]) #Isolate numbers
master.dat = case.data_all[match(name.isolates, case.data_all$Isolates),]
master.seq.names = master.dat$Isolates


######################################################################
## Compute matrices 
######################################################################


## Time differences between isolates
time_mat = abs(outer(master.dat$Collection.Year,master.dat$Collection.Year,"-"))
colnames(time_mat) = master.seq.names
rownames(time_mat) = master.seq.names
diag(time_mat)<-NA


## Genetic distances between isolates 
dist.mat<-cophenetic.phylo(tree)
seq.names<-row.names(dist.mat)
b<-match(as.character(master.dat$Isolates), seq.names)
b<-b[which(is.na(b)==F)]

gene_mat = matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))[b,b]
colnames(gene_mat) = master.seq.names
rownames(gene_mat) = master.seq.names
diag(gene_mat)<-NA

MRCA_mat = (gene_mat - time_mat)/2 
MRCA_mat[which(MRCA_mat == 0)] = 1E-6



## Matrix same region (1: pair from the same region, 0: pair NOT from the same region)
master.dat_tmp = master.dat
geo_mat_region_world<-matrix(0,length(master.dat_tmp$Country),length(master.dat_tmp$Country))
master.dat_tmp$Region = as.factor(master.dat_tmp$Region)
for (i in levels(master.dat_tmp$Region)){
  n = which(master.dat_tmp$Region == i)
  geo_mat_region_world[n,n] = 1
}
colnames(geo_mat_region_world) = master.seq.names
rownames(geo_mat_region_world) = master.seq.names
diag(geo_mat_region_world)<-NA


######################################################################





#####################################################################################
### First computation: Average proportion of pairs with MRCA<2y per region
#####################################################################################
## MRCA to define a transmission chain in years
threshold_MRCA = 2
## Minimum of sequences to have per region to peform this analysis
threshold_seqs = 2

## Names of regions
regions = levels(as.factor(master.dat$Region))[-1]

## Create list to store the proportions
boot.out = rep(NA, nboot*nsim)

## Sequences on which to perform the computation
A = which(is.na(master.dat$Region) == F)
nseq = length(A)

## Corresponding matrices to consider
geo_mat = geo_mat_region_world[A,A]
time_mat2 = time_mat[A,A]<=1
MRCA_mat2 = MRCA_mat[A,A]

## Bootstrap to create the confidence intervals
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  MRCA_mat = MRCA.mats[,,j]
  MRCA_mat2 = MRCA_mat[A,A]
  for (i in (1:(nboot))){
    tmp = sample(nseq, nseq, replace = T) ## resample sequences at each bootstrap iteration
    
    geo_mat_tmp = geo_mat[tmp, tmp]
    time_mat2_tmp = time_mat2[tmp,tmp]
    MRCA_mat2_tmp = MRCA_mat2[tmp,tmp]
    
    mat_tmp = (MRCA_mat2_tmp<=threshold_MRCA)*geo_mat_tmp*time_mat2_tmp
    mat_tmp[which(mat_tmp == 0)] = NA
    mat_tmp2 = geo_mat_tmp*time_mat2_tmp
    mat_tmp2[which(mat_tmp2 == 0)] = NA
    
    boot.out[(j-1)*nboot + i] = sum(mat_tmp, na.rm = T)/sum(mat_tmp2, na.rm = T)
  }
}

## Compute the confidence intervals
boot.ci1 = quantile(boot.out, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m = quantile(boot.out, probs = c(0.5), na.rm = T)

## Store parameters and results
res = list('boot.ci.m' = boot.ci.m,
           'boot.ci1' = boot.ci1,
           'nboot' = nboot,
           'MRCA' = '2y')

## Save results
saveRDS(file = 'Proportion_MRCA/Proportion_pairs_2y_regions_worldwide.rds', res)
#####################################################################################





#####################################################################################
## Main computation: proportion of pairs with MRCA<2y, 
##             location by population size                    
#####################################################################################
## Data on population size, for each region
pop_size = read.csv2('Data/Pop_size_regions.csv')
a = match(master.dat$Region, pop_size$Regions)
master.dat$pop_size_regon = pop_size$Pop_size[a]

## Number of bootstrap iterations
nboot = 200

## Windows of population sizes (on a log scale)
Pmax = c(seq(5.5, 6.8, 0.1), 7.7)
win = 0.4
Pmin = Pmax-win
Pmin[length(Pmin)] = Pmax[length(Pmax)-1]
Pmin[1] = 4
Pmax = 10^Pmax
Pmin = 10^Pmin
pmid1<-(Pmin+Pmax)/2

## Number of sequences to resample per bin
n_seq_sampling = 500

## MRCA to define a transmission chain in years
threshold_MRCA = 2

## Create a matrices to store results
mean_pop_size = matrix(NA, 4, length(pmid1))
boot.out_all = NULL
boot.out = matrix(NA, length(pmid1), nboot)

## Compute the proportions, for each population size bin, for each bootstrap iteration
for(i in 1:length(pmid1)){
  a = which(master.dat$Collection.Year >= 1900 & 
              (master.dat$pop_size_regon >= Pmin[i] & master.dat$pop_size_regon < Pmax[i]))
  nseq = length(a)

  geo_mat = geo_mat_region_world[a,a]
  time_mat2 = time_mat[a,a]<=1
  MRCA_mat2 = MRCA_mat[a,a]
  
  mean_pop_size[1,i] = mean(master.dat$pop_size_regon[a])
  
  ##Bootstrap to create the confidence intervals
  for (j in 1:nboot){
    print(paste0('nboot : ', j, '/', nboot+1))
    tmp = sample(nseq, nseq, replace = T)
    
    geo_mat_tmp = geo_mat[tmp, tmp]
    time_mat2_tmp = time_mat2[tmp,tmp]
    MRCA_mat2_tmp = MRCA_mat2[tmp,tmp]
    
    mat_tmp = (MRCA_mat2_tmp<=threshold_MRCA)*geo_mat_tmp*time_mat2_tmp
    mat_tmp[which(mat_tmp == 0)] = NA
    mat_tmp2 = geo_mat_tmp*time_mat2_tmp
    mat_tmp2[which(mat_tmp2 == 0)] = NA
    
    boot.out[i,j] = sum(mat_tmp, na.rm = T)/sum(mat_tmp2, na.rm = T)
  }
}

## Compute the confidence intervals
boot.ci1_worldwide = apply(boot.out, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_worldwide = apply(boot.out, 1, quantile, probs = c(0.5), na.rm = T)
boot.out_all[[1]] = boot.out

## Store results in list
boot.ci1_all = list(boot.ci1_worldwide)
boot.ci.m_all = list(boot.ci.m_worldwide)

## Store paramters
boot.out_all$mean_pop_size = mean_pop_size
boot.out_all$nboot = nboot
boot.out_all$pmid1 = pmid1

## Save results
saveRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_m.rds', boot.ci1_all)
saveRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_ci.rds', boot.ci.m_all)
saveRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_bootstraps.rds', boot.out_all)
#####################################################################################



#####################################################################################
## Plot results
#####################################################################################
grDevices::windows(width = 10, height = 5)
par(mai = c(0.8, 0.8, 0.8, 0.1), mfrow = c(1,2))

boot.ci1_all = readRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_m.rds')
boot.ci.m_all = readRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_ci.rds')
boot.out_all = readRDS(file = 'Proportion_MRCA/Proportion_pairs_population_2y_bootstraps.rds')

nboot = boot.out_all$nboot
mean_pop_size = boot.out_all$mean_pop_size
pmid1 = boot.out_all$pmid1

estimates = matrix(NA, 1, nboot)
model_fit_per_bootstrap = NULL
for(i in 1:1){
  model_fit_per_bootstrap[[i]] = matrix(NA, length(seq(1E5, 1E8, length.out = 1000)), nboot)
  for(j in 1:nboot){
    y = 1/as.vector(boot.out_all[[i]][which(is.na(mean_pop_size[i,]) == F),j])
    y[which(y > 100000)] = NA
    x = rep(mean_pop_size[i,which(is.na(mean_pop_size[i,]) == F)], nboot)
    data = data.frame('x' = log(x), 'y' = log(y))
    mod = lm(y~x, data = data)
    
    estimates[i,j] = mod$coefficients[2]
    newdata = data.frame('x' = log(seq(1E5, 1E8, length.out = 1000)))
    model_fit_per_bootstrap[[i]][,j] = exp(predict.lm(mod, newdata = newdata))
  }
}


## Plot proportions
plot(NULL, pch = 20, ylim = c(0,0.4), xlim = c(2E5, 2E7), bty = 'n', log = 'x',
     xlab = 'Population size', ylab = 'Transmission chains, MRCA<2y', 
     main = '',
     yaxt = 'n')
axis(2, las = 2)

i=1
fit_ci = apply(model_fit_per_bootstrap[[i]], 1, quantile, probs = c(0.025, 0.975), na.rm = T)
fit_m = apply(model_fit_per_bootstrap[[i]], 1, quantile, probs = c(0.5), na.rm = T)

points(mean_pop_size[i,], boot.ci.m_all[[i]], pch = 20, col = 'black')
lines(seq(1E5, 1E8, length.out = 1000), 1/fit_m, lty = 2, col = 'black')
polygon(x = c(seq(1E5, 1E8, length.out = 1000), rev(seq(1E5, 1E8, length.out = 1000))),
        y = c(1/fit_ci[1,],rev(1/fit_ci[2,])),
        col = adjustcolor('black', alpha.f = 0.2), border = F)

## Plot transmission chains
plot(NULL, pch = 20, ylim = c(3, 200), xlim = c(2E5, 2E7), bty = 'n', log = 'xy',
     xlab = 'Population size', ylab = 'Transmission chains, MRCA<2y', 
     main = '',
     yaxt = 'n')
axis(2, las = 2)

i=1
fit_ci = apply(model_fit_per_bootstrap[[i]], 1, quantile, probs = c(0.025, 0.975), na.rm = T)
fit_m = apply(model_fit_per_bootstrap[[i]], 1, quantile, probs = c(0.5), na.rm = T)

points(mean_pop_size[i,], 1/boot.ci.m_all[[i]], pch = 20, col = 'black')
lines(seq(1E5, 1E8, length.out = 1000), fit_m, lty = 2, col = 'black')
polygon(x = c(seq(1E5, 1E8, length.out = 1000), rev(seq(1E5, 1E8, length.out = 1000))),
        y = c(fit_ci[1,],rev(fit_ci[2,])),
        col = adjustcolor('black', alpha.f = 0.2), border = F)

#####################################################################################


