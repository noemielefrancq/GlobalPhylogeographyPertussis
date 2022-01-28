#######################################################
## Figure 1I: Relation between genetic distance and geographical distance
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


## Load Tree
trees= readRDS('Data/Beast_10_posterior_trees_26012022.rds')


## Create the global master dataframes
case.data_all = read.csv(file = 'Data/Metadata_analyses_26012022.csv', sep = ';')
colnames(case.data_all)[2] = 'Isolates'
tree_tmp = trees[[1]]
tree_tmp = drop.tip(tree_tmp, tip = 'Tohama-I')
name.isolates = sapply(tree_tmp$tip.label, function(x)strsplit(x,split="_")[[1]][[1]]) #Isolate numbers
master.dat = case.data_all[match(name.isolates, case.data_all$Isolates),]
master.seq.names = master.dat$Isolates


######################################################################
## Compute matrices 
######################################################################

##  Time differences between isolates
time_mat = abs(outer(master.dat$Collection.Year,master.dat$Collection.Year,"-"))
colnames(time_mat) = master.seq.names
rownames(time_mat) = master.seq.names
diag(time_mat)<-NA

## Genetic distances between isolates 
# Store genetic distances into a list, 
# which corresponds to each tree of the posterior distribution
nsim = 10
ntrees = length(trees)
MRCA.mats<-array(NaN,c(length(master.seq.names),length(master.seq.names),nsim))

for (ii in 1:nsim){
  select.tree<-ii
  tree_tmp = trees[[select.tree]]
  tree_tmp = drop.tip(tree_tmp, tip = 'Tohama-I_1954')
  dist.mat<-cophenetic.phylo(tree_tmp)
  seq.names<-row.names(dist.mat)
  b<-match(as.character(master.dat$Isolates), seq.names)
  b<-b[which(is.na(b)==F)]
  
  gene_mat = matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))[b,b]
  colnames(gene_mat) = master.seq.names
  rownames(gene_mat) = master.seq.names
  diag(gene_mat)<-NA
  
  MRCA_mat = (gene_mat - time_mat)/2
  MRCA_mat[which(MRCA_mat == 0)] = 1E-6
  
  MRCA.mats[,,ii]<- MRCA_mat
  print(paste0(ii, ' / ', nsim))
}
rm(trees)

## Matrix same countries (1:pair from the same country, 0: pair NOT from the same country)
geo_mat_country<-matrix(0,length(master.dat$Country),length(master.dat$Country))
master.dat$Country = as.factor(master.dat$Country)
for (i in levels(master.dat$Country)){
  n = which(master.dat$Country == i)
  geo_mat_country[n,n] = 1
}
colnames(geo_mat_country) = master.seq.names
rownames(geo_mat_country) = master.seq.names
diag(geo_mat_country)<-NA

## Geographical distances between isolates 
get_dist = function(x1,x2,y1,y2){ 
  ## function to compute geographical distances 
  ## between 2 points of coordinates 
  ## (x1,y1) and (x2, y2)
  R=6371000 ## radius of Earth
  x1 = deg2rad(x1)
  x2 = deg2rad(x2)
  y1 = deg2rad(y1)
  y2 = deg2rad(y2)
  delta.x = x2 - x1
  delta.y = y2 - y1
  a = sin(delta.x/2.0)**2+ cos(x1)*cos(x2)*sin(delta.y/2.0)**2
  c=2*atan2(sqrt(a),sqrt(1-a))
  
  self.meters=R*c 
  km=self.meters/1000.0 # output distance in kilometers
  
  return(km)
}
geo_mat_km<-matrix(0,length(master.seq.names),length(master.seq.names))
master.dat$Longitude = as.numeric(gsub(',', '.', as.character(master.dat$Longitude)))
master.dat$Latitude = as.numeric(gsub(',', '.', as.character(master.dat$Latitude)))
for (i in 1:length(master.seq.names)){
  geo_mat_km[i,] = get_dist(rep(master.dat$Latitude[which(master.dat$Isolates == master.seq.names[i])], length(master.seq.names)), master.dat$Latitude, rep(master.dat$Longitude[which(master.dat$Isolates == master.seq.names[i])], length(master.seq.names)), master.dat$Longitude)
}
diag(geo_mat_km)<-NA
######################################################################





######################################################################
## Main computation: compute geographical distance for pairs of various genetic distance 
######################################################################
## Number of bootstrap iterations
nboot = 10
sampling_size_global = 1000 # number of sequences to sample at each iteration

## Genetic distances to consider (in years)
Pmax = c(seq(0.25,5,0.25), seq(6,40,2)) #

## Create matrix to store results
boot.out = list()
boot.out$'Global_weighted_country'=  matrix(NA, nrow = length(Pmax), ncol = nboot*nsim)

## Weight of each sequence, to allow for even sampling by location
a = c(which(master.dat$Precision_loc == 'State' | master.dat$Precision_loc == 'City' | master.dat$Precision_loc == 'Region'))
w = 1/(freq(master.dat$Region[a])$n/sum(freq(master.dat$Region[a])$n)) 
l = levels(as.factor(master.dat$Region[a]))
w_vec = w[match(master.dat$Region[a], l)]

## Select the geo matrix to work on
nseq = length(a)
geo_mat_country_tmp = geo_mat_country[a,a]
geo_mat_country_tmp[which(geo_mat_country_tmp == 0)] = NA
geo_mat = geo_mat_km[a,a]*geo_mat_country_tmp

## Boostrap loop, to compute each time the mean distance between pairs of isolates, for each MRCA window
for(k in 1:nsim){
  print(paste0('K : ', k, '/', nsim))
  for(j in 1:length(Pmax)){
    MRCA_mat2 = MRCA.mats[,,k] # Choose the genetic distances to consider, which corresponds to one tree of the posterior distribution
    MRCA_mat2 = MRCA_mat2[a,a] 
    for(i in 1:nboot){
      tmp = sample(nseq, sampling_size_global, prob = w_vec, replace = T) # Resample isolates
    
      MRCA_mat2.tmp = (MRCA_mat2<=Pmax[j]) # Select only pairs that verify the MRCA distance
      rr.out = mean(geo_mat[tmp,tmp]*MRCA_mat2.tmp[tmp,tmp], na.rm = T) # Compute mean geo distance
      boot.out$'Global_weighted_country'[j,(k-1)*nboot + i] = rr.out# Store result
    }
  }
}
######################################################################

######################################################################
## Store parameters in the list
######################################################################
boot.out$Pmax = Pmax 
boot.out$nboot = nboot
boot.out$sampling_size_global = sampling_size_global
######################################################################



######################################################################
## Save the results
######################################################################
saveRDS(boot.out, 'Distance_MRCA/Mean_distance_MRCA_global.rds')
######################################################################



######################################################################
## Plot results
######################################################################
mean.and.ci <-function(v){
  ## Function to compute mean and confidence intervas, from the bootstrap distributions
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

## Plot
grDevices::windows(width = 7, height = 7)
par(mai = c(1,1,1,1), mfrow = c(1,1))
a = apply(boot.out$"Global_weighted_country", MARGIN = 1, function(x)mean.and.ci(x))
Pmax = boot.out$Pmax
plot(Pmax, a[1,], xlim = c(0,6), ylim = c(0.25,100), log = 'y', type = 'l', lwd = 2, bty = 'n', yaxt = 'n',
     xlab = 'Evolutionary time (years)',
     ylab = 'Spatial distance (km)')
axis(2, las = 2)
polygon(x = c(Pmax, rev(Pmax)), y = c(a[2,], rev(a[3,])), border = F, col = adjustcolor('grey', alpha.f = 0.5))
######################################################################
