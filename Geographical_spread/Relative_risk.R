########################################################################################
## Relative risk of pairs of sequences to come from the same location versus different locations
########################################################################################
library(ape)
library(stringr)
library(lubridate)
library(REdaS)
library(questionr)

##Load the data--------------------------------------------------------------------------------------
##Tree
trees= readRDS('Data/Beast_10_posterior_trees_26012022.rds')

##Create the global master dataframes and lists------------------------------------------------------
case.data_all = read.csv(file = 'Data/Metadata_analyses_26012022.csv', sep = ';')
colnames(case.data_all)[2] = 'Isolates'
tree_tmp = trees[[1]]
tree_tmp = drop.tip(tree_tmp, tip = 'Tohama-I')
name.isolates = sapply(tree_tmp$tip.label, function(x)strsplit(x,split="_")[[1]][[1]]) #Isolate numbers
master.dat = case.data_all[match(name.isolates, case.data_all$Isolates),]
master.seq.names = master.dat$Isolates

##Compute matrices ----------------------------------------------------------------------------------

########################  TIME  ####################################
#Compute time distances
time_mat = abs(outer(master.dat$Collection.Year,master.dat$Collection.Year,"-"))
colnames(time_mat) = master.seq.names
rownames(time_mat) = master.seq.names
diag(time_mat)<-NA
####################################################################

######################## GENETIC ###################################
##Store genetic distance matrices into a huge matrix################
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
  
  MRCA_mat = (gene_mat - time_mat)/2 #make the genetic matrix "independant" of the difference in the isolation dates 
  MRCA_mat[which(MRCA_mat == 0)] = 1E-6
  
  MRCA.mats[,,ii]<- MRCA_mat
  print(paste0(ii, ' / ', nsim))
}
rm(trees)
######################################################################

######################## GEOGRAPHY ###################################
#Matrix region, world
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

#Matrix States in the US
a = length(which(master.dat$Country=='US'))
master.dat_tmp = master.dat[which(master.dat$Country=='US'),]
geo_mat_states_US<-matrix(0,a,a)
master.dat_tmp$Region = as.factor(master.dat_tmp$Region)
for (i in levels(master.dat_tmp$Region)){
  n = which(master.dat_tmp$Region == i)
  geo_mat_states_US[n,n] = 1
}
colnames(geo_mat_states_US) = master.seq.names[which(master.dat$Country=='US')]
rownames(geo_mat_states_US) = master.seq.names[which(master.dat$Country=='US')]
diag(geo_mat_states_US)<-NA

#Matrix countries
geo_mat_country<-matrix(0,length(master.dat$Country),length(master.dat$Country))
master.dat$Country = as.factor(master.dat$Country)
for (i in levels(master.dat$Country)){
  n = which(master.dat$Country == i)
  geo_mat_country[n,n] = 1
}
colnames(geo_mat_country) = master.seq.names
rownames(geo_mat_country) = master.seq.names
diag(geo_mat_country)<-NA

#Matrix continent
geo_mat_continent<-matrix(0,length(master.dat$Continent),length(master.dat$Continent))
master.dat$Continent = as.factor(master.dat$Continent)
for (i in levels(master.dat$Continent)){
  n = which(master.dat$Continent == i)
  geo_mat_continent[n,n] = 1
}
colnames(geo_mat_continent) = master.seq.names
rownames(geo_mat_continent) = master.seq.names
diag(geo_mat_continent)<-NA

######################################################################



######################################################################
## Compute the relative risks                                        #   
######################################################################
### Parameters for the RR plots, in "continuous" MRCA intervals
# Time to MRCA from 0 to 10 years
Pmax = seq(1, 20, 0.1)
windows = 1.5
Pmin = Pmax - windows
Pmin[which(Pmin<0)]<-0
pmid1<-(Pmin+Pmax)/2

nboot = 10

## Function to compute the RR
ratio_bootstrap <- function(x, geo_mat, time_mat2, MRCA_mat){
  geo_mat.tmp = geo_mat[x,x]
  time_mat.tmp = time_mat2[x,x]
  MRCA_mat.tmp2 = MRCA_mat[x,x] 
  
  tmp = geo_mat.tmp * time_mat.tmp
  tmp[which(tmp == 0)] = NA
  
  tmp2 = (1-geo_mat.tmp) * time_mat.tmp
  tmp2[which(tmp2 == 0)] = NA
  
  a1 = cumsum(hist(tmp*MRCA_mat.tmp2, breaks = c(0,Pmax,1E10), plot = F)$counts)
  a2 = cumsum(hist(tmp*MRCA_mat.tmp2, breaks = c(0,Pmin,1E10), plot = F)$counts)
  a = a1 - a2
  
  b = sum(tmp,na.rm=T)
  
  c1 = cumsum(hist(tmp2 * MRCA_mat.tmp2, breaks = c(0,Pmax,1E20), plot = F)$counts)
  c2 = cumsum(hist(tmp2 * MRCA_mat.tmp2, breaks = c(0,Pmin,1E20), plot = F)$counts)
  c = c1 - c2
  
  d = sum(tmp2,na.rm=T)
  
  rr = (a/b)/(c/d) 
  
  return(rr[-length(rr)])
}

###########################################################################################################################################################
## 1 : Continuous Plot RR in Europe by country
###########################################################################################################################################################
# Create relative risk matrix
boot.out = matrix(NA, length(pmid1), nsim*nboot)

## Data on which we are computing the RR
A = c(which(master.dat$Dataset == 'Pasteur_Brisse+EU'), which(master.dat$Dataset == 'EU_dataset'))
nseq = length(A)
##Bootstrap to create the ci
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_mat = geo_mat_country[A,A]
  time_mat2 = time_mat[A,A]<=1
  MRCA_mat = MRCA.mats[,,j]
  MRCA_mat2 = MRCA_mat[A,A]
  for (i in (1:(nboot))){
    tmp = sample(nseq, nseq, replace = T)
    rr = ratio_bootstrap(tmp, geo_mat, time_mat2, MRCA_mat2)
    boot.out[,(j-1)*nboot + i] = rr
  }
}

##Compute the ci over the 2 to nboot +1 columns 
boot.ci1_EU = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_EU = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.5), na.rm = T)

list  = list('pmid' = pmid1, 
             'median' = boot.ci.m_EU,
             'ci' = boot.ci1_EU)
saveRDS(list, 'Relative_risk/RR_Pertussis_country_europe.rds')

###########################################################################################################################################################
## 2 : Continuous Plot RR in USA by state
###########################################################################################################################################################
# Create relative risk matrix
boot.out = matrix(NA, length(pmid1), nsim*nboot)

## Data on which we are computing the RR
A = which(master.dat$Country == 'US' & master.dat$Collection.Year >= 2010)
nseq = length(A)
B = which(master.dat[which(master.dat$Country == "US"),]$Collection.Year>= 2010)
##Bootstrap to create the ci
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_mat = geo_mat_states_US[B,B]
  time_mat2 = time_mat[A,A]<=1
  MRCA_mat = MRCA.mats[,,j]
  MRCA_mat2 = MRCA_mat[A,A]
  for (i in (1:(nboot))){
    tmp = sample(nseq, nseq, replace = T)
    rr = ratio_bootstrap(tmp, geo_mat, time_mat2, MRCA_mat2)
    boot.out[,(j-1)*nboot + i] = rr
  }
}

##Compute the ci over the 2 to nboot +1 columns
boot.ci1_US = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_US = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.5), na.rm = T)

list  = list('pmid' = pmid1, 
             'median' = boot.ci.m_US,
             'ci' = boot.ci1_US)
saveRDS(list, 'Relative_risk/RR_Pertussis_US_states.rds')

###########################################################################################################################################################
## 3 : Plot RR between continents (average across different continents)
###########################################################################################################################################################
# Create relative risk matrix
boot.out = matrix(NA, length(pmid1), nsim*nboot)

## Data on which we are computing the RR
A = c(which(master.dat$Collection.Year>=2010))
nseq = length(A)
nb_seq_per_continent = 300
continents = c(levels(as.factor(master.dat$Continent)))[-1] ## Remove Africa:  not enough sequences

##Bootstrap to create the ci
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_mat = geo_mat_continent[A,A]
  time_mat2 = time_mat[A,A]<=1
  MRCA_mat = MRCA.mats[,,j]
  MRCA_mat2 = MRCA_mat[A,A]
  for (i in (1:(nboot))){
    ## resample isolates, making sure to have equal representation of all continents
    tmp = NULL
    for(k in 1:length(continents)){
      b = which(master.dat$Continent[which(master.dat$Collection.Year>=2010)] == continents[k])
      tmp = c(tmp, sample(b, nb_seq_per_continent, replace = T))
    }
    # tmp = sample(nseq, nseq, replace = T) ## old
    rr = ratio_bootstrap(tmp, geo_mat, time_mat2, MRCA_mat2)
    boot.out[,(j-1)*nboot + i] = rr
  }
}

##Compute the ci over the 2 to nboot +1 columns 
boot.ci1_Continents = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_Continents = apply(boot.out[,(2:(nboot+1))], 1, quantile, probs = c(0.5), na.rm = T)

list  = list('pmid' = pmid1, 
             'median' = boot.ci.m_Continents,
             'ci' = boot.ci1_Continents)

saveRDS(list, 'Relative_risk/RR_Pertussis_average_continents.rds')
#####################################################################################






###########################################################################################################################################################
## Plot RR
###########################################################################################################################################################
US_states = readRDS('Relative_risk/RR_Pertussis_US_states.rds')
EU_countries = readRDS('Relative_risk/RR_Pertussis_country_europe.rds')
Continents = readRDS('Relative_risk/RR_Pertussis_average_continents.rds')

grDevices::windows(width = 10, height=5)
par(mfrow = c(1, 2),oma = c(0.1, 0.1, 0.1, 0.1), mai = c(1,1,1,1))

titles = c("Within countries in Europe / states in the US", "Within continents")

boot.ci.m.compiled = list(EU_countries$median, US_states$median, Continents$median)
boot.ci1_compiled = list(EU_countries$ci, US_states$ci, Continents$ci)

boot.ci.m.compiled[[4]][which(boot.ci.m.compiled[[4]] >= Inf)] = 1000
boot.ci1_compiled[[4]][which(boot.ci1_compiled[[4]] >= Inf)] = 1000
pmid1 = FR_regions$pmid

color = NULL
color[1] = 'royalblue2'
color[2] = 'red3'
color[3] = 'grey40'

plot(NULL, xlim = c(0,20), ylim=c(0.25, 1000),log="y",
     cex=1,xlab="Time to MRCA (years)",ylab="Relative risk", bty = 'n', yaxt = 'n',
     main = titles[1], cex.axis=1, cex.main = 1)
abline(h = 1, col = 'black', lty = 2, lwd = 1.2)
for (i in (1:2)){
  lines(pmid1[c(-190, -191)],boot.ci.m.compiled[[i]][c(-190, -191)], lwd = 1.5, col = adjustcolor(color[i], alpha.f = 0.7))
  axis(2, las = 2, at = c(0.25, 1, 10, 100, 1000), labels = c(0.25, 1, 10, 100, '>1000'))
  boot.ci1_compiled[[i]][which(boot.ci1_compiled[[i]] == 0)] = 0.001
  polygon(x = c(pmid1[c(-190, -191)], rev(pmid1[c(-190, -191)])), y = c(boot.ci1_compiled[[i]][1,c(-190, -191)], rev(boot.ci1_compiled[[i]][2,c(-190, -191)])), 
          col = adjustcolor(color[i], alpha.f = 0.3), border = F)
}

i=3
plot(NULL, xlim = c(0,20), ylim=c(0.25, 1000),log="y",
     cex=1,xlab="Time to MRCA (years)",ylab="Relative risk", bty = 'n', yaxt = 'n',
     main = titles[2], cex.axis=1, cex.main = 1)
abline(h = 1, col = 'black', lty = 2, lwd = 1.2)
lines(pmid1,boot.ci.m.compiled[[i]], lwd = 1.5, col = adjustcolor(color[i], alpha.f = 0.7))
axis(2, las = 2, at = c(0.25, 1, 10, 100, 1000), labels = c(0.25, 1, 10, 100, '>1000'))
boot.ci1_compiled[[i]][which(boot.ci1_compiled[[i]] == 0)] = 0.001
polygon(x = c(pmid1, rev(pmid1)), y = c(boot.ci1_compiled[[i]][1,], rev(boot.ci1_compiled[[i]][2,])), 
        col = adjustcolor(color[i], alpha.f = 0.5), border = F)

###########################################################################################################################################################

