########################################################################################
## Compute the proportion of pairs by MRCA bin
########################################################################################
library(ape)
library(stringr)
library(lubridate)
library(REdaS)
library(questionr)

##Load the data--------------------------------------------------------------------------------------
##Tree
tree = read.nexus('Data/Beast_tree_MCC_26012022.tree')

##Create the global master dataframes and lists------------------------------------------------------
case.data_all = read.csv(file = 'Data/Metadata_analyses_26012022.csv', sep = ';')
colnames(case.data_all)[2] = 'Isolates'
tree_tmp = trees[[1]]
tree_tmp = drop.tip(tree_tmp, tip = 'Tohama-I')
name.isolates = sapply(tree_tmp$tip.label, function(x)strsplit(x,split="_")[[1]][[1]]) #Isolate numbers
master.dat = case.data_all[match(name.isolates, case.data_all$Isolates),]
master.seq.names = master.dat$Isolates

##Compute matrices ----------------------------------------------------------------------------------

########################  TIME  ###################################
#Compute time distances
time_mat = abs(outer(master.dat$Collection.Year,master.dat$Collection.Year,"-"))
colnames(time_mat) = master.seq.names
rownames(time_mat) = master.seq.names
diag(time_mat)<-NA
######################################################################


######################## GENETIC ###################################
#Genetic distance matrix
dist.mat<-cophenetic.phylo(tree)
seq.names<-row.names(dist.mat)
b<-match(as.character(master.dat$Isolates), seq.names)
b<-b[which(is.na(b)==F)]

gene_mat = matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))[b,b]
colnames(gene_mat) = master.seq.names
rownames(gene_mat) = master.seq.names
diag(gene_mat)<-NA

MRCA_mat = (gene_mat - time_mat)/2 #make the genetic matrix "independant" of the difference in the isolation dates 
MRCA_mat[which(MRCA_mat == 0)] = 1E-6
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
## Compute the proportion of pairs by MRCA bin                       #   
######################################################################
breaks = c(0, 2, 5, 10, 20, 1E10)
nb_cat = length(breaks)-1
loc = rep(c('BSame region', 'CSame country', 'DSame state', 'ESame continent'), each = nb_cat)
cat = c(rep(rev(LETTERS[1:nb_cat]), times = 4))
freq = rep(NA, nb_cat*4)

## Regions, averaged across different countries
a = which((master.dat$Precision_loc == "City" | master.dat$Precision_loc == "Province" | 
            master.dat$Precision_loc == "Region") & 
            master.dat$Collection.Year >= 2000)
master.dat_tmp = master.dat[a,]
c = levels(as.factor(master.dat_tmp$Country))
geo_mat = geo_mat_region_world[a,a]
time_mat2 = time_mat[a,a]<=1
MRCA_mat2 = MRCA_mat[a,a]
res_prop_region_per_country = matrix(NA, nrow = length(c), ncol = nb_cat) ## For supplementary
for(i in 1:length(c)){
  tmp = which(master.dat_tmp$Country == c[i])
  if(length(tmp) > 0){
    mat_tmp = MRCA_mat2[tmp,tmp]*geo_mat[tmp,tmp]*time_mat2[tmp,tmp]
    mat_tmp[which(mat_tmp == 0)] = NA
    res_prop_region_per_country[i,] = hist(mat_tmp, breaks = breaks, plot = F)$count
    res_prop_region_per_country[i,] = res_prop_region_per_country[i,]/(sum(res_prop_region_per_country[i,]))
  }
}
freq[(1:nb_cat)+(0*nb_cat)] = apply(res_prop_region_per_country, MARGIN = 2, function(x)mean(x, na.rm = T))*1000

## EU countries
a = c(which(master.dat$Dataset == 'Pasteur_Brisse+EU'), which(master.dat$Dataset == 'EU_dataset'))
geo_mat = geo_mat_country[a,a]
time_mat2 = time_mat[a,a]<=1
MRCA_mat2 = MRCA_mat[a,a]
mat_tmp = MRCA_mat2*geo_mat*time_mat2
mat_tmp[which(mat_tmp == 0)] = NA
freq[(1:nb_cat)+(1*nb_cat)] = hist(mat_tmp, breaks = breaks, plot = F)$count

## US states
a = which(master.dat$Country == 'US' & master.dat$Collection.Year >= 2010)
b = which(master.dat[which(master.dat$Country == "US"),]$Collection.Year>= 2010)
geo_mat = geo_mat_states_US[b,b]
time_mat2 = time_mat[a,a]<=1
MRCA_mat2 = MRCA_mat[a,a]
mat_tmp = MRCA_mat2*geo_mat*time_mat2
mat_tmp[which(mat_tmp == 0)] = NA
freq[(1:nb_cat)+(2*nb_cat)] = hist(mat_tmp, breaks = breaks, plot = F)$count

## Continents, averaged across different continents
a = which(master.dat$Collection.Year >= 2010)
master.dat_tmp = master.dat[a,]
c = levels(as.factor(master.dat_tmp$Continent))
geo_mat = geo_mat_continent[a,a]
time_mat2 = time_mat[a,a]<=1
MRCA_mat2 = MRCA_mat[a,a]
res_prop_contient_per_continent = matrix(NA, nrow = length(c), ncol = nb_cat)
for(i in 1:length(c)){
  tmp = which(master.dat_tmp$Continent == c[i])
  if(length(tmp) > 0){
    mat_tmp = MRCA_mat2[tmp,tmp]*geo_mat[tmp,tmp]*time_mat2[tmp,tmp]
    mat_tmp[which(mat_tmp == 0)] = NA
    res_prop_contient_per_continent[i,] = hist(mat_tmp, breaks = breaks, plot = F)$count
    res_prop_contient_per_continent[i,] = res_prop_contient_per_continent[i,]/(sum(res_prop_contient_per_continent[i,]))
  }
}
freq[(1:nb_cat)+(3*nb_cat)] = apply(res_prop_contient_per_continent[c(2,3,5),], MARGIN = 2, function(x)mean(x, na.rm = T))*1000

######################################################################
## Save results
######################################################################
Data = data.frame(loc, cat, freq)
saveRDS(Data, 'Proportion_MRCA/Proportion_pairs_MRCA_per_different_locs.rds')

######################################################################
## Plot results
######################################################################
grDevices::windows(width = 5, height=5) 

Data = readRDS('Proportion_MRCA/Proportion_pairs_MRCA_per_different_locs.rds')
g = ggplot(Data, aes(x = loc, y = freq, fill = cat)) +
      geom_bar(position = "fill", stat="identity")+
      scale_fill_brewer(palette="Blues", name = 'MRCA', 
                        labels = rev(c(paste0('<', breaks[2], 'y'), 
                                   paste0(breaks[2:(length(breaks)-2)], 'y - ', breaks[3:(length(breaks)-1)], 'y'),
                                   paste0('>', breaks[length(breaks)-1], 'y'))))+
      scale_x_discrete(labels=c('BSame region' = 'Same region', 
                                'CSame country' = 'Same country', 
                                'DSame state' = 'Same state', 
                                'ESame continent' = 'Same continent'))+
      labs(x = 'Location', y = 'Proportion')+
      theme_classic()+
      theme(plot.title = element_text (face = 'bold',size = 12, hjust = 0.5),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.key.size = unit(0.25,"cm"),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_text(size = 12, angle = 45, hjust=1),
              axis.title.y = element_text(size = 12),
              axis.text.y = element_text(size = 12, hjust=0.5),
              strip.text.x = element_text(size = 12, colour = "black", angle = 0),
              strip.text.y = element_text(size = 12, colour = "black", angle = 0),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
g
######################################################################