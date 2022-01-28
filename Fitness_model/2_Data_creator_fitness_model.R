#######################################################
## Code to prepare data for the fitness model
#######################################################
### Author: Noemie Lefrancq
### Last modification: 28/01/2022
#######################################################


## Load necessary packages
library(questionr)
library(treeio)
library(ape)
library(rstan)
library(gridExtra)
library(grid)
library(gtable)
library(ggpubr)
library(binom)
library(stringr)

## Load metadata ###############################################
data = readRDS('Estimate_relative_fitness/data_tips_and_nodes_clades.rds') 
data_antigens = read.csv('Data/Metadata_analyses_26012022.csv', sep = ';')

## Load tree
tree.beast = read.beast('Data/Beast_tree_MCC_26012022.tree')
tree = tree.beast@phylo

## Get distance between nodes
mat.dist.nodes = dist.nodes(tree)
diag(mat.dist.nodes) = NA
nb_tips = length(tree$tip.label)

## Check the continent labels of each node
## For continents, do not consider nodes that are more than 5 years away from any sequence
threshold = 5 ## threshold of minimum distance to at least one tip
node_vect = (nb_tips+1):(nb_tips*2-1)
list_to_remove = NULL
list_nb = NULL
for (i in node_vect){
  a = which(mat.dist.nodes[i,1:nb_tips] < threshold)
  list_nb = c(list_nb, length(a))
  if (length(a) == 0){list_to_remove = c(list_to_remove, i)}
}

## Take out these values
data = data[-list_to_remove,]

## Common values
data$date = 2019-round(as.numeric(data$height))
max_date = 2030
min_date = min(data$date)
min_date = 1940
seasons = seq(min_date, max_date, 1)

## List of clades
clades_prn_pos = levels(as.factor(data$clade))[c(1,2,3)]
clades_prn_neg = as.numeric(levels(as.factor(data$clade))[-c(1,2,3)])[order(as.numeric(levels(as.factor(data$clade))[-c(1,2,3)]))]

### Data ACV introduction - to cut ACV PERIOD
vaccine_introduction =  c(1999, 1999, 1997, 2012, 2007, 1997, 2005, 2005, 2004, 1996, 2002, 2019, 1995, 2005, 1998, 1998, 2019, 2000, 1997, 1981)
vaccine_introduction_country = c('AU', 'BE', 'CA', 'CN', 'CZ', 'DK', 'ES', 'FI', 'FR', 'IE', 'IL', 'IR', 'IT', 'NL', 'NO', 'SE', 'TN', 'UK', 'US', 'JP')
cut_at_ACV = F

### Data ACV introduction - to cut WCV PERIOD
vaccine_introduction =  c(1999, 1999, 1997, 2012, 2007, 1997, 2005, 2005, 2004, 1996, 2002, 2019, 1995, 2005, 1998, 1998, 2019, 2000, 1997, 1981)
vaccine_introduction_country = c('AU', 'BE', 'CA', 'CN', 'CZ', 'DK', 'ES', 'FI', 'FR', 'IE', 'IL', 'IR', 'IT', 'NL', 'NO', 'SE', 'TN', 'UK', 'US', 'JP')
cut_at_WCV = F

####################################################################################################
## Divide those clades into country specific clades
####################################################################################################
###### Remove country state when posterior probability is <0.9 
data$country = data$Countries
a = which(as.numeric(data$Countries.prob) < 0.9)
data$country[a] = NA
data = data[which(!is.na(data$country)),]

data$clade = as.factor(data$clade)
t = table(data$clade)
a = which(t < 5)
a = as.data.frame(a)
m = match(data$clade, row.names(a))
data = data[which(is.na(m)),]

clades = levels(as.factor(as.character(data$clade)))
################################################
## Compute frequencies on a window of time
################################################
data$date = floor(as.numeric(data$date))
compute_clade_numbers_country_1 = function(country){
  data_tmp = data
  if(cut_at_ACV == T){
    data_tmp = data[which(data$date <= vaccine_introduction[which(vaccine_introduction_country == country)]),]
  }
  if(cut_at_WCV == T){
    data_tmp = data[which(data$date >= vaccine_introduction[which(vaccine_introduction_country == country)]),]
  }

  list_x_prn_neg = as.list(rep(NA, length(clades_prn_neg)))
  for (i in 1:length(list_x_prn_neg)) list_x_prn_neg[[i]] = rep(NA, length(seasons))
  list_x_prn_pos = as.list(rep(NA, length(clades_prn_pos)))
  for (i in 1:length(list_x_prn_pos)) list_x_prn_pos[[i]] = rep(NA, length(seasons))
  number_total = rep(NA, length(seasons)-2)
  
  for (i in 1:(length(seasons))){
    sum = 0
    # print(i)
    for (j in 1:length(list_x_prn_neg)) {
      list_x_prn_neg[[j]][i] =  round(length(which((data_tmp$date >= seasons[i] & 
                                                      data_tmp$date <= seasons[i]) &
                                                     (data_tmp$country == country | is.na(data_tmp$country) == T) &
                                                     data_tmp$clade == clades_prn_neg[j])), digits = 0)
      sum = sum + list_x_prn_neg[[j]][i] 
    }
    for (k in 1:length(list_x_prn_pos)) {
      list_x_prn_pos[[k]][i] =   round(length(which((data_tmp$date >= seasons[i] & 
                                                       data_tmp$date <= seasons[i]) & 
                                                      (data_tmp$country == country | is.na(data_tmp$country) == T) &
                                                      data_tmp$clade == clades_prn_pos[k])), digits = 0)
      sum = sum + list_x_prn_pos[[k]][i]
    }
    number_total[i] = round(sum, digits = 0)
  }
  list_clades_neg = clades_prn_neg
  return(list(list_x_prn_neg, list_x_prn_pos, list_clades_neg, number_total))
}

list_AU = compute_clade_numbers_country_1('AU')
list_BE = compute_clade_numbers_country_1('BE')
list_CA = compute_clade_numbers_country_1('CA')
list_CN = compute_clade_numbers_country_1('CN')
list_CZ = compute_clade_numbers_country_1('CZ')
list_DK = compute_clade_numbers_country_1('DK')
list_ES = compute_clade_numbers_country_1('ES')
list_FI = compute_clade_numbers_country_1('FI')
list_FR = compute_clade_numbers_country_1('FR')
list_IE = compute_clade_numbers_country_1('IE')
list_IL = compute_clade_numbers_country_1('IL')
list_IR = compute_clade_numbers_country_1('IR')
list_IT = compute_clade_numbers_country_1('IT')
list_JP = compute_clade_numbers_country_1('JP')
list_NL = compute_clade_numbers_country_1('NL')
list_NO = compute_clade_numbers_country_1('NO')
list_SE = compute_clade_numbers_country_1('SE')
list_TN = compute_clade_numbers_country_1('TN')
list_UK = compute_clade_numbers_country_1('UK')
list_US = compute_clade_numbers_country_1('US')

################################################
## Years (seasons) used
################################################
seasons_plot = seasons
nb_years = length(seasons)

################################################
## Add genotype info
################################################
compute_data_obs_country = function(list_x_prn_neg, list_x_prn_pos, list_clades_neg){
  data_obs = NULL
  ## Set prn - data
  for (i in 1:length(list_clades_neg)){
    a =  which(data$clade == as.character(list_clades_neg[i]))
    seqs = data$label[a]
    seqs = seqs[which(is.na(seqs) == F)]
    b = match(seqs, data_antigens$ï..Isolate)
    c = data_antigens[b,]
    fim = freq(as.factor(as.character(c$fim3)))
    fim = row.names(fim)[which(fim$`val%`>50)]
    if(length(fim) == 0) fim = NA
    ptxP = freq(as.factor(as.character(c$ptxP)))
    ptxp = row.names(ptxP)[which(ptxP$`val%`>50)]
    if(length(ptxp) == 0) ptxp = NA
    data_obs[[i]] = data.frame('Clade' = rep(list_clades_neg[i], length(seasons_plot)),
                               'count' = list_x_prn_neg[[i]],
                               'Year' = seasons_plot,
                               'Prn_exp' = rep(0, length(seasons_plot)), 
                               'ptxp_allele' = rep(ptxp, length(seasons_plot)), 
                               'Fim_3_allele' = rep(fim, length(seasons_plot)))
  }
  
  ## Set prn+ data 
  for (i in 1: length(clades_prn_pos)){
    if (clades_prn_pos[i] == '0a'){
      fim = 'fim3-1'
      ptxp = 'ptxP1'
    }
    if (clades_prn_pos[i] == '0b'){
      fim = 'fim3-1'
      ptxp = 'ptxP3'
    }
    if (clades_prn_pos[i] == '0c'){
      fim = 'fim3-2'
      ptxp = 'ptxP3'
    }
    data_obs[[length(list_clades_neg) + i]] = data.frame('Clade' = rep(clades_prn_pos[i], length(seasons_plot)),
                                                         'count' = list_x_prn_pos[[i]],
                                                         'Year' = seasons_plot,
                                                         'Prn_exp' = rep(1, length(seasons_plot)), 
                                                         'ptxp_allele' = rep(ptxp, length(seasons_plot)), 
                                                         'Fim_3_allele' = rep(fim, length(seasons_plot)))
  }
  return(data_obs)
}

data_obs_AU = compute_data_obs_country(list_AU[[1]], list_AU[[2]], list_AU[[3]])
data_obs_BE = compute_data_obs_country(list_BE[[1]], list_BE[[2]], list_BE[[3]])
data_obs_CA = compute_data_obs_country(list_CA[[1]], list_CA[[2]], list_CA[[3]])
data_obs_CN = compute_data_obs_country(list_CN[[1]], list_CN[[2]], list_CN[[3]])
data_obs_CZ = compute_data_obs_country(list_CZ[[1]], list_CZ[[2]], list_CZ[[3]])
data_obs_DK = compute_data_obs_country(list_DK[[1]], list_DK[[2]], list_DK[[3]])
data_obs_ES = compute_data_obs_country(list_ES[[1]], list_ES[[2]], list_ES[[3]])
data_obs_FI = compute_data_obs_country(list_FI[[1]], list_FI[[2]], list_FI[[3]])
data_obs_FR = compute_data_obs_country(list_FR[[1]], list_FR[[2]], list_FR[[3]])
data_obs_IE = compute_data_obs_country(list_IE[[1]], list_IE[[2]], list_IE[[3]])
data_obs_IL = compute_data_obs_country(list_IL[[1]], list_IL[[2]], list_IL[[3]])
data_obs_IR = compute_data_obs_country(list_IR[[1]], list_IR[[2]], list_IR[[3]])
data_obs_IT = compute_data_obs_country(list_IT[[1]], list_IT[[2]], list_IT[[3]])
data_obs_JP = compute_data_obs_country(list_JP[[1]], list_JP[[2]], list_JP[[3]])
data_obs_NL = compute_data_obs_country(list_NL[[1]], list_NL[[2]], list_NL[[3]])
data_obs_NO = compute_data_obs_country(list_NO[[1]], list_NO[[2]], list_NO[[3]])
data_obs_SE = compute_data_obs_country(list_SE[[1]], list_SE[[2]], list_SE[[3]])
data_obs_TN = compute_data_obs_country(list_TN[[1]], list_TN[[2]], list_TN[[3]])
data_obs_UK = compute_data_obs_country(list_UK[[1]], list_UK[[2]], list_UK[[3]])
data_obs_US = compute_data_obs_country(list_US[[1]], list_US[[2]], list_US[[3]])

################################################
## Merge clades per genotype, per country
################################################
genotypes_clades = matrix(NA, length(data_obs_FR), 3)
for(i in 1:nrow(genotypes_clades)){
  genotypes_clades[i,1] = as.character(data_obs_FR[[i]][1,c(4)])
  genotypes_clades[i,2] = as.character(data_obs_FR[[i]][1,c(5)])
  genotypes_clades[i,3] = as.character(data_obs_FR[[i]][1,c(6)])
}

merge_clades_per_genotype = function(d){
  data_merged = NULL
  
  ## ancestral: ptxP1 - fim3-1 - prn+
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(1, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "1" &
              genotypes_clades[,2] == "ptxP1" &
              genotypes_clades[,3] == "fim3-1")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn+', 'ptxP1', 'fim3-1')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[1]] = tmp
  
  ## ptxP1 - fim3-1 - prn-
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(2, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "0" &
              genotypes_clades[,2] == "ptxP1" &
              genotypes_clades[,3] == "fim3-1")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn-', 'ptxP1', 'fim3-1')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[2]] = tmp
  
  ## ptxP3 - fim3-1 - prn+
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(3, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "1" &
              genotypes_clades[,2] == "ptxP3" &
              genotypes_clades[,3] == "fim3-1")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn+', 'ptxP3', 'fim3-1')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[3]] = tmp
  
  ## ptxP3 - fim3-1 - prn-
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(4, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "0" &
              genotypes_clades[,2] == "ptxP3" &
              genotypes_clades[,3] == "fim3-1")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn-', 'ptxP3', 'fim3-1')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[4]] = tmp
  
  ## ptxP3 - fim3-2 - prn+
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(5, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "1" &
              genotypes_clades[,2] == "ptxP3" &
              genotypes_clades[,3] == "fim3-2")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn+', 'ptxP3', 'fim3-2')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[5]] = tmp
  
  ## ptxP3 - fim3-2 - prn-
  tmp = matrix(0, nb_years, 6)
  colnames(tmp) = c('Clade', 'count', 'Year',  'Prn_exp', 'ptxp_allele', 'Fim_3_allele')
  tmp[,1] = rep(6, nb_years)
  tmp[,3] = d[[1]][,3]
  a = which(genotypes_clades[,1] == "0" &
              genotypes_clades[,2] == "ptxP3" &
              genotypes_clades[,3] == "fim3-2")
  for(i in 1:nb_years){
    tmp[i,4:6] = c('prn-', 'ptxP3', 'fim3-2')
  }
  for(i in a){
    tmp[,2] = as.numeric(tmp[,2])+as.numeric(d[[i]][,2])
  }
  data_merged[[6]] = tmp
  
  return(data_merged)
}
nb_years = nrow(data_obs_FR[[1]])
data_obs_AU_geno = merge_clades_per_genotype(data_obs_AU)
data_obs_BE_geno = merge_clades_per_genotype(data_obs_BE)
data_obs_CA_geno = merge_clades_per_genotype(data_obs_CA)
data_obs_CN_geno = merge_clades_per_genotype(data_obs_CN)
data_obs_CZ_geno = merge_clades_per_genotype(data_obs_CZ)
data_obs_DK_geno = merge_clades_per_genotype(data_obs_DK)
data_obs_ES_geno = merge_clades_per_genotype(data_obs_ES)
data_obs_FI_geno = merge_clades_per_genotype(data_obs_FI)
data_obs_FR_geno = merge_clades_per_genotype(data_obs_FR)
data_obs_IE_geno = merge_clades_per_genotype(data_obs_IE)
data_obs_IL_geno = merge_clades_per_genotype(data_obs_IL)
data_obs_IR_geno = merge_clades_per_genotype(data_obs_IR)
data_obs_IT_geno = merge_clades_per_genotype(data_obs_IT)
data_obs_JP_geno = merge_clades_per_genotype(data_obs_JP)
data_obs_NL_geno = merge_clades_per_genotype(data_obs_NL)
data_obs_NO_geno = merge_clades_per_genotype(data_obs_NO)
data_obs_SE_geno = merge_clades_per_genotype(data_obs_SE)
data_obs_TN_geno = merge_clades_per_genotype(data_obs_TN)
data_obs_UK_geno = merge_clades_per_genotype(data_obs_UK)
data_obs_US_geno = merge_clades_per_genotype(data_obs_US)

############################################################################
## Compute total numbers function
############################################################################
recompute_total_number_per_year = function(data_tmp){
  res = rep(0, length(seasons_plot))
  for (i in 1:length(data_tmp)){
    res = res + data_tmp[[i]]$count
  }
  return(res)
}
################################################################################







################################################################################
## Prepare data in the right fotmat for fitnel model
################################################################################
nb_clades = length(data_obs_FR_geno)
nb_years = nrow(data_obs_FR[[1]])
nb_alleles = 3
nb_parameters = nb_alleles*2
nb_countries = 20
threshold = 0
ref_clade = 3
subsample = NA
################################################################################

clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_ref = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
clade_number_array_freq_non_zeros = array(0, dim = c(nb_clades, nb_years, nb_countries))
clade_number_array_freq_non_zeros_ref = array(0, dim = c(nb_clades-1, nb_years, nb_countries))
clade_covariates_array = array(0, dim = c(nb_clades, nb_parameters, nb_countries))
clade_number_array_freq_first_non_zeros = rep(0, nb_countries)
total_number = matrix(0, nrow = nb_countries, ncol = length(seasons_plot))

### ACV introduction as main vaccination
vaccine_introduction = rep(NA, nb_countries)
y = min_date + (1:nb_years)*2
vaccine_introduction[1] = 1999 ## Australia
vaccine_introduction[2] = 1999 ## Belgium
vaccine_introduction[3] = 1997 ## Canada
vaccine_introduction[4] = 2012 ## China
vaccine_introduction[5] = 2007 ## Czech
vaccine_introduction[6] = 1997 ## Denmark
vaccine_introduction[7] = 2005 ## Spain
vaccine_introduction[8] = 2005 ## Finland
vaccine_introduction[9] = 2004 ## France
vaccine_introduction[10] = 1996 ## Ireland
vaccine_introduction[11] = 2002 ## Israel
vaccine_introduction[12] = 1995 ## Italy
vaccine_introduction[13] = 2005 ## NL
vaccine_introduction[14] = 1998 ## Norway
vaccine_introduction[15] = 1998 ## Sweden
vaccine_introduction[16] = 2000 ## UK
vaccine_introduction[17] = 1997 ## US
if(nb_countries >= 18) {vaccine_introduction[18] = 1981} ## Japan
if(nb_countries >= 19) {vaccine_introduction[19] = 2050} ## Iran -- Not implemented (yet)
if(nb_countries >= 20) {vaccine_introduction[20] = 2050} ## Tunisia -- Not implemented (yet)
vaccine_introduction = vaccine_introduction-min_date

## Booster ACV introduction in the pop (with or without WCV as primo vaccination)
booster_introduction = rep(NA, nb_countries)
y = min_date + (1:nb_years)*2
booster_introduction[1] = 1997 ## Australia
booster_introduction[2] = 1999 ## Belgium
booster_introduction[3] = 1997 ## Canada
booster_introduction[4] = 2007 ## China
booster_introduction[5] = 2007 ## Czech
booster_introduction[6] = 1997 ## Denmark
booster_introduction[7] = 1999 ## Spain
booster_introduction[8] = 2003 ## Finland
booster_introduction[9] = 1998 ## France
booster_introduction[10] = 1996 ## Ireland
booster_introduction[11] = 2002 ## Israel
booster_introduction[12] = 1995 ## Italy
booster_introduction[13] = 2001 ## NL
booster_introduction[14] = 1998 ## Norway
booster_introduction[15] = 1996 ## Sweden
booster_introduction[16] = 2000 ## UK
booster_introduction[17] = 1992 ## US
if(nb_countries >= 18) {booster_introduction[18] = 1981} ## Japan
if(nb_countries >= 19) {booster_introduction[19] = 2050} ## Iran -- Not implemented (yet)
if(nb_countries >= 20) {booster_introduction[20] = 2050} ## Tunisia -- Not implemented (yet)
booster_introduction = booster_introduction-min_date

## Year WCV implementation
WCV_introduction = rep(NA, nb_countries)
y = min_date + (1:nb_years)*2
WCV_introduction[1] = 1942 ## Australia
WCV_introduction[2] = 1950 ## Belgium
WCV_introduction[3] = 1943 ## Canada
WCV_introduction[4] = 1960 ## China
WCV_introduction[5] = 1958 ## Czech
WCV_introduction[6] = 1961 ## Denmark
WCV_introduction[7] = 1965 ## Spain
WCV_introduction[8] = 1952 ## Finland
WCV_introduction[9] = 1959 ## France
WCV_introduction[10] = 1952 ## Ireland
WCV_introduction[11] = 1957 ## Israel
WCV_introduction[12] = 1961 ## Italy
WCV_introduction[13] = 1953 ## NL
WCV_introduction[14] = 1952 ## Norway
WCV_introduction[15] = 1953 ## Sweden
WCV_introduction[16] = 1957 ## UK
WCV_introduction[17] = 1940 ## US
if(nb_countries >= 18) {WCV_introduction[18] = 1947} ## Japan
if(nb_countries >= 19) {WCV_introduction[19] = 1950} ## Iran 
if(nb_countries >= 20) {WCV_introduction[20] = 1978} ## Tunisia
WCV_introduction = WCV_introduction-min_date

## Year to reach 80% coverage in vaccine
coverage80p = rep(NA, nb_countries)
y = min_date + (1:nb_years)*2
coverage80p[1] = 1988 ## Australia
coverage80p[2] = 1981 ## Belgium
coverage80p[3] = 1988 ## Canada
coverage80p[4] = 1989 ## China
coverage80p[5] = 1994 ## Czech
coverage80p[6] = 1981 ## Denmark
coverage80p[7] = 1985 ## Spain
coverage80p[8] = 1981 ## Finland
coverage80p[9] = 1983 ## France
coverage80p[10] = 1999 ## Ireland
coverage80p[11] = 1981 ## Israel
coverage80p[12] = 1991 ## Italy
coverage80p[13] = 1981 ## NL
coverage80p[14] = 1984 ## Norway
coverage80p[15] = 1981 ## Sweden
coverage80p[16] = 1991 ## UK
coverage80p[17] = 1981 ## US
if(nb_countries >= 18) {coverage80p[18] = 1983} ## Japan
if(nb_countries >= 19) {coverage80p[19] = 1988} ## Iran 
if(nb_countries >= 20) {coverage80p[20] = 1987} ## Tunisia 
coverage80p = coverage80p-min_date

## Year to reach 90% coverage in vaccine
coverage90p = rep(NA, nb_countries)
y = min_date + (1:nb_years)*2
coverage90p[1] = 1990 ## Australia
coverage90p[2] = 1981 ## Belgium
coverage90p[3] = 1994 ## Canada
coverage90p[4] = 1989 ## China
coverage90p[5] = 1994 ## Czech
coverage90p[6] = 1986 ## Denmark
coverage90p[7] = 1985 ## Spain
coverage90p[8] = 1981 ## Finland
coverage90p[9] = 1984 ## France
coverage90p[10] = 2006 ## Ireland
coverage90p[11] = 1982 ## Israel
coverage90p[12] = 1992 ## Italy
coverage90p[13] = 1981 ## NL
coverage90p[14] = 1984 ## Norway
coverage90p[15] = 1981 ## Sweden
coverage90p[16] = 1993 ## UK
coverage90p[17] = 1981 ## US
if(nb_countries >= 18) {coverage90p[18] = 1991} ## Japan
if(nb_countries >= 19) {coverage90p[19] = 1990} ## Iran -- ACV Not implemented (yet)
if(nb_countries >= 20) {coverage90p[20] = 1988} ## Tunisia -- ACV Not implemented (yet)
coverage90p = coverage90p-min_date

for (kkk in 1:nb_countries){
  ## Total numbers
  if (kkk == 1) data_obs = data_obs_AU
  if (kkk == 2) data_obs = data_obs_BE
  if (kkk == 3) data_obs = data_obs_CA
  if (kkk == 4) data_obs = data_obs_CN
  if (kkk == 5) data_obs = data_obs_CZ
  if (kkk == 6) data_obs = data_obs_DK
  if (kkk == 7) data_obs = data_obs_ES
  if (kkk == 8) data_obs = data_obs_FI
  if (kkk == 9) data_obs = data_obs_FR
  if (kkk == 10) data_obs = data_obs_IE
  if (kkk == 11) data_obs = data_obs_IL
  if (kkk == 12) data_obs = data_obs_IT
  if (kkk == 13) data_obs = data_obs_NL
  if (kkk == 14) data_obs = data_obs_NO
  if (kkk == 15) data_obs = data_obs_SE
  if (kkk == 16) data_obs = data_obs_UK
  if (kkk == 17) data_obs = data_obs_US
  if (kkk == 18) data_obs = data_obs_JP
  if (kkk == 19) data_obs = data_obs_IR
  if (kkk == 20) data_obs = data_obs_TN
  total_number[kkk,] = recompute_total_number_per_year(data_obs)
  if(is.na(subsample) == F) {    
    sampling_prop = total_number[kkk,]/subsample
    total_number[kkk,which(total_number[kkk,]>subsample)] = subsample
  }
  
  ## Other params, sich as coutns etc
  if (kkk == 1) data_obs = data_obs_AU_geno
  if (kkk == 2) data_obs = data_obs_BE_geno
  if (kkk == 3) data_obs = data_obs_CA_geno
  if (kkk == 4) data_obs = data_obs_CN_geno
  if (kkk == 5) data_obs = data_obs_CZ_geno
  if (kkk == 6) data_obs = data_obs_DK_geno
  if (kkk == 7) data_obs = data_obs_ES_geno
  if (kkk == 8) data_obs = data_obs_FI_geno
  if (kkk == 9) data_obs = data_obs_FR_geno
  if (kkk == 10) data_obs = data_obs_IE_geno
  if (kkk == 11) data_obs = data_obs_IL_geno
  if (kkk == 12) data_obs = data_obs_IT_geno
  if (kkk == 13) data_obs = data_obs_NL_geno
  if (kkk == 14) data_obs = data_obs_NO_geno
  if (kkk == 15) data_obs = data_obs_SE_geno
  if (kkk == 16) data_obs = data_obs_UK_geno
  if (kkk == 17) data_obs = data_obs_US_geno
  if (kkk == 18) data_obs = data_obs_JP_geno
  if (kkk == 19) data_obs = data_obs_IR_geno
  if (kkk == 20) data_obs = data_obs_TN_geno
  
  clade_number = lapply(data_obs, function(x)x[,2])
  clade_number_mat = clade_number[[1]]
  clade_covariates = lapply(data_obs, function(x)x[1,4])
  clade_covariates_clean1 = clade_covariates
  clade_covariates = lapply(data_obs, function(x)x[1,5])
  clade_covariates_clean2 = clade_covariates
  clade_covariates = lapply(data_obs, function(x)x[1,6])
  clade_covariates_clean3 = clade_covariates
  
  tmp = c(0,0,0,0,0,0)
  if(clade_covariates_clean1[[1]] == 'prn-')  tmp [1] = 1 ## prn deficient strains
  if(clade_covariates_clean1[[1]] == 'prn+')  tmp [2] = 1 ## prn WT strains
  if(clade_covariates_clean2[[1]] == 'ptxP3')  tmp [3] = 1 ## ptxP3 strains
  if(clade_covariates_clean2[[1]] == 'ptxP1')  tmp [4] = 1 ## ptxP1 strains
  if(clade_covariates_clean3[[1]] == 'fim3-2')  tmp [5] = 1 ## fim3-2 strains
  if(clade_covariates_clean3[[1]] == 'fim3-1')  tmp [6] = 1 ## fim3-1 strains
  
  clade_covariates_mat = tmp
  for (i in 2:nb_clades){
    tmp = c(0,0,0,0,0,0)
    if(clade_covariates_clean1[[i]] == 'prn-')  tmp [1] = 1 ## prn deficient strains
    if(clade_covariates_clean1[[i]] == 'prn+')  tmp [2] = 1 ## prn WT strains
    if(clade_covariates_clean2[[i]] == "ptxP3")  tmp [3] = 1 ## ptxP3 strains
    if(clade_covariates_clean2[[i]] == 'ptxP1')  tmp [4] = 1 ## ptxP1 strains
    if(clade_covariates_clean3[[i]] == 'fim3-2')  tmp [5] = 1 ## fim3-2 strains
    if(clade_covariates_clean3[[i]] == 'fim3-1')  tmp [6] = 1 ## fim3-1 strains
    clade_covariates_mat = rbind(clade_covariates_mat, tmp)
    clade_number_mat = rbind(clade_number_mat, clade_number[[i]])
  }
  
  clade_number_array[,,kkk] = apply(clade_number_mat, MARGIN = 2, as.numeric)
  clade_number_mat_freq = apply(clade_number_mat, MARGIN = 2, as.numeric)
  
  if(is.na(subsample) == F) {
    a = which(sampling_prop>1)
    if(length(a)>0){
      for(i in 1:length(a)){
        clade_number_mat_freq[,a[i]] = round(clade_number_mat_freq[,a[i]]/sampling_prop[a[i]])
      }
      clade_number_array[,,kkk] = clade_number_mat_freq
      total_number[kkk,] =  apply(clade_number_mat_freq, MARGIN = 2, sum) ## some values will not be 'subsample' number, but very very close (like 101 for the US)
    }
  }
  
  for(i in 1:nrow(clade_number_mat_freq)){
    clade_number_mat_freq[i,] = as.numeric(clade_number_mat_freq[i,])/t(total_number)[,kkk]
  }
  clade_number_mat_freq[which(is.na(clade_number_mat_freq))] = 0
  clade_number_array_freq[,,kkk] = clade_number_mat_freq
  clade_number_array_freq_non_zeros[,,kkk] = clade_number_mat_freq
  for(i in 1:nb_years){
    if(sum(clade_number_array[,i,kkk]) > threshold){ ## previous version
      clade_number_array_freq_non_zeros[,i,kkk] = 1
    }
    if(sum(clade_number_array[,i,kkk]) <= threshold){
      clade_number_array_freq_non_zeros[,i,kkk] = 0
    }
    count = 1
    for(j in 1:(nb_clades)){
      if(j != ref_clade){
        if((clade_number_array[count,i,kkk])+ (clade_number_array[ref_clade,i,kkk])> threshold){ ## previous version
          clade_number_array_freq_non_zeros_ref[count,i,kkk] = 1
        }
        if((clade_number_array[count,i,kkk])+ (clade_number_array[ref_clade,i,kkk])<= threshold){
          clade_number_array_freq_non_zeros_ref[count,i,kkk] = 0
        }
        count = count+1
      }
    }
  }
  zero = T
  c = 1
  while(zero == T){
    if(sum(clade_number_array_freq_non_zeros[,c,kkk])==0) {c = c+1}
    if(sum(clade_number_array_freq_non_zeros[,c,kkk])>0) {zero = F}
  }
  clade_number_array_freq_first_non_zeros[kkk] = c
  print(which(rowSums(apply(clade_number_mat, MARGIN = 2, as.numeric))>0))
  clade_covariates_array[,,kkk] = clade_covariates_mat
}

## compute freq with respect to a ref clade
for(kkk in 1:nb_countries){
  tmp = 1
  for(i in 1:nb_clades){
    if(i != ref_clade){
      clade_number_array_freq_ref[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[3,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref[which(is.na(clade_number_array_freq_ref) == T)] = 0
clade_number_array_ref = clade_number_array[ref_clade,,]
clade_number_array_paired = clade_number_array[-ref_clade,,]


## Gather data in one list
data.MCMC = list(nb_genotypes = nb_clades,
                 nb_years = nb_years,
                 nb_countries = nb_countries,
                 
                 data_genotype_non_ref = clade_number_array_paired,
                 data_genotype_ref = clade_number_array_ref,
                 data_total_number = t(total_number),
                 non_zero_country_year = clade_number_array_freq_non_zeros[1,,],
                 non_zero_country_year_genotype = clade_number_array_freq_non_zeros_ref,
                 number_zeros_country_year = length(which(clade_number_array_freq_non_zeros[1:(nb_clades-1),,] == 0)),
                 number_zeros_country_year_genotype = length(which(clade_number_array_freq_non_zeros_ref == 0)),
                 
                 vaccine_introduction = vaccine_introduction,
                 booster_introduction = booster_introduction,
                 WCV_introduction = WCV_introduction,
                 coverage90p = coverage90p,
                 coverage80p = coverage80p)

################################################################################
## Save data
################################################################################
saveRDS(data.MCMC, 'Estimate_relative_fitness/Data_model_allcountries_refgeno3.rds')
################################################################################

