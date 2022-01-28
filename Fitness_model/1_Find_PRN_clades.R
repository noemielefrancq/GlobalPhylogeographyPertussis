#######################################################
## Code to find the PRN-deficient clades 
## in a phylogenetic tree
#######################################################
### Author: Noemie Lefrancq
### Last modification: 28/01/2022
#######################################################


## Load necessary packages
library(ape)
library(treeio)

#######################################################
### NEED A TREE ORDERED BY 'INCRESING NODE ORDER' (!) 
### FIGTREE IS THE SIMPLEST
##########################################################

## Load tree
tree.beast = read.beast('Data/Beast_tree_MCC_26012022.tree')


## Load data
data_tip = read.csv(file = 'Data/Metadata_analyses_26012022.csv', sep = ';')
data_tip$ï..Isolate[which(data_tip$ï..Isolate == 'Tohama')] = 'Tohama-I'


## Clean prn labels
data_tip$prn.expression = as.character(data_tip$prn.expression)
data_tip$prn.expression[which(data_tip$prn.expression == 'prn-')] = 0
data_tip$prn.expression[which(data_tip$prn.expression == 'prn+')] = 1
to_drop = which(data_tip$prn.expression == 'Unkown' | is.na(data_tip$prn.expression) == T)
data_tip$prn.expression[which(data_tip$prn.expression == 'Unkown' | is.na(data_tip$prn.expression) == T)] = NA
data_tip$prn.expression = as.numeric(data_tip$prn.expression)


## Drop unknwon PRNs
prn_to_drop = data_tip$ï..Isolate[to_drop] ## tips with unkwown PRN labels
tree.beast = treeio::drop.tip(tree.beast, tip = as.character(prn_to_drop[which(is.na(prn_to_drop) == F)])) ## drop tips with unkwown PRN labels 


## Tree in phylo class
tree = tree.beast@phylo


## Get tip lables
a = match(tree$tip.label, data_tip$ï..Isolate)
data_tip = data_tip[a,]
prn_labels = data_tip$prn.expression


## Write new beast tree
write.beast(tree.beast, 'Estimate_relative_fitness/Beast_tree_MCC_26012022_known_prn_expression.tree')
######################################################################################################################################




######################################################################################################################################
## Find PRN- clades
######################################################################################################################################
clade_number = rep(1, length(prn_labels))
k = 2
ii = 1
while (ii < length(prn_labels)){
  tips_to_check = tree$tip.label[ii:(ii+1)]
  prn_to_check = prn_labels[ii:(ii+1)]
  j=1
  while (sum(prn_to_check) == 0 & is.monophyletic(tree,tips_to_check) == T){
    tips_to_check = c(tips_to_check, tree$tip.label[ii+1+j])
    prn_to_check = c(prn_to_check, prn_labels[ii+1+j])
    j = j+1
  }
  if (j == 1) {
    if (prn_to_check[1] == 0) {
      clade_number[ii] = k
      k = k+1
    }
  }
  if (j > 1) {
    if (sum(prn_to_check[-length(prn_to_check)]) == 0){
    j = j-1
    clade_number[ii:(ii+j)] = k
    k = k+1
    }
  }
  ii = ii+j
  print(ii)
}
d = data.frame('prn' = prn_labels, 'clade' = clade_number)


######################################################################################################################################
## Collapse as much as possible the clades into monophyletic, prn- ones
######################################################################################################################################
clade_number_save = clade_number
clade_number = clade_number_save
n_true = 1
while(n_true != 0){
  clade_number_tmp = rep(1, length(prn_labels))
  k = 2
  ii = 2
  n_true = 0
  while (ii <= length(levels(as.factor(clade_number)))){
    a = which(clade_number == ii)
    b = which(clade_number == ii+1)
    xx = c(a,b)
    tips_to_check = tree$tip.label[c(a,b)]
    prn_to_check = prn_labels[c(a,b)]
    j=1
    while ((sum(prn_to_check) == 0) & is.monophyletic(tree,tips_to_check) == T & length(b) >0){
      n_true = n_true + 1
      b = which(clade_number == ii+1+j)
      tips_to_check = c(tips_to_check, tree$tip.label[b])
      prn_to_check = c(prn_to_check, prn_labels[b])
      xx = c(xx,b)
      j=j+1
    }
    if (j == 1) {
      if (sum(prn_to_check[1:length(a)]) == 0){
        clade_number_tmp[a] = k
        k = k+1
      } 
    }
    if (j > 1) {
      if (sum(prn_to_check[1:(length(xx)- length(b))]) == 0){
        clade_number_tmp[head(xx, length(xx)- length(b))] = k
        k = k+1
      }
    }
    ii = ii+j
    print(paste0(ii, ' / ', length(levels(as.factor(clade_number)))))
  }
  print(paste0('NTRUE = ', n_true))
  clade_number = clade_number_tmp
}
d = data.frame('prn' = prn_labels, 'clade' = clade_number)
tips_names = tree$tip.label
d$tips_names = tree$tip.label


######################################################################################################################################
## Produce a file with information about internal nodes (genotypes and co)
######################################################################################################################################
library(treeio)
library(tibble)
x <- as_tibble(tree.beast)
d <- tibble(label = tree$tip.label)
y <- full_join(x, d, by = 'label')

######################
## Distinguish strains based on genotype 
######################
## fim31 & ptxp1 (it is the orginal geneotype, so label all the treedata as such, then make some changes)
######################
d = y
d$clade = rep("0a", nrow(d)) ## original  clade
d$prn = rep(1, nrow(d))

######################
## fim3-1 & ptxp3
######################
a = which(clade_number == 1)
tips = tips_names[a]
data_tip_tmp = data_tip[which(data_tip$fim3 == 'fim3-1' & 
                                                            data_tip$ptxP == 'ptxP3'),]
b = match(tips, data_tip_tmp$ï..Isolate)
MRCA = getMRCA(tree.beast@phylo, tips[which(!is.na(b))])
e = as.data.frame(offspring(y, MRCA))
m = match(e$node, d$node)
d$clade[m] = rep('0b', length(m))

######################
## fim3-2 & ptxp3
######################
a = which(clade_number == 1)
tips = tips_names[a]
data_tip_tmp = data_tip[which(data_tip$fim3 == 'fim3-2' & 
                                                            data_tip$ptxP == 'ptxP3'),]
b = match(tips, data_tip_tmp$ï..Isolate)
MRCA = getMRCA(tree.beast@phylo, tips[which(!is.na(b))])
e = as.data.frame(offspring(y, MRCA))
m = match(e$node, d$node)
d$clade[m] = rep('0c', length(m))

######################
## Other clades
######################
for (ii in 2:max(clade_number)){ ## other clades
  a = which(clade_number == ii)
  tips = tips_names[a]
  if (length(tips) > 1){
    MRCA = getMRCA(tree.beast@phylo, tips)
    e = as.data.frame(offspring(y, MRCA))
    m = match(e$node, d$node)
    d$clade[m] = rep(ii, length(m))
    d$prn[m] = rep(0, length(m))
  }
  if (length(tips) == 1){
    d$clade[which(d$label == tips)] = ii
    d$prn[which(d$label == tips)] = 0
  }
}

## add date
d$date = 2019-as.numeric(d$height)

######################################################################################################################################
## Write output: the list of clades
######################################################################################################################################
saveRDS(d, 'Estimate_relative_fitness/data_tips_and_nodes_clades.rds')
write.csv(data_tip, 'Estimate_relative_fitness/data_tip_with_clades.csv')

######################################################################################################################################
