#######################################################
##       Script used to filter vcf files 
#######################################################
### Author: Noemie Lefrancq
### Last modification: 27/01/2020
#######################################################

## Load libraries
library(seqinr)
library(stringr)
library(vcfR)

## Parameters to select the SNPs:
q_second = 100 #Phred quality score 
n_reads_F_R = 2 # minimum number of reads per strand
DP_min = 5  # Minimum read depth
p_reads_min1 = 0.40 #Percentage of reads to replace by IUPAC code (REF or potential SNP) 
p_reads_min2 = 0.80 #Percentage of reads to have to be a SNP
p_ref = 1 

## Function used to replace uncertain SNP by IUPAC code
replace <- function(y){
  y = str_replace(y, "A,C,G,T", "N")
  y = str_replace(y, "A,C,T,G", "N") 
  y = str_replace(y, "A,T,G,C", "N")
  y = str_replace(y, "A,T,C,G", "N")
  y = str_replace(y, "A,G,C,T", "N")
  y = str_replace(y, "A,G,T,C", "N")
  y = str_replace(y, "C,A,G,T", "N")
  y = str_replace(y, "C,A,T,G", "N")
  y = str_replace(y, "C,T,G,A", "N")
  y = str_replace(y, "C,T,A,G", "N")
  y = str_replace(y, "C,G,T,A", "N")
  y = str_replace(y, "C,G,A,T", "N")
  y = str_replace(y, "T,A,C,G", "N")
  y = str_replace(y, "T,A,G,C", "N")
  y = str_replace(y, "T,C,A,T", "N")
  y = str_replace(y, "T,C,T,A", "N")
  y = str_replace(y, "T,G,C,A", "N")
  y = str_replace(y, "T,G,A,C", "N")
  y = str_replace(y, "G,A,C,T", "N")
  y = str_replace(y, "G,A,T,C", "N")
  y = str_replace(y, "G,C,T,A", "N")
  y = str_replace(y, "G,C,A,T", "N")
  y = str_replace(y, "G,T,A,C", "N")
  y = str_replace(y, "G,T,C,A", "N")
  
  y = str_replace(y, "C,G,T", "B")
  y = str_replace(y, "C,T,G", "B")
  y = str_replace(y, "G,C,T", "B")
  y = str_replace(y, "G,T,C", "B")
  y = str_replace(y, "T,G,C", "B")
  y = str_replace(y, "T,C,G", "B")
  
  y = str_replace(y, "A,G,T", "D")
  y = str_replace(y, "A,T,G", "D")
  y = str_replace(y, "T,A,G", "D")
  y = str_replace(y, "T,G,A", "D")
  y = str_replace(y, "G,T,A", "D")
  y = str_replace(y, "G,A,T", "D")
  
  y = str_replace(y, "A,C,T", "H")
  y = str_replace(y, "A,T,C", "H")
  y = str_replace(y, "T,A,C", "H")
  y = str_replace(y, "T,C,A", "H")
  y = str_replace(y, "C,A,T", "H")
  y = str_replace(y, "C,T,A", "H")
  
  y = str_replace(y, "A,C,G", "V")
  y = str_replace(y, "A,G,C", "V")
  y = str_replace(y, "G,C,A", "V")
  y = str_replace(y, "G,A,C", "V")
  y = str_replace(y, "C,A,G", "V")
  y = str_replace(y, "C,G,A", "V")
  
  y = str_replace(y, "A,G", "R")
  y = str_replace(y, "G,A", "R")
  y = str_replace(y, "C,T", "Y")
  y = str_replace(y, "T,C", "Y")
  y = str_replace(y, "G,C", "S")
  y = str_replace(y, "C,G", "S")
  y = str_replace(y, "A,T", "W")
  y = str_replace(y, "T,A", "W")
  y = str_replace(y, "G,T", "K")
  y = str_replace(y, "T,G", "K")
  y = str_replace(y, "A,C", "M")
  y = str_replace(y, "C,A", "M")
}
## Function used to take out phage regions
positions_phage_region_N = function(x){
  pos_rm = which(x$POS>=37926 & x$POS<=46978) #phage region 1
  pos_rm = c(pos_rm, which(x$POS>=494583 & x$POS<=531881)) #phage region 2
  pos_rm = c(pos_rm, which(x$POS>=3573119 & x$POS<=3592359)) #phage region 3
  pos_rm = c(pos_rm, which(x$POS>=3592730 & x$POS<=3602410)) #phage region 4
  x = x[pos_rm,]
  return(x)
}

## Main function used to filter SNPs
write_filtered_fasta = function(data2, names_year, names_files, names_isolates, place, directory_write, directory_read){
  
  #Compute the lists of SNPs which respect the criterias (see before)
  #Add to the dataframe the DP4 values: number of reads for the REF and for the ALT
  DP4 = lapply(data2, function(x)data.frame("DP"=x$DP, "DP4"=x$`x$gt$gt_SB`))
  for (i in (1:length(DP4))){
    DP4_tmp = matrix(0, nrow = nrow(DP4[[i]]), ncol = 4)
    for (j in (1:length(DP4[[i]]$DP4))){
      DP4_tmp[j,] = as.numeric(unlist(strsplit(as.character(DP4[[i]]$DP4[j]), ',')))
    }
    data2[[i]]$DP4_REF = DP4_tmp[,1]+DP4_tmp[,2] #number of reads for the REF
    data2[[i]]$DP4_ALT = DP4_tmp[,3]+DP4_tmp[,4] #number of reads for the ALT
    data2[[i]]$DP4_ALT_F = DP4_tmp[,3] #number of reads for the REF
    data2[[i]]$DP4_ALT_R = DP4_tmp[,4] #number of reads for the ALT
    data2[[i]]$DP4_R =  data2[[i]]$DP4_ALT_F + DP4_tmp[,1]#number of R reads
    data2[[i]]$DP4_F = data2[[i]]$DP4_ALT_R+ DP4_tmp[,2] #number of F reads
    data2[[i]]$sum = data2[[i]]$DP4_REF + data2[[i]]$DP4_ALT #total number of reads
  }
  
  #distinguish NON REF snps (N) from potential SNPs 
  dataX = data2
  for (i in (1: length(data2))){
    dataX[[i]]= data2[[i]][which(data2[[i]]$ALT == '<NON_REF>'),]
    data2[[i]]= data2[[i]][which(data2[[i]]$ALT != '<NON_REF>' ),]
  }
  
  for (i in (1: length(data2))){
    a = unlist(lapply(data2[[i]]$ALT, function(x)str_split(x, ',')[[1]][1]))
    for (j in which(a == '*')){
      data2[[i]]$ALT[j] = str_sub(data2[[i]]$ALT[j], start=3)
    }
    data2[[i]]$ALT = unlist(lapply(data2[[i]]$ALT, function(x)str_split(x, ',<NON_REF>')[[1]][1]))
    a = unlist(lapply(data2[[i]]$ALT, function(x)str_split(x, ',')[[1]][length(str_split(x, ',')[[1]])]))
    for (j in which(a == '*')){
      data2[[i]]$ALT[j] = str_sub(data2[[i]]$ALT[j], end=-3)
    }
    a = unlist(lapply(data2[[i]]$ALT, function(x)str_split(x, ',')[[1]][2]))
    for (j in which(a == '*')){
      data2[[i]]$ALT[j] = paste0(str_split(data2[[i]]$ALT[j], ',')[[1]][1], ',', str_split(data2[[i]]$ALT[j], ',')[[1]][3]) 
    }
    a = unlist(lapply(data2[[i]]$ALT, function(x)str_split(x, ',')[[1]][3]))
    for (j in which(a == '*')){
      data2[[i]]$ALT[j] = paste0(str_split(data2[[i]]$ALT[j], ',')[[1]][1], ',', str_split(data2[[i]]$ALT[j], ',')[[1]][2], ',', str_split(data2[[i]]$ALT[j], ',')[[1]][4]) 
    }
  }

  #SELECTION OF THE SNPS: AT LEAST 80% (p_reads) OF THE READS AND AT LEAST 10(n_reads) READS
  select_N <- function(x){
    return(x[which(((x$DP4_ALT)/x$sum < p_reads_min1) | (x$sum<DP_min) | (x$DP4_R<n_reads_F_R) | (x$DP4_F<n_reads_F_R) | (x$QUAL<q_second)),])
  }
  data3 = lapply(data2, select_N) # to be replace by th reference
  
  #SELECTION OF THE SNPS: AT LEAST 80% (p_reads) OF THE READS AND AT LEAST 10(n_reads) READS
  select_IUPAC <- function(x){
    return(x[which((x$DP4_ALT)/x$sum >= p_reads_min1 & (x$DP4_ALT)/x$sum < p_reads_min2 & x$DP>=DP_min & x$DP4_R>=n_reads_F_R & x$DP4_F>=n_reads_F_R & x$QUAL>=q_second),])
  }
  data3_IUPAC = lapply(data2, select_IUPAC) # to be replace by th reference
  
  for (i in (1: length(data3_IUPAC))){
    if (nrow(data3_IUPAC[[i]])> 0){
      data3_IUPAC[[i]]$ALT = paste0(data3_IUPAC[[i]]$REF, ',', data3_IUPAC[[i]]$ALT)
      data3_IUPAC[[i]]$ALT = replace(data3_IUPAC[[i]]$ALT ) 
    }
  }
  
  ## Store de positions to change
  ## Phage regions 
  data4 = lapply(data2, positions_phage_region_N)
  
  dataMULTI = data2
  for (i in (1: length(data2))){
    a = lapply(data2[[i]]$ALT, function(x)str_split(x, ',')[[1]])
    b = lapply(a, function(x)length(x))
    b = which(unlist(b) > 1)
    dataMULTI[[i]] = data2[[i]][b,]
    data2[[i]]$ALT = replace(data2[[i]]$ALT )
  }
  
  #read tohama
  tohama = getSequence(read.fasta("TohamaNC_0029292.fasta"))
  tohama = toupper(tohama[[1]])
  
  ############# Handle fastas
  setwd(directory_read)
  if (length(place) != 0){
    files = files[-place]
    names_files = names_files[-place]
    data3 = data3[-place]
    data3_IUPAC = data3_IUPAC[-place]
    data4 = data4[-place]
    dataX = dataX[-place]
    data2 = data2[-place]
  }
  files_for_coverage = lapply(files, function(x)paste0(str_sub(x, end=-10)))
  detach(package:vcfR)
  library(seqinr)
  for (i in (1: length(names_files))){
    f = tohama
    c = read.csv(paste0(directory_read,'coverage/',files_for_coverage[[i]], "_trim_coverage_each_position.txt"))
    c = rbind(str_split(colnames(c), 'X')[[1]][[2]], c)
    f[data2[[i]]$POS] = data2[[i]]$ALT
    f[data3_IUPAC[[i]]$POS] = data3_IUPAC[[i]]$ALT 
    f[dataX[[i]]$POS] = 'N' ## replace unsure SNPs by reference
    if (nrow(data3[[i]]) >0 ){f[data3[[i]]$POS] = 'N' }## replace unsure SNPs by reference
    f[data4[[i]]$POS] = 'N' ## replace uncovered regions by 'N'
    f[which(as.numeric(c[,1])<=DP_min)] = '-' ## replace uncovered regions by 'N''
    write.fasta(f, names = names_year[[i]],  file.out = paste0(directory_write,names_year[[i]], extension), open='w')
    print(paste0("writing sequence ", i, "/", length(names_files), '  length sequence = ', length(f), 'bp'))
  }
}
