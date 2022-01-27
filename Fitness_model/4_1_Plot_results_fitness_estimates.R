library(ggplot2)
library(rstan)

setwd('E:/Stages/Project_Pertussis/fitness_genotypes/')

## Load fit
fit = readRDS(file = 'Estimate_relative_fitness/Output_WCV_booster_ref3_fit_all.rds')

## Chains
Chains=rstan::extract(fit$fit)

################################################################################
## Functions
################################################################################
mean.and.ci <-function(v){ 
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}
# define the summary function for plot in ggplot2
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
f2 <- function(x) {
  r <- quantile(x, probs = c(0.025,0.5,0.975))
  names(r) <- c("ymin","y","ymax")
  r
}
################################################################################

################################################################################
## Compute data to plot
################################################################################

################################################################################
## Growth rate (fitness) weighted by proportions
################################################################################
vaccination = fit$data$booster_introduction
nb_chains = length(Chains$lp__)
nb_countries = fit$data$nb_countries
nb_years = fit$data$nb_years
length_vect_fitness_post_vacc = fit$data$number_R_post_vacc
length_vect_fitness_pre_vacc = fit$data$number_R_pre_vacc
nb_genotypes = fit$data$nb_genotypes
ref_genotype = 3
ref_genotype_header = 'Genotype 3'
Genotype = c('Genotype 1', 'Genotype 2', "Genotype 4", "Genotype 5", "Genotype 6")

df_overall_fitness = data.frame('Values' = NA,
                                'Time' = NA,
                                'Genotype' = NA)

vaccinationACV = fit$data$yearIntroduction
vaccinationWCV = fit$data$yearF0

for(k in 1:nb_countries){
  ## Compute reference fitness first:
  ## pre vacc
  mat_tmp = matrix(NA, ncol = 6, nrow = nb_chains*length_vect_fitness_pre_vacc)
  vacc_ACV_k = vaccinationACV[k]
  vacc_WCV_k = vaccinationWCV[k]
  if(vacc_ACV_k>=nb_years) vacc_ACV_k = nb_years-1;
  if(vacc_WCV_k<=0) vacc_WCV_k = 1;
  if(vacc_WCV_k == (vacc_ACV_k-1)) vacc_WCV_k = vacc_ACV_k-2;
  if(vacc_WCV_k >= (vacc_ACV_k)) vacc_WCV_k = vacc_ACV_k-2;
  for(gref in 1:(nb_genotypes-1)){
    for(i in 1:length_vect_fitness_pre_vacc){
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_pre_vacc[,gref,i])
      mean_freq = Chains$pred_absolute_freq[,k,gref, vacc_WCV_k:(vacc_ACV_k-1)]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
    }
  }
  for(i in 1:length_vect_fitness_pre_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = rep(0, nb_chains)-rep(0, nb_chains)
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_WCV_k:(vacc_ACV_k-1)]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(-1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  ## post vacc
  mat_tmp = matrix(NA, ncol = 6, nrow = nb_chains*length_vect_fitness_post_vacc)
  for(gref in 1:(nb_genotypes-1)){
    for(i in 1:length_vect_fitness_post_vacc){
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] =  rep(0, nb_chains)-as.vector(Chains$fitness_genotypes_post_vacc[,gref,i])
      mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_ACV_k:nb_years]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
    }
  }
  for(i in 1:length_vect_fitness_post_vacc){
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] =  rep(0, nb_chains)-rep(0, nb_chains)
    mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_ACV_k:nb_years]
    mean_freq = apply(mean_freq, MARGIN = 1, mean)
    mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
  }
  df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                            'Time' = rep(1, each = nb_chains),
                                            'Genotype' = rep(ref_genotype_header, nb_chains*length_vect_fitness_post_vacc)))
  for(g in 1:(nb_genotypes-1)){
    ## pre vacc
    mat_tmp = matrix(NA, ncol = 6, nrow = nb_chains*length_vect_fitness_pre_vacc)
    for(gref in 1:(nb_genotypes-1)){
      for(i in 1:length_vect_fitness_pre_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_pre_vacc[,g,i])-as.vector(Chains$fitness_genotypes_pre_vacc[,gref,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_WCV_k:(vacc_ACV_k-1)]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    for(i in 1:length_vect_fitness_pre_vacc){
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_pre_vacc[,g,i])-rep(0, nb_chains)
      mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_WCV_k:(vacc_ACV_k-1)]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
    }
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                              'Time' = rep(-1, each = nb_chains),
                                              'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
    
    ## post vacc
    mat_tmp = matrix(NA, ncol = 6, nrow = nb_chains*length_vect_fitness_post_vacc)
    for(gref in 1:(nb_genotypes-1)){
      for(i in 1:length_vect_fitness_post_vacc){
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = as.vector(Chains$fitness_genotypes_post_vacc[,g,i])-as.vector(Chains$fitness_genotypes_post_vacc[,gref,i])
        mean_freq = Chains$pred_absolute_freq[,k,gref,vacc_ACV_k:nb_years]
        mean_freq = apply(mean_freq, MARGIN = 1, mean)
        mat_tmp[(1:nb_chains)+((i-1)*nb_chains),gref] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),gref] * mean_freq
      }
    }
    for(i in 1:length_vect_fitness_post_vacc){
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = as.vector(Chains$fitness_genotypes_post_vacc[,g,i])-rep(0, nb_chains)
      mean_freq = Chains$pred_absolute_freq[,k,nb_genotypes,vacc_ACV_k:nb_years]
      mean_freq = apply(mean_freq, MARGIN = 1, mean)
      mat_tmp[(1:nb_chains)+((i-1)*nb_chains),nb_genotypes] = mat_tmp[(1:nb_chains)+((1-1)*nb_chains),nb_genotypes] * mean_freq
    }
    df_overall_fitness = rbind(df_overall_fitness, data.frame('Values' = exp(apply(mat_tmp, MARGIN = 1, sum)),
                                              'Time' = rep(1, each = nb_chains),
                                              'Genotype' = rep(Genotype[g], nb_chains*length_vect_fitness_post_vacc)))
  }
}
df_overall_fitness = df_overall_fitness[-1,] ## Remove the NA from the beginning


################################################################################
## Compute the effect of prn shift
################################################################################
df_effect_prn_shift = data.frame('Values' = NA,
                                       'Time' = NA,
                                       'Genotype' = NA)

## ACV
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 2')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 1')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 1', nb_chains)))

a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 4')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 3')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 3', nb_chains)))

a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 6')
b = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 5')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 5', nb_chains)))
df_effect_prn_shift = df_effect_prn_shift[-1,] ## Remove the NA from the beginning

## Average
df_effect_prn_shift_averages = data.frame('Values' = NA,
                                           'Time' = NA)
a = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 2'))
b = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 4'))
c = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 6'))
prnneg = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

a = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 1'))
b = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 3'))
c = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 5'))
prnpos = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

df_effect_prn_shift_averages = rbind(df_effect_prn_shift_averages, data.frame('Values' = prnneg/prnpos,
                                                            'Time' = rep(1, each = length(prnneg))))

## WCV
a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 2')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 1')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(-1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 1', nb_chains)))

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 4')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 3')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(-1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 3', nb_chains)))

a = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 6')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 5')
df_effect_prn_shift = rbind(df_effect_prn_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(-1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 5', nb_chains)))

## Compute average across background genotypes
a = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 2'))
b = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 4'))
c = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 6'))
prnneg = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

a = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 1'))
b = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 3'))
c = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 5'))
prnpos = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

df_effect_prn_shift_averages = rbind(df_effect_prn_shift_averages, data.frame('Values' = prnneg/prnpos,
                                                                              'Time' = rep(-1, each = length(prnneg))))
################################################################################

################################################################################
## Compute effect of vaccine shift, for each genotype
################################################################################
df_effect_vaccine_shift = data.frame('Values' = NA,
                                       'Time' = NA,
                                       'Genotype' = NA)
## G1
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 1')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 1')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                        'Time' = rep(1, each = nb_chains),
                                                                        'Genotype' = rep('Genotype 1', nb_chains)))

## G2
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 2')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 2')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                    'Time' = rep(1, each = nb_chains),
                                                                    'Genotype' = rep('Genotype 2', nb_chains)))

## G3
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 3')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 3')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                    'Time' = rep(1, each = nb_chains),
                                                                    'Genotype' = rep('Genotype 3', nb_chains)))

## G4
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 4')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 4')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                    'Time' = rep(1, each = nb_chains),
                                                                    'Genotype' = rep('Genotype 4', nb_chains)))

## G5
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 5')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 5')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                    'Time' = rep(1, each = nb_chains),
                                                                    'Genotype' = rep('Genotype 5', nb_chains)))

## G6
a = which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 6')
b = which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 6')
df_effect_vaccine_shift = rbind(df_effect_vaccine_shift, data.frame('Values' = df_overall_fitness$Values[a]/df_overall_fitness$Values[b],
                                                                    'Time' = rep(1, each = nb_chains),
                                                                    'Genotype' = rep('Genotype 6', nb_chains)))

df_effect_vaccine_shift = df_effect_vaccine_shift[-1,] ## Remove the NA from the beginning

## Fitness of PRN- strains, ACV vs WCV
a = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 2'))
b = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 4'))
c = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 6'))
fitnessPRNneg = (df_effect_vaccine_shift$Values[a]+df_effect_vaccine_shift$Values[b]+df_effect_vaccine_shift$Values[c])/3
mean.and.ci(fitnessPRNneg)

## Fitness of PRN+ strains, ACV vs WCV
a = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 1'))
b = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 3'))
c = c(which(df_effect_vaccine_shift$Genotype == 'Genotype 5'))
fitnessPRNpos = (df_effect_vaccine_shift$Values[a]+df_effect_vaccine_shift$Values[b]+df_effect_vaccine_shift$Values[c])/3
mean.and.ci(fitnessPRNpos)

## Average
df_effect_vaccine_shift_averages = data.frame('Values' = NA,
                                              'Prn' = NA)
a = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 2'))
b = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 4'))
c = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 6'))
prnnegACV = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

a = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 2'))
b = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 4'))
c = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 6'))
prnnegWCV = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

df_effect_vaccine_shift_averages = rbind(df_effect_vaccine_shift_averages, data.frame('Values' = prnnegACV/prnnegWCV,
                                                                              'Prn' = rep("prnneg", each = length(prnnegACV))))

a = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 1'))
b = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 3'))
c = c(which(df_overall_fitness$Time == '1' & df_overall_fitness$Genotype == 'Genotype 5'))
prnposACV = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

a = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 1'))
b = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 3'))
c = c(which(df_overall_fitness$Time == '-1' & df_overall_fitness$Genotype == 'Genotype 5'))
prnposWCV = (df_overall_fitness$Values[a]+df_overall_fitness$Values[b]+df_overall_fitness$Values[c])/3

df_effect_vaccine_shift_averages = rbind(df_effect_vaccine_shift_averages, data.frame('Values' = prnposACV/prnposWCV,
                                                                                      'Prn' = rep("prnpos", each = length(prnposACV))))
# mean.and.ci(prnposACV/prnposWCV)
################################################################################





################################################################################
## Plot results
################################################################################

################################################################################
## Overall fitness, per clade
################################################################################
grDevices::windows(width = 15, height = 5)

Overall_fitness = ggplot(df_overall_fitness, aes(x = Time, y=Values, color = Genotype)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f2, lwd = 0.2, alpha=1, position = position_dodge2(width = 0.2), geom = "pointrange")+ 
  theme_classic()+
  facet_grid(.~Genotype)+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  scale_x_continuous(limits = c(-2,2), 
                     breaks = c(-1,1),
                     labels = c('WCV', 'ACV'))+
  scale_y_continuous(trans = 'log', limits = c(0.6,1.6),
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10), 
                     labels = c('<0.1', 0.25,0.5, 0.75, 1.0, 1.5, 2, 3, 10))+
  labs(title = "", x = '', y = 'Relative growth rate')+
  scale_color_manual(name = "Clades", values = c("royalblue", "firebrick","royalblue", "firebrick", "royalblue","firebrick"))+                                                  
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 10, hjust=0.5),
        axis.text.y = element_text(size = 10,hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

Overall_fitness
################################################################################

################################################################################
## Effect prn, overall, by clade, WCV
################################################################################
grDevices::windows(width = 4, height = 5)

df_env = data.frame(x = c(seq(0,4,1), rev(seq(0,4,1))),
                    y = c(rep(mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == -1),]$Values)[2], length(seq(0,4,1))),
                          rep(mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == -1),]$Values)[3], length(seq(0,4,1)))))

Effect_PRN_WCV = ggplot(df_effect_prn_shift[which(df_effect_prn_shift$Time == -1),], aes(x = Genotype, y = Values)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f, width = 0.4, alpha=0.5, lwd = 0.2, geom="boxplot", col = 'Black')+ 
  theme_classic()+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  scale_y_continuous(trans = 'log', limits = c(0.65,1.7), 
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2,5,10), 
                     labels = c('<0.1',0.25,0.5,0.75, 1, 1.5, 2,5,10))+
  labs(title = "Effect of PRN deficiency in WCV", x = 'Genotype', y = 'Relative growth rate')+                                               
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, hjust=0.5),
        axis.text.y = element_text(size = 10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == -1),]$Values)[1]) +
  ggplot2::geom_polygon(data = df_env,
                        ggplot2::aes(x=x, y=y),
                        fill="black",  alpha=0.1)

Effect_PRN_WCV
################################################################################

################################################################################
## Effect prn, overall, by clade, ACV
################################################################################
grDevices::windows(width = 4, height = 5)

df_env = data.frame(x = c(seq(0,4,1), rev(seq(0,4,1))),
                    y = c(rep(mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == 1),]$Values)[2], length(seq(0,4,1))),
                          rep(mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == 1),]$Values)[3], length(seq(0,4,1)))))

Effect_PRN_ACV = ggplot(df_effect_prn_shift[which(df_effect_prn_shift$Time == 1),], aes(x = Genotype, y = Values)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f, width = 0.4, alpha=0.5, lwd = 0.2, geom="boxplot", col = 'Black')+ 
  theme_classic()+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  scale_y_continuous(trans = 'log', limits = c(0.65,1.7), 
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2,5,10), 
                     labels = c('<0.1',0.25,0.5,0.75, 1, 1.5, 2,5,10))+
  labs(title = "Effect of PRN deficiency in ACV", x = 'Genotype', y = 'Relative growth rate')+                                               
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, hjust=0.5),
        axis.text.y = element_text(size = 10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = mean.and.ci(df_effect_prn_shift_averages[which(df_effect_prn_shift_averages$Time == 1),]$Values)[1]) +
  ggplot2::geom_polygon(data = df_env,
                        ggplot2::aes(x=x, y=y),
                        fill="black",  alpha=0.1)

Effect_PRN_ACV
################################################################################

################################################################################
## Plot effect vaccine switch
################################################################################
grDevices::windows(width = 6, height = 5)

df_env_PRNdef = data.frame(x = c(seq(0,7,1), rev(seq(0,7,1))),
                    y = c(rep(mean.and.ci(df_effect_vaccine_shift_averages[which(df_effect_vaccine_shift_averages$Prn == 'prnneg'),]$Values)[2], length(seq(0,7,1))),
                          rep(mean.and.ci(df_effect_vaccine_shift_averages[which(df_effect_vaccine_shift_averages$Prn == 'prnneg'),]$Values)[3], length(seq(0,7,1)))))
df_env_PRNpos = data.frame(x = c(seq(0,7,1), rev(seq(0,7,1))),
                           y = c(rep(mean.and.ci(df_effect_vaccine_shift_averages[which(df_effect_vaccine_shift_averages$Prn == 'prnpos'),]$Values)[2], length(seq(0,7,1))),
                                 rep(mean.and.ci(df_effect_vaccine_shift_averages[which(df_effect_vaccine_shift_averages$Prn == 'prnpos'),]$Values)[3], length(seq(0,7,1)))))

Vaccine_switch_effect = ggplot(df_effect_vaccine_shift, aes(x = Genotype, y = Values)) + 
  geom_abline(slope = 0, intercept = log(1), linetype = "dashed", colour = 'grey60') +
  stat_summary(fun.data = f, width = 0.4, alpha=0.5, lwd = 0.2, geom="boxplot", col = c("royalblue", "firebrick","royalblue", "firebrick", "royalblue","firebrick"))+ 
  theme_classic()+
  geom_vline(xintercept = 0, linetype = "longdash", color = 'grey')+
  scale_y_continuous(trans = 'log', limits = c(0.65,1.7), 
                     breaks = c(0.1, 0.25,0.5, 0.75, 1.0, 1.5, 2,4,10), 
                     labels = c('<0.1',0.25,0.5,0.75, 1, 1.5, 2,4,10))+
  labs(title = "Effect of vaccine switch", x = 'Genotype', y = 'Relative growth rate')+                                                
  theme(plot.title = element_text (face = 'bold',size = 10,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 10, hjust=0.5),
        axis.text.y = element_text(size = 10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") 

Vaccine_switch_effect
################################################################################


