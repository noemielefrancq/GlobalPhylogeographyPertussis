// ##################################################################
// ##              Estimation of relative fitness                  ##
// ##                         Stan code                            ##
// ##                                                              ##
// ## Author: Noemie Lefrancq                                      ##
// ## Email: ncmjl2@cam.ac.uk                                      ##
// ## Last update: 28/01/2022                                      ##
// ##                                                              ##
// ##                                                              ##  
// ##################################################################


functions {
	real[] logistic_model(real y0, real fitness_clades_pre_vacc, real fitness_clades_post_vacc, int nb_years, int vaccine_introduction, int ACV_introduction){
		real res[nb_years];
		real tmp;
    tmp = y0;
	  if(tmp < 1E-30) {tmp = 1E-30;}
	  if(tmp >= 1-1E-30) {tmp = 1-1E-30;}
	  for(i in 1:nb_years){
	      res[i] = tmp;
	      if(i <= vaccine_introduction && i<= ACV_introduction){
	        res[i] = tmp;
	      }
	      if(i > vaccine_introduction && i <=  ACV_introduction){
	        if(fitness_clades_pre_vacc<=0){ // 1 - usual logistic
	          res[i] = 1 - 1/(1+(1/(1-tmp)-1)*exp(fitness_clades_pre_vacc*(1)));
	        }
	        if(fitness_clades_pre_vacc > 0){ // usual logistic
	          res[i] = 1/(1+(1/tmp-1)*exp(-fitness_clades_pre_vacc*(1)));
	        }
  	    }
  	    if(i > vaccine_introduction && i > ACV_introduction){
  	      if(fitness_clades_post_vacc<=0){ // 1 - usual logistic
	          res[i] = 1 - 1/(1+(1/(1-tmp)-1)*exp(fitness_clades_post_vacc*(1)));
	        }
	        if(fitness_clades_post_vacc > 0){ // usual logistic
	          res[i] = 1/(1+(1/tmp-1)*exp(-fitness_clades_post_vacc*(1)));
	        }
  	    }
    	  tmp = res[i];
	      if(res[i] <= 1E-30) {res[i] = 1E-30;}
	      if(res[i] >= 1) {res[i] = 1-1E-30;}
	  }
  	return (res);
	}
}

data {
  // Shape data
	int nb_genotypes; // number of clades
	int nb_years; // length of the data 
	int nb_countries; // number of countries 
  
  // Data
	int data_genotype_non_ref[nb_genotypes-1, nb_years, nb_countries]; // counts
	int data_genotype_ref[nb_years, nb_countries]; // counts 
	int data_total_number[nb_years, nb_countries]; // reference for each year // real because of the divison in lp for the freq
	int non_zero_country_year[nb_years, nb_countries];
	int non_zero_country_year_genotype[nb_genotypes-1, nb_years, nb_countries];
	int number_zeros_country_year;
	int number_zeros_country_year_genotype;
	
  // ACV implementation
	int yearIntroduction[nb_countries]; // dates of introdution for each country
	
	// F0 estimation
	int yearF0[nb_countries];
	
	// Fitness vector
  int R_every_pre_vacc;
  int R_every_post_vacc;
  int number_R_pre_vacc;
  int number_R_post_vacc;
}

parameters {
	real<lower=-5, upper=5> fitness_genotypes_pre_vacc[nb_genotypes-1, number_R_pre_vacc];
	real<lower=-5, upper=5> fitness_genotypes_post_vacc[nb_genotypes-1, number_R_post_vacc];
	
	simplex[nb_genotypes] f[nb_countries];
	real<lower=0.5, upper=1-1E-30> alpha;
	
  real<lower=0.25, upper=10> overdispersion_inv;
}

transformed parameters{
  real<lower=-5, upper=5> fitness_genotypes_vector_pre_vacc[nb_genotypes-1, number_R_pre_vacc];
	real<lower=-5, upper=5> fitness_genotypes_vector_post_vacc[nb_genotypes-1, number_R_post_vacc];
	
	simplex[nb_genotypes] f0[nb_countries];
		
	real<lower=0, upper=1> pred_paired_freq[nb_countries,nb_genotypes-1, nb_years] ; // pred frequencies
	real pred_paired_freq_tmp[nb_countries,nb_genotypes, nb_years] ; // pred numbers
	real pred_absolute_freq[nb_countries,nb_genotypes, nb_years] ; // pred numbers
	
	real pred_number_non_ref[nb_countries,nb_genotypes-1, nb_years] ; // pred numbers non ref
	real pred_number_ref[nb_countries, nb_years] ; // pred numbers ref
	
  real<lower=0, upper=10> overdispersion;
	
	// overdispersion
	overdispersion = 1/sqrt(overdispersion_inv);

  // compute strating frequencies
	for(j in 1:nb_countries){
	  f0[j,1] = alpha + f[j,1]*(1-alpha);
	  for(i in 2:nb_genotypes){
	    f0[j,i] = f[j,i]*(1-alpha);
	  }
	}
	
	// construct fitness vector
	fitness_genotypes_vector_pre_vacc = fitness_genotypes_pre_vacc;
	fitness_genotypes_vector_post_vacc = fitness_genotypes_post_vacc;
  
  // compute predicted numbers
	for (k in 1:nb_countries){
	  // compute pred fequencies
	  for (j in 1:(nb_genotypes-1)){
	    pred_paired_freq[k,j,] = logistic_model(f0[k,j]/(f0[k,j]+f0[k,nb_genotypes]), fitness_genotypes_vector_pre_vacc[j,1], fitness_genotypes_vector_post_vacc[j,1], nb_years, yearF0[k], yearIntroduction[k]);
	  }
    for(l in 1:nb_years){
  	    for(j in 1:(nb_genotypes-1)){
            pred_paired_freq_tmp[k,j,l] = (pred_paired_freq[k,j,l])*10000/(1-pred_paired_freq[k,j,l]);
            if(pred_paired_freq_tmp[k,j,l] < 0) pred_paired_freq_tmp[k,j,l] = 0;
    	  }
    	  pred_paired_freq_tmp[k,nb_genotypes,l] = 10000;
    	  for(j in 1:nb_genotypes){
    	    pred_absolute_freq[k,j,l] = 0;
    	    pred_absolute_freq[k,j,l] = pred_paired_freq_tmp[k,j,l]/sum(pred_paired_freq_tmp[k,,l]);
    	    if(pred_absolute_freq[k,j,l]<0) pred_absolute_freq[k,j,l] = 0;
    	  }
    	  
    	  for(j in 1:(nb_genotypes-1)){
    	    pred_number_non_ref[k,j,l] = pred_paired_freq[k,j,l]*(data_genotype_non_ref[j,l,k] + data_genotype_ref[l,k]); 
    	    if (pred_number_non_ref[k,j,l] <1E-30) {
    		  	pred_number_non_ref[k,j,l] = 1E-30;
    		  }
	      }
  	    // Compute ref clade numbers
    	  pred_number_ref[k,l] = 0;
    	  pred_number_ref[k,l] = data_total_number[l,k] - sum(pred_number_non_ref[k,,l]);
    	  if (pred_number_ref[k,l] <1E-30) {
    		  pred_number_ref[k,l] = 1E-30;
    		}
    }
	}
}


model {
  overdispersion_inv ~ normal(0, 1);
  
  for (i in 1:(nb_genotypes-1)){
    fitness_genotypes_post_vacc[i,1] ~ cauchy(0, 0.15);
    fitness_genotypes_pre_vacc[i,1] ~ cauchy(fitness_genotypes_post_vacc[i,1], 0.15);
  }
  
  alpha ~ normal(1, 0.1);
  
	for (k in 1:nb_countries){
    for(l in 1:nb_years){
  	    for (j in 1:(nb_genotypes-1)){
  	      if(non_zero_country_year_genotype[j,l,k]) {
  	 		    data_genotype_non_ref[j,l,k] ~ neg_binomial_2(pred_number_non_ref[k,j,l], pow(pred_number_non_ref[k,j,l], overdispersion));
  	      }
	  	  }
    }
	}
}

generated quantities {
  vector[nb_countries*(nb_genotypes-1)*nb_years-number_zeros_country_year_genotype] log_lik;
  int n;
  n = 1;
  for (k in 1:nb_countries){
	  for(l in 1:nb_years){
	      for (j in 1:(nb_genotypes-1)){
	        if(non_zero_country_year_genotype[j,l,k]) {
    	      log_lik[n] = 1E-30;
  	 		    log_lik[n] = neg_binomial_2_lpmf(data_genotype_non_ref[j,l,k]| pred_number_non_ref[k,j,l], pow(pred_number_non_ref[k,j,l], overdispersion));
  	 		    n=n+1;
	        }
	      }
  	}
	}
}

