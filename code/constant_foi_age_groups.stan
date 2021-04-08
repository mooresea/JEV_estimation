// Calculate FOI from age-specific incidence
data {
  int<lower=0> num_ages; //number of age classes
  int<lower=0> age_range;
  //int<lower=0> cases; //total number of cases
  
  int cases_age[num_ages]; //cases in each age class
  simplex[num_ages] age_distr; //pct of population in each age
  int lower_bound[num_ages]; //lower bound of age group
  int upper_bound[num_ages];  //upper bound of age group
 
  }

parameters{
  real<lower=0> lambdaA; //foi
}

model {
  vector[num_ages] prob_age; //probabilities of cases being in each  age class
  vector[num_ages] lambda_age; //
  int min_age;
  int max_age;
  //real lambda_temp;
  
  //lambdaA ~ beta(2,30);
  lambdaA ~ normal(0,10) T[0.0,];
  
  //lambda_age = gamma_age .* age_distr;
  for (iter in 1:num_ages){
    min_age = lower_bound[iter];
    max_age = upper_bound[iter]+1;
    prob_age[iter] = exp(-lambdaA*min_age)-exp(-lambdaA*max_age); 
  }
  //divide by total probability because a small fraction may be infected after max age (particularly if FOI is low)
  //prob_age = prob_age / sum(prob_age); 
  
  //Account for age distribution of population
  lambda_age = prob_age .* age_distr;
  lambda_age = lambda_age / sum(lambda_age);

  //statistical model of our observations
  cases_age ~ multinomial(lambda_age);  

}

generated quantities{
  vector[num_ages] prob_age; //probabilities of cases being in each  age class
  vector[num_ages] lambda_age; //
  int min_age;
  int max_age;
  int exp_cases_age[num_ages];
  
  //lambda_age = gamma_age .* age_distr;
  for (iter in 1:num_ages){
    min_age = lower_bound[iter];
    max_age = upper_bound[iter]+1;
    prob_age[iter] = exp(-lambdaA*min_age)-exp(-lambdaA*max_age); 
  }

  //Account for age distribution of population
  lambda_age = prob_age .* age_distr;
  lambda_age = lambda_age / sum(lambda_age);
  
  exp_cases_age = multinomial_rng(lambda_age,sum(cases_age));

}
