data {
  int<lower=1> N; // number of observations
  int<lower=1> N_test; // number of test observations
  int<lower=1> B2; // number of environments (to cluster)
  int<lower=1> B1; // number of genotypes
  vector[N] y; // response variable
  int g1[N]; // genotype variable
  vector[N_test] y_test; // test dependent variable
  int g1_test[N_test]; // test genotype variable
  vector<lower=0>[B2] alpha; // dirichlet prior parameter
  vector[N] var2_1; // continuous environmental variable for observations
  vector[N_test] var2_1_test; // continuous environmental variable for test observations
}

parameters {
  simplex[B2] Pi; // cluster membership probabilities
  ordered[B2] mu_g2; // environment effect parameters
  vector[B1] b1; // genotype effect parameters
  vector[B2] b2; // environment effect parameters
  real<lower=0> sigma; // standard deviation for phenotype
  real<lower=0> sigma_var2_1; // standard deviation for environment
  real<lower=0> sigma_b2; // standard deviation for environment effect
  real<lower=0> sigma_b1; // standard deviation for genotype effect
  //vector<lower = 0>[N] g2; // environment parameter for observations
  //vector[N_test] g2_test; // environment parameter for test observations
}

transformed parameters{
  vector[B2] Pinew = Pi/sum(Pi);
}

model {
  int g2[N]; // environment parameter for observations
  int g2_test[N_test]; // environment parameter for test observations

  Pi ~ dirichlet(alpha); // prior on cluster membership probabilities
  mu_g2 ~ normal(0, 100); // prior on environment effect parameters
  b1 ~ normal(0, sigma_b1); // prior on genotype effect parameters
  b2 ~ normal(0, sigma_b2); // prior on environment effect parameters
  sigma ~ cauchy(0, 10);
  sigma_var2_1 ~ cauchy(0, 10);
  sigma_b1 ~ cauchy(0, 10);
  sigma_b2 ~ cauchy(0, 10);

  for (i in 1:N) {
    g2 ~ categorical(Pinew);
    y[i] ~ normal(b1[g1[i]] + b2[g2[i]], sigma);
    var2_1[i] ~ normal(mu_g2[g2[i]], sigma_var2_1);
  }

  for(i in 1:N_test) {
    y_test[i] ~ normal(b1[g1_test[i]] + b2[g2_test[i]], sigma);
    g2_test[i] ~ categorical(Pinew);
    var2_1_test[i] ~ normal(mu_g2[g2_test[i]], sigma_var2_1);
  }
}
