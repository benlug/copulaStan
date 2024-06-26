functions {
  #include clayton_cop_loglik.stan
}

data {
  int<lower=0> N;                  // number of observations
  vector[N] y1;                    // variable 1
  vector[N] y2;                    // variable 2
  int<lower=1, upper=2> dist1;     // distribution type for variable 1
  int<lower=1, upper=2> dist2;     // distribution type for variable 2
}

parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real<lower=0> theta;  // dependence parameter theta for the Clayton copula
}

model {
  vector[N] x1;
  vector[N] x2;

  sigma ~ cauchy(0, 2);
  theta ~ exponential(1); // prior for theta - lognormal good alternative

  // Marginal distributions
  if (dist1 == 1) {
    y1 ~ normal(mu[1], sigma[1]);
    for (n in 1:N) {
      x1[n] = normal_cdf(y1[n], mu[1], sigma[1]);
    }
  } else if (dist1 == 2) {
    y1 ~ lognormal(mu[1], sigma[1]);
    for (n in 1:N) {
      x1[n] = lognormal_cdf(y1[n], mu[1], sigma[1]);
    }
  }

  if (dist2 == 1) {
    y2 ~ normal(mu[2], sigma[2]);
    for (n in 1:N) {
      x2[n] = normal_cdf(y2[n], mu[2], sigma[2]);
    }
  } else if (dist2 == 2) {
    y2 ~ lognormal(mu[2], sigma[2]);
    for (n in 1:N) {
      x2[n] = lognormal_cdf(y2[n], mu[2], sigma[2]);
    }
  }

  for (n in 1:N) {
    target += log(clayton_copula_density(x1[n], x2[n], theta));
  }
}
