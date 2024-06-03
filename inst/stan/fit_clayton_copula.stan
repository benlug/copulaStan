functions {
  real clayton_copula_log_density(real u, real v, real theta) {
    real log_c = log1p(theta);
    real term1 = (-theta - 1) * (log(u) + log(v));
    real term2 = -(2 * theta + 1) / theta * log(pow(u, -theta) + pow(v, -theta) - 1);

    return log_c + term1 + term2;
  }
}

data {
  int<lower=0> N;                  // number of observations
  vector<lower=0, upper=1>[N] x1;  // variable 1
  vector<lower=0, upper=1>[N] x2;  // variable 2
}

parameters {
  real<lower=0> theta;  // dependence parameter theta for the Clayton copula
}

model {
  theta ~ exponential(1); // prior for theta - lognormal good alternative

  for (n in 1:N) {
    target += clayton_copula_log_density(x1[n], x2[n], theta);
  }
}

generated quantities {
  matrix[2, 2] Omega;
  Omega[1, 1] = 1;
  Omega[1, 2] = 0;
  Omega[2, 1] = 0;
  Omega[2, 2] = 1;
}
