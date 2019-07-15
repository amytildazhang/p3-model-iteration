// Fits a mixture of two normals, where membership in a mixture is informed by the feature matrix X
// Feature matrix X is assuemd to have fixed effects

data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;  // response values
  matrix[n, p] X; // matrix of genetic information
}

parameters {
  vector[p] B;
  ordered[2] mu;
  real<lower=0> sigma[2];
  //real<lower=0, upper=1> theta;
}

transformed parameters {
    vector<lower=0, upper=1>[n] theta = Phi(X * B);

}

model {
  mu ~ student_t(3, 0, 3);
  sigma ~ normal(0, 3); // 

  // genetic info coeffs
  B[1] ~ student_t(3, 0, sqrt(10.0));
  B[2:p] ~ student_t(3, 0, sqrt(7));  // for now, no shrinking bc using PCA/PLS

  for (i in 1:n) {
    target += log_mix(theta[i],
     normal_lpdf(y[i] | mu[1], sigma[1]),
     normal_lpdf(y[i] | mu[2], sigma[2]));
  }
}
