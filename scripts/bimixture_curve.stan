// Fits a mixture of two normals, where membership in a mixture is informed by the feature matrix X
// Feature matrix X is assuemd to have fixed effects

data {
  int<lower=0> n;  // num cell lines
  int<lower=0> p;  // num parameters
  int<lower=0> d;  // num dosages
  int<lower=0> ns; // num spline basis functions

  matrix[n, d] y;  // response values
  matrix[n, p] X; // matrix of genetic information
  matrix[d, ns] fd;    // basis values for dosage response
}

parameters {
  vector[p] B;
  real<lower=0> sigma[2];
  vector[ns] S1;
  vector[ns] S2;
  ordered[2] mu;
 // vector[n] c_eff;

}

transformed parameters {
    vector<lower=0, upper=1>[n] theta = Phi(X * B);
    vector[d] mfunc1 = fd * S1 + mu[1]; // issue: 
    vector[d] mfunc2 = fd * S2 + mu[2];
}

model {
 // c_eff ~ student_t(3, 0, 4);
  mu ~ student_t(3, 0, 10);

  sigma[1] ~ cauchy(0, 10);
  S1 ~ student_t(3, 0, 8);
  
  sigma[2] ~ cauchy(0, 8);
  S2 ~ student_t(3, 0, 10);

  // genetic info coeffs
  B[1] ~ student_t(3, 0, sqrt(10.0));
  B[2:p] ~ student_t(3, 0, 7);  // for now, no shrinking bc using PCA/PLS

  for (i in 1:n) {
    target += log_mix(theta[i],
     normal_lpdf(y[i,] | mfunc1, sigma[1]),
     normal_lpdf(y[i,] | mfunc2, sigma[2]));
  }
}

