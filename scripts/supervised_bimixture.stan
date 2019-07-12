// Fits a mixture of two normals, where membership in a mixture is informed by the feature matrix X
// Feature matrix X is assuemd to have fixed effects

data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;  // response values
  matrix[n, p] X; // matrix of genetic information
  matrix[n, 2] m_X; //matrix of class indicators
}

parameters {
  vector[p] B;
  ordered[2] mu;
  real<lower=0> sigma;
  //real<lower=0, upper=1> theta;
}

model {
  mu ~ student_t(3, 0, 3);
  sigma ~ normal(0, 3); // 

  // genetic info coeffs
  B ~ student_t(3, 0, sqrt(7));  // for now, no shrinking bc using PCA/PLS

  y ~ normal(m_X * mu + X * B, sigma);
}
