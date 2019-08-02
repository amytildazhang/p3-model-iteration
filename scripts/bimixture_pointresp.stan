// Fits a mixture of two normals, where membership in a mixture is informed by the feature matrix X
// Feature matrix X is assumed to have fixed effects



data {
// initialize required data for the model

  int<lower=0> n;       // number of samples. syntax means: integer named 'n' with value >= 0
  int<lower=0> p;	// number of features
  vector[n] y; 		// response values, column vector of dimension n. by default, any real number
  matrix[n, p] X; 	// model matrix X with dimensions n x p
}

parameters {
// initialize parameters to estimate
  vector[p] B;		// coefficients
  ordered[2] mu;	// means for the two normal mixtures, ordered to prevent label switching:
                        // https://mc-stan.org/docs/2_18/stan-users-guide/mixture-modeling-chapter.html
  real<lower=0> sigma[2];  // std dev parameters for the two normal mixtures, must be positive
}

transformed parameters {
// initialize & define any parameters which are functions of others
// typically faster to do any matrix calculation in this block rather than in the model block
// REF 

// probability for each data point that it comes from the first distribution
// defined essentially as a binomial GLM
    vector<lower=0, upper=1>[n] theta = Phi(X * B); 	

}

model {
// define the model density by placing priors on all parameters in the parameters block


  // means for the normal mixtures
  // fairly narrow prior, because the response is scaled (for consistency)
  // place a student's t distribution instead of normal because normal tends to over-shrink
  mu ~ student_t(3, 0, 3);   // documentation: ref 4a

  // std devs for the normal mixtures
  // again fairly narrow prior because of the scale of the response
  sigma ~ normal(0, 3);      // since sigma > 0, this becomes a half-normal on the standard deviation 

  // broader prior for the probability of being in either mixture
  // the features are not standardised int eh way the response is
 // not a very 'uninformative' prior -- see ref 3
  B ~ student_t(3, 0, sqrt(10));  // 

  // specifies how the parameters are combined in the joint density
  // Stan tries to maximize the log likelihood of the target
  for (i in 1:n) {
    // mixture model is defined as recommended, see ref 4b
    target += log_mix(theta[i],
     normal_lpdf(y[i] | mu[1], sigma[1]),
     normal_lpdf(y[i] | mu[2], sigma[2]));
  }
}


// References
// 1. STAN User reference
// 2. Gelman, Andrew. "Prior distributions for variance parameters in hierarchical models." 
//    Bayesian analysis 1.3 (2006): 515-534.  
// 3. Gelman, Andrew, Daniel Simpson, and Michael Betancourt. 
//    "The prior can often only be understood in the context of the likelihood."
//     Entropy 19.10 (2017): 555.
// 4. Specific STAN reference pages:
//      a. https://mc-stan.org/docs/2_18/functions-reference/student-t-distribution.html
//      b. https://mc-stan.org/docs/2_19/stan-users-guide/vectorizing-mixtures.html


