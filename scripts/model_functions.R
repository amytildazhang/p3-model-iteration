
#'
#' Draws posterior predictive samples from a bimodal normal mixture model, 
#' specifically for those fit using the \code{bimixture_pointresp.stan} parameterization.
#' --------
#'
#' The bimodal normal mixture model must be defined so that the normal densities 
#' are represented by 2-dimensional vector \code{mu} for their means and 2-dimensional vector
#' \code{sigma} for their standard deviations. Any coefficients which are used to determine
#' the mixing coefficient are represented by the P-dimensional vector \code{B}. 
#' 
#' Parameters
#' ----------
#'
#' @param stanfit A \code{stanfit} object returned from fitting an \code{rstan} model.
#' @param sampnames A character vector of distinguishing names for observations
#' @param X The \code{NxP} design matrix of observations to draw posterior predictive samples for, 
#'         where \code{N} is the number of observations and \code{P} is the number of features. 
#'          Should match the format of the design matrix used in running the \code{rstan} model.
#' @param resp  A numeric vector of response values corresponding to the observations in \code{X}.
#'
#' Returns
#' -------
#'
#' @return Two concatenated numeric matrices corresponding to posterior predictive samples 
#'          (columns beginning with \code{yhat_}) and the log conditional predictive ordinate (CPO_i) 
#'          (columns beginning with \code{lcpoi_}). Each row represents a draw from the 
#'          posterior predictive density of \code{stanfit}. Columns are named based off 
#'          of the sample names given by \code{sampnames}. 
#'          The number of rows (draws) corresponds to the number of samples kept in 
#'          Hamiltonian Monte Carlo, i.e. the number of samples from calling \code{rstan::extract}.
#'
#' References
#' ----------
#' 1. Gelman, Andrew, Xiao-Li Meng, and Hal Stern. "Posterior predictive assessment of model fitness 
#'    via realized discrepancies." Statistica sinica (1996): 733-760.
#' 2. Gelfand A.E.,  Dey D.K.,  Chang H.. Model determination using predictive distributions with
#'    implementation via sampling-based methods (with discussion), 1992 Department of Statistics,
#'    Stanford University, Tech. Rep., Stanford, California, 462
#'

draw_post_bimodal <- function(stanfit, sampnames, X, resp = NULL) {
    # draw posterior samples for all parameters within the stanfit object
    samps <- rstan::extract(stanfit)

     # calculate posterior probability of sample in X belonging to the first (lower) mode 
     # 'mixing coefficient'
    prob_mix <- pnorm(X %*% t(samps$B))
  
    # draw posterior predictive samples from the two modes of the bimodal mixture
    mixtures <- map(1:2, function(a) {
        rnorm(nrow(samps$mu), mean = samps$mu[, a], sd = samps$sigma[, a])
    })
    
   # use mixing coefficient to weight the samples drawn from the two modes
   # and obtain posterior predictive samples for each observation in X 
    yhats <- sapply(1:nrow(prob_mix), function(i) {
        prob_mix[i, ] * mixtures[[1]] + (1 - prob_mix[i, ]) * mixtures[[2]]
    })

    # Similar to above, gather posterior samples for the means of the two modes 
    # to get posterior distribution for each observation in X
    # and obtain the likelihood of true data defined in `resp` under this posterior distribution.
    logcpoi <- sapply(1:nrow(prob_mix), function(i) {
        mn <- prob_mix[i, ] * samps$mu[,1] + (1 - prob_mix[i, ]) * samps$mu[,2]
        sd <- prob_mix[i, ] * samps$sigma[,1] + (1 - prob_mix[i, ]) * samps$sigma[,2]
        dnorm(resp[i], mean = mn, sd = sd, log = T)
    })
    
    # label columns
    colnames(yhats) <- paste0("yhat_", sampnames)
    colnames(logcpoi) <- paste0("lcpoi_", sampnames)

    # return as one matrix for ease
    return(cbind(yhats, logcpoi))
}


#'
#'
#' Draws posterior predictive samples from a normal linear model fit
#' using the \code{R} package \code{rstanarm}.
#' --------
#'
#' Although \code{rstanarm} provides its own posterior prediction function, for unknown
#' reasons, its results differ strikingly from those provided by other packages, hence
#' this implementation. 
#'
#' 
#' Parameters
#' ----------
#'
#' @param armfit A \code{armfit} object returned from fitting an \code{rstanarm} linear model.
#' @param sampnames A character vector of distinguishing names for observations
#' @param X The \code{NxP} design matrix of observations to draw posterior predictive samples for, 
#'         where \code{N} is the number of observations and \code{P} is the number of features. 
#'          Should match the format of the design matrix used in running the \code{rstan} model.
#' @param resp  A numeric vector of response values corresponding to the observations in \code{X}.
#'
#' Returns
#' -------
#'
#' @return Two concatenated numeric matrices corresponding to posterior predictive samples 
#'          (columns beginning with \code{yhat_}) and the log conditional predictive ordinate (CPO_i) 
#'          (columns beginning with \code{lcpoi_}). Each row represents a draw from the 
#'          posterior predictive density of \code{stanfit}. Columns are named based off 
#'          of the sample names given by \code{sampnames}. 
#'          The number of rows (draws) corresponds to the number of samples kept in 
#'          Hamiltonian Monte Carlo, i.e. the number of samples from calling \code{rstan::extract}.
#'
#' References
#' ----------
#' 1. Gelman, Andrew, Xiao-Li Meng, and Hal Stern. "Posterior predictive assessment of model fitness 
#'    via realized discrepancies." Statistica sinica (1996): 733-760.
#' 2. Gelfand A.E.,  Dey D.K.,  Chang H.. Model determination using predictive distributions with
#'    implementation via sampling-based methods (with discussion), 1992 Department of Statistics,
#'    Stanford University, Tech. Rep., Stanford, California, 462
#'
draw_post_linear <- function(armfit, sampnames, X,  resp = NULL) {
    # draw posterior samples for all parameters within the stanfit object
    samps <- rstan::extract(armfit$stanfit)
  
    # calculate posterior mean for each observation in X
    coeffs <- cbind(samps$alpha, samps$beta)
    mu <- (X %*% t(coeffs))
    
    # treat 1-sample case separately as objects become vectors rather than matrices
    if (nrow(X) > 1) {
        # For each observation in X, draw a predictive sample 'yhat' from the posterior mean
        yhats <- sapply(1:nrow(X), function(i) {
            rnorm(ncol(mu), mean = mu[i, ], sd = samps$aux)
        })
        
        # For each observation in X, calculate log likelihood of true resposnse 
        # under the posterior mean
        logcpoi <- sapply(1:nrow(X), function(i) {
            dnorm(resp[i], mean = mu[i, ], sd = samps$aux, log = T)
        })
        
    } else {
        yhats <- cbind(rnorm(length(mu), mean = mu, sd = samps$aux))
        logcpoi <- cbind(as.vector(dnorm(resp, mean = mu, sd = samps$aux, log = T)))
    }    
  
     # label columns
    colnames(yhats) <- paste0("yhat_", sampnames)
    colnames(logcpoi) <- paste0("lcpoi_", sampnames)

     # return as one matrix
     return(cbind(yhats, logcpoi))
}
