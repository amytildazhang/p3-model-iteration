#
# TODO: comment / repeat for each function... since this is where some of the key
# computations will be performed, it's worth explaining things a bit more verbosely
# here / including any useful references.
# 
# For documenting functions, you can follow http://r-pkgs.had.co.nz/man.html#man-workflow
# I added the section headers (overview, etc.) for convenience..
#

#'
#' SHORT DESC
#'
#' Overview
#' --------
#'
#' LONGER DESC
#' 
#' Parameters
#' ----------
#'
#' @param stanfit <type> <description>
#' @param sampnames 
#' @param X
#' @param resp
#'
#' Returns
#' -------
#'
#' @return <desc>
#'
#' References
#' ----------
#' 1. xx
#' 2. yy
#'

# TODO: use underscores for all variables/function names 
draw_post.bimodal <- function(stanfit, sampnames, X, resp = NULL) {
    # bivariate mixture model
    samps <- rstan::extract(stanfit)

    # TODO: comment each basic step inside function definitions.. this is mostly for 
    # my benefit since a lot of this
    # is going to be new to me... :)
    prob_mix <- pnorm(X %*% t(samps$B))
    mixtures <- map(1:2, function(a) {
        rnorm(nrow(samps$mu), mean = samps$mu[, a], sd = samps$sigma[, a])
    })
    
    yhats <- sapply(1:nrow(prob_mix), function(i) {
        prob_mix[i, ] * mixtures[[1]] + (1 - prob_mix[i, ]) * mixtures[[2]]
    })

    logcpoi <- sapply(1:nrow(prob_mix), function(i) {
        mn <- prob_mix[i, ] * samps$mu[,1] + (1 - prob_mix[i, ]) * samps$mu[,2]
        sd <- prob_mix[i, ] * samps$sigma[,1] + (1 - prob_mix[i, ]) * samps$sigma[,2]
        dnorm(resp[i], mean = mn, sd = sd, log = T)
    })
    
    colnames(yhats) <- paste0("yhat_", sampnames)
    colnames(logcpoi) <- paste0("lcpoi_", sampnames)

    return(cbind(yhats, logcpoi))
}


# draw post predictive samples
draw_post.linear <- function(armfit, sampnames, X,  resp = NULL) {
    # bivariate mixture model
    samps <- rstan::extract(armfit$stanfit)
    coeffs <- cbind(samps$alpha, samps$beta)
    mu <- (X %*% t(coeffs))
    
    if (nrow(X) > 1) {
        yhats <- sapply(1:nrow(X), function(i) {
            rnorm(ncol(mu), mean = mu[i, ], sd = samps$aux)
        })
        
        logcpoi <- sapply(1:nrow(X), function(i) {
            dnorm(resp[i], mean = mu[i, ], sd = samps$aux, log = T)
        })
        
    } else {
        yhats <- cbind(rnorm(length(mu), mean = mu, sd = samps$aux))
        logcpoi <- cbind(as.vector(dnorm(resp, mean = mu, sd = samps$aux, log = T)))
    }    
    colnames(yhats) <- paste0("yhat_", sampnames)
    colnames(logcpoi) <- paste0("lcpoi_", sampnames)

    return(cbind(yhats, logcpoi))
}
