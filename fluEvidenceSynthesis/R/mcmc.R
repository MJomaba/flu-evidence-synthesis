#' Adaptive MCMC algorithm
#'
#' MCMC which adapts its proposal distribution for faster convergence following:
#' Sherlock, C., Fearnhead, P. and Roberts, G.O. The Random Walk Metrolopois: Linking Theory and Practice Through a Case Study. Statistical Science 25, no.2 (2010): 172-190.
#'
#' @param lprior A function returning the log prior probability of the parameters 
#' @param llikelihood A function returning the log likelihood of the parameters given the data
#' @param nburn Number of iterations of burn in
#' @param initial Vector with starting parameter values
#' @param nbatch Number of batches to run (number of samples to return)
#' @param blen Length of each batch
#' @param ... Extra parameters passed to the log likelihood function
#' 
#' @return Returns a list with the accepted samples and the corresponding llikelihood values
#'
#' @seealso \code{\link{adaptive.mcmc}} Used internally by this function.
adaptive.mcmc <- function( lprior, llikelihood, nburn, initial, nbatch, blen = 1, ...  )
{
    adaptive.mcmc.cpp( lprior, function(pars) llikelihood( pars, ... ),
                  nburn, initial, nbatch, blen )
}
