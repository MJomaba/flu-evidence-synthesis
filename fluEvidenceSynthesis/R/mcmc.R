adaptive.mcmc <- function( lprior, llikelihood, nburn, initial, nbatch, blen = 1, ...  )
{
    adaptive.mcmc.cpp( lprior, function(pars) llikelihood( pars, ... ),
                  nburn, initial, nbatch, blen )
}
