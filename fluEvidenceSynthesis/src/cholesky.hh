#ifndef FLU_PROPOSAL_HH
#define FLU_PROPOSAL_HH

#include"rcppwrap.hh"
#include<RcppEigen.h>

// [[Rcpp::Export]]
Eigen::MatrixXd cholesky_factorization(
        const Eigen::MatrixXd &A);

#endif

 
