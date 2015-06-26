#ifndef FLU_PROPOSAL_HH
#define FLU_PROPOSAL_HH

#include<Rcpp.h>
#include<RcppEigen.h>

// [[Rcpp::Export]]
Eigen::MatrixXd cholesky_factorization(
        const Eigen::MatrixXd &A);

#endif

 
