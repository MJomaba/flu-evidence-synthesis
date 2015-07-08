#include "rcppwrap.hh"
#include<RcppEigen.h>

#include "proposal.hh"

/**
 * This file contains some thin wrappers around C++ code. Mostly used
 * for testing the C++ code from R.
 *
 * The wrappers are needed, because in the C++ code we often want to pass
 * by const reference, which is impossible in a function called from R
 */

// [[Rcpp::export]]
Eigen::VectorXd updateMeans( Eigen::VectorXd means,
        Eigen::VectorXd v, size_t n )
{
    return flu::proposal::updateMeans( means, v, n );
}

// [[Rcpp::export]]
Eigen::MatrixXd updateCovariance( Eigen::MatrixXd cov, 
        Eigen::VectorXd v, Eigen::VectorXd means, size_t n )
{
    return flu::proposal::updateCovariance( cov, v, means, n );
}
