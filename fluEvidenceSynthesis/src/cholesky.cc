#include "cholesky.hh"

// [[Rcpp::export]]
Eigen::MatrixXd cholesky( SEXP A ) {
    return cholesky_factorization( Rcpp::as<Eigen::MatrixXd>(A) );
}

Eigen::MatrixXd cholesky_factorization(
        const Eigen::MatrixXd &A)
{
    assert( A.rows() == A.cols() );

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(A.rows(), A.cols() );

    for(size_t i=0;i<A.rows();i++)
    {
        for(size_t j=i;j<A.cols();j++)
        {
            double sum_L2=A(i,j);
            for(size_t k=0;k<i;k++)
                sum_L2-=res(i,k)*res(j,k);
            if(i==j)
                res(i,i)=sqrt(sum_L2);
            else
                res(j,i)=sum_L2/res(i,i);
        }
    }
    return res;
}
