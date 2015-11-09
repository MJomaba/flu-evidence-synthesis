#ifndef INFERENCE_HH
#define INFERENCE_HH

#include "rcppwrap.h"

namespace flu {
    struct mcmc_result_inference_t
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
            batch;
        Eigen::VectorXd llikelihoods;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 
            Eigen::RowMajor>
            contact_ids;
    };
}
#endif
