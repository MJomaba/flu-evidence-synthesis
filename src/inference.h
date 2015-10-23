#ifndef INFERENCE_HH
#define INFERENCE_HH

namespace flu {
    mcmc_result_t binomial_inference( std::vector<size_t> age_sizes, 
            Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
            Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
            flu::vaccine::vaccine_t vaccine_calendar,
            flu::contacts::contacts_t polymod_data,
            Eigen::VectorXd initial, 
            int nburn = 10000, 
            int nbatch = 1000, int blen = 1 );
}
#endif
