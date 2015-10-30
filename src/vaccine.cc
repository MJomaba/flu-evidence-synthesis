#include "rcppwrap.h"
#include "vaccine.h"

#include "model.h"
#include "data.h"
#include "contacts.h"

/*
namespace flu {
    namespace vaccine {
    };
};*/

//' Calculate number of influenza cases given a vaccination strategy
//'
//' @param age_sizes A vector with the population size by each age {1,2,..}
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param polymod_data Contact data for different age groups
//' @param sample The parameters needed to run the ODE model (typically one of the posterior sample created when running the inference)
//' @return A data frame with the total number of influenza cases in that year
//'
// [[Rcpp::export]]
std::vector<double> vaccinationScenario( std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar,
        flu::contacts::contacts_t polymod_data,
        flu::state_t sample ) {

    auto age_data = flu::data::group_age_data( age_sizes );
    auto pop_vec = flu::data::separate_into_risk_groups( age_data );

    //auto vac_cal = Rcpp::as< flu::vaccine::vaccine_t >(vaccine_calendar);

    auto nag = vaccine_calendar.efficacy_year.size();

    //translate into an initial infected population
    auto init_inf = Eigen::VectorXd::Constant( 
            nag, pow(10,sample.parameters.init_pop) );


    flu::data::age_data_t ages;
    ages.age_sizes = age_sizes;
    ages.age_group_sizes = age_data;

    auto contact_matrix = flu::contacts::to_symmetric_matrix( 
            flu::contacts::shuffle_by_id( polymod_data, 
                sample.contact_ids ), ages );

    auto result_simu = flu::one_year_SEIR_with_vaccination(pop_vec, init_inf, sample.time_latent, sample.time_infectious, sample.parameters.susceptibility, contact_matrix, sample.parameters.transmissibility, vaccine_calendar)
        .cases;

    auto final_sizes = std::vector<double>( result_simu.cols(), 0.0 ); 

    for(int i=0;i<result_simu.rows();i++)
    {
        for(int j=0; j<result_simu.cols(); j++)
        {
            final_sizes[j]+=result_simu(i,j);
        }
    }

    return final_sizes;
}

