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
//' @description Superseded by \code{vaccination_scenario}
//'
//' @param age_sizes A vector with the population size by each age {1,2,..}
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param polymod_data Contact data for different age groups
//' @param contact_ids IDs (row numbers) of the contact data used when modelling this scenario 
//' @param parameters The parameters to use
//' 
//' @keywords internal
//'
//' @seealso \code{\link{vaccination_scenario}}
//' 
//' @return A data frame with the total number of influenza cases in that year
//'
// [[Rcpp::export]]
std::vector<double> vaccinationScenario( std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data, std::vector<size_t> contact_ids,
        Eigen::VectorXd parameters ) {
    ::Rf_warning("\'vaccinationScenario\' is deprecated\nUse \'vaccination_scenario\' instead.\nSee help(Deprecated).");

    std::vector<size_t> age_group_limits = {1,5,15,25,45,65};
    auto age_data = flu::data::group_age_data( age_sizes, age_group_limits );

    auto nag = age_data.size();

    Eigen::MatrixXd risk_proportions = Eigen::MatrixXd( 
            2, nag );
    risk_proportions << 
        0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
        0, 0, 0, 0, 0, 0, 0;

    auto pop_vec = flu::data::separate_into_risk_groups( age_data, 
            risk_proportions );

    //auto vac_cal = Rcpp::as< flu::vaccine::vaccine_t >(vaccine_calendar);


    //translate into an initial infected population
    auto init_inf = Eigen::VectorXd::Constant( 
            nag, pow(10,parameters[8] ) );


    flu::data::age_data_t ages;
    ages.age_sizes = age_sizes;
    ages.age_group_sizes = age_data;

    auto contact_matrix = flu::contacts::to_symmetric_matrix( 
            flu::contacts::shuffle_by_id( 
                flu::contacts::table_to_contacts( polymod_data,
                    age_group_limits ), 
                contact_ids ), ages );

    auto time_latent = 0.8;
    auto time_infectious = 1.8;
    Eigen::VectorXd susc( 7 );
    susc << parameters[5], parameters[5], parameters[5], 
         parameters[6], parameters[6], parameters[6], 
         parameters[7]; 




    auto result_simu = flu::one_year_SEIR_with_vaccination(pop_vec, init_inf, time_latent, time_infectious, susc, contact_matrix, parameters[4], 
            vaccine_calendar)
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

