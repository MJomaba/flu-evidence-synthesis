#include "rcppwrap.hh"
#include "vaccine.hh"

#include "model.hh"
#include "data.hh"
#include "contacts.hh"

/*
namespace flu {
    namespace vaccine {
    };
};*/

//' Calculate number of influenza cases given a vaccination strategy
// [[Rcpp::export]]
std::vector<double> vaccinationScenario( std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar, flu::state_t sample,
        flu::contacts::contacts_t polymod_uk ) {

    auto age_data = flu::data::group_age_data( age_sizes );
    auto pop_vec = flu::data::separate_into_risk_groups( age_data );

    //auto vac_cal = Rcpp::as< flu::vaccine::vaccine_t >(vaccine_calendar);

    //translate into an initial infected population
    double init_inf[NAG];
    for(size_t i=0;i<NAG;i++)
        init_inf[i]=pow(10,sample.parameters.init_pop);


    flu::data::age_data_t ages;
    ages.age_sizes = age_sizes;
    ages.age_group_sizes = age_data;

    auto contact_matrix = flu::contacts::to_symmetric_matrix( 
            flu::contacts::shuffle_by_id( polymod_uk, 
                sample.contact_ids ), ages );

    auto result_simu = flu::one_year_SEIR_with_vaccination(pop_vec, init_inf, sample.time_latent, sample.time_infectious, sample.parameters.susceptibility, contact_matrix, sample.parameters.transmissibility, vaccine_calendar)
        .cases;

    auto final_sizes = std::vector<double>( result_simu.cols(), 0.0 ); 

    for(size_t i=0;i<result_simu.rows();i++)
    {
        for(size_t j=0; j<result_simu.cols(); j++)
        {
            final_sizes[j]+=result_simu(i,j);
        }
    }

    return final_sizes;
}

