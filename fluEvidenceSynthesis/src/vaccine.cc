#include "rcppwrap.hh"
#include "vaccine.hh"

#include "model.hh"
#include "data.hh"


/*namespace flu {
    namespace vaccine {
    };
};*/

// [[Rcpp::export]]
bool vaccinationScenario( std::vector<size_t> age_sizes, 
        vaccine_t vaccine_calendar, state_t sample ) {

    double result_simu[7644];
    double FinalSize[21];

    auto pop_vec = flu::data::separate_into_risk_groups( 
                flu::data::group_age_data( age_sizes ) );

    Rcpp::Rcout << vaccine_calendar.efficacy_year[0] << std::endl;
    
    //auto vac_cal = Rcpp::as< flu::vaccine::vaccine_t >(vaccine_calendar);

    //translate into an initial infected population
    double init_inf[NAG];
    for(size_t i=0;i<NAG;i++)
        init_inf[i]=pow(10,state.parameters.init_pop);


    one_year_SEIR_with_vaccination(result_simu, pop_vec, prop_init_inf, state.time_latent, state.time_infectious, state.parameters.susceptibility, contact_mat, state.parameters.transmissibility, vaccine_calendar);
    for(size_t j=0; j<21; j++)
    {
        FinalSize[j]=0.0;
    }
    for(size_t i=0;i<364;i++)
    {
        for(size_t j=0; j<21; j++)
        {
            FinalSize[j]+=result_simu[21*i+j];
        }
    }

    return true;
}

