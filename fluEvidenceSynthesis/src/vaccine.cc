#include "vaccine.hh"

#include "model.hh"
#include "data.hh"

namespace flu {
    namespace vaccine {
    };
};

//using namespace flu::vaccine;
// Outside of namespace so that R(cpp) can access it
// [[Rcpp::export]]
bool vaccination_scenario( std::vector<size_t> age_sizes, SEXP vaccine_calendar ) {

    double result_simu[7644];
    double FinalSize[21];

    auto pop_vec = flu::data::separate_into_risk_groups( 
                flu::data::group_age_data( age_sizes ) );
    
    auto vac_cal = Rcpp::as< flu::vaccine::vaccine_t >(vaccine_calendar);

    /*one_year_SEIR_with_vaccination(result_simu, pop_vec, prop_init_inf, state.time_latent, state.time_infectious, state.parameters.susceptibility, contact_mat, state.parameters.transmissibility, vaccine_scenarios[scen]);
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
    }*/

    return true;
}

#include <Rcpp.h>
namespace Rcpp {
    template <> flu::vaccine::vaccine_t as( SEXP rVac )
    {
        flu::vaccine::vaccine_t vac_cal;
        auto rListVac = as<List>(rVac);
        NumericVector eff = rListVac["efficacy"];
        for (size_t i = 0; i < vac_cal.efficacy_year.size(); ++i)
            vac_cal.efficacy_year[i] = eff[i];

        Rcout << eff[0] << std::endl;

        NumericMatrix cal = rListVac["calendar"];

        for(size_t i=0;i<21;++i) {
            for(size_t j=0;j<123;j++)
            {
                vac_cal.calendar[j*21+i] = cal(j,i);
            }
        }
        return vac_cal;
    }
}
