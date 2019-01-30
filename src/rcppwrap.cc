#include "rcppwrap.h"
#include<RcppEigen.h>

template <> flu::vaccine::vaccine_t Rcpp::as( SEXP rVac )
{
    flu::vaccine::vaccine_t vac_cal;
    auto rListVac = Rcpp::as<List>(rVac);
    Rcpp::NumericVector eff = rListVac["efficacy"];

    auto nag = eff.size();

    if (eff.size() == vac_cal.calendar.cols())
    {
        if( eff.size()%3==0 )
            nag /= 3; // assume 3 risk groups
        if( eff.size()%2==0 )
            nag /= 2; // assume 2 risk groups
    }
    vac_cal.efficacy = Eigen::VectorXd::Zero( 3*nag );
    for (int i = 0; i < nag; ++i)
        if (nag == eff.size()) {
            vac_cal.efficacy[i] = eff[i];
            vac_cal.efficacy[i+nag] = eff[i];
            vac_cal.efficacy[i+2*nag] = eff[i];
        }
        else 
            vac_cal.efficacy[i] = eff[i];

    vac_cal.calendar = Rcpp::as<Eigen::MatrixXd>(rListVac["calendar"]);

    // HARDCODED risk groups!
    if (vac_cal.calendar.cols() < 3*nag) {
        auto dim = vac_cal.calendar.cols();
        vac_cal.calendar.conservativeResize( vac_cal.calendar.rows(), 
                3*nag );
        for (size_t j=dim; j<vac_cal.calendar.cols(); ++j)
        {
            for (size_t i = 0; i < vac_cal.calendar.rows(); ++i)
                vac_cal.calendar(i,j)=0;
        }
    }

    if (rListVac.containsElementNamed("dates")) 
    {
        // Ideally would work with Datetime, but automatic conversion of
        // date vector to datetime vector seems to fail. Even though
        // R does not seem (double check) separate datetime type
        auto rDates = Rcpp::as<Rcpp::DateVector>( rListVac["dates"] )
            .getDates();

        for( auto rDate : rDates )
        {
            // TODO: use rapi date_to_ptime
            vac_cal.dates.push_back(
                    boost::posix_time::ptime(
                        boost::gregorian::date( rDate.getYear(),
                            rDate.getMonth(), rDate.getDay() ),
                        boost::posix_time::time_duration(12, 0, 0 ) 
                        )
                    );
        }
    }
    return vac_cal;
}

template <> flu::parameter_set Rcpp::as( SEXP rPar ) 
{
    flu::parameter_set pars;
    auto rListPar = Rcpp::as<List>(rPar);
    pars.epsilon = Rcpp::as<Eigen::VectorXd>( rListPar["epsilon"] );
    pars.psi = rListPar["psi"];
    pars.transmissibility = rListPar["transmissibility"];
    pars.susceptibility = Rcpp::as<Eigen::VectorXd>( 
            rListPar["susceptibility"] );
    pars.init_pop = rListPar["init_pop"];
    return pars;
}

template <> SEXP Rcpp::wrap( const flu::parameter_set &parameters )
{
    List rPar;
    rPar["epsilon"] = wrap(parameters.epsilon);
    rPar["psi"] = wrap(parameters.psi);
    rPar["transmissibility"] = wrap(parameters.transmissibility);
    rPar["susceptibility"] = wrap(parameters.susceptibility);
    rPar["init_pop"] = wrap(parameters.init_pop);
    return Rcpp::wrap(rPar);
}

/// Keeps the current state of the model/mcmc
/*struct state_t 
{
    parameter_set parameters;
    double time_infectious, time_latent;

    std::vector<double> positivity_ij = std::vector<double>(260);
    std::vector<size_t> contact_ids;
};*/
template <> flu::state_t Rcpp::as( SEXP rState )
{
    flu::state_t state;
    auto rListState = Rcpp::as<List>(rState);

    state.parameters = Rcpp::as<parameter_set>( rListState["parameters"] );
    state.time_infectious = rListState["time_infectious"];
    state.time_latent = rListState["time_latent"];
    state.contact_ids = Rcpp::as<std::vector<size_t> >( 
            rListState["contact_ids"] );
    if (rListState.containsElementNamed("likelihood"))
        state.likelihood = rListState["likelihood"];

    return state;
}

template <> SEXP Rcpp::wrap( const flu::state_t &sample )
{
    List rState;
    rState["parameters"] = wrap( sample.parameters );
    rState["time_infectious"] = wrap( sample.time_infectious );
    rState["time_latent"] = wrap( sample.time_latent );
    rState["contact_ids"] = wrap( sample.contact_ids );
    rState["likelihood"] = wrap( sample.likelihood );
    return Rcpp::wrap(rState);
}

template <> SEXP Rcpp::wrap( const flu::mcmc_result_inference_t &mcmcResult )
{
    Rcpp::List rState;
    rState["batch"] = Rcpp::wrap( mcmcResult.batch );
    rState["llikelihoods"] = Rcpp::wrap( mcmcResult.llikelihoods );
    rState["contact.ids"] = Rcpp::wrap( mcmcResult.contact_ids );
    return rState;
}
