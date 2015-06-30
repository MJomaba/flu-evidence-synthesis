#include "rcppwrap.hh"

template <> flu::vaccine::vaccine_t Rcpp::as( SEXP rVac )
{
    flu::vaccine::vaccine_t vac_cal;
    auto rListVac = Rcpp::as<List>(rVac);
    Rcpp::NumericVector eff = rListVac["efficacy"];
    for (size_t i = 0; i < vac_cal.efficacy_year.size(); ++i)
        vac_cal.efficacy_year[i] = eff[i];

    Rcpp::DataFrame cal = Rcpp::as<DataFrame>(rListVac["calendar"]);
    //Rcpp::Rcout << cal(0,0) << std::endl;

    for(size_t i=0;i<21;++i) {
        for(size_t j=0;j<123;j++)
        {
            NumericVector column = cal[i];
            vac_cal.calendar[j*21+i] = column[j];
        }
    }
    return vac_cal;
}

/* 
/// Parameters of the model
struct parameter_set {
    std::vector<double> epsilon = std::vector<double>(5);
    double psi;
    double transmissibility;
    std::vector<double> susceptibility = std::vector<double>(7);
    double init_pop;
};*/
template <> parameter_set Rcpp::as( SEXP rPar ) 
{
    flu::parameter_set pars;
    auto rListPar = Rcpp::as<List>(rPar);
    pars.epsilon = Rcpp::as<std::vector<double> >( rListPar["epsilon"] );
    pars.psi = rListPar["psi"];
    pars.transmissibility = rListPar["transmissibility"];
    pars.susceptibility = Rcpp::as<std::vector<double> >( 
            rListPar["susceptibility"] );
    pars.init_pop = rListPar["init_pop"];
    return pars;
}

/// Keeps the current state of the model/mcmc
/*struct state_t 
{
    parameter_set parameters;
    double time_infectious, time_latent;

    std::vector<double> positivity_ij = std::vector<double>(260);
    std::vector<size_t> contact_ids;
};*/
template <> state_t Rcpp::as( SEXP rState )
{
    flu::state_t state;
    auto rListState = Rcpp::as<List>(rState);

    state.parameters = Rcpp::as<parameter_set>( rListState["parameters"] );
    state.time_infectious = rListState["time_infectious"];
    state.time_latent = rListState["time_latent"];
    state.contact_ids = Rcpp::as<std::vector<size_t> >( 
            rListState["contact_ids"] );

    return state;
}
