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
