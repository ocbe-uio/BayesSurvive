#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List settingInterval_cpp(const arma::vec y, const arma::vec delta_, const arma::vec s_, const unsigned int J_)
{
    // set a finite partition of the time axis to define the indicator matrices for risk sets and failure sets
    // in order to calculate the increment in the cumulative baseline hazard in each interval
    // also in order to construct the grouped data likelihood
    
    arma::mat ind_d_ = arma::zeros<arma::mat>(y.n_elem, J_);
    arma::mat ind_r_ = ind_d_;

    arma::uvec case0yleq;
    arma::uvec case0ygeq;
    arma::uvec case1yleq;
    arma::uvec case1ygeq;

    double smax = max(s_);
    case0yleq = arma::find(delta_ == 0. && y <= smax);
    case0ygeq = arma::find(delta_ == 0. && y > smax);
    case1yleq = arma::find(delta_ == 1. && y <= smax);
    case1ygeq = arma::find(delta_ == 1. && y > smax);

    int cen_j;
    for (unsigned int i = 0; i < case1yleq.n_elem; ++i)
    {
        cen_j = min(arma::find(s_ >= y(case1yleq(i))));
        ind_d_(case1yleq(i), cen_j) = 1.;
        ind_r_.submat(case1yleq(i), 0, case1yleq(i), cen_j).fill(1.);
        // Rcout << cen_j << ",";
    }

    for (unsigned int i = 0; i < case0yleq.n_elem; ++i)
    {
        cen_j = min(arma::find(s_ >= y(case0yleq(i))));
        ind_r_.submat(case0yleq(i), 0, case0yleq(i), cen_j).fill(1.);
    }

    if (case0ygeq.n_elem + case1ygeq.n_elem > 0)
    {
        arma::uvec union_case01ygeq;
        union_case01ygeq = arma::unique(arma::join_cols(case0ygeq, case1ygeq));
        ind_r_.rows(union_case01ygeq).fill(1.);
    }

    arma::mat ind_r_d_ = ind_r_ - ind_d_;
    arma::vec d_ = arma::sum(ind_d_.t(), 1);
    
    return Rcpp::List::create(
        Rcpp::Named("s") = s_,
        Rcpp::Named("J") = J_,
        Rcpp::Named("ind.r") = ind_r_,
        Rcpp::Named("ind.d") = ind_d_,
        Rcpp::Named("d") = d_,
        Rcpp::Named("ind.r_d") = ind_r_d_
      );
}
