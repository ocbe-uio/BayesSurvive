  #include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat list_to_matrix(Rcpp::List r_list) {
  // Determine the number of columns (length of the list)
  int n_cols = r_list.size();

  // Determine the number of rows (length of the first vector in the list)
  int n_rows = Rcpp::as<arma::vec>(r_list[0]).n_elem;

  // Initialize the matrix
  arma::mat result(n_rows, n_cols);

  // Fill the matrix
  for (int j = 0; j < n_cols; ++j) {
    arma::vec col_vec = Rcpp::as<arma::vec>(r_list[j]);
    result.col(j) = col_vec;
  }

  return result;
}
