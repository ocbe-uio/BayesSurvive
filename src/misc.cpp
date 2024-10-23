  #include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat list_to_matrix(Rcpp::List r_list) {
  // Converts a *list of vectors* into a matrix
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

//[[Rcpp::export]]
arma::cube list_to_cube(Rcpp::List r_list) {
  // Converts a *list of matrices* into a cube

  // Determine the number of slices (length of the list)
  int n_slices = r_list.size();

  // Determine the number of rows (length of the first matrix in the list)
  int n_rows = Rcpp::as<arma::mat>(r_list[0]).n_rows;
  int n_cols = Rcpp::as<arma::mat>(r_list[0]).n_cols;

  // Initialize the cube
  arma::cube result(n_rows, n_cols, n_slices);

  // Fill the cube
  for (int g = 0; g < n_slices; ++g) {
    arma::mat mat = Rcpp::as<arma::mat>(r_list[g]);
    result.slice(g) = mat;
  }

  return result;
}
