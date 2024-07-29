#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat construct_G_MRF(arma::mat G, arma::vec b, uint S, int p, bool MRF_2b) {
  arma::mat G_MRF(G.n_rows, G.n_cols);
  if (MRF_2b) {
    // TODO: develop test for this case. Off-by-one minefield!
    // two different values for b in MRF prior for subgraphs G_ss and G_rs
    arma::uvec p_seq = arma::regspace<arma::uvec>(0, p - 1);
    for (uint g = 0; g < S; g++) {
      // b1 * G_ss
      arma::uvec g_seq(p, arma::fill::value(g * p));
      arma::uvec idx = g_seq + p_seq;
      G_MRF.submat(idx, idx) *= b(0);
    }
    for (uint g = 0; g < S - 1; g++) {
      // b2 * G_rs
      for (uint r = g; r < S - 1; r++) {
        arma::uvec g_seq(p, arma::fill::value(g * p));
        arma::uvec r_seq(p, arma::fill::value(r * p));
        arma::uvec idx = g_seq + p_seq;
        arma::uvec idx2 = r_seq + p_seq;
        G_MRF.submat(idx, idx2) *= b(1);
        G_MRF.submat(idx2, idx) *= b(1);
      }
    }
  } else {
    // one value for b in MRF prior for all subgraphs
    G_MRF = b[0] * G;
  }
  return G_MRF;
}

// [[Rcpp::export]]
Rcpp::List func_MCMC_graph_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const uint S,
  const std::string method,
  const bool MRF_2b
) {
  // Extracting data
  Rcpp::List n = sobj["n"];
  uint p = Rcpp::as<uint>(sobj["p"]);
  Rcpp::List SSig = sobj["SSig"];

  double pii = Rcpp::as<double>(hyperpar["pi.G"]);
  double a = Rcpp::as<double>(hyperpar["a"]);
  arma::vec b(2);
  double lambda = Rcpp::as<double>(hyperpar["lambda"]);
  arma::mat V0 = Rcpp::as<arma::mat>(hyperpar["V0"]);
  arma::mat V1 = Rcpp::as<arma::mat>(hyperpar["V1"]);

  arma::mat G = Rcpp::as<arma::mat>(ini["G.ini"]);
  Rcpp::List V = Rcpp::as<Rcpp::List>(ini["V.ini"]);
  Rcpp::List Sig = Rcpp::as<Rcpp::List>(ini["Sig.ini"]);
  Rcpp::List C = Rcpp::as<Rcpp::List>(ini["C.ini"]);

  Rcpp::List gamma_ini_list = Rcpp::as<Rcpp::List>(ini["gamma.ini"]);
  arma::vec gamma_ini = Rcpp::as<arma::vec>(gamma_ini_list[0]);

  if (MRF_2b) {
    // two different values for b in MRF prior for subgraphs G_ss and G_rs
    b = Rcpp::as<arma::vec>(hyperpar["b"]);
  } else {
    b = {hyperpar["b"], hyperpar["b"]};
  }
  arma::mat G_MRF = construct_G_MRF(G, b, S, p, MRF_2b);

  // Update of precision matrix and graph within each subgroup
  // (analogous to SSSL algorithm for concentration (precision) graph models in
  // function 'BayesGGM_SSVS_FixedV0V1' (Wang, 2015))
  for (uint g = 0; g < S; g++) { // loop through subgroups
    arma::vec V_g = V[g];
    arma::vec C_g = C[g];

    arma::uvec p_seq = arma::regspace<arma::uvec>(0, p - 1);
    arma::uvec g_seq(p, arma::fill::value(g * p));
    arma::uvec idx = g_seq + p_seq;
    arma::mat G_g = G.submat(idx, idx);
    uint n_g = Rcpp::as<uint>(n[g]);
    arma::mat S_g = SSig[g];
    arma::mat Sig_g = Sig[g];

    // TODO: code i loop through genes

    V[g] = V_g;
    C[g] = C_g;
    G.submat(idx, idx) = G_g;
    Sig[g] = Sig_g;
  }

  // Assembling output
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("Sig.ini") = Sig,
    Rcpp::Named("G.ini") = G,
    Rcpp::Named("V.ini") = V,
    Rcpp::Named("C.ini") = C
  );
  return out;
}
