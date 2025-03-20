#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat construct_G_MRF(arma::mat G_MRF, arma::vec b, unsigned int S, unsigned int p, bool MRF_2b) {
  // arma::mat G_MRF = (G.n_rows, G.n_cols);
  if (MRF_2b) {
    // TODO: develop test for this case. Off-by-one minefield!
    // two different values for b in MRF prior for subgraphs G_ss and G_rs
    arma::uvec p_seq = arma::regspace<arma::uvec>(0, p - 1);
    for (unsigned int g = 0; g < S; g++) {
      // b1 * G_ss
      // arma::uvec g_seq(p, arma::fill::value(g * p));
      // arma::uvec idx = g_seq + p_seq;
      arma::uvec idx = p_seq + g * p;
      G_MRF.submat(idx, idx) *= b(0);
    }
    for (unsigned int g = 0; g < S - 1; g++) {
      // b2 * G_rs
      for (unsigned int r = g + 1; r < S; r++) {
        // arma::uvec g_seq(p, arma::fill::value(g * p));
        // arma::uvec r_seq(p, arma::fill::value(r * p));
        // arma::uvec idx = g_seq + p_seq;
        // arma::uvec idx2 = r_seq + p_seq;
        arma::uvec idx = p_seq + g * p;
        arma::uvec idx2 = p_seq + r * p;
        G_MRF.submat(idx, idx2) *= b(1);
        G_MRF.submat(idx2, idx) *= b(1);
      }
    }
  } else {
    // one value for b in MRF prior for all subgraphs
    G_MRF = b(0) * G_MRF;
  }
  return G_MRF;
}

arma::vec randMvNormal(const arma::vec &m, const arma::mat &Sigma) {
  unsigned int d = m.n_elem;
  //check
  if(Sigma.n_rows != d || Sigma.n_cols != d ) {
    Rcpp::stop("Dimension not matching in the multivariate normal sampler.");
  }

  arma::mat A;
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec res(d);

  if( arma::chol(A, Sigma) ) {
    res =  A.t() * Rcpp::as<arma::vec>(Rcpp::rnorm(d)) ;
  } else {
    if( eig_sym(eigval, eigvec, Sigma) ) {
      res = eigvec * arma::diagmat(arma::sqrt(eigval)) * Rcpp::as<arma::vec>(Rcpp::rnorm(d));
    } else {
      // res.fill(0.);
      Rcpp::stop("randMvNorm failing because of singular Sigma matrix.");
    }
  }
  
  return res + m;
}

// [[Rcpp::export]]
Rcpp::List func_MCMC_graph_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const unsigned int S,
  const std::string method,
  const bool MRF_2b
) {
  // Extracting data
  Rcpp::List n = sobj["n"];
  unsigned int p = Rcpp::as<unsigned int>(sobj["p"]);
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

  // Rcpp::List gamma_ini_list = Rcpp::as<Rcpp::List>(ini["gamma.ini"]);
  // arma::vec gamma_ini = Rcpp::as<arma::vec>(gamma_ini_list[0]); // FIXME: use list_to_matrix() and change with S <== Fixed by George as follows
  arma::vec gamma_ini(p * S); // vectorize/unlist ini["gamma.ini"]
  std::size_t position = 0;
  for (unsigned int g = 0; g < S; g++) {
    arma::vec component = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(ini["gamma.ini"])[g]);
    gamma_ini.subvec(position, position + component.n_elem - 1) = component;
    position += component.n_elem;
  }

  if (MRF_2b) {
    // two different values for b in MRF prior for subgraphs G_ss and G_rs
    b = Rcpp::as<arma::vec>(hyperpar["b"]);
  } else {
    double b_val = Rcpp::as<double>(hyperpar["b"]);
    b.fill(b_val);
  }
  arma::mat G_MRF = G; 
  G_MRF = construct_G_MRF(G, b, S, p, MRF_2b);

  // Update of precision matrix and graph within each subgroup
  // (analogous to SSSL algorithm for concentration (precision) graph models in
  // function 'BayesGGM_SSVS_FixedV0V1' (Wang, 2015))
  for (unsigned int g = 0; g < S; g++) { // loop through subgroups
    arma::mat V_g = V[g];
    arma::mat C_g = C[g];

    // arma::uvec p_seq = arma::regspace<arma::uvec>(0, p - 1);
    // arma::uvec g_seq(p, arma::fill::value(g * p));
    // arma::uvec idx = g_seq + p_seq;
    arma::uvec idx = arma::regspace<arma::uvec>(g * p, (g + 1) * p - 1);
    arma::mat G_g = G.submat(idx, idx);
    unsigned int n_g = Rcpp::as<unsigned int>(n[g]);
    arma::mat S_g = SSig[g];
    arma::mat Sig_g = Sig[g];

    for (unsigned int i = 0; i < p; i++) {
      arma::uvec ind_noi = arma::regspace<arma::uvec>(0, p - 1);
      ind_noi.shed_row(i);
      arma::vec v_temp = V_g.submat(ind_noi, arma::uvec({i}));
      // v_temp.elem(arma::find(v_temp < 1.0e-10)).fill(1.0e-10);

      arma::mat Sig11 = Sig_g.submat(ind_noi, ind_noi);
      arma::vec Sig12 = Sig_g.submat(ind_noi, arma::uvec({i}));

      // Omega_11^(-1)
      arma::mat invC11 = Sig11 - Sig12 * Sig12.t() / Sig_g(i, i);

      // C^(-1)
      arma::mat Ci = (S_g(i, i) + lambda) * invC11 + arma::diagmat(1. / v_temp);

      Ci = 0.5 * (Ci + Ci.t());
      // arma::mat Ci_chol = arma::chol(Ci);

      // arma::vec mu_i = -arma::solve(Ci_chol, arma::solve(Ci_chol.t(), S_g.submat(ind_noi, arma::uvec({i}))));
      // // arma::mat beta = mu_i + arma::solve(Ci_chol, arma::randn<arma::vec>(p - 1));
      // arma::vec beta = mu_i + arma::solve(Ci_chol, Rcpp::as<arma::vec>(Rcpp::rnorm(p - 1)));

      // arma::vec mu_i = -arma::inv_sympd(Ci, arma::inv_opts::allow_approx) * S_g.submat(ind_noi, arma::uvec({i})); 
      // // using inverse directly instead of chol() & solve() by Madjar
      arma::mat invCi;
      arma::vec mu_i, beta;
      if ( arma::inv_sympd(invCi, Ci) ) {
        mu_i = - invCi * S_g.submat(ind_noi, arma::uvec({i}));
        beta = randMvNormal(mu_i, Ci); 
      } else {
        // arma::inv(invCi, Ci + 0.1 * arma::eye(p - 1, p - 1), arma::inv_opts::allow_approx);
        // arma::inv(invCi, Ci);
        arma::inv(invCi, Ci, arma::inv_opts::allow_approx);
        mu_i = - invCi * S_g.submat(ind_noi, arma::uvec({i}));
        beta = mu_i;
        // beta = randMvNormal(mu_i, Ci); 
        // if (beta.has_nan()) beta = mu_i;
      }

      // Update of last column in Omega_gg
      C_g.submat(ind_noi, arma::uvec({i})) = beta;
      C_g.submat(arma::uvec({i}), ind_noi) = beta.t();

      double a_gam = 0.5 * n_g + 1;
      double b_gam = (S_g(i, i) + lambda) * 0.5;
      // double gam = arma::as_scalar(arma::randg(1, arma::distr_param(a_gam, 1 / b_gam)));
      double gam = R::rgamma( a_gam, 1. / b_gam );

      double c = arma::as_scalar(beta.t() * invC11 * beta);
      C_g(i, i) = gam + c;

      // Below updating covariance matrix according to one-column change of precision matrix

      arma::mat invC11beta = invC11 * beta;

      Sig_g.submat(ind_noi, ind_noi) = invC11 + invC11beta * invC11beta.t() / gam;
      Sig_g.submat(ind_noi, arma::uvec({i})) = -invC11beta / gam;
      Sig_g.submat(arma::uvec({i}), ind_noi) = arma::trans(-invC11beta / gam);
      Sig_g(i, i) = 1 / gam;

      // Update of variance matrix V_g and subgraph G_g

      if (i < p) {
        arma::vec beta2 = arma::zeros<arma::vec>(p);
        beta2.elem(ind_noi) = beta;

        for (unsigned int j = i + 1; j < p; j++) {
          // G where g_ss,ij=1 or 0:
          arma::mat G_prop1 = G_MRF;
          arma::mat G_prop0 = G_MRF;
          arma::uword gpj = g * p + j;
          arma::uword gpi = g * p + i;
          G_prop1(gpj, gpi) = b(0);
          G_prop1(gpi, gpj) = b(0);
          G_prop0(gpj, gpi) = 0;
          G_prop0(gpi, gpj) = 0;

          double v0 = V0(j, i);
          double v1 = V1(j, i);

          double w1 = -0.5 * std::log(v1) - 0.5 * std::pow(beta2(j), 2) / v1 + std::log(pii);
          double w2 = -0.5 * std::log(v0) - 0.5 * std::pow(beta2(j), 2) / v0 + std::log(1 - pii);

          double wa = arma::as_scalar(w1 + (a * arma::sum(gamma_ini) + gamma_ini.t() * G_prop1 * gamma_ini));
          double wb = arma::as_scalar(w2 + (a * arma::sum(gamma_ini) + gamma_ini.t() * G_prop0 * gamma_ini));

          double w_max = std::max(wa, wb);

          double w = std::exp(wa - w_max) / (std::exp(wa - w_max) + std::exp(wb - w_max));

          bool z = R::runif(0.0, 1.0) < w;
          double v = z ? v1 : v0;
          V_g(j, i) = v;
          V_g(i, j) = v;

          G_g(j, i) = static_cast<double>(z);
          G_g(i, j) = static_cast<double>(z);
        }
      }
    }
    V[g] = V_g;
    C[g] = C_g;
    G.submat(idx, idx) = G_g;
    Sig[g] = Sig_g;
  }

  if (method == "CoxBVSSL") {
    Rcpp::stop("This is not yet implemented with argument method == 'CoxBVSSL'.");
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
