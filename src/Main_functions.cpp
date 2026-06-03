#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int samp_com_for_birth_cpp(NumericMatrix x) {
  int nrow = x.nrow(); int ncol = x.ncol();
  NumericVector sums(nrow);
  NumericVector probabilities(nrow);
  
  double csize = 0.0;
  
  for (int i = 0; i < nrow; i++) {
    double total = 0.0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
      csize += x(i, j);
    }
    sums[i] = total;
  }
  
  for (int i = 0; i < nrow; i++) {
    probabilities[i] = sums[i] / csize;
  }
  
  IntegerVector vec(nrow);
  for (int i = 0; i < nrow; i++) {
    vec[i] = i;
  }
  
  IntegerVector out_vec = Rcpp::sample(vec, 1, false, probabilities);
  return out_vec[0];
}

// [[Rcpp::export]]
NumericVector calculate_species_probs_cpp(NumericMatrix x, int index, List params) {
  int ncol = x.ncol();
  NumericVector freqs(ncol);
  double csize = 0.0;
  
  for (int i = 0; i < ncol; i++) {
    csize += x(index, i);
  }
  
  for (int i = 0; i < ncol; i++) {
    freqs[i] = x(index, i) / csize;
  }
  
  NumericVector fd = as<NumericVector>(params["fd"]);
  NumericMatrix s = as<NumericMatrix>(params["s"]);
  
  NumericVector coeff(ncol);
  NumericVector fds(ncol);
  NumericVector ss(ncol);
  
  for (int i = 0; i < ncol; i++) {
    double fd_clause = fd[i] * (freqs[i] - (1.0/ncol));
    fds[i] = fd_clause;
    double sel_clause = log(s(index, i));
    ss[i] = sel_clause;
    coeff[i] = exp(fd_clause + sel_clause);
  }
  
  NumericVector probs(ncol);
  double denom = 0.0;
  
  for (int i = 0; i < ncol; i++) {
    denom += coeff[i]*freqs[i];
  }
  
  for (int i = 0; i < ncol; i++) {
    probs[i] = (coeff[i]*freqs[i])/denom;
  }
  
  return probs;
  
}


// [[Rcpp::export]]
int samp_species_for_birth_cpp(NumericMatrix x, NumericVector probs) {
  int ncol = x.ncol();
  IntegerVector vec(ncol);
  
  for (int i = 0; i < ncol; i++) {
    vec[i] = i;
  }
  
  IntegerVector out_vec = Rcpp::sample(vec, 1, false, probs);
  return out_vec[0];
} 


// [[Rcpp::export]]
int samp_com_for_death_cpp(NumericMatrix x,
                       int species_index,
                       int com_index,  
                       List params) {
  int nrow = x.nrow();
  IntegerVector vec(nrow);
  
  for (int i = 0; i < nrow; i++) {
    vec[i] = i;
  }
  
  Rcpp::List migmats = as<Rcpp::List>(params["mig"]);
  NumericMatrix mat = as<NumericMatrix>(migmats[species_index]);
  NumericVector probs = mat(_,com_index);
  
  IntegerVector out_vec = Rcpp::sample(vec, 1, false, probs);
  return out_vec[0];
} 

// [[Rcpp::export]]
int samp_species_for_death_cpp(NumericMatrix x,
                               int com_index) {
  int ncol = x.ncol();
  NumericVector vec = x(com_index, _);
  
  double n_ind = 0.0;
  
  for (int i = 0; i < ncol; i++) {
    n_ind += vec[i];
  }
  
  NumericVector probs(ncol);
  
  for (int i = 0; i < ncol; i++) {
    probs[i] = vec[i]/n_ind;
  }
  
  IntegerVector ind_vec(ncol);
  
  for (int i = 0; i < ncol; i++) {
    ind_vec[i] = i;
  }
  
  IntegerVector out_vec = Rcpp::sample(ind_vec, 1, false, probs);
  return out_vec[0];
}


// [[Rcpp::export]]
NumericMatrix update_meta_cpp(NumericMatrix x,
                              int com_index,
                              int b_species_index,
                              int d_species_index) {
  x(com_index, d_species_index) = x(com_index, d_species_index) - 1;
  x(com_index, b_species_index) = x(com_index, b_species_index) + 1;
  
  return x;
}


// [[Rcpp::export]]
NumericMatrix birth_death_process_cpp(NumericMatrix x,
                                      Rcpp::List params) {
  int b_com = samp_com_for_birth_cpp(x);
  NumericVector s_probs = calculate_species_probs_cpp(x, b_com, params);
  int s_birth = samp_species_for_birth_cpp(x, s_probs);
  int d_com = samp_com_for_death_cpp(x, s_birth, b_com, params);
  int s_death = samp_species_for_death_cpp(x, d_com);
  x = update_meta_cpp(x, d_com, s_birth, s_death);
  return x;
}


// [[Rcpp::export]]
int sample_pos_spec(NumericVector x) {
  int vec_size = x.size();
  IntegerVector pos_spec;
  
  double pos_spec_sum = 0.0;
  for (int i = 0; i < vec_size; i++) {
    if (x[i] > 1) {
      pos_spec.push_back(i);
      pos_spec_sum += x[i];
    }
  }
  
  int n_pos_spec = pos_spec.size();

  NumericVector pos_spec_probs(n_pos_spec);
  
  for (int i=0; i < n_pos_spec; i++) {
    pos_spec_probs[i] = x[pos_spec[i]]/pos_spec_sum;
  }
  
  IntegerVector sel_spec = Rcpp::sample(pos_spec, 1, false, pos_spec_probs);
  return sel_spec[0];
}


// [[Rcpp::export]]
Rcpp::List speciate_cpp(NumericMatrix x,
                  Rcpp::List params,
                  bool change_params = false,
                  double prop_new = 0.00,
                  double diff_sd = 0.05) {
  int n_coms = x.nrow();
  int n_spec = x.ncol();
  IntegerVector coms(n_coms);
  NumericVector fd = as<NumericVector>(params["fd"]);
  NumericMatrix s = as<NumericMatrix>(params["s"]);
  List migmats = as<List>(params["mig"]);

  for (int i = 0; i < n_coms; i++) {
    coms[i] = i;
  }
  
  Rcpp::NumericVector colsums(n_spec);
  
  for (int i = 0; i < n_spec; i++) {
    double tmpsum = 0.0;
    for (int j = 0; j < n_coms; j++) {
      tmpsum += x(j, i);
    }
    colsums[i] = tmpsum;
  }
  
  for (int i = 0; i < n_spec; i++) {
    if (colsums[i] == 0) {
      IntegerVector com_ind = Rcpp::sample(coms, 1, false);

      double com_size = 0.0;
      
      for (int k = 0; k < n_spec; k++) {
        com_size += x(com_ind[0], k);
      }
      
      int n_to_add = std::ceil(prop_new*com_size);
      int n_new = 0;
      
      while (n_new < n_to_add) {
        int sel_spec = sample_pos_spec(x(com_ind[0], _));
        x(com_ind[0], sel_spec) = x(com_ind[0], sel_spec) - 1;
        n_new += 1;
      }
      
      x(com_ind[0], i) = n_new;
      
      if (change_params) {
        for (int j = 0; j < n_coms; j++) {
          double dif = R::rnorm(0, diff_sd);
          s(j, i) = s(j, i) + dif;
          if (s(j, i) <= 0) {
            s(j, i) = 0.00000001;
          }
        }
      
        double dif_fd = R::rnorm(0, diff_sd);
        fd[i] = fd[i] + dif_fd;
        if (fd[i] >= 0.0) {
          fd[i] = 0.0;
        }
      }
    }
  }
  
  List Lparams = List::create(Named("s") = s , _["fd"] = fd, _["mig"] = migmats);
  List L = List::create(x, Lparams);
  return L;
  
}


// [[Rcpp::export]]
NumericVector calc_col_freqs(NumericMatrix x) {

  int ncol = x.ncol(); 
  double J = sum(x);
  
  NumericVector out(ncol);
  
  for (int j = 0; j < ncol; j++) {
    double colsum = sum(x(_, j));
    out[j] = colsum/J;
  }
  
  return out;
  
}

// [[Rcpp::export]]
Rcpp::List moran_deme_cpp(NumericMatrix x, int t, Rcpp::List params,
                          bool change_params = false,
                          double prop_new = 0.0,
                          double diff_sd = 0.05) {
  
  int ncol = x.ncol();
  int nrow = x.nrow();
  int J = sum(x);
  
  List params_out = List::create(x, params);
  
  NumericMatrix outfreqs(t, ncol);
  
  outfreqs(0, _) = calc_col_freqs(x);
  
  int ncycles = t*J;
  int Gen = 1;
  
  for (int i = 1; i < ncycles; i++) {
    
    x = birth_death_process_cpp(x, params);
    
    if (i % J == 0) {
      speciate_cpp(x, params, change_params, prop_new, diff_sd);
      outfreqs(Gen, _) = calc_col_freqs(x);
      Gen += 1;
    }
  }
  
  List L = List::create(Named("metacommunity") = x, _["metafreq"] = outfreqs, _["input"] = params_out);
  return L;
}


