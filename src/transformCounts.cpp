#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".apply_transformation_difference")]]
S4 apply_transformation_difference(NumericMatrix mat) {
    int n_features = mat.nrow();
    int n_samples  = mat.ncol();
    int n_pairs    = n_features * (n_features - 1) / 2;
  
    // Prepare output vectors
    std::vector<int> i_vec, j_vec;
    std::vector<double> x_vec;
  
    int pair_idx = 0;
    for( int i = 0; i < n_features - 1; ++i ){
        for( int j = i + 1; j < n_features; ++j ){
            for( int k = 0; k < n_samples; ++k ){
                double diff = mat(i, k) - mat(j, k);
                if( diff != 0.0 ){
                    i_vec.push_back(pair_idx);
                    j_vec.push_back(k);
                    x_vec.push_back(diff);
                }
            }
            ++pair_idx;
        }
    }
  
    // Create dgTMatrix (triplet sparse matrix)
    S4 tmat("dgTMatrix");
    tmat.slot("i") = IntegerVector(i_vec.begin(), i_vec.end());
    tmat.slot("j") = IntegerVector(j_vec.begin(), j_vec.end());
    tmat.slot("x") = NumericVector(x_vec.begin(), x_vec.end());
    tmat.slot("Dim") = IntegerVector::create(n_pairs, n_samples);
    tmat.slot("Dimnames") = List::create(R_NilValue, colnames(mat));
  
    // Convert dgTMatrix to dgCMatrix
    Function as("as");
    S4 cmat = as(tmat, "CsparseMatrix");
  
    // Add attribute
    cmat.attr("mia") = "difference";
    return cmat;
}

// [[Rcpp::export(name = ".apply_transformation_division")]]
S4 apply_transformation_division(NumericMatrix mat, double pseudocount = 1e-6) {
    int n_features = mat.nrow();
    int n_samples  = mat.ncol();
    int n_pairs    = n_features * (n_features - 1) / 2;
    
    // Prepare output vectors
    std::vector<int> i_vec, j_vec;
    std::vector<double> x_vec;
  
    int pair_idx = 0;
    for( int i = 0; i < n_features - 1; ++i ){
        for( int j = i + 1; j < n_features; ++j ){
            for( int k = 0; k < n_samples; ++k ){
                double denom = mat(j, k);
                if( denom != 0.0 && std::isfinite(denom) ){
                    double ratio = mat(i, k) / denom;
                    if( ratio != 0.0 && std::isfinite(ratio) ){
                        i_vec.push_back(pair_idx);
                        j_vec.push_back(k);
                        x_vec.push_back(ratio);
                    }
                }
            }
            ++pair_idx;
        }
    }
    
    // Create dgTMatrix
    S4 tmat("dgTMatrix");
    tmat.slot("i") = IntegerVector(i_vec.begin(), i_vec.end());
    tmat.slot("j") = IntegerVector(j_vec.begin(), j_vec.end());
    tmat.slot("x") = NumericVector(x_vec.begin(), x_vec.end());
    tmat.slot("Dim") = IntegerVector::create(n_pairs, n_samples);
    tmat.slot("Dimnames") = List::create(R_NilValue, colnames(mat));
  
    // Convert dgTMatrix to dgCMatrix
    Function as("as");
    S4 cmat = as(tmat, "CsparseMatrix");
    
    // Add attribute
    cmat.attr("mia") = "division";
    return cmat;
}
