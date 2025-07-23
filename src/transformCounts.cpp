#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".apply_transformation_difference")]]
S4 apply_transformation_difference(NumericMatrix mat) {
    int n_features = mat.nrow();
    int n_samples  = mat.ncol();
    int n_pairs    = n_features * (n_features - 1) / 2;
  
    // Check if over 1000 features
    if( n_features > 1000 ){
        Rcpp::warning("The input matrix has over 1000 features, which may "
                      "cause memory or performance issues. Consider subsetting "
                      "or using fewer features.");
    }
  
    // Generate default rownames if not present
    CharacterVector original_rownames = rownames(mat);
    bool has_valid_rownames = !Rf_isNull(original_rownames) && 
         original_rownames.size() == n_features;
  
    if( !has_valid_rownames ){
        Rcpp::warning("No rownames found in the matrix. "
                      "Generated labels like Feature1 will be used.");
    }
  
    // Prepare output vectors
    std::vector<int> i_vec, j_vec;
    std::vector<double> x_vec;
    CharacterVector rownames(n_pairs);
  
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
            // Rownames
            std::string name_i = has_valid_rownames
                ? as<std::string>(original_rownames[i])
                : "Feature" + std::to_string(i + 1);
          
            std::string name_j = has_valid_rownames
                ? as<std::string>(original_rownames[j])
                : "Feature" + std::to_string(j + 1);
          
            rownames[pair_idx] = "diff_" + name_i + "-" + name_j;
            ++pair_idx;
        }
    }
  
    // Create dgTMatrix (triplet sparse matrix)
    S4 tmat("dgTMatrix");
    tmat.slot("i") = IntegerVector(i_vec.begin(), i_vec.end());
    tmat.slot("j") = IntegerVector(j_vec.begin(), j_vec.end());
    tmat.slot("x") = NumericVector(x_vec.begin(), x_vec.end());
    tmat.slot("Dim") = IntegerVector::create(n_pairs, n_samples);
    tmat.slot("Dimnames") = List::create(rownames, colnames(mat));
  
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
  
    // Check if over 1000 features
    if( n_features > 1000 ){
        Rcpp::warning("The input matrix has over 1000 features, which may "
                      "cause memory or performance issues. Consider subsetting "
                      "or using fewer features.");
    }
    
    // Add pseudocount if needed
    if( pseudocount > 0 ){
        for( int i = 0; i < n_features; ++i ){
            for( int j = 0; j < n_samples; ++j ){
                if( mat(i, j) == 0.0 ){
                    Rcpp::warning("Zero values detected in input matrix. "
                                  "Pseudocount will be added.");
                    i = n_features;
                    break;
                }
            }
        }
        for( int i = 0; i < n_features; ++i )
            for( int j = 0; j < n_samples; ++j )
                mat(i, j) += pseudocount;
    }
    
    // Generate default rownames if not present
    CharacterVector original_rownames = rownames(mat);
    bool has_valid_rownames = !Rf_isNull(original_rownames) && 
        original_rownames.size() == n_features;
    
        if( !has_valid_rownames ){
            Rcpp::warning("No rownames found in the matrix. "
                          "Generated labels like Feature1 will be used.");
        }
    
    // Prepare output vectors
    std::vector<int> i_vec, j_vec;
    std::vector<double> x_vec;
    CharacterVector rownames(n_pairs);
  
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
            // Rownames
            std::string name_i = has_valid_rownames
                ? as<std::string>(original_rownames[i])
                : "Feature" + std::to_string(i + 1);
          
            std::string name_j = has_valid_rownames
                ? as<std::string>(original_rownames[j])
                : "Feature" + std::to_string(j + 1);
          
            rownames[pair_idx] = "div_" + name_i + "/" + name_j;
            ++pair_idx;
        }
    }
    
    // Create dgTMatrix
    S4 tmat("dgTMatrix");
    tmat.slot("i") = IntegerVector(i_vec.begin(), i_vec.end());
    tmat.slot("j") = IntegerVector(j_vec.begin(), j_vec.end());
    tmat.slot("x") = NumericVector(x_vec.begin(), x_vec.end());
    tmat.slot("Dim") = IntegerVector::create(n_pairs, n_samples);
    tmat.slot("Dimnames") = List::create(rownames, colnames(mat));
  
    // Convert dgTMatrix to dgCMatrix
    Function as("as");
    S4 cmat = as(tmat, "CsparseMatrix");
    
    // Add attribute
    cmat.attr("mia") = "division";
    return cmat;
}
