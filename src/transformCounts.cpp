#include <Rcpp.h>
using namespace Rcpp;

// Helper for difference
void calculate_difference(
        const NumericMatrix& mat, int i, int j, int k, int pair_idx,
        std::vector<int>& i_vec, std::vector<int>& j_vec,
        std::vector<double>& x_vec) {
    double diff = mat(k, i) - mat(k, j);
    if( diff != 0.0 ){
        i_vec.push_back(k);
        j_vec.push_back(pair_idx);
        x_vec.push_back(diff);
    }
}

// Helper for division
void calculate_division(
        const NumericMatrix& mat, int i, int j, int k, int pair_idx,
        std::vector<int>& i_vec, std::vector<int>& j_vec,
        std::vector<double>& x_vec) {
    double denom = mat(k, j);
    if( denom != 0.0 && std::isfinite(denom) ){
        double ratio = mat(k, i) / denom;
        if( ratio != 0.0 && std::isfinite(ratio) ){
            i_vec.push_back(k);
            j_vec.push_back(pair_idx);
            x_vec.push_back(ratio);
        }
    }
}

// [[Rcpp::export(name = ".apply_transformation_difference_or_division")]]
S4 apply_transformation_difference_or_division(
        NumericMatrix mat, std::string method = "difference") {
    // samples = rows, features = cols
    int n_samples  = mat.nrow();
    int n_features = mat.ncol();
    int n_pairs    = n_features * (n_features - 1) / 2;
    CharacterVector original_feature_names = colnames(mat);

    // Prepare output vectors
    std::vector<int> i_vec, j_vec;
    std::vector<double> x_vec;
    CharacterVector colnames(n_pairs);

    int pair_idx = 0;
    for( int i = 0; i < n_features - 1; ++i ){
        for( int j = i + 1; j < n_features; ++j ){
            for( int k = 0; k < n_samples; ++k ){
                if( method == "difference" ){
                    calculate_difference(
                        mat, i, j, k, pair_idx, i_vec, j_vec, x_vec);
                } else if( method == "division" ){
                    calculate_division(
                        mat, i, j, k, pair_idx, i_vec, j_vec, x_vec);
                } else{
                    stop("Unknown method: use 'difference' or 'division'");
                }
            }
            // Generate colname for feature pair
            std::string name_i = as<std::string>(original_feature_names[i]);
            std::string name_j = as<std::string>(original_feature_names[j]);
            colnames[pair_idx] = name_i + "-" + name_j;
            ++pair_idx;
        }
    }

    // Create dgTMatrix (triplet sparse matrix)
    S4 tmat("dgTMatrix");
    // Samples as rows
    tmat.slot("i") = IntegerVector(i_vec.begin(), i_vec.end());
    // Pairs as columns
    tmat.slot("j") = IntegerVector(j_vec.begin(), j_vec.end());
    tmat.slot("x") = NumericVector(x_vec.begin(), x_vec.end());
    tmat.slot("Dim") = IntegerVector::create(n_samples, n_pairs);
    tmat.slot("Dimnames") = List::create(rownames(mat), colnames);

    // Convert dgTMatrix to dgCMatrix
    Function as("as");
    S4 cmat = as(tmat, "CsparseMatrix");

    return cmat;
}
