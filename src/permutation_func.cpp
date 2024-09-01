#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// function for adjusted p value from BH method
arma::vec adjustPValues(const arma::vec& p_values) {
  Environment stats("package:stats");  // Access the 'stats' package
  Function p_adjust = stats["p.adjust"];  // Get the p.adjust function

  // Call the p.adjust function from R and return the result
  NumericVector res = p_adjust(wrap(p_values), "BH");

  // Convert the Rcpp::NumericVector result to an arma::vec
  arma::vec adjusted_p_values = Rcpp::as<arma::vec>(res);

  return adjusted_p_values;
}


// function for calculating geometric mean
arma::vec geometricMean(const arma::mat& mat) {
  int nRows = mat.n_rows;
  int nCols = mat.n_cols;
  arma::vec result(nRows);

  for (int i = 0; i < nRows; i++) {
    double product = 1;
    for (int j = 0; j < nCols; j++) {
      product *= mat(i, j);  // Multiply elements in the row
    }
    result[i] = pow(product, 1.0 / nCols);  // nth root of the product
  }

  return result;
}

// function for calculating geomatric mean
arma::vec callNnls(const arma::mat& x, const arma::mat& y) {
  Environment stats("package:nnls");
  Function nnls = stats["nnls"];
  NumericMatrix x_rcpp = wrap(x);
  NumericMatrix y_rcpp = wrap(y);
  List result = nnls(Named("A") = x_rcpp, Named("b") = y_rcpp);
  NumericVector coefficients_rcpp = result["x"];
  arma::vec coefficients = as<arma::vec>(coefficients_rcpp);
  return coefficients;
}


// pre-function for local and regional permutation
List Local_Regional_permuteIndices(const arma::mat& permutationMatrix,
                                   const arma::mat& permut_col,
                                   const arma::mat& cellPropMatrix,
                                   const arma::mat& spotGeneMatrix,
                                   const arma::vec& LigandVectorIndex,
                                   const arma::vec& ReceptorVectorIndex,
                                   const arma::mat& null_expression,
                                   unsigned int nE) {

  arma::vec permutedIndices = permutationMatrix.col(nE - 1);

  arma::vec col_indices = permut_col.col(nE - 1);

  arma::uvec uPermutedIndices = conv_to<uvec>::from(permutedIndices - 1);
  arma::uvec ucol_indices = conv_to<uvec>::from(col_indices - 1);
  //std::cout << "ucol_indices:\n" << ucol_indices << std::endl;

  arma::mat selectpre = cellPropMatrix.rows(uPermutedIndices);
  arma::mat selectCellProp = selectpre.cols(ucol_indices);

  arma::mat avgLigand(selectCellProp.n_cols, LigandVectorIndex.size());
  arma::mat avgReceptor(selectCellProp.n_cols, ReceptorVectorIndex.size());

  arma::rowvec col_sums = sum(selectCellProp, 0);
  arma::vec col_sums_transposed = col_sums.t();
  double n = selectCellProp.n_rows;

  for (size_t l = 0; l < LigandVectorIndex.size(); ++l) {
    // Extract expression levels for the specified gene across selected spots
    arma::vec exp_L_vector = spotGeneMatrix.row( LigandVectorIndex[l] - 1 ).t();
    exp_L_vector = exp_L_vector.elem(uPermutedIndices);
    arma::mat exp_L_matrix = exp_L_vector.as_col();
    arma::vec exp = callNnls(selectCellProp, exp_L_matrix );


    double sum_y1 = arma::sum(exp_L_vector);
    double mean_y1 = sum_y1 / exp_L_vector.n_elem;

    //double threshold1 = mean_y1 / exp.n_elem;
    arma::vec threshold1 = mean_y1*(col_sums_transposed / n);
    arma::uvec indices = arma::find(exp > 0); // Find indices of elements > 0
    double count_nonzero = indices.n_elem; // Number of non-zero elements
    double threshold2 = count_nonzero > 0 ? 1.0 / count_nonzero : 0; // Check to prevent division by zero
    arma::mat swept = selectCellProp.each_row() % exp.t();
    arma::vec col_means = arma::mean(swept, 0).t();
    arma::uvec condition_indices = arma::find(col_means < threshold1 || (exp / arma::sum(exp)) <  threshold2);
    arma::mat swept1 = selectCellProp.each_row() % exp.t();
    arma::vec col_means1 = arma::mean(swept1, 0).t();

    avgLigand.col(l) = exp;




  }

  for (size_t r = 0; r < ReceptorVectorIndex.size(); ++r) {
    arma::vec exp_R_vector = spotGeneMatrix.row(ReceptorVectorIndex[r] - 1).t();
    exp_R_vector = exp_R_vector.elem(uPermutedIndices);
    arma::mat exp_R_matrix = exp_R_vector.as_col();
    arma::vec exp = callNnls(selectCellProp, exp_R_matrix );

    double sum_y1 = arma::sum(exp_R_vector);
    double mean_y1 = sum_y1 / exp_R_vector.n_elem;

    //double threshold1 = mean_y1 / exp.n_elem;
    arma::vec threshold1 = mean_y1*(col_sums_transposed / n);
    arma::uvec indices = arma::find(exp > 0); // Find indices of elements > 0
    double count_nonzero = indices.n_elem; // Number of non-zero elements
    double threshold2 = count_nonzero > 0 ? 1.0 / count_nonzero : 0; // Check to prevent division by zero
    arma::mat swept = selectCellProp.each_row() % exp.t();
    arma::vec col_means = arma::mean(swept, 0).t();
    arma::uvec condition_indices = arma::find(col_means < threshold1 || (exp / arma::sum(exp)) <  threshold2);
    arma::mat swept2 = selectCellProp.each_row() % exp.t();
    arma::vec col_means2 = arma::mean(swept2, 0).t();



    avgReceptor.col(r) = exp;

  }

  //std::cout << "avgLigand:\n" << avgLigand << std::endl;
  //std::cout << "avgReceptor:\n" << avgReceptor << std::endl;

  arma::vec finalAvgLigand = geometricMean(avgLigand);
  arma::vec finalAvgReceptor = geometricMean(avgReceptor);


  //std::cout << "Column Sums:\n" << col_sums_transposed << std::endl;


  //std::cout << "n:\n" << n << std::endl;
  arma::vec norm_col_sums = col_sums_transposed / n;

  //std::cout << "Normalized Column Sums:\n" << norm_col_sums << std::endl;


  arma::mat outer_product = norm_col_sums * norm_col_sums.t(); // .t() for transpose

  //std::cout << "Outer Product:\n" << outer_product << std::endl;


  arma::mat interact = finalAvgLigand * finalAvgReceptor.t();

  //std::cout << "interact:\n" << interact << std::endl;

  arma::mat transformed_interact = interact / (0.5 + interact);

  //std::cout << "Transformed Interact Matrix:\n" << transformed_interact << std::endl;

  arma::mat per_result = transformed_interact % outer_product ;


  //std::cout << "Element-wise Product of Transformed and Outer Product:\n" << per_result << std::endl;

  arma::mat null_expression1 = null_expression;

  //std::cout << "Transformed Interact Matrix:\n" << null_expression1 << std::endl;

  // Compute the outer product
  arma::umat result = (per_result > null_expression1 );



  // Create a list to return the results
  List results;
  results["interact"] = result;
  return  results;
}



// function for local and regional permutation
// [[Rcpp::export]]
arma::mat Local_Regional_Permutations(const arma::mat& permutationMatrix,
                                       const arma::mat& permut_col,
                                       const arma::mat& cellPropMatrix,
                                       const arma::mat& spotGeneMatrix,
                                       const arma::vec& LigandVectorIndex,
                                       const arma::vec& ReceptorVectorIndex,
                                       const arma::mat& null_expression,
                                       int nBoot) {
  // Initialize a matrix to hold the sum of all interactions
  arma::mat sumInteraction;
  bool firstIteration = true;

  for (int nE = 1; nE <= nBoot; ++nE) {
    List result = Local_Regional_permuteIndices(permutationMatrix,permut_col,  cellPropMatrix, spotGeneMatrix,
                                                LigandVectorIndex, ReceptorVectorIndex,
                                                null_expression, nE);
    arma::mat interact = as<mat>(result["interact"]);

    if (firstIteration) {
      sumInteraction = interact;
      firstIteration = false;
    } else {
      sumInteraction += interact; // Sum the interaction matrices
    }
  }

  // Divide the sumInteraction matrix by the number of bootstraps (nBoot) to get the average
  arma::mat averageInteraction = sumInteraction / nBoot;
  return averageInteraction;
}



















// pre-function for global permutation
List Global_permuteIndices(const arma::mat& permutationMatrix,
                                   const arma::mat& permut_null_regionMatrix,
                                   const arma::mat& permut_col,
                                   const arma::mat& cellPropMatrix,
                                   const arma::mat& spotGeneMatrix,
                                   const arma::vec& LigandVectorIndex,
                                   const arma::vec& ReceptorVectorIndex,
                                   const arma::mat& null_expression,
                                   unsigned int nE) {

  arma::vec permutedIndices = permutationMatrix.col(nE - 1);

  arma::vec col_indices = permut_col.col(nE - 1);

  arma::vec permut_null_regionIndices = permut_null_regionMatrix.col(0);
  arma::uvec uPermutedIndices = conv_to<uvec>::from(permutedIndices - 1);
  arma::uvec upermut_null_regionIndices = conv_to<uvec>::from(permut_null_regionIndices - 1);
  arma::uvec ucol_indices = conv_to<uvec>::from(col_indices - 1);
  //std::cout << "ucol_indices:\n" << ucol_indices << std::endl;

  arma::mat selectCellProp = cellPropMatrix.rows(uPermutedIndices);

  arma::mat avgLigand(selectCellProp.n_cols, LigandVectorIndex.size());
  arma::mat avgReceptor(selectCellProp.n_cols, ReceptorVectorIndex.size());

  arma::rowvec col_sums = sum(selectCellProp, 0);
  arma::vec col_sums_transposed = col_sums.t();
  double n = selectCellProp.n_rows;

  for (size_t l = 0; l < LigandVectorIndex.size(); ++l) {
    // Extract expression levels for the specified gene across selected spots
    arma::vec exp_L_vector = spotGeneMatrix.row( LigandVectorIndex[l] - 1 ).t();
    exp_L_vector = exp_L_vector.elem(upermut_null_regionIndices);
    arma::mat exp_L_matrix = exp_L_vector.as_col();
    arma::vec exp = callNnls(selectCellProp, exp_L_matrix );
    avgLigand.col(l) = exp;
  }

  for (size_t r = 0; r < ReceptorVectorIndex.size(); ++r) {
    arma::vec exp_R_vector = spotGeneMatrix.row(ReceptorVectorIndex[r] - 1).t();
    exp_R_vector = exp_R_vector.elem(upermut_null_regionIndices);
    arma::mat exp_R_matrix = exp_R_vector.as_col();
    arma::vec exp = callNnls(selectCellProp, exp_R_matrix );

    avgReceptor.col(r) = exp;

  }


  arma::vec finalAvgLigand = geometricMean(avgLigand);
  arma::vec finalAvgReceptor = geometricMean(avgReceptor);


  arma::vec norm_col_sums = col_sums_transposed / n;
  //std::cout << "Normalized Column Sums:\n" << norm_col_sums << std::endl;


  arma::mat outer_product = norm_col_sums * norm_col_sums.t(); // .t() for transpose

  //std::cout << "Outer Product:\n" << outer_product << std::endl;


  arma::mat interact = finalAvgLigand * finalAvgReceptor.t();

  //std::cout << "interact:\n" << interact << std::endl;

  arma::mat transformed_interact = interact / (0.5 + interact);

  //std::cout << "Transformed Interact Matrix:\n" << transformed_interact << std::endl;

  arma::mat per_result = transformed_interact % outer_product ;

  //std::cout << "Element-wise Product of Transformed and Outer Product:\n" << per_result << std::endl;

  arma::mat null_expression1 = null_expression;


  //std::cout << "Transformed Interact Matrix:\n" << null_expression1 << std::endl;

  // Compute the outer product
  arma::umat result = (per_result > null_expression1 );


  // Create a list to return the results
  List results;
  results["interact"] = result;
  return  results;
}

// function for global permutation
// [[Rcpp::export]]
arma::mat Global_Permutations(const arma::mat& permutationMatrix,
                                       const arma::mat& permut_null_regionMatrix,
                                       const arma::mat& permut_col,
                                       const arma::mat& cellPropMatrix,
                                       const arma::mat& spotGeneMatrix,
                                       const arma::vec& LigandVectorIndex,
                                       const arma::vec& ReceptorVectorIndex,
                                       const arma::mat& null_expression,
                                       int nBoot) {
  // Initialize a matrix to hold the sum of all interactions
  arma::mat sumInteraction;
  bool firstIteration = true;

  for (int nE = 1; nE <= nBoot; ++nE) {
    List result = Global_permuteIndices(permutationMatrix,permut_null_regionMatrix,permut_col,  cellPropMatrix, spotGeneMatrix,
                                                LigandVectorIndex, ReceptorVectorIndex,
                                                null_expression, nE);
    arma::mat interact = as<mat>(result["interact"]);

    if (firstIteration) {
      sumInteraction = interact;
      firstIteration = false;
    } else {
      sumInteraction += interact; // Sum the interaction matrices
    }
  }

  // Divide the sumInteraction matrix by the number of bootstraps (nBoot) to get the average
  arma::mat averageInteraction = sumInteraction / nBoot;
  return averageInteraction;
}













