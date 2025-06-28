#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
using namespace Rcpp;

// Fast median computation using partial_sort_copy (preserves input data)
inline double fast_median(const NumericVector& data, int start, int end) {
  int n = end - start;
  if (n == 0) return NA_REAL;
  
  int mid = n / 2;
  
  if (n % 2 == 1) {
    // Odd number: need (mid+1) smallest elements, return the middle one
    std::vector<double> temp(mid + 1);
    std::partial_sort_copy(data.begin() + start, data.begin() + end,
                          temp.begin(), temp.end());
    return temp[mid];
  } else {
    // Even number: need (mid+1) smallest elements, return average of last two
    std::vector<double> temp(mid + 1);
    std::partial_sort_copy(data.begin() + start, data.begin() + end,
                          temp.begin(), temp.end());
    return (temp[mid-1] + temp[mid]) / 2.0;
  }
}

// Efficient bootstrap sampling with pre-allocated memory
inline void bootstrap_sample_inplace(const NumericVector& source, 
                                   NumericVector& target,
                                   int n,
                                   std::uniform_int_distribution<int>& dist,
                                   std::mt19937& gen) {
  for (int i = 0; i < n; ++i) {
    target[i] = source[dist(gen)];
  }
}

//' Optimized bootstrap analysis for broad depth testing
//' 
//' @param ref_sup_depths Depth values for superficial reference group
//' @param ref_deep_depths Depth values for deep reference group
//' @param interest_depths Depth values for types of interest
//' @param coefficient Threshold coefficient (default 0.25)
//' @param n_bootstrap Number of bootstrap iterations (default 1000)
//' @param conf_int Confidence interval percentage (default 95)
//' @param threshold_type Type of threshold ("variable" or "fixed")
//' @param seed Random seed for reproducibility (optional)
//' @return List with bootstrap results
// [[Rcpp::export]]
List perform_broad_depth_bootstrap(NumericVector ref_sup_depths,
                                 NumericVector ref_deep_depths,
                                 NumericVector interest_depths,
                                 double coefficient = 0.25,
                                 int n_bootstrap = 1000,
                                 double conf_int = 95.0,
                                 std::string threshold_type = "variable",
                                 Rcpp::Nullable<int> seed = R_NilValue) {
  
  // Input validation
  if (ref_sup_depths.size() == 0 || ref_deep_depths.size() == 0 || interest_depths.size() == 0) {
    stop("Empty input vectors provided");
  }
  
  if (conf_int <= 0 || conf_int >= 100) {
    stop("Confidence interval must be between 0 and 100");
  }
  
  // Remove NA values using Rcpp sugar (creates new vectors)
  NumericVector sup_data = na_omit(ref_sup_depths);
  NumericVector deep_data = na_omit(ref_deep_depths);
  NumericVector interest_data = na_omit(interest_depths);
  
  int n_sup = sup_data.size();
  int n_deep = deep_data.size();
  int n_interest = interest_data.size();
  
  // Pre-allocate work vectors for bootstrap sampling
  NumericVector work_sup(n_sup);
  NumericVector work_deep(n_deep);
  NumericVector work_interest(n_interest);
  
  // Copy data to work vectors for initial median calculation
  std::copy(sup_data.begin(), sup_data.end(), work_sup.begin());
  std::copy(deep_data.begin(), deep_data.end(), work_deep.begin());
  std::copy(interest_data.begin(), interest_data.end(), work_interest.begin());
  
  // Calculate observed statistics
  double observed_median_sup = fast_median(work_sup, 0, n_sup);
  double observed_median_deep = fast_median(work_deep, 0, n_deep);
  double observed_median_interest = fast_median(work_interest, 0, n_interest);
  
  double observed_sup_d = std::abs(observed_median_interest - observed_median_sup);
  double observed_deep_d = std::abs(observed_median_interest - observed_median_deep);
  double observed_distance_diff = observed_deep_d - observed_sup_d;
  double observed_delta_thres_base = std::abs(observed_median_deep - observed_median_sup);
  double observed_delta_thres = observed_delta_thres_base * coefficient;
  
  // Pre-allocate bootstrap result vectors
  NumericVector bootstrap_sup_d(n_bootstrap);
  NumericVector bootstrap_deep_d(n_bootstrap);
  NumericVector bootstrap_delta_thres_base(n_bootstrap);
  NumericVector bootstrap_distance_diff(n_bootstrap);
  NumericVector bootstrap_test_stats(n_bootstrap);
  
  // Calculate observed test statistic
  double observed_test_stat = std::abs(observed_distance_diff) - observed_delta_thres;
  
  // Set up random number generator
  std::mt19937 gen;
  if (seed.isNotNull()) {
    gen.seed(as<int>(seed));
  } else {
    std::random_device rd;
    gen.seed(rd());
  }
  
  // Set up distributions for independent bootstrap sampling from each group
  std::uniform_int_distribution<int> dist_sup(0, n_sup - 1);
  std::uniform_int_distribution<int> dist_deep(0, n_deep - 1);
  std::uniform_int_distribution<int> dist_interest(0, n_interest - 1);
  
  // Counter for bootstrap test statistics > 0 (one-sided p-value)
  int count_positive = 0;
  
  // Perform bootstrap iterations with independent resampling
  for (int i = 0; i < n_bootstrap; ++i) {
    // Independently resample each group with replacement
    for (int j = 0; j < n_sup; ++j) {
      work_sup[j] = sup_data[dist_sup(gen)];
    }
    for (int j = 0; j < n_deep; ++j) {
      work_deep[j] = deep_data[dist_deep(gen)];
    }
    for (int j = 0; j < n_interest; ++j) {
      work_interest[j] = interest_data[dist_interest(gen)];
    }
    
    // Calculate bootstrap medians
    double boot_median_sup = fast_median(work_sup, 0, n_sup);
    double boot_median_deep = fast_median(work_deep, 0, n_deep);
    double boot_median_interest = fast_median(work_interest, 0, n_interest);
    
    // Calculate bootstrap distances
    double boot_sup_d = std::abs(boot_median_interest - boot_median_sup);
    double boot_deep_d = std::abs(boot_median_interest - boot_median_deep);
    double boot_distance_diff = boot_deep_d - boot_sup_d;
    
    // Calculate threshold (always use variable threshold for independent bootstrap)
    double boot_delta_thres_base_val = std::abs(boot_median_deep - boot_median_sup);
    double boot_delta_thres = boot_delta_thres_base_val * coefficient;
    
    // Store bootstrap values
    bootstrap_sup_d[i] = boot_sup_d;
    bootstrap_deep_d[i] = boot_deep_d;
    bootstrap_delta_thres_base[i] = boot_delta_thres_base_val;
    bootstrap_distance_diff[i] = boot_distance_diff;
    
    // Calculate bootstrap test statistic: |distance_diff| - coef * |ref_diff|
    double boot_test_stat = std::abs(boot_distance_diff) - boot_delta_thres;
    bootstrap_test_stats[i] = boot_test_stat;
    
    // Count positive test statistics for one-sided p-value
    if (boot_test_stat > 0) {
      count_positive++;
    }
  }
  
  // Calculate one-sided p-value for test statistic > 0
  double p_value_exceeds_threshold = 1.0 - static_cast<double>(count_positive) / n_bootstrap;
  
  // Calculate confidence intervals
  // Clone vectors for sorting (to preserve original bootstrap results)
  NumericVector sorted_distance_diff = clone(bootstrap_distance_diff);
  NumericVector sorted_delta_thres = clone(bootstrap_delta_thres_base);
  NumericVector sorted_test_stats = clone(bootstrap_test_stats);
  
  std::sort(sorted_distance_diff.begin(), sorted_distance_diff.end());
  std::sort(sorted_delta_thres.begin(), sorted_delta_thres.end());
  std::sort(sorted_test_stats.begin(), sorted_test_stats.end());
  
  double alpha = (100.0 - conf_int) / 100.0;
  int lower_idx = static_cast<int>(n_bootstrap * alpha / 2.0);
  int upper_idx = static_cast<int>(n_bootstrap * (1.0 - alpha / 2.0));
  int median_idx = static_cast<int>(n_bootstrap * 0.5);
  
  // Ensure indices are within bounds
  lower_idx = std::max(0, std::min(lower_idx, n_bootstrap - 1));
  upper_idx = std::max(0, std::min(upper_idx, n_bootstrap - 1));
  median_idx = std::max(0, std::min(median_idx, n_bootstrap - 1));
  
  List result;
  result["observed_sup_d"] = observed_sup_d;
  result["observed_deep_d"] = observed_deep_d;
  result["observed_distance_diff"] = observed_distance_diff;
  result["observed_delta_thres_base"] = observed_delta_thres_base;
  result["observed_delta_thres"] = observed_delta_thres;
  result["observed_test_statistic"] = observed_test_stat;
  result["bootstrap_distance_diff_median"] = sorted_distance_diff[median_idx];
  result["bootstrap_distance_diff_lower"] = sorted_distance_diff[lower_idx];
  result["bootstrap_distance_diff_upper"] = sorted_distance_diff[upper_idx];
  result["bootstrap_delta_thres_base_median"] = sorted_delta_thres[median_idx];
  result["bootstrap_delta_thres_base_lower"] = sorted_delta_thres[lower_idx];
  result["bootstrap_delta_thres_base_upper"] = sorted_delta_thres[upper_idx];
  result["bootstrap_test_stat_median"] = sorted_test_stats[median_idx];
  result["bootstrap_test_stat_lower"] = sorted_test_stats[lower_idx];
  result["bootstrap_test_stat_upper"] = sorted_test_stats[upper_idx];
  result["p_value_exceeds_threshold"] = p_value_exceeds_threshold;
  result["bootstrap_sup_d"] = bootstrap_sup_d;
  result["bootstrap_deep_d"] = bootstrap_deep_d;
  result["bootstrap_delta_thres_base"] = bootstrap_delta_thres_base;
  result["bootstrap_distance_diff"] = bootstrap_distance_diff;
  result["bootstrap_test_stats"] = bootstrap_test_stats;
  
  return result;
}