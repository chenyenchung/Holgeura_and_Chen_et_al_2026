#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
using namespace Rcpp;

// Fast median computation using nth_element (O(n) complexity)
inline double fast_median(NumericVector& data, int start, int end) {
  int n = end - start;
  if (n == 0) return NA_REAL;
  
  int mid = n / 2;
  
  if (n % 2 == 1) {
    // Odd number of elements
    std::nth_element(data.begin() + start, data.begin() + start + mid, data.begin() + end);
    return data[start + mid];
  } else {
    // Even number of elements
    std::nth_element(data.begin() + start, data.begin() + start + mid, data.begin() + end);
    double val1 = data[start + mid];
    
    // Find the previous element (mid-1) - it's the max of the lower half
    std::nth_element(data.begin() + start, data.begin() + start + mid - 1, data.begin() + start + mid);
    double val2 = data[start + mid - 1];
    
    return (val1 + val2) / 2.0;
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
  
  // Pool all groups for bootstrap under null hypothesis
  std::vector<double> pooled_all;
  pooled_all.reserve(n_sup + n_deep + n_interest);
  pooled_all.insert(pooled_all.end(), sup_data.begin(), sup_data.end());
  pooled_all.insert(pooled_all.end(), deep_data.begin(), deep_data.end());
  pooled_all.insert(pooled_all.end(), interest_data.begin(), interest_data.end());
  
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
  
  // Set up distribution for bootstrap sampling from pooled data
  std::uniform_int_distribution<int> dist_pooled(0, static_cast<int>(pooled_all.size()) - 1);
  
  // Counter for bootstrap test statistics >= observed
  int count_extreme = 0;
  
  // Perform bootstrap iterations under null hypothesis
  for (int i = 0; i < n_bootstrap; ++i) {
    // Bootstrap under null: resample from pooled data to original group sizes
    std::shuffle(pooled_all.begin(), pooled_all.end(), gen);
    
    // Split shuffled data back into three groups of original sizes
    std::vector<double> boot_sup_vec(pooled_all.begin(), pooled_all.begin() + n_sup);
    std::vector<double> boot_deep_vec(pooled_all.begin() + n_sup, pooled_all.begin() + n_sup + n_deep);
    std::vector<double> boot_interest_vec(pooled_all.begin() + n_sup + n_deep, pooled_all.end());
    
    // Copy to NumericVector for median calculation
    std::copy(boot_sup_vec.begin(), boot_sup_vec.end(), work_sup.begin());
    std::copy(boot_deep_vec.begin(), boot_deep_vec.end(), work_deep.begin());
    std::copy(boot_interest_vec.begin(), boot_interest_vec.end(), work_interest.begin());
    
    // Calculate bootstrap medians
    double boot_median_sup = fast_median(work_sup, 0, n_sup);
    double boot_median_deep = fast_median(work_deep, 0, n_deep);
    double boot_median_interest = fast_median(work_interest, 0, n_interest);
    
    // Calculate bootstrap distances
    double boot_sup_d = std::abs(boot_median_interest - boot_median_sup);
    double boot_deep_d = std::abs(boot_median_interest - boot_median_deep);
    double boot_distance_diff = boot_deep_d - boot_sup_d;
    
    // Calculate threshold
    double boot_delta_thres_base_val = std::abs(boot_median_deep - boot_median_sup);
    double boot_delta_thres;
    
    if (threshold_type == "variable") {
      boot_delta_thres = boot_delta_thres_base_val * coefficient;
    } else {
      boot_delta_thres = observed_delta_thres;
      boot_delta_thres_base_val = observed_delta_thres_base;
    }
    
    // Store bootstrap values
    bootstrap_sup_d[i] = boot_sup_d;
    bootstrap_deep_d[i] = boot_deep_d;
    bootstrap_delta_thres_base[i] = boot_delta_thres_base_val;
    bootstrap_distance_diff[i] = boot_distance_diff;
    
    // Calculate bootstrap test statistic and compare to observed
    double boot_test_stat = std::abs(boot_distance_diff) - boot_delta_thres;
    if (boot_test_stat >= observed_test_stat) {
      count_extreme++;
    }
  }
  
  // Calculate bootstrap p-value
  double p_value_exceeds_threshold = static_cast<double>(count_extreme) / n_bootstrap;
  
  // Calculate confidence intervals
  // Clone vectors for sorting (to preserve original bootstrap results)
  NumericVector sorted_distance_diff = clone(bootstrap_distance_diff);
  NumericVector sorted_delta_thres = clone(bootstrap_delta_thres_base);
  
  std::sort(sorted_distance_diff.begin(), sorted_distance_diff.end());
  std::sort(sorted_delta_thres.begin(), sorted_delta_thres.end());
  
  double alpha = (100.0 - conf_int) / 100.0;
  int lower_idx = static_cast<int>(n_bootstrap * alpha / 2.0);
  int upper_idx = static_cast<int>(n_bootstrap * (1.0 - alpha / 2.0));
  int median_idx = static_cast<int>(n_bootstrap * 0.5);
  
  // Ensure indices are within bounds
  lower_idx = std::max(0, std::min(lower_idx, n_bootstrap - 1));
  upper_idx = std::max(0, std::min(upper_idx, n_bootstrap - 1));
  median_idx = std::max(0, std::min(median_idx, n_bootstrap - 1));
  
  return List::create(
    Named("observed_sup_d") = observed_sup_d,
    Named("observed_deep_d") = observed_deep_d,
    Named("observed_distance_diff") = observed_distance_diff,
    Named("observed_delta_thres_base") = observed_delta_thres_base,
    Named("observed_delta_thres") = observed_delta_thres,
    Named("observed_test_statistic") = observed_test_stat,
    Named("bootstrap_distance_diff_median") = sorted_distance_diff[median_idx],
    Named("bootstrap_distance_diff_lower") = sorted_distance_diff[lower_idx],
    Named("bootstrap_distance_diff_upper") = sorted_distance_diff[upper_idx],
    Named("bootstrap_delta_thres_base_median") = sorted_delta_thres[median_idx],
    Named("bootstrap_delta_thres_base_lower") = sorted_delta_thres[lower_idx],
    Named("bootstrap_delta_thres_base_upper") = sorted_delta_thres[upper_idx],
    Named("p_value_exceeds_threshold") = p_value_exceeds_threshold,
    Named("bootstrap_sup_d") = bootstrap_sup_d,
    Named("bootstrap_deep_d") = bootstrap_deep_d,
    Named("bootstrap_delta_thres_base") = bootstrap_delta_thres_base,
    Named("bootstrap_distance_diff") = bootstrap_distance_diff
  );
}