#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
using namespace Rcpp;

// Fast quantile computation for sorted vectors
double quantile_sorted(const std::vector<double>& sorted_data, double p) {
  if (sorted_data.empty()) return 0.0;
  
  double index = p * (sorted_data.size() - 1);
  int lower = static_cast<int>(std::floor(index));
  int upper = static_cast<int>(std::ceil(index));
  
  if (lower == upper) {
    return sorted_data[lower];
  }
  
  double weight = index - lower;
  return sorted_data[lower] * (1 - weight) + sorted_data[upper] * weight;
}

// Compute quantiles for a vector
std::vector<double> compute_quantiles(std::vector<double> data, const std::vector<double>& probs) {
  std::sort(data.begin(), data.end());
  std::vector<double> quantiles;
  quantiles.reserve(probs.size());
  
  for (double p : probs) {
    quantiles.push_back(quantile_sorted(data, p));
  }
  
  return quantiles;
}

// Compute quantile Wasserstein distance between two samples
double compute_wasserstein(std::vector<double> sample1, 
                          std::vector<double> sample2, 
                          const std::vector<double>& probs) {
  auto q1 = compute_quantiles(sample1, probs);
  auto q2 = compute_quantiles(sample2, probs);
  
  double sum_abs_diff = 0.0;
  for (size_t i = 0; i < q1.size(); ++i) {
    double diff = q1[i] - q2[i];
    sum_abs_diff += std::abs(diff);
  }
  
  return sum_abs_diff / q1.size();
}

//' Fast bootstrap for Quantile Wasserstein test
//' 
//' @param sample1 First sample
//' @param sample2 Second sample  
//' @param n_quantiles Number of quantiles (default 100)
//' @param n_bootstrap Number of bootstrap iterations (default 1000)
//' @param conf_int Confidence interval percentage (default 95)
//' @return List with observed statistic and bootstrap statistics
// [[Rcpp::export]]
List quantile_wasserstein_bootstrap(NumericVector sample1, 
                                   NumericVector sample2,
                                   int n_quantiles = 100,
                                   int n_bootstrap = 1000,
                                   double conf_int = 95.0) {
  
  // Input validation
  if (sample1.size() < 10 || sample2.size() < 10) {
    stop("Samples must have at least 10 observations each");
  }
  
  if (conf_int <= 0 || conf_int >= 100) {
    stop("Confidence interval must be between 0 and 100");
  }
  
  // Convert to std::vector
  std::vector<double> s1(sample1.begin(), sample1.end());
  std::vector<double> s2(sample2.begin(), sample2.end());
  
  // Generate quantile probabilities (excluding 0 and 1)
  std::vector<double> probs;
  probs.reserve(n_quantiles);
  for (int i = 1; i <= n_quantiles; ++i) {
    probs.push_back(static_cast<double>(i) / (n_quantiles + 1));
  }
  
  // Compute observed test statistic
  double observed_stat = compute_wasserstein(s1, s2, probs);
  
  // Pool samples for bootstrap
  std::vector<double> pooled;
  pooled.reserve(s1.size() + s2.size());
  pooled.insert(pooled.end(), s1.begin(), s1.end());
  pooled.insert(pooled.end(), s2.begin(), s2.end());
  
  size_t n1 = s1.size();
  size_t n2 = s2.size();
  
  // Bootstrap null distribution
  std::vector<double> null_stats;
  null_stats.reserve(n_bootstrap);
  
  // Set up random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  
  for (int b = 0; b < n_bootstrap; ++b) {
    // Shuffle pooled samples
    std::vector<double> shuffled = pooled;
    std::shuffle(shuffled.begin(), shuffled.end(), gen);
    
    // Split into two bootstrap samples
    std::vector<double> boot_s1(shuffled.begin(), shuffled.begin() + n1);
    std::vector<double> boot_s2(shuffled.begin() + n1, shuffled.end());
    
    // Compute bootstrap statistic
    double boot_stat = compute_wasserstein(boot_s1, boot_s2, probs);
    null_stats.push_back(boot_stat);
  }
  
  
  // Compute confidence intervals
  std::sort(null_stats.begin(), null_stats.end());
  
  double offset = (100.0 - conf_int) / 100.0 / 2.0;
  double lower_perc = offset;
  double upper_perc = 1.0 - offset;
  double median_perc = 0.5;
  
  double boot_median = quantile_sorted(null_stats, median_perc);
  double boot_lower = quantile_sorted(null_stats, lower_perc);
  double boot_upper = quantile_sorted(null_stats, upper_perc);
  
  return List::create(
    Named("observed_statistic") = observed_stat,
    Named("bootstrap_median") = boot_median,
    Named("bootstrap_lower") = boot_lower,
    Named("bootstrap_upper") = boot_upper,
    Named("null_distribution") = wrap(null_stats),
    Named("conf_level") = conf_int
  );
}

//' Fast pairwise bootstrap tests
//' 
//' @param data_list List of numeric vectors (samples)
//' @param names Character vector of sample names
//' @param n_quantiles Number of quantiles (default 100)
//' @param n_bootstrap Number of bootstrap iterations (default 1000)
//' @param conf_int Confidence interval percentage (default 95)
//' @return DataFrame with test results
// [[Rcpp::export]]
DataFrame pairwise_wasserstein_tests(List data_list,
                                    CharacterVector names,
                                    int n_quantiles = 100,
                                    int n_bootstrap = 1000,
                                    double conf_int = 95.0) {
  
  int n_samples = data_list.size();
  if (n_samples < 2) {
    stop("Need at least 2 samples for pairwise tests");
  }
  
  std::vector<std::string> group1_vec, group2_vec;
  std::vector<double> median_vec, lower_vec, upper_vec;
  
  // Perform all pairwise comparisons
  for (int i = 0; i < n_samples - 1; ++i) {
    for (int j = i + 1; j < n_samples; ++j) {
      NumericVector sample1 = as<NumericVector>(data_list[i]);
      NumericVector sample2 = as<NumericVector>(data_list[j]);
      
      List result = quantile_wasserstein_bootstrap(sample1, sample2, 
                                                  n_quantiles, n_bootstrap, conf_int);
      
      group1_vec.push_back(as<std::string>(names[i]));
      group2_vec.push_back(as<std::string>(names[j]));
      median_vec.push_back(as<double>(result["bootstrap_median"]));
      lower_vec.push_back(as<double>(result["bootstrap_lower"]));
      upper_vec.push_back(as<double>(result["bootstrap_upper"]));
    }
  }
  
  return DataFrame::create(
    Named("group1") = wrap(group1_vec),
    Named("group2") = wrap(group2_vec),
    Named("test_statistic_median") = wrap(median_vec),
    Named("test_statistic_lower") = wrap(lower_vec),
    Named("test_statistic_upper") = wrap(upper_vec)
  );
}