# Statistical Analysis Module

### 1. Distribution Comparison Tests

To examine if the depths of synapses differ among neurons with different
developmental origins or serving different functions, we used three
complementary methods to robustly compare meaningful differences with the
enormous number of synapses.

- Quantile Wasserstein Distance with Bootstrap Confidence Intervals
  - Measures how different two synaptic distributions are overall
- Subsampled Kolmogorov-Smirnov Test
  - Takes random samples to prevent over-sensitivity from huge sample sizes
  - Resample to assess instability from sampling
- Deep-superficial innervation test
  - Measures if each group (developmental or functional) prefers the
  superficial or deep part of the neuropil
  - Using reference cell types to mark the superficial or deep and test
  if the examined preference is more prominent than the reference types
  with bootstrapping
    - ME: Dm types (superficial) vs Pm types (deep)
    - LOP: T4/T5a,b (superficial) vs T4/T5c,d (deep) 
    - LO: Tm1,2,4 (superficial) vs Tm20,5a-f (deep)

### 2. Functional Enrichment Analysis
- Tests if birth time and function are linked with Fisher's exact test
- Corrects for multiple comparisons using FDR

## Configuration Parameters

Several parameters can be tuned in the analysis, and the defaults used in
this analysis are:

```yaml
# Bootstrap settings (Wasserstein distance)
n_bootstrap: 1000          # More iterations = more accurate but slower
conf_int: 95              # Confidence level (%)

# Quantile settings
n_quantiles: 1000         # Resolution of distribution comparison (updated default)

# K-S test settings
ks_subsample_size: 1000   # Points sampled per group
ks_n_iterations: 1000     # Number of subsampling iterations
ks_correction_method: "fdr" # Multiple testing correction method (fdr/bonferroni)

# Broad depth analysis settings
broad_depth_coefficient: 0.5    # Threshold coefficient for variable threshold
broad_depth_n_bootstrap: 1000    # Bootstrap iterations for broad depth analysis
broad_depth_conf_int: 95         # Confidence level for broad depth analysis

# Data filtering
sparse_limit: 100         # Minimum neurons per group for testing
```

## Technical Requirements

- **R packages**: data.table, ggplot2, Rcpp, openxlsx, dplyr, patchwork
- **C++ compiler**: For optimized statistical computations
- **Memory**: ~16GB RAM for full analysis with bootstrap iterations
- **Compute time**: 1-4 hours depending on data size and bootstrap iterations
