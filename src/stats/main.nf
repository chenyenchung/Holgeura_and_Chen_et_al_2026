#!/usr/bin/env nextflow
nextflow.preview.output = true

// Default parameters
params.presetf = 'data/viz_preset.csv'
params.metaf = 'data/viz_meta.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.utilsf = 'src/utils.r'
params.depth_stats_cppf = 'src/stats/bin/depth_stats.cpp'
params.combine_scriptf = 'src/stats/bin/combine_results.r'
params.n_quantiles = 100
params.n_bootstrap = 1000
params.conf_int = 95
params.ks_subsample_size = 1000
params.ks_n_iterations = 1000
params.ks_correction_methods = 'bonferroni'
params.subsample = 10000
params.sparse_limit = 100

process DepthStatsAnalysis {
  cpus '1'
  memory '8GB'
  time '1h'
  module 'r/gcc/4.4.0'

  input:
  tuple val(np), path(syn), val(stype), val(preset)
  path ann
  path meta
  path presetf
  path utils
  path cppsrc
  val n_quantiles
  val n_bootstrap
  val conf_int
  val slimit
  val ks_subsample_size
  val ks_n_iterations
  val ks_correction_methods

  output:
  tuple val("${np}"), val("${preset}"), val("${stype}"), 
        path('*.csv')

  script:
  """
  depth_stats.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --use_preset ${preset} \
    --ann ${ann} \
    --meta ${meta} \
    --preset ${presetf} \
    --utils ${utils} \
    --cppsrc ${cppsrc} \
    --sparse_limit ${slimit} \
    --n_quantiles ${n_quantiles} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int} \
    --ks_subsample_size ${ks_subsample_size} \
    --ks_n_iterations ${ks_n_iterations} \
    --ks_correction_methods ${ks_correction_methods}
  """
}

process VisualizeStatSummary {
  cpus '1'
  memory '8GB'
  time '30m'
  module 'r/gcc/4.4.0'

  input:
  tuple val(np), val(preset), val(stype), path(results_csv)
  path meta
  path utils

  output:
  tuple val("${np}"), val("${preset}"), val("${stype}"), path('*.pdf'), optional: true

  script:
  def output_prefix = results_csv.baseName.replaceAll('_depth_stats', '')
  """
  visualize_stat_summary.r \
    --input_file ${results_csv} \
    --output_prefix ${output_prefix} \
    --utils ${utils}
  """
}

process CombineResults {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  path csvs
  path combine_script

  output:
  path 'combined_depth_stats_results.csv', emit: csv
  path 'statistical_summary.txt', emit: summary

  script:
  """
  combine_results.r
  """
}

workflow {
  main:
  // Define analysis parameters
  def NP = ['ME_L', 'ME_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R']
  def STYPE = ['pre', 'post']
  def MAT_PREFIX = 'int/idv_mat/'
  def SUBSAMPLE_TO = params.subsample
  def SPARSE_LIMIT = params.sparse_limit

  // Create input channel
  cond_ch = channel
    .fromList(NP)
    .map { it ->
      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
      return [it, file(synp)]
    }
    .combine(channel.fromList(STYPE))
    .combine(
      channel.fromPath(file(params.presetf))
        .splitCsv(header:true)
        .map { row -> row.preset }
        .take( 4 )
    )
    
    
  // Perform statistical analysis
  analysis_ch = DepthStatsAnalysis(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.presetf),
    file(params.utilsf),
    file(params.depth_stats_cppf),
    params.n_quantiles,
    params.n_bootstrap,
    params.conf_int,
    SPARSE_LIMIT,
    params.ks_subsample_size,
    params.ks_n_iterations,
    params.ks_correction_methods
  )

  // Generate visualizations
  viz_ch = VisualizeStatSummary(
    analysis_ch,
    file(params.metaf),
    file(params.utilsf)
  )

  // Combine all results
  combined_ch = CombineResults(
    analysis_ch.map { it -> it[3] }.collect(),
    file(params.combine_scriptf)
  )

  publish:
  stats_plots = viz_ch
  summary = combined_ch.summary
  csv = combined_ch.csv
}

output {
  stats_plots {
    path { input ->
      def np = input[0]
      def preset = input[1]
      def stype = input[2]
      return "stats/${preset}/${np}_${stype}"
    }
  }
  summary { path "stats/" }
  csv { path "stats/" }
}
