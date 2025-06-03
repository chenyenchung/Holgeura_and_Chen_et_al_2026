#!/usr/bin/env nextflow
nextflow.preview.output = true

params.presetf = 'data/viz_preset.csv'
params.metaf = 'data/viz_meta.csv'
params.annf = 'data/visual_neurons_20250602.csv'

process Visualize {
  cpus '1'
  memory '8GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  tuple val(np), path(syn), val(stype), val(preset), val(den)
  path ann
  path meta
  path presetf
  val subsample
  val slimit

  output:
  path '*.pdf'

  script:
  """
  visualize.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --use_preset ${preset} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --preset ${presetf} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

workflow {
  def NP = ['ME_L', 'ME_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R']
  def STYPE = ['pre', 'post']
  def PRESET = ['temporal', 'subsystem_ann', 'new_type']
  def DEN_ALGO = ['asis', 'pertype']
  def MAT_PREFIX = 'int/idv_mat/'
  def SUBSAMPLE_TO = 10000
  def SPARSE_LIMIT = 100
  cond_ch = channel
    .fromList(NP)
    .map { it ->
      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
      return [it, file(synp)]
    }
    .combine(channel.fromList(STYPE))
    .combine(channel.fromList(PRESET))
    .combine(channel.fromList(DEN_ALGO))
  out_ch = Visualize(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.presetf),
    SUBSAMPLE_TO,
    SPARSE_LIMIT,
  )

  publish:
  viz = out_ch
}

output {
  viz {
    path '20250603'
  }
}
