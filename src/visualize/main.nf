#!/usr/bin/env nextflow
nextflow.preview.output = true

params.presetf = 'data/viz_preset.csv'
params.metaf = 'data/viz_meta.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.tsf = 'data/P15_tf.csv'
params.camf = 'data/P15_CAM.csv'
params.tslutf = 'data/to_selector.csv'

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
  path utils
  val subsample
  val slimit

  output:
  tuple val("${np}"), val("${preset}"), path('*.pdf'), optional: true

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
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

process VisualizeSelector {
  cpus '1'
  memory '8GB'
  time '30m'
  module 'r/gcc/4.4.0'

  input:
  tuple val(np), path(syn), val(stype), val(den), path(ts), path(tslut), val(gselection)
  path ann
  path meta
  path utils
  val subsample
  val slimit

  output:
  tuple val("${np}"), val("${gselection}"), path('*.pdf'), optional: true

  script:
  """
  v_selector.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${ts} \
    --tslut ${tslut} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

workflow {
  main:
  def NP = ['ME_L', 'ME_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R']
  def STYPE = ['pre', 'post']
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
    .combine(
      channel.fromPath(file(params.presetf))
        .splitCsv(header:true)
        .map { row -> row.preset }
    )
    .combine(channel.fromList(DEN_ALGO))
  utils_file = file('src/utils.r')
  out_ch = Visualize(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.presetf),
    utils_file,
    SUBSAMPLE_TO,
    SPARSE_LIMIT,
  )
  

  scond_ch = channel
    .fromList(NP)
    .map { it ->
      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
      return [it, file(synp)]
    }
    .combine(channel.fromList(STYPE))
    .combine(channel.fromList(DEN_ALGO))
    
  def TS_BIN = [file(params.tsf), file(params.camf)]
  def TS_LUT = [file(params.tslutf), file(params.tslutf)]
  def TS_LABEL = ['P15_TF', 'P15_CAM']
  g_ch = channel.fromList(TS_BIN)
    .merge(channel.fromList(TS_LUT))
    .merge(channel.fromList(TS_LABEL))
  sel_ch = VisualizeSelector(
    scond_ch.combine(g_ch),
    file(params.annf),
    file(params.metaf),
    utils_file,
    SUBSAMPLE_TO,
    SPARSE_LIMIT
  )

  publish:
  annov = out_ch
  selv = sel_ch
}

output {
  annov {
    path { input -> 
      def subset_dict = [
        'known': 'Previously_known',
        'new': 'All_confidently_annotated_highlighting_new',
        'all': 'All_confidently_annotated',
        'putative': 'Putative_OPC_types'
      ]
      def color_dict = [
        'temporal': 'Color_by_temporal_origin',
        'subsystem': 'Color_by_function_subsystem',
        'type': 'Color_by_cell_type'
      ]
      def side = input[0].split('_')[1]
      def ctype = input[1].split('_')[0]
      def subset = input[1].split('_')[1]
      def clabel = color_dict[ctype]
      def slabel = subset_dict[subset]
      return "${clabel}/${slabel}_${side}"
    }
  }
  selv {
    path { input ->
      def side = input[0].split('_')[1]
      def gsel = input[1]
      return "${gsel}/${side}"
    }
  }
}
