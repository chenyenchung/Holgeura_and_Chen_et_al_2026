#!/usr/env/bin nextflow
nextflow.preview.output = true

process AnnotateNeuropil {
  cpus '1'
  memory '64GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  path syn
  path ann
  path type
  path utils

  output:
  path "ol_passed_filter.csv.gz"

  script:
  """
  filter_synapses.R \
    --syn ${syn} \
    --ann ${ann} \
    --type ${type} \
    --utils ${utils}
  """
}

process SplitNeuropil {
  cpus '1'
  memory '8GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  path syn
  path utils

  output:
  path "*.csv.gz"

  script:
  """
  split_neuropil.R --syn ${syn} --utils ${utils}
  """
}

process RotateNeuropil {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  tuple path(syn), path(utils)

  output:
  path "*_rotated.csv.gz"

  script:
  """
  rotate.R --syn ${syn} --utils ${utils}
  """
}


workflow {
  main:
  SYN_FEATHER = 'data/flywire_synapses_783.feather'
  SYN_ANN = 'data/connections_princeton_no_threshold.csv.gz'
  TYPE_ANN = 'data/consolidated_cell_types.csv.gz'
  utils_file = file('src/utils.r')
  syn_ch = AnnotateNeuropil(file(SYN_FEATHER), file(SYN_ANN), file(TYPE_ANN), utils_file)
  npl_ch = SplitNeuropil(syn_ch, utils_file).flatten()
    | map { it -> [it, utils_file] }
    | RotateNeuropil
 
  publish:
  idv = npl_ch
}

output {
  idv {
    path "idv_mat"
  }
}
