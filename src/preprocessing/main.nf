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

  output:
  path "ol_passed_filter.csv.gz"

  script:
  """
  filter_synapses.R \
    --syn ${syn} \
    --ann ${ann} \
    --type ${type}
  """
}

process SplitNeuropil {
  cpus '1'
  memory '8GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  path syn

  output:
  path "*.csv.gz"

  script:
  """
  split_neuropil.R --syn ${syn}
  """
}

process RotateNeuropil {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/gcc/4.4.0'

  input:
  path syn

  output:
  path "*_rotated.csv.gz"

  script:
  """
  rotate.R --syn ${syn}
  """
}


workflow {
  main:
  SYN_FEATHER = 'data/flywire_synapses_783.feather'
  SYN_ANN = 'data/connections_princeton_no_threshold.csv.gz'
  TYPE_ANN = 'data/consolidated_cell_types.csv.gz'
  syn_ch = AnnotateNeuropil(file(SYN_FEATHER), file(SYN_ANN), file(TYPE_ANN))
  npl_ch = SplitNeuropil(syn_ch).flatten()
    | RotateNeuropil
 
  publish:
  idv = npl_ch
}

output {
  idv {
    path "idv_mat"
  }
}
