# Demo Preprocessing Inputs

This directory contains a reduced raw-input dataset for running the full
pipeline quickly while keeping the same input schema expected by
`src/preprocessing/main.nf`.

## Files

- `fafb_v783_princeton_synapse_table.csv.gz`
- `connections_princeton_no_threshold.csv.gz`
- `consolidated_cell_types.csv.gz`
- `manifest.json`

## Subsampling Targets

- `ME_L`: 50,000 synapses
- `ME_R`: 50,000 synapses
- `LO_L`: 25,000 synapses
- `LO_R`: 25,000 synapses
- `LOP_L`: 25,000 synapses
- `LOP_R`: 25,000 synapses

## Reproducibility

Generated with:

```bash
python src/preprocessing/bin/make_demo_data.py --outdir data/demo_data
```

with random seed `1` (see `manifest.json` for full metadata).
