# Preprocessing Module

#### 1. AnnotateNeuropil (`bin/filter_synapses.R`)

**Inputs**: 
- Synapse coordinates from FlyWire dataset (`flywire_synapses_783.feather`)
- List of verified synapses (proofread annotations)
- Cell type classifications
- Utility functions for data processing

**Processing**:
- Keeps only synapses in our regions of interest:
  - ME_L/ME_R (medulla left/right)
  - LO_L/LO_R (lobula left/right)  
  - LOP_L/LOP_R (lobula plate left/right)
- Removes weak connections (4 or fewer synapses between neurons)
- Adds cell type labels for both sides of each connection

**Output**: `ol_passed_filter.csv.gz` - Filtered and annotated synapse data

#### 2. SplitNeuropil (`bin/split_neuropil.R`)

Divides data into separate files by brain region

**Output**: Six files, one per region (e.g., `ME_L.csv.gz`, `LO_R.csv.gz`)

#### 3. RotateNeuropil (`bin/rotate.R`)

**What it does**: Standardizes spatial coordinates across brain regions

**How it works**:
- Combines all synapse locations in a region
- Finds the main directions of variation using PCA
- Rotates all coordinates to align with these main axes
- Adds new columns: `pre_rx`, `pre_ry`, `pre_rz` (rotated pre-synaptic coordinates)
  and `post_rx`, `post_ry`, `post_rz` (rotated post-synaptic coordinates)

**Output**: Rotated coordinate files (`*_rotated.csv.gz`)

### Resource Requirements

- **AnnotateNeuropil**: 64GB RAM (handles millions of synapses)
- **SplitNeuropil**: 8GB RAM
- **RotateNeuropil**: 4GB RAM

## Demo Dataset Generation

To generate a smaller raw input set that still runs the full preprocessing
workflow and produces standard `*_rotated.csv.gz` outputs:

```bash
python src/preprocessing/bin/make_demo_data.py \
  --syn data/fafb_v783_princeton_synapse_table.csv.gz \
  --ann data/connections_princeton_no_threshold.csv.gz \
  --types data/consolidated_cell_types.csv.gz \
  --outdir data/demo_data \
  --seed 1 \
  --me-count 50000 \
  --lo-count 25000 \
  --lop-count 25000
```

Then run preprocessing against demo inputs:

```bash
nextflow run src/preprocessing/main.nf \
  --synf data/demo_data/fafb_v783_princeton_synapse_table.csv.gz \
  --annf data/demo_data/connections_princeton_no_threshold.csv.gz \
  --typef data/demo_data/consolidated_cell_types.csv.gz
```

## Data Source

The initial synapse data comes from the FlyWire connectome project:
- **Download**: [Zenodo Repository](https://zenodo.org/records/10676866)
- **File**: `flywire_synapses_783.feather`
- **Reference**: Dorkenwald, S., Matsliah, A., Sterling, A.R. et al. Neuronal wiring diagram of an adult brain. *Nature* **634**, 124–138 (2024).
