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

## Data Source

The initial synapse data comes from the FlyWire connectome project:
- **Download**: [Zenodo Repository](https://zenodo.org/records/10676866)
- **File**: `flywire_synapses_783.feather`
- **Reference**: Dorkenwald, S., Matsliah, A., Sterling, A.R. et al. Neuronal wiring diagram of an adult brain. *Nature* **634**, 124–138 (2024).
