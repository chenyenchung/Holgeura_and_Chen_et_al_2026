## Files

### P15_CAM.csv

Cell Adhesion Molecule (CAM) expression matrix

- **Structure**: Binary matrix (TRUE/FALSE values)
- **Rows**: CAM genes (e.g., "18w", "alpha-Spec", "Apc2", "baz")
- **Columns**: Cell type IDs (numeric identifiers like 153, 200, 128, etc.)
- **Source**: Özel, M.N., Simon, F., Jafari, S. et al. Neuronal diversity and convergence in a visual system developmental atlas. *Nature* **589**, 88–95 (2021). https://doi.org/10.1038/s41586-020-2879-3
- **Used by**:
  - Workflow: `src/visualize/main.nf`
  - Scripts: Passed as parameter `params.camf` to visualization processes

### P15_tf.csv

Transcription Factor (TF) expression matrix

- **Structure**: Binary matrix (TRUE/FALSE values)
- **Rows**: TF genes (e.g., "ab", "Abd-B", "ac", "achi")
- **Columns**: Cell type IDs (numeric identifiers)
- **Source**: Özel, M.N., Simon, F., Jafari, S. et al. Neuronal diversity and convergence in a visual system developmental atlas. *Nature* **589**, 88–95 (2021). https://doi.org/10.1038/s41586-020-2879-3
- **Used by**:
  - Workflow: `src/visualize/main.nf`
  - Scripts: `src/visualize/bin/v_selector.r` (selector gene visualization)
  - Passed as parameter `params.tsf`

### TypeToTypeDistances.csv

Cell type similarity/distance matrix

- **Structure**: Square matrix of pairwise distances
- **Values**: Range from 0 (identical) to 1 (most dissimilar)
- **Rows/Columns**: Visual neuron cell types (Am1, C2, C3, CT1, Dm1, etc.)
- **Source**: Matsliah, A., Yu, Sc., Kruk, K. et al. Neuronal parts list and wiring diagram for a visual system. *Nature* **634**, 166–180 (2024). https://doi.org/10.1038/s41586-024-07981-1
- **Used by**:
  - Workflow: `src/visualize/main.nf`
  - Scripts: `src/visualize/bin/similarity_tree.r` (generates hierarchical clustering trees)
  - Passed as parameter `params.distances`

### visual_neurons_anno.csv

Visual neuron annotations

- **Structure**: Multi-column metadata table
- **Key columns**:
  - `cell_type`: Cell type name
  - `temporal_label`: Temporal origin classification
  - `subsystem`: Functional subsystem (e.g., Motion)
  - `func`: Functional annotation
  - `Confident_annotation`: Confidence flag (Y/N)
  - Additional columns for various annotations and clusters
- **Used by**:
  - Workflows: 
    - `src/visualize/main.nf` (parameter `params.annf`)
    - `src/stats/main.nf` (parameter `params.annf`)
  - Scripts:
    - `src/stats/bin/depth_stats.r` (depth statistics)
    - `src/stats/bin/functional_enrichment.r` (functional enrichment analysis)
    - `src/visualize/bin/v_selector.r` (selector visualization)
    - `src/visualize/bin/similarity_tree.r` (similarity analysis)
    - `src/visualize/bin/visualize.r` (main visualization)

### viz_meta.csv

Visualization spatial metadata defining the field of view and axes in visualization.

- **Structure**: Neuropil-specific visualization parameters
- **Rows**: Different neuropils (ME_L, ME_R, LOP_L, LOP_R)
- **Key columns**:
  - `x_axis`, `y_axis`: Axis definitions
  - `axis_1_func`, `axis_2_func`: Scaling functions
  - `min1`, `max1`, `min2`, `max2`: Spatial boundaries
  - `outlayout`: Layout orientation (landscape/portrait)
- **Used by**:
  - Workflows:
    - `src/stats/main.nf` (parameter `params.metaf`)
    - `src/visualize/main.nf` (parameter `params.metaf`)
  - Scripts:
    - `src/stats/bin/depth_stats.r` (spatial statistics)
    - `src/visualize/bin/v_selector.r` (spatial rendering)
    - `src/visualize/bin/visualize.r` (main visualization rendering)

### viz_preset.csv

Visualization color and filtering presets to generate plots of different comparisons used in this work.

- **Structure**: Preset definitions for different visualization modes
- **Key columns**:
  - `preset`: Preset name
  - `palette`: Color palette function
  - `color_guide`: Legend title
  - `filter_func`: Data filtering function
  - `color_by`: Column to use for coloring
  - `do_highlight`: Whether to highlight specific elements
  - `known_only`: Whether to filter out newly annotated types
- **Used by**:
  - Workflows:
    - `src/stats/main.nf` (parameter `params.vizpresetf`)
    - `src/visualize/main.nf` (parameter `params.vizpresetf`)
  - Scripts:
    - `src/stats/bin/depth_stats.r` (applies color presets to statistics)
    - `src/visualize/bin/visualize.r` (applies presets for visualization)
