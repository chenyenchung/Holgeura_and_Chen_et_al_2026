# FlyWire Visual System Connectome Analysis

> Temporal and Notch identity determine layer targeting and synapse location of medulla neurons
I. Holguera, Y-C. Chen, Y-C-D. Chen, F. Simon, A.G. Gaffney, J.D. Rodas, S. Córdoba, C. Desplan
You can find the manuscript associated with the analysis in [preprint](https://doi.org/10.1101/2025.01.06.631439).

This repository analyzes how connectivity in the fly visual system is associated with their developmental origin using FlyWire connectome.

## Overview

The analysis pipeline consists of three main modules:

### 1. Data Preparation (src/preprocessing/)

- Filters out unreliable connections with fewer than 5 synapses
- Labels each connection with the cell types involved
- Splits data by neuropil for efficient processing
- Rotates 3D coordinates to a Fischbach-esque orientation

### 2. Statistical Analysis (src/stats/)

- Compare how the depths of synapses vary with developmental origin and function
- Tests if neurons born at similar times tend to have similar functions

### 3. Visualization (src/visualize/)

- Make scatter plots showing where synapses are located in a Fischbach-esque perspective
- Color connections by cell type, birth time, or function
- Create trees showing which cell types are most similar in their connectivity patterns
- Generates donut charts showing the partners of neuronal types of interest
- Highlights synapses from neuron expressing specific genes

## System Requirement

To reproduce the analysis, you will need:

### Dependencies

#### System Dependencies
- Nextflow (workflow manager)
- R (version 4.0 or higher)
- C++ compiler (for statistical modules)
- Access to FlyWire dataset files

#### R Packages
- data.table
- ggplot2
- Rcpp
- openxlsx
- dplyr
- patchwork
- cowplot
- ggdendro

### Running the Pipeline
```bash
# Run specific modules
nextflow run src/preprocessing/main.nf
nextflow run src/stats/main.nf
nextflow run src/visualize/main.nf
```

## Data Sources

### Not included in the repository

- **Brain connectivity data**: FlyWire connectome dataset (Dorkenwald et al., *Nature* 2024)
- **Cell type classifications**: Visual neuron annotations from FlyWire consortium

### Included in `data/`

- **Similarity measurements**: Connectivity dissimilarity matrix from Matsliah et al., *Nature* 2024
- **Binarized gene expression**: Transcription factor and cell adhesion molecule patterns (Özel et al., *Nature* 2020)

For detailed information about each module, see the individual README files in the `src/` subdirectories.
