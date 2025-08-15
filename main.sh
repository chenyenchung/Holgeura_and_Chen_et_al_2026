#!/bin/bash
set -e
nextflow run src/preprocessing/main.nf -N ycc520@nyu.edu
nextflow run src/stats/main.nf -N ycc520@nyu.edu
nextflow run src/visualize/main.nf -N ycc520@nyu.edu
