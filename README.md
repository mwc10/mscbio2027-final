# Final Project for MSCBIO2027 Bioimaging Spring 2025

The final project for MSCBIO2027 was:

- Read [“Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications”](https://www.nature.com/articles/s41467-024-54145-w)
- Reimplement their R code in Python to generate Figures 4 and 5
- Critique the paper and their results
- Apply other spatial NMF methods to their data to identify something new

This repository holds the code from our group!

## Overview

Each of the folders has code and results for different aspects of the project. All code should be able to run from a cloned repository, but we are dealing with Python here.  

### Data Wrangling
These scripts download and process the raw counts and annotation data from the paper. They convert from the R data into `AnnData` files.

### Regression
This is the morphological regressor used to estimate the annotation of all ST spots.

### Deconv
A re-implementation of the NB deconvolutional model used to generate tumor and stroma pseudobulks.

### Figures
Code and data necessary to recreate figures 4 and 5.

### Spatial NMF
The extension of their data to other spatial transcriptomic methods.

## How to Run 

1) Install (`uv`)[https://docs.astral.sh/uv/]
2) Use `uv` to install dependencies and create the virtual environment `.venv`:
```bash
uv sync
```
3) Follow the instruction in convert-data.ipynb to download and convert the raw counts and annotation data from Wang et al.
4) Have fun!
