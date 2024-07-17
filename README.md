# Code for "Quantifying higher-order epistasis: beware the chimera"

This repository contains the code for the analyses in the manuscript "Quantifying higher-order epistasis: beware the chimera". In these analyses, we compare the multiplicative and the "chimeric" epistasis measures for quantifying higher-order interactions between mutations, drug, or other perturbations.

Required packages:
* Numpy
* Matplotlib
* Pandas
* itertools
* Scipy
* Seaborn

## Yeast epistasis analysis

The notebook `Kuzmin_data_analysis.ipynb` in the directory `yeast_screens` contains all the analyses to reproduce our results on three-way interactions in yeast. 

> WARNING: refactoring of this code is in progress.

## Drug interaction analysis

The notebook `drug_interaction_reanalysis_additive.ipynb` in the directory `drug_data` contains all the analyses to reproduce our results on higher-order interactions between antibiotics. 

> WARNING: refactoring of this code is in progress.


## Protein epistasis analysis

The notebook `protein_epistasis.ipynb` contains code to reproduce the analyses in the manuscript. We have included all data files in the folder `protein_data`, one file for each protein analyzed. We obtained these files from the original publications; we append the name of the protein to the beginning of the name of the original data file.
