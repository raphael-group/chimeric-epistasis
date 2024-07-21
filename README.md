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

## Citations

The manuscript is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.17.603976v1). If you use our code for your work, please cite our paper.

```
@article {Chitra2024.07.17.603976,
	author = {Chitra, Uthsav and Arnold, Brian J and Raphael, Benjamin},
	title = {Quantifying higher-order epistasis: beware the chimera},
	elocation-id = {2024.07.17.603976},
	year = {2024},
	doi = {10.1101/2024.07.17.603976},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/07/19/2024.07.17.603976},
	eprint = {https://www.biorxiv.org/content/early/2024/07/19/2024.07.17.603976.full.pdf},
	journal = {bioRxiv}
}

```
