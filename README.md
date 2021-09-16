# brainscaling
A collection of code and data for "Scaling principles of white matter connectivity in the human and non-human primate brain" (Ardesch et al., 2021)

## Contents

This repository contains the following data and code:
* volumetric and connectivity data of the primate species included in the study (`data` folder)
* code for the main analyses in MATLAB
* code for phylogenetic generalized least squares regression in R (`R` folder)
* precomputed results for phylogenetic generalized least squares regression (`pgls` folder)
* code and data for supplementary analyses (`suppl` folder)

The main scripts used for the analyses presented in the main text are:
* `analysisVolumetrics.m`
* `analysisDensity.m`
* `analysisConnectionLength.m`
* `analysisNetworkMetrics.m`
* `analysisAsymmetry.m`

## Dependencies

Some external toolboxes and packages are necessary to run the scripts.


### MATLAB:
* [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) 

### R:
* [`nlme`](https://CRAN.R-project.org/package=nlme)
* [`phytools`](https://CRAN.R-project.org/package=phytools)
* [`ape`](https://CRAN.R-project.org/package=ape)
* [`geiger`](https://CRAN.R-project.org/package=geiger)
* [`caper`](https://CRAN.R-project.org/package=caper)
* [`treeio`](https://guangchuangyu.github.io/software/treeio/)
* [`ggtree`](https://guangchuangyu.github.io/software/ggtree/)
