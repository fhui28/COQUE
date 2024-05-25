# COQUE -- Simultaneous coefficient clustering and sparsity for multivariate mixed models

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11296754.svg)](https://doi.org/10.5281/zenodo.11296754)

<!-- badges: end -->

`COQUE` (Composite Quadratic Estimator) is an R package associated with the manuscript "Simultaneous coefficient clustering and sparsity for multivariate mixed models" by [Hui](https://francishui.netlify.app/), [Dang](https://sites.google.com/view/khuedungdang/home?authuser=0), and [Maestrini](https://sites.google.com/view/lucamaestrini), which is currently in review.

# Installation

Currently, `COQUE` is available and can be installed from github with the help of `devtools` package using:

```         
devtools::install_github("fhui28/COQUE")
```

Alternatively, or if the above does not work, you may download a (hopefully!) stable release of `COQUE` by choosing the latest release on the right hand side of this Github webpage.

# Getting started

1.  For general usage, please check the help file for the main `coque` function after installation.

2.  For code to reproduce the simulation study in the manuscript, please check head to the `manuscript/simulations` folder. Here, there are template `R` scripts, for example:

    -   The script `setting1_binary_n200.R` contains code for simulating multivariate longitudinal binary data with $N = 200$ clusters from a multivariate GLMM, and thenapplying COQUE along with several other methods. Note the script is set up to be run as an array job on a HPC cluster;
    -   The script `setting1_binary_analyzeresults.R` is to be run *after* `setting1_binary_n100.R` and `setting1_binary_n200.R` are completed and the corresponding .`RData` files saved. This script is used to process the results and obtain performance metrics such as root-mean-squared error (RMSE) and $F_1$ score.

3.  Finally, the folder `manuscript/alternative_methods` **can be safely ignored for general usage of COQUE.** The folder contains `R` scripts for implementing some of the competing methods used in the simulation study of the manuscript e.g., fitting stacked (potentially penalized) univariate mixed models using [rpql](https://cran.r-project.org/web/packages/rpql/index.html) and [glmmPen](https://cran.r-project.org/web/packages/glmmPen/index.html), among others.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
