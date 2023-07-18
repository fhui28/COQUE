# COQUE -- Simultaneous coefficient clustering and sparsity for multivariate mixed models

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8092711.svg)](https://doi.org/10.5281/zenodo.8092711) -->

<!-- badges: end -->

<!-- `application_GBR`, which contains scripts applying the proposed method to the Great Barrier Reef data in the manuscript. Note the original data sources is from [Pichler et al., 2007](http://www.frdc.com.au/Archived-Reports/FRDC%20Projects/2003-021-DLD.pdf). We provide an example dataset that possess the same structure as the actual data used for analysis in the manuscript, but the responses have been masked to alter their original values. Also contained in this folder are a set of html files of `plotly` objects showing the regularization path (for the fit to the actual dataset), and hence the clustering/variable selection across the species for each covariate, as a function of the tuning parameter. -->

This repository code contain template code associated with the manuscript "Simultaneous coefficient clustering and sparsity for multivariate mixed models" by [Hui](https://francishui.netlify.app/), [Dang](https://sites.google.com/view/khuedungdang/home?authuser=0), and [Maestrini](https://sites.google.com/view/lucamaestrini).

# Getting started

There are currently three directories in this repository:

-   `code`, which contains `R` scripts implementing the proposed COQUE (composite quadratic estimator) method along with other scripts for simulating multivariate longitudinal data and fitting stacked penalized univariate mixed models using [rpql](https://cran.r-project.org/web/packages/rpql/index.html) and [glmmPen](https://cran.r-project.org/web/packages/glmmPen/index.html). Many of the functions in these scripts contain pseudo-help files;

-   `simulations`, which contains template scripts to implement the simulation study in the manuscript. **Users are recommended to start here by examining one of `setting1_xxx_n200.R` scripts to understand how to apply COQUE, among other methods.**

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem (Francis tends to make a lot of mistakes in my code, so some may be easy amendments!), then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
