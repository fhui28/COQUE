# COQUE -- Simultaneous coefficient clustering and sparsity for multivariate mixed models

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

[![DOI](https://zenodo.org/badge/667732583.svg)](https://zenodo.org/doi/10.5281/zenodo.8174027)

<!-- badges: end -->

`COQUE` is an R package associated with the manuscript "Simultaneous coefficient clustering and sparsity for multivariate mixed models" by [Hui](https://francishui.netlify.app/), [Dang](https://sites.google.com/view/khuedungdang/home?authuser=0), and [Maestrini](https://sites.google.com/view/lucamaestrini).


# Installation

Currently, `COQUE` is available and can be installed from github with the help of `devtools` package using:

```         
devtools::install_github("fhui28/COQUE")
```

Alternatively, or if the above does not work, you may download a (supposedly) stable release of `COQUE` by choosing the latest release on the right hand side of this Github webpage, and install it manually on your machine.


# Getting started

Currently, there are three directories in this repository:

-   `code`, which contains `R` scripts implementing the proposed COQUE (composite quadratic estimator) method along with other scripts for simulating multivariate longitudinal data, and fitting stacked (potentially penalized) univariate mixed models using [rpql](https://cran.r-project.org/web/packages/rpql/index.html) and [glmmPen](https://cran.r-project.org/web/packages/glmmPen/index.html) among others. The latter are found in the subfolder `alternative_methods`. Many of the functions in these scripts contain pseudo-help files;

-   `simulations`, which contains template scripts to implement the simulation study in the manuscript. **Users are recommended to start here by examining one of `setting1_xxx_n200.R` scripts to understand how to apply COQUE, among other methods.**

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
