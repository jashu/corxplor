# corxplor: Interactive tools for exploring correlation matrices

`corxplor` has two main objectives:

- to facilitate the exploration of large correlation matrices by providing a data structure that can be quickly and interactively queried to view different subsets of bivariate relationships.

- to provide resampling procedures to generate bootstrapped confidence intervals for correlation coefficients and to control for false discoveries when calculating p-values for multiple correlations.

## Installation instructions

To install `corxplor` in R, you first need to install the `devtools` package if you haven’t already:

```
	install.packages("devtools")
```

Once `devtools` is installed, use the following command to install `corxplor` on your system:

```
	devtools::install_github("jashu/corxplor")
```
