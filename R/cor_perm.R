#' Permutation Tests for Multiple Correlations
#'
#' Calculates the empirical probabilities of obtaining correlation coefficients
#' equal or greater in squared magnitude than the ones observed, given the
#' null hypothesis that the true correlations all equal 0.
#'
#' To calculate permutation-based P values for all bivariate relationships
#' between the columns of \eqn{X} and the columns of \eqn{Y}, the observations
#' in each column of \eqn{X} are randomly reshuffled (which simulates how the
#' data would be distributed if the null hypothesis were true), and all
#' bivariate correlations between the reshuffled \eqn{X} and \eqn{Y} are
#' recalculated, with the largest squared coefficient recorded. This procedure
#' is repeated \code{n_perm} times to construct an empirical sampling
#' distribution of the largest correlation that one obtains by chance when the
#' null hypothesis is true and one computes the same number of correlations. The
#' correlations calculated from the real data are then squared and compared to
#' this distribution to determine an empirical adjusted P value. If \eqn{n}
#' permutations were carried out and \eqn{r} of the \eqn{n} largest squared
#' coefficients are equal to or greater than a given squared coefficient from
#' the actual correlation matrix, then the empirical adjusted P value for that
#' coefficient is given by \eqn{(r + 1)/(n + 1)}, i.e., the frequency of values
#' in the null distribution that are at least as large as what was observed,
#' with an added constant of 1 to prevent a probability of 0 in the event that
#' the observed value lies entirely outside the null distribution.
#'
#' @seealso \code{\link{cor_test}}, \code{\link{cor_list}},
#'
#' @param x a numeric matrix or data frame.
#'
#' @param y NULL (default) or a matrix or data frame with compatible dimensions
#'  to x. The default is equivalent to y = x (but more efficient).
#'
#' @param use an optional character string giving a method for computing
#'  covariances in the presence of missing values. This must be (an abbreviation
#'  of) one of the strings "everything", "all.obs", "complete.obs",
#'  "na.or.complete", or "pairwise.complete.obs" (default).
#'
#' @param method a character string indicating which correlation coefficient is
#'  to be computed. One of "pearson" (default), "kendall", or "spearman": can be
#'  abbreviated.
#'
#' @param n_perm The number of permutations.
#'
#' @param seed Integer used to seed the random number generator.
#'
#' @param n_cores Integer indicating number of processes to be used in parallel
#' processing. Set to 1 fewer than the number of available CPUs by default.
#'
#' @return Object of class "cor_perm", containing a \code{\link{cor_list}}
#' object with an additional vector of empirically adjusted p-values that
#' correct for multiple testing.
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' cor_perm(iris[,-5], n_perm = 1000)
#'
#' @export

cor_perm <- function(x, y = NULL, use = "pairwise", method = "pearson",
                     n_perm = 10000, seed = 42, n_cores = NULL, ...){
  # construct cor_list. do this first so that cor_list can check that `x` and
  # `y` have named columns
  output <- cor_list(x, y, use, method)
  # do not proceed if cor_list returned empty object
  if(length(output$x) == 0){
    warning("There are no complete cases. Permutations not performed.")
    return(output)
  }
  set.seed(seed)
  x_perm <- replicate(n_perm, purrr::map_df(x, sample), simplify = FALSE)
  .perm_cor <- function(X, Y, use, method){
    XY <- suppressWarnings(cor_list(X, Y, use, method))
    max(XY$coef^2)
  }
  if(is.null(y)) y <- x
  if(is.null(n_cores)) n_cores <- parallel::detectCores(logical = FALSE) - 1L
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, "cor_list")
  null_r2 <- parallel::parSapply(cl, x_perm, .perm_cor, Y = y,
                                 use = use, method = method)
  parallel::stopCluster(cl)
  get_p <- function(observed, distribution){
    (sum(distribution >= observed) + 1) / (n_perm + 1)
  }
  output$p <- sapply(output$coef^2, get_p, null_r2)
  structure(output, class = c("cor_perm", "cor_list"), n_perm = n_perm)
}
