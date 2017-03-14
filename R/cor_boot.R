#' Bootstrapped Confidence Intervals for Correlations
#'
#' Calculates the bias-corrected and accelerated (BCa) bootstrap confidence
#' interval for correlation coefficients.
#'
#' @seealso \code{\link[boot]{boot.ci}}, \code{\link[stats]{cor}},
#'  \code{\link{cor_list}},
#'
#' @param x a numeric matrix or data frame.
#'
#' @param y NULL (default) or a matrix or data frame with compatible dimensions
#'  to x. The default is equivalent to y = x (but more efficient).
#'
#' @param use an optional character string giving a method for handling missing
#'  values. This must be (an abbreviation of) one of the strings "everything",
#'  "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs"
#'  (default). See \code{\link[stats]{cor}} for explanation of these options.
#'
#' @param method a character string indicating which correlation coefficient is
#'  to be computed. One of "pearson" (default), "kendall", or "spearman": can be
#'  abbreviated.
#'
#' @param n_rep The number of bootstrap replicates.
#'
#' @param conf A scalar or vector containing the confidence level(s) of the
#' required interval(s). 95\% intervals used by default.
#'
#' @param seed Integer used to seed the random number generator.
#'
#' @param n_cores Integer indicating number of processes to use in parallel
#' processing. Set to 2 by default. To maximize speed, set to the number of
#' available CPUs.
#'
#' @return Object of class "cor_boot", a \code{\link{cor_list}} object with the
#'  following additional vectors:
#' \describe{
#'  \item{lower}{lower endpoint of the confidence interval}
#'  \item{upper}{upper endpoint of the confidence interval}
#' }
#' and the attribute "CI" indicating the size of the confidence interval being
#' returned.
#' @examples
#' # Correlations with 95% CI for the numeric variables from the iris data set
#' cor_boot(iris[,-5], n_rep = 1000)
#'
#' # with 99% CI
#' cor_boot(iris[,-5], n_rep = 1000, conf = 0.99)
#'
#' @export

cor_boot <- function(x, y = NULL, use = "pairwise", method = "pearson",
                     n_rep = 10000, conf = 0.95, seed = 42, n_cores = 2){
  # construct cor_list. do this first so that cor_list can check that `x` and
  # `y` have named columns
  output <- cor_list(x, y, use, method)
  # do not proceed if cor_list returned empty object
  if(length(output$x) == 0){
    warning("There are no complete cases. Bootstrapping not performed.")
    return(output)
  }
  # calculation of BCa confidence intervals requires n_rep >= sample size
  # in addition, require at least 1000 replicates
  min_rep <- nrow(x); if(min_rep < 1000) min_rep <- 1000
  if(n_rep < min_rep){
    if(min_rep > 1000){
      message(paste("Calculation of BCa intervals requires the number of",
                    "\nbootstrap replicates to equal or exceed sample size."))
    } else {
      message(paste("At least 1000 bootstrap replicates recommended for",
                    "\nthe accurate calculation of BCa confidence intervals."))
    }
    message(paste("\n`n_rep` has been automatically increased from ", n_rep,
                  " to ", min_rep, ".\n", sep = ""))
    n_rep <- min_rep
  }
  # identify which columns in `data` belong to `x`; this is needed to work
  # around the requirement by `boot` that data be passed as a single argument,
  # enabling `boot_cor` to restore `data` to separate `x` & `y` args.
  x_cols <- 1:ncol(x)
  # use tibble data structure to prevent possible coercion to a vector when
  # .boot_cor dereferences data back into x and y parameters
  data <- tibble::as_tibble(x)
  # if y exists, merge with x so that it can be passed to boot as single frame
  if(!is.null(y)) data <- dplyr::bind_cols(data, y)
  # `statistic` function passed to `boot` returns just the coefficient vector
  # from the cor_list object
  .boot_cor <- function(data, rows, ..., x_cols){
    x <- data[rows, x_cols]
    y <- data[rows, -x_cols]; if(ncol(y) == 0) y <- NULL
    suppressWarnings(cor_list(x, y, ...)$coef)
  }
  # ============================================================================
  # parallel bootstrap
  # ----------------------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores)
  set.seed(seed)
  boot_out <- boot::boot(data = data,
                         statistic = .boot_cor,
                         R = n_rep,
                         x_cols = x_cols,
                         use = use,
                         method = method,
                         parallel = "snow",
                         ncpus = n_cores,
                         cl = cl)
  parallel::stopCluster(cl)
  # ============================================================================

  # ============================================================================
  # obtain confidence intervals
  # ----------------------------------------------------------------------------
  .get_ci <- function(i){
    boot_ci <-  boot::boot.ci(boot_out, conf = conf, type = "bca", index = i)
    boot_ci$bca[4:5]
  }
  ci <- lapply(seq_along(boot_out$t0), .get_ci)
  # ============================================================================

  # join confidence intervals to cor_list output
  output$lower <- sapply(ci, function(x) x[1])
  output$upper <- sapply(ci, function(x) x[2])
  structure(output, class = c("cor_boot", "cor_list"), CI = conf)
}
