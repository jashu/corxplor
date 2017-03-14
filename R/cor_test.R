#' Test for Multiple Correlations Between Paired Samples
#'
#' Test for multiple associations between paired samples, using one of Pearson's
#' product moment correlation coefficient, Kendall's tau or Spearman's rho.
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' iris_cors <- cor_test(iris[,-5])
#' summary(iris_cors)
#'
#' # Calculate a bootstrap confidence interval
#' iris_cors <- cor_test(iris[,-5], boot_ci = TRUE, n_rep = 1000)
#' summary(iris_cors)
#'
#' # Obtain unadjusted p-values with no correction for false discovery
#' iris_cors <- cor_test(iris[,-5], p_adjust = "none")
#' summary(iris_cors)
#'
#' # Use permutation test to adjust p-value for family-wise error
#' iris_cors <- cor_test(iris[,-5], p_adjust = "permute", n_perm = 1000)
#' summary(iris_cors)
#'
#' # Use Bonferroni correction to adjust p-value for family-wise error
#' iris_cors <- cor_test(iris[,-5], p_adjust = "bonferroni")
#' summary(iris_cors)
#'
#' #' Calculate spearman's rho
#' iris_cors <- cor_test(iris[,-5], method = "spearman")
#' summary(iris_cors)
#'
#' @inheritParams cor_list
#'
#' @param boot_ci Logical value indicating whether or not to generate a
#'  bootstrapped confidence interval for the correlation coefficient. Defaults
#'  to \code{FALSE}.
#'
#' @param p_adjust Character string naming the correction method to adjust for
#'  multiple correlation tests. In addition to the standard available
#'  \code{\link[stats]{p.adjust.methods}}, one may also choose \code{"permute"},
#'  in which case \code{\link{cor_perm}} will be used to empirically determine
#'  and adjust for the family-wise error rate. Default method is to control
#'  the false discovery rate (\code{"fdr"}) with the Benjamini–Hochberg
#'  (\code{"BH"}) procedure.
#'
#' @param ... Optional named arguments accepted by \code{\link{cor_boot}},
#'  \code{\link{cor_perm}}, and/or \code{\link[stats]{cor.test}}
#'
#' @return a \code{\link{cor_test}} object
#'
#' @seealso \code{\link[stats]{cor.test}}, \code{\link{cor_list}},
#' \code{\link{summary.cor_list}},  \code{\link{cor_boot}},
#' \code{\link{cor_perm}}, \code{\link[stats]{p.adjust}}
#'
#' @export

cor_test <- function(x, y = NULL, use = "pairwise", method = "pearson",
                 boot_ci = FALSE, p_adjust = "fdr", ...){
  p_adjust <- match.arg(p_adjust, c("permute", stats::p.adjust.methods))
  use <- match.arg(use, c("all.obs", "complete.obs", "pairwise.complete.obs",
                       "everything"))
  if(use == "all.obs" && (any(is.na(x)) || any(is.na(y)))){
    stop("missing observations")
  }
  if(use == "complete.obs"){
    x_cols <- 1:ncol(x)
    temp <- dplyr::bind_cols(x, y)
    temp <- stats::na.omit(temp)
    if(nrow(temp) == 0) stop("no cases with complete observations")
    if(nrow(temp) < nrow(x))
      warning(paste(nrow(x) - nrow(temp), "incomplete cases deleted"))
    x <- temp[, x_cols]
    y <- temp[, -x_cols]
    if(length(y) == 0) y <- NULL
  }
  if(is.matrix(x)) x <- as.data.frame(x)
  if(is.null(y))  y <- x
  if(is.matrix(y)) y <- as.data.frame(y)
  if(use == "everything"){
    x_to_drop <- apply(x, 2, function(a) any(is.na(a)))
    y_to_drop <- apply(y, 2, function(a) any(is.na(a)))
    if(all(x_to_drop) || all(y_to_drop)){
      stop("no variables with complete observations")
    }
    if(sum(x_to_drop) > 0 || sum(y_to_drop) > 0){
      warning(paste("The following variables with incomplete observations ",
                    "were deleted:\n\t",
                    paste0(unique(c(names(x)[x_to_drop], names(y)[y_to_drop])),
                           collapse = "\n\t"), sep = ""))
    }
    x <- x[, !x_to_drop]
    y <- y[, !y_to_drop]
  }
  # construct cor_list.
  out <- if(boot_ci){
    suppressWarnings(cor_boot(x, y, use, method, ...))
  } else {
    suppressWarnings(cor_list(x, y, use, method))
  }
  # do not proceed if cor_list returned empty object
  if(length(out$x) == 0){
    stop("There are no complete cases. Testing not performed.")
  }
  # create list of all bivariate comparisons
  all_pairs <- purrr::cross2(x, y)
  pair_names <- purrr::cross2(names(x), names(y))
  # remove self-correlations from list
  diagonal <- sapply(pair_names, function(a) identical(a[[1]], a[[2]]))
  all_pairs <- all_pairs[!diagonal]
  pair_names <- pair_names[!diagonal]
  # if applicable, reduce pairs to pairwise complete obs
  if(use == "pairwise.complete.obs"){
    all_pairs <- lapply(all_pairs, function(a){
      x <- a[[1]]; y <- a[[2]]
      pairwise_complete <- !(is.na(x) | is.na(y))
      list(x = x[pairwise_complete], y = y[pairwise_complete])
    })
  }
  # compute n and join to out
  df_left <- as.data.frame(out)
  df_right <- dplyr::data_frame(
    x = sapply(pair_names, function(x) x[[1]]),
    y = sapply(pair_names, function(y) y[[2]]),
    n = sapply(all_pairs, function(a) length(a[[1]]))
  )
  df <- dplyr::left_join(df_left, df_right, by = c("x", "y"))
  out$n <- df$n
  # remove duplicates from list
  x_names <- df_right$x
  y_names <- df_right$y
  duplicate <- vector("logical", length(pair_names))
  for(i in seq_along(duplicate)){
    if(i > 1){
      duplicate[i] <- any(sapply(1:(i-1), function(a)
        identical(list(y_names[i], x_names[i]), pair_names[[a]])))
    }
  }
  all_pairs <- all_pairs[!duplicate]
  pair_names <- pair_names[!duplicate]
  # obtain n for each correlation and remove pairs with less than 5 cases
  # do not proceed if there are no cases with at least 5 observations
  n <- sapply(all_pairs, function(a) length(a[[1]]))
  low_n <- n < 5
  if(all(low_n)) stop("insufficient number of complete observation pairs")
  if(any(low_n)) warning(paste(
    "not enough complete pairs of observations to test these relationships:\n\t",
    paste0(
      paste0(unlist(pair_names[low_n]), collapse = " <-> "), collapse = "\n\t"),
    sep = ""))
  all_pairs <- all_pairs[!low_n]
  pair_names <- pair_names[!low_n]
  x_names <- sapply(pair_names, function(x) x[[1]])
  y_names <- sapply(pair_names, function(y) y[[2]])
  if(p_adjust == "permute"){
    perm <- cor_perm(x, y, use, method, ...)
    out$p <- perm$p
    attr(out, "class") <- c("cor_perm", attr(out, "class"))
    attr(out, "n_perm") <- attr(perm, "n_perm")
  } else {
    p <- sapply(all_pairs, function(p, method, ...)
      suppressWarnings(
        stats::cor.test(x = p[[1]], y = p[[2]], method = method, ...))$p.value,
        method = method)
    p <- stats::p.adjust(p, method = p_adjust)
    df_left <- as.data.frame(out)
    df_right <- dplyr::bind_rows(
      dplyr::data_frame(x = x_names, y = y_names, p = p),
      dplyr::data_frame(x = y_names, y = x_names, p = p)
    )
    df <- dplyr::left_join(df_left, df_right, by = c("x", "y"))
    out$p <- df$p
  }
  attr(out, "class") <- c("cor_test", attr(out, "class"))
  attr(out, "p_adjust") <- p_adjust
  out
}
