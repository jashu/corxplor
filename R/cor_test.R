#' Test for Multiple Correlations Between Paired Samples
#'
#' Test for multiple associations between paired samples, using one of Pearson's
#' product moment correlation coefficient, Kendall's tau or Spearman's rho.
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' cor_test(iris[,-5])
#'
#' # Calculate a bootstrap confidence interval
#' cor_test(iris[,-5], boot_ci = TRUE, n_rep = 1000)
#'
#' # Obtain unadjusted p-values with no correction for false discovery
#' cor_test(iris[,-5], p_adjust = "none")
#'
#' # Use permutation test to adjust p-value for family-wise error
#' cor_test(iris[,-5], p_adjust = "permute", n_perm = 1000)
#'
#' # Use Bonferroni correction to adjust p-value for family-wise error
#' cor_test(iris[,-5], p_adjust = "bonferroni")
#'
#' #' Calculate spearman's rho
#' cor_test(iris[,-5], method = "spearman")
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
#' \code{\link{summarise.cor_list}},  \code{\link{cor_boot}},
#' \code{\link{cor_perm}}, \code{\link[stats]{p.adjust}}
#'
#' @importFrom magrittr `%>%`
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

  if(is.matrix(y)) y <- as.data.frame(y)
  if(use == "everything"){
    x_to_drop <- apply(x, 2, function(a) any(is.na(a)))
    y_to_drop <- NULL
    if(!is.null(y)) y_to_drop <- apply(y, 2, function(a) any(is.na(a)))
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
    if(!is.null(y)) y <- y[, !y_to_drop]
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
  if(is.null(y)) y <- x
  all_pairs <- purrr::cross2(x, y)
  # remove self-correlations from list
  diagonal <- purrr::map_lgl(all_pairs, ~identical(.x[[1]], .x[[2]]))
  all_pairs <- all_pairs[!diagonal]
  pair_names <- list(out$x, out$y) %>% purrr::transpose() %>%
    purrr::map(purrr::as_vector)
  # if applicable, reduce pairs to pairwise complete obs
  if(use == "pairwise.complete.obs"){
    all_pairs <- lapply(all_pairs, function(a){
      x <- a[[1]]; y <- a[[2]]
      pairwise_complete <- !(is.na(x) | is.na(y))
      list(x = x[pairwise_complete], y = y[pairwise_complete])
    })
  }
  out$n <- purrr::map_int(all_pairs, ~ length(.x[[1]]))
  duplicate <- vector("logical", length(pair_names))
  for(i in seq_along(duplicate)){
    if(i > 1){
      duplicate[i] <- any(
        purrr::map_lgl(
          1:(i-1), ~ any(identical(pair_names[[i]], rev(pair_names[[.x]])))
        )
      )
    }
  }
  all_pairs <- all_pairs[!duplicate]
  pair_names <- pair_names[!duplicate]
  # obtain n for each correlation and remove pairs with less than 5 cases
  # do not proceed if there are no cases with at least 5 observations
  low_n <- out$n[!duplicate] < 5
  if(all(low_n)) stop("insufficient number of complete observation pairs")
  if(any(low_n)) warning(paste(
    "not enough complete pairs of observations to test these relationships:\n\t",
    paste0(
      paste0(unlist(pair_names[low_n]), collapse = " <-> "), collapse = "\n\t"),
    sep = ""))
  all_pairs <- all_pairs[!low_n]
  pair_names <- pair_names[!low_n]
  x_names <- purrr::map_chr(pair_names, ~.x[1])
  y_names <- purrr::map_chr(pair_names, ~.x[2])
  if(p_adjust == "permute"){
    perm <- cor_perm(x, y, use, method, ...)
    out$p <- perm$p
    attr(out, "class") <- c("cor_perm", attr(out, "class"))
    attr(out, "n_perm") <- attr(perm, "n_perm")
  } else {
    p <- purrr::map_dbl(all_pairs, function(p, method, ...)
      suppressWarnings(
        stats::cor.test(x = p[[1]], y = p[[2]], method = method, ...))$p.value,
        method = method)
    p <- stats::p.adjust(p, method = p_adjust)
    df_left <- as.data.frame(out)
    df_right <- dplyr::bind_rows(
      dplyr::tibble(x = x_names, y = y_names, p = p),
      dplyr::tibble(x = y_names, y = x_names, p = p)
    )
    df <- dplyr::left_join(df_left, df_right, by = c("x", "y"))
    out$p <- df$p
  }
  attr(out, "class") <- c("cor_test", attr(out, "class"))
  attr(out, "p_adjust") <- p_adjust
  out
}
