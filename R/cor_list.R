#' Correlation List
#'
#' Return the correlations between the columns of a matrix or data frame as
#' a \code{"cor_list"} object.
#'
#' The \code{\link[stats]{cor}} function from the R stats package returns
#' bivariate correlations as a matrix of values, which can be difficult to parse
#' when there are many correlations. \code{cor_list} converts the correlation
#' matrix into a list of all row-column pairs, removing the diagonal entries,
#' and stores the type of correlation coefficient (Pearson's \emph{r}, Kendall's
#' \emph{tau}, or Spearman's \emph{rho}) as an attribute in a \code{cor_list}
#' object with specialized \code{print} and \code{summarize} methods. The
#' \code{print} method automatically outputs unique bivariate relationships,
#' sorted from strongest to weakest. The \code{\link{summarise.cor_list}} method
#' enables the use of \code{\link[tidyselect]{select_helpers}} expressions to return
#' specific bivariate relationships that you want to look at.
#'
#' \code{cor_list} differs from R's native \code{\link[stats]{cor}} in the
#' following ways:
#' \enumerate{
#'   \item \code{x} (and, if given, \code{y}) must be a matrix or data frame
#'     with named columns; unnamed vectors are not allowed. (You can still
#'     pass data for a single variable as \code{x} or \code{y}, but it must
#'     be in the form of a one-column matrix or data frame rather than a
#'     numeric vector.)
#'   \item The default argument for \code{use} is pairwise-complete observations
#'     instead of "everything".
#'   \item The \code{cor_list} object stores the type of correlation
#'     coefficient as an attribute.
#'   \item The \code{cor_list} object has special summarize methods for
#'     selectively viewing relationships of interest. See
#'     \code{\link{summarize.cor_list}}.
#'     }
#'
#' @seealso \code{\link{cor_boot}}, \code{\link{summarise.cor_list}}
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
#' @return A \code{cor_list} object with the following three vectors:
#'  \describe{
#'    \item{\code{x}}{variables that form the rows of the correlation matrix}
#'    \item{\code{y}}{variables that form the columns of the correlation matrix}
#'    \item{\code{coef}}{the correlation coefficient between corresponding
#'      elements of \code{x} and \code{y}}
#'    }
#' and the attribute "coef" indicating the statistic that is being returned
#' (Pearson's \emph{r}, Kendall's \emph{tau}, or Spearman's \emph{rho}).
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' iris_cors <- cor_list(iris[,-5])
#'
#' # Print returns all unique bivariate relationships sorted by strength
#' iris_cors
#'
#' # Look at all correlations with Sepal.Length
#' summarize(iris_cors, x = Sepal.Length)
#'
#' # Look at all correlations between sepal measurements and petal measurements
#' summarize(iris_cors, x = starts_with("Sepal"), y = starts_with("Petal"))
#'
#' # Look at all correlations with Sepal.Length, excluding Sepal.Width
#' summarize(iris_cors, x = Sepal.Length, y = -Sepal.Width)
#'
#' # Look at all correlations in their original order
#' summarize(iris_cors, sort = FALSE)
#'
#' @export
cor_list <- function(x, y = NULL, use = "pairwise", method = "pearson"){
  # require x (and y if it exists) to be matrix or data frame
  if(is.vector(x) || is.vector(y)){
    stop(paste("attempt to pass data as a vector.",
               "\nConvert single vectors to one-column matrix or data frame."))
  }
  # require x (and y if it exists) to have named columns
  if(is.null(colnames(x)) || (!is.null(y) && is.null(colnames(y))))
    stop("data columns must be named")
  # check method argument and store name of statistic
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  stat <- "r"
  if(method != "pearson"){
    if(method == "kendall") stat <- "tau" else stat <- "rho"
  }
  cor_matrix <- stats::cor(x, y = y, use = use, method = method)
  # if cor function is given only one bivariate relationship between two
  # vectors, it returns a scalar; in this case, reformat as matrix with dimnames
  if(class(cor_matrix) == "numeric"){
    x_name <- colnames(x)[1]
    y_name <- colnames(y)
    if(is.null(y_name)) y_name <- colnames(x)[2]
    if(is.null(y_name)) y_name <- "y"
    cor_matrix <- matrix(c(1, rep(output, 2), 1), nrow = 2, ncol = 2,
                         dimnames = list(c(x_name, y_name), c(x_name, y_name)))
  }
  # convert correlation matrix to list
  x <- rep(rownames(cor_matrix), each = ncol(cor_matrix))
  y <- rep(colnames(cor_matrix), times = nrow(cor_matrix))
  coef <- as.vector(cor_matrix)
  output <- tibble::tibble(x = x, y = y, coef = coef)
  # remove self-correlations
  output <- dplyr::filter(output, x != y)
  # if there are missing correlations, issue warning and then remove
  if(any(is.na(coef))){
    warning("One or more correlations could not be computed.")
    output <- dplyr::filter(output, !is.na(coef))
  }
  structure(output, class = "cor_list", coef = stat)
}
