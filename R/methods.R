#' Summarize Correlation List
#'
#' Summary method for the \code{cor_list} class.
#'
#' The \code{summary} method for the \code{cor_list} class allows for
#' interactive exploration of correlation results. Use
#' \code{\link[dplyr]{select_helpers}} expressions to return specific
#' bivariate relationships that you want to look at.
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' iris_cors <- cor_list(iris[,-5])
#'
#' # Look at all correlations with Sepal.Length
#' summary(iris_cors, x = Sepal.Length)
#'
#' # Look at all correlations between sepal measurements and petal measurements
#' summary(iris_cors, x = starts_with("Sepal"), y = starts_with("Petal"))
#'
#' # Look at all correlations with Sepal.Length, excluding Sepal.Width
#' summary(iris_cors, x = Sepal.Length, y = -Sepal.Width)
#'
#' # Look at all correlations in their original order
#' summary(iris_cors, sort = FALSE)
#'
#' @param object An object of class \code{\link{cor_list}}.
#'
#' @param x Variables to include/exclude from one side of the correlation.
#' You can use the same specifications as in \code{\link[dplyr]{select}}. If
#' missing, defaults to all variables.
#'
#' @param y Variables to include/exclude from the other side of the
#' correlation. You can use the same specifications as in
#' \code{\link[dplyr]{select}}. If missing, defaults to all variables.
#'
#' @param sort Logical value indicating whether the output should be sorted. If
#' \code{TRUE} (the default), output will be sorted by the absolute value of
#' the correlation coefficients (in descending order). If \code{FALSE}, output
#' will appear in the same order as the columns of the matrix or data frame
#' passed to \code{\link{cor}}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @seealso \code{\link{cor_boot}}, \code{\link{cor_list}},
#'   \code{\link[dplyr]{select}}, \code{\link[dplyr]{select_helpers}},
#'
#' @importFrom dplyr everything
#' @export
summary.cor_list <- function (object,
                              x = everything(),
                              y = everything(),
                              sort = TRUE,
                              ...){
  out <- as.data.frame(object)
  x_names <- unique(out$x); y_names <- unique(out$y)
  x <- lazyeval::lazy(x); y <- lazyeval::lazy(y)
  x <- dplyr::select_vars_(x_names, x); y <- dplyr::select_vars_(y_names, y)
  out <- out[out$x %in% x & out$y %in% y,]
  paths <- paste(out$x, out$y, sep = " <-> ")
  reverse_paths <- paste(out$y, out$x, sep = " <-> ")
  out <- cbind(correlates = paths, out, stringsAsFactors = FALSE)
  duplicate <- vector("logical", length(paths))
  for(i in seq_along(paths)){
    if(i > 1){
      duplicate[i] <- reverse_paths[i] %in% paths[1:(i-1)]
    }
  }
  out <- out[!duplicate,]
  out <- dplyr::select(out, -x, -y)
  if(sort) out <- dplyr::arrange(out, dplyr::desc(abs(out[,2])))
  CI <- if("cor_boot" %in% class(object)) attr(object, "CI") else NULL
  n_perm <- if("cor_perm" %in% class(object)) attr(object, "n_perm") else NULL
  p_adjust <- if("cor_test" %in% class(object)) attr(object, "p_adjust") else NULL
  structure(out, class = "cor_list_summary",
            CI = CI, n_perm = n_perm, p_adjust = p_adjust)
}

#' @export
as.data.frame.cor_list <- function (x, ...){
  output <- data.frame(x = x[[1]],
                       y = x[[2]],
                       coef = x[[3]],
                       stringsAsFactors = FALSE)
  if(!is.null(x$n)){
    output$n <- x$n
  }
  if(!is.null(x$lower)){
    output$lower <- x$lower
    output$upper <- x$upper
  }
  if(!is.null(x$p)){
    output$p <- x$p
  }
  coef_name <- attr(x, "coef")
  if(!is.null(coef_name)) names(output)[3] <- coef_name
  output
}

#' @export
as.matrix.cor_list <- function (x, ...){
  x_names <- unique(x$x); y_names <- unique(x$y)
  cor_mat <- matrix(1, nrow = length(x_names), ncol = length(y_names))
  rownames(cor_mat) <- x_names
  colnames(cor_mat) <- y_names
  for(i in seq_along(x$x)){
    cor_mat[x$x[i], x$y[i]] <- x$coef[i]
  }
  return(cor_mat)
}

#' @export
print.cor_list <- function(x, ...){
  temp <- summary(x)
  temp <- purrr::map_if(temp, is.numeric, round, digits = 2)
  output <-  data.frame(correlates = temp[[1]],
                        coef = temp[[2]],
                        stringsAsFactors = FALSE)
  names(output)[2] <- names(temp)[2]
  if(nrow(output) > 50){
    print(utils::head(output, 25, addrownums = FALSE))
    if(nrow(output) > 25){
      cat("............................................\n",
          "\t(", nrow(output) - 50, " correlations omitted)\n",
          "............................................\n",sep="")
      print(utils::tail(output, 25, addrownums = FALSE))
    }
  } else print(output, row.names = FALSE)
}

#' @export
print.cor_boot <- function(x, ...){
  temp <- summary(x)
  temp <- purrr::map_if(temp, is.numeric, round, digits = 2)
  output <-  data.frame(correlates = temp[[1]],
                        coef = temp[[2]],
                        stringsAsFactors = FALSE)
  names(output)[2] <- names(temp)[2]
  output$conf.int <- paste("[", temp$lower, ", ", temp$upper, "]", sep = "")
  names(output)[3] <- paste(attr(x, "CI") * 100, "% CI", sep = "")
  if(nrow(output) > 50){
    print(utils::head(output, 25, addrownums = FALSE))
    if(nrow(output) > 25){
      cat("............................................\n",
          "\t(", nrow(output) - 50, " correlations omitted)\n",
          "............................................\n",sep="")
      print(utils::tail(output, 25, addrownums = FALSE))
    }
  } else print(output, row.names = FALSE)
}

#' @export
print.cor_perm <- function(x, ...){
  temp <- summary(x)
  output <-  data.frame(correlates = temp[[1]],
                        coef = round(temp[[2]],2),
                        p = round(temp$p, 3),
                        stringsAsFactors = FALSE)
  names(output)[2] <- names(temp)[2]
  output$p <- round(output$p, 3)
  min_p <- round(1 / (attr(x, "n_perm") + 1), 3)
  low_p <- output$p <= min_p
  output$p <- as.character(output$p)
  min_p <- as.character(min_p)
  output$p <- gsub("^0.", " .", output$p)
  min_p <- gsub("^0.", "<.", min_p)
  output$p[low_p] <- min_p
  if(nrow(output) > 50){
    print(utils::head(output, 25, addrownums = FALSE))
    if(nrow(output) > 25){
      cat("............................................\n",
          "\t(", nrow(output) - 50, " correlations omitted)\n",
          "............................................\n",sep="")
      print(utils::tail(output, 25, addrownums = FALSE))
    }
  } else print(output, row.names = FALSE)
}


#' @export
print.cor_list_summary <- function(x, ...){
  if(!is.null(x$p)){
    x$p <- round(x$p, 3)
    min_p <- round(1 / (attr(x, "n_perm") + 1), 3)
    if(length(min_p) == 0 || min_p == 0) min_p <- .001
    low_p <- x$p <= min_p
    x$p <- as.character(x$p)
    min_p <- as.character(min_p)
    x$p <- gsub("^0.", " .", x$p)
    min_p <- gsub("^0.", "<.", min_p)
    x$p[low_p] <- min_p
  }
  temp <- purrr::map_if(x, is.numeric, round, digits = 2)
  out <- data.frame(correlates = temp[[1]], coef = temp[[2]],
                    stringsAsFactors = FALSE)
  names(out)[2] <- names(x)[2]
  if(!is.null(x$lower)){
    out$conf.int <- paste("[", temp$lower, ", ",
                             temp$upper, "]", sep = "")
    names(out)[3] <- paste(attr(x, "CI") * 100,
                              "% CI", sep = "")
  }
  if(!is.null(x$n)) out$n <- x$n
  if(!is.null(attr(x, "p_adjust"))) out$`p*` <- x$p else{
    if(!is.null(x$p)) out$p <- x$p
  }
  print(out, row.names = FALSE)
  if(!is.null(attr(x, "p_adjust"))){
    cat(paste("\n\  *p-value adjustment:", attr(x, "p_adjust")))
  }
}

#' @export
print.cor_test <- function(x, ...){
  print(summarize(x))
}
