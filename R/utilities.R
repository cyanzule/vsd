#' Title
#'
#' @param formula
#' @param model
#'
#' @return
#'
#' @examples
getStrata <- function(formula, model) {
  if(inherits(formula, "coxph")) {
    # discards any columns not starting with strata()
    # assumedly it's only one but...
    columns <- which(grepl("^strata\\(", colnames(model)))
    if(length(columns) > 1) {
      return(strata(model[, columns]))
    } else if (length(columns) == 1) {
      return(model[, columns])
    }
  } else {
    # discard the left side (which is always the Surv object)
    if(ncol(model) > 1) {
      return(strata(model[, -1]))
    }
  }
}
