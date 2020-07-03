#' Title
#'
#' @param formula TBD
#' @param model TBD
#'
#' @return
getStrata <- function(formula, model) {
  if (inherits(formula, "coxph")) {
    # discards any columns not starting with strata()
    # assumedly it's only one but...
    columns <- which(grepl("^strata\\(", colnames(model)))
    if (length(columns) > 1) {
      return(strata(model[, columns]))
    } else if (length(columns) == 1) {
      return(model[, columns])
    }
  } else {
    # the whole right side of the formula IS the strata
    # discard the left side (which is always the Surv object)
    if (ncol(model) >= 2) {
      if (ncol(model) == 2 && !is.factor(model[, 2])) {
        return(as.factor(model[, 2]))
      }
      return(strata(model[,-1]))
    }
  }
}


#' Title
#'
#' @param plots TBD
#'
#' @return
unlistPlots <- function(plots) {
  result <- list()

  for (name in names(plots)) {
    item <- plots[[name]]


    if (inherits(item, c("ggplot", "ggsurvplot"))) {
      result[[name]] <- list(item)
    } else {
      result[[name]] <- item
    }
  }

  return (unlist(result, recursive = FALSE))
}
