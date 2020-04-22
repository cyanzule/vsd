#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted
#'  survival data analysis data into a graphically pleasing output, while
#'  also explaining any data processing that has occured during the way.
#' Or it will one day, for now it only goes `plot()`.
#' @param surv The input Survival object, from Surv() (left side of the input function on survfit.formula)
#' @param terms Factor that "buckets" the individuals (right side of the function on survfit.formula)
#' @export
#' @examples
#' vsd(Surv(...)): Does a Keplein-Meyer survival fit, also displays a smooth hazard line if possible
vsd <- function(surv, terms = NULL,
                weights = NULL, subset = NULL,
                xlab = "", ylab = "Survival probability",
                legend = NULL) {
  require(survival)
  if (is.Surv(surv)) {
    if (is.null((terms))) {
      state <- par()
      par(mfrow = c(2, 1))
      fit <- survfit.formula(surv ~ 1, weights = weights, subset = subset)
      plot(fit, mark.time = TRUE, xlab = xlab, ylab = ylab)
      title("Kaplan-Meier Survival Curve")

      # smooth hazard line
      m_surv <- as.matrix(surv)
      require(muhaz)
      plot(muhaz(m_surv[, 1], m_surv[, 2]))
      title("Hazard Plot")
      par(state)
      return(fit)
    } else {
      factors <- table(terms)
      colors <- seq(1, nrow(factors))

      fit <- survfit.formula(as.formula(surv ~ terms),
        weights = weights, subset = subset
      )
      plot(fit, mark.time = TRUE, xlab = xlab, ylab = ylab, col = colors)
      if (is.null(legend)) {
        legend <- rownames(factors)
      }
      legend("bottomleft", legend, col = colors, lty = rep(1, nrow(factors)))
      title("Kaplan-Meier Survival Curves")
      return(fit)
    }
  } else {
    warning("Not a necognized survival model (use `Surv`)")
  }
}