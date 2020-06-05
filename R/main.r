#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted survival data analysis data into a graphically pleasing output, while also explaining any data processing that has occured during the way.
#' Or it will one day, for now it only goes \code{"plot()"}.
#' @import survival
#' @param surv The input Survival object, from \code{"Surv()"}
#' @param fit The estimated survival curve for the model, defaults to a single Kaplain-Meyer estimation of \code{"surv"}
#' @param \dots Optional arguments to be passed to graphic generation
#' @export
#' @examples
#' # non-models are cohersed into a survfit with default arguments
#' vsd(with(aml, Surv(time, status)))
#' vsd(coxph(Surv(time, status) ~ sex + rx + adhere, data = colon))
#'
#' # example of models (with a repeated coxph)
#' vsd(survfit(Surv(futime, fustat) ~ rx, data = ovarian))
#' vsd(survfit(coxph(Surv(time, status) ~ sex + rx + adhere, data = colon)), colon)
#' vsd(flexsurv::flexsurvreg(Surv(rectime, censrec) ~ group, data = flexsurv::bc, dist = "gengamma"))
vsd <- function(fit, data = NULL, ...) {
  plots <- list()
  model <- NULL

  # constructs default-option survival models from common survival-related objects
  fit.original <- fit
  if(inherits(fit.original, 'Surv')) {
    # Surv object
    # TODO: more than just right-censored survival
    data <- as.data.frame(as.matrix(fit))

    formula <- Surv(time, status) ~ 1
    fit <- survfit(formula, data)
    model <- model.frame(formula, data)
    strata <- NULL
    subset <- NULL
  } else if(inherits(fit.original, "coxph")) {
    # coxph object
    if(is.null(data)) {
      data <- eval(fit.original$call$data)
      if(is.null(data)) {
        stop("Original data structure couldn't be extracted, supply it to function call instead")
      }
    }

    formula <- fit.original$formula
    fit <- survfit(fit.original, data)
    model <- model.frame(formula, data)
    strata <- .getStrata(formula, model)
    subset <- NULL # TODO: is this true?
  }

  # survival fit curve
  if(inherits(fit, 'survfit')) {
    # retrieving mid-objects
    if(is.null(data)) {
      data <- eval(fit$call$data)
      if(is.null(data)) {
        stop("Original data structure couldn't be extracted, supply it to function call instead")
      }
    }

    if(is.null(model)) {
      formula <- fit$call$formula

      # todo: use grepl with "^coxph(" instead
      if(inherits(eval(formula, data), "coxph")) {
        formula <- eval(formula, data)
        model <- model.frame(formula$formula, data)
      } else {
        model <- model.frame(formula, data)
      }

      strata <- .getStrata(formula, model)
    }

    subset <- fit$subset
    surv <- as.data.frame(as.matrix(model[, 1]))
    surv.type <- fit$type

    #### PLOT$FIT
    requireNamespace("survminer", quietly = TRUE)

    plots$fit <- survminer::ggsurvplot(fit, data, surv.median.line = "hv", ...)

    #### PLOT$FOREST (for coxph)
    # BUG: ggforest doesn't like strata...
    # TODO: remake
    requireNamespace("survminer", quietly = TRUE)

    if(is.null(strata)) {
      if(inherits(formula, "coxph")) {
        plots$forest <- survminer::ggforest(formula, data)
      }
    }


    #### PLOT$HAZARD
    requireNamespace("muhaz", quietly = TRUE)

    if(is.null(strata)) {
      # make a simple muhaz graphic
      hazard <- muhaz::muhaz(surv$time, surv$status)
      hazard.df <- data.frame(x = hazard$est.grid, y = hazard$haz.est)

      plots$hazard <- ggplot2::ggplot(hazard.df, ggplot2::aes(x, y, color=I(2)))
    } else {
      # TODO: do several muhaz's
      # append each table, with a schema factor (which increments)
      # that way you hand-craft the melted input table
      # then it's as easy as doing color=schema
      hazard.df <- data.frame(x = numeric(), y = numeric(), strata = numeric())

      for (i in levels(strata)) {
        # creates a sub-table with each muhaz graphic, appends the corresponding strata
        hazard <- muhaz::muhaz(surv$time, surv$status, strata == i)
        hazard.df.level <- data.frame(x = hazard$est.grid, y  = hazard$haz.est)
        hazard.df.level$strata <- rep(i, nrow(hazard.df.level))
        hazard.df <- rbind(hazard.df, hazard.df.level)
      }

      hazard.df$Strata <- factor(hazard.df$strata, levels(strata))

      # ggplot is smart enough to render all subgraphs as long as the data values are in the tall format
      # which they already are by their construction order
      plots$hazard <- ggplot2::ggplot(hazard.df, ggplot2::aes(x, y, color=Strata))
    }

    plots$hazard <- plots$hazard + ggplot2::geom_line(size = 1) + ggplot2::labs(x = "Time", y = "Hazard Rate", title = "Hazard Plot") + ggpubr::theme_pubr()

  } else if(inherits(fit, "flexsurvreg")) {
    requireNamespace("flexsurv", quietly = TRUE)

    plots$fit <- survminer::ggflexsurvplot(fit, model.frame(fit), ...)
  } else {
    stop("Unknown survival model structure " + names(fit) + " in fit")
  }

  return(plots)
}


.getStrata <- function(formula, model) {
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


# eval(survfit(...)$call$formula)


# todo next time:
## do parametric estimation (as in the book)
##
## check surCurve: plots survival models from the survival package and curves of multistate model from the mstate package
## check survMisc
