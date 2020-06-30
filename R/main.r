#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted survival
#' data analysis data into a graphically pleasing output, while also explaining
#' any data processing that has occured during the way. Or it will one day, for
#' now it only goes \code{"plot()"}.
#'
#' @param fit The estimated survival curve for the model, defaults to a single
#' @param data Dataframe from where the model fetches its variable, if left blank will be extracted from the model, if possible
#' @param ... Kaplain-Meyer estimation of \code{"surv"}
#' @import survival
#' @import ggplot2
#' @export
#' @return A list of graphs, relevant to the model
#' @examples
#' # non-models are cohersed into a survfit with default arguments
#' \dontrun{
#' vsd(with(aml, Surv(time, status)))
#' }
#'
#' # survival fit model
#' vsd(survfit(Surv(futime, fustat) ~ rx, data = ovarian))
#'
#' # coxph (with and without strata)
#' vsd(coxph(Surv(time, status) ~ sex + rx + adhere, data = colon))
#' vsd(survfit(coxph(Surv(time, status) ~ sex + strata(rx) + adhere, data = colon)), colon)
#'
#' # parametric models
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
    fit <- survfit(formula = formula, data = data)
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

    formula <- fit.original
    fit <- survfit(formula = fit.original, data = data)
    model <- model.frame(formula, data)
    strata <- getStrata(formula, model)
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

      strata <- getStrata(formula, model)
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

    if(inherits(formula, "coxph")) {
      if(!is.null(strata)) {
        # strategy: remake formula coxph with strata removed
        # (using call and grep, 'optional:(+\w*)?strata\(.+\)' with '')
        # then do several ggforests, using the strata to filter the *data* off

        plots$forest <- list()

        fit.expression <- deparse(formula$call)
        fit.expression <- str2expression(gsub('\\+?\\s?(strata\\(.+\\)) ', '', fit.expression))

        fit.strataless <- eval(fit.expression)

        plots$forest <- survminer::ggforest(fit.strataless, data, main = "Hazard ratio (all datapoints)")

        plots$forest.strata

        for (i in levels(strata)) {
          # does a forest for each strata, separatedly!
          subdata <- data[strata == i,]
          fit.expression[[1]]$data <- subdata
          plots$forest.strata[[i]] <- survminer::ggforest(eval(fit.expression), subdata, main = paste("Hazard ratio (", i, ")", sep = ""))
        }
      } else {
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
      # make several separate hazard maps
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

    plots$fit <- survminer::ggflexsurvplot(fit, data, ...)
  } else {
    stop("Unknown survival model structure " + names(fit) + " in fit")
  }

  return(plots)
}


