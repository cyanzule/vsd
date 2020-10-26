#' @seealso vsd
#' @include options.R
#' @include util.R
NULL

#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted survival
#' data analysis data into a graphically pleasing output.
#'
#' More info will come here soon.
#'
#' @param fit The model, or Surv object, to fit the survival data
#' @param data Dataframe from where the model fetches its variables, if left
#'   blank will be extracted from the model, if possible
#' @param arguments Collection of list of arguments, indexed by the specific
#'   type of graph they should be passed to, has priority over \dots
#' @param interactive Allows to explore the generated graphs before returning
#'   (use with \code{'plotly'} for best results)
#' @param .include Graph types to output if relevant, defaults to all possible
#' @param ... Miscellaneous arguments, passed to ALL graphs
#'
#' @import survival
#' @import ggplot2
#' @importFrom stats model.frame
#' @export
#' @return A list of graphs, relevant to the model
#' @examples
#' # non-models are cohersed into a survfit object with default arguments
#' vsd(with(aml, Surv(time, status)))
#' vsd(Surv(time, status) ~ ph.ecog, data=lung)
#'
#' # survival fit model
#' vsd(survfit(Surv(futime, fustat) ~ rx, data = ovarian))
#'
#' # coxph (with and without strata)
#' vsd(coxph(Surv(time, status) ~ sex + rx + adhere, data = colon))
#' vsd(survfit(coxph(Surv(time, status) ~ sex + strata(rx) + adhere, data = colon)))
#'
#' # parametric models
#' vsd(flexsurv::flexsurvreg(Surv(rectime, censrec) ~ group, data = flexsurv::bc, dist = 'gengamma'))
vsd <- function(fit,
                data = NULL,
                arguments = list(),
                interactive = FALSE,
                .include = c("fit", "parametric", "forest", "residuals", "hazard"),
                ...) {
    UseMethod("vsd")
}

#' @describeIn vsd Wraps \code{Surv ~ (...)} in a survfit object (Kaplan-Meier
#' model)
#' @export
vsd.formula <- function(fit,
                         data = NULL,
                         arguments = list(),
                         interactive = FALSE,
                         .include = c("fit", "hazard"),
                         ...) {
  # (Assumedly) '~' call (TODO: fail first?)
  if (is.null(data)) {
    # stop('Data structure required with fit object of type call.')
  }
  new_model <- survfit(fit, data)
  new_model$call$formula <- eval(fit, data)

  vsd(new_model, data, arguments, interactive, .include, ...)
}

#' @describeIn vsd \code{Surv(...)} in a survfit object (Kaplan-Meier
#' model)
#' @export
vsd.Surv <- function(fit,
                     data = NULL,
                     arguments = list(),
                     interactive = FALSE,
                     .include = c("fit", "hazard"),
                     ...) {
  # Surv object TODO: more than just right-censored survival
  data <- as.data.frame(as.matrix(fit))
  # warning("Fetched `data`: ", data, call. = FALSE, immediate. = TRUE)
  new_model <- survfit(Surv(time, status) ~ 1, data)
  # warning("New model: ", new_model, call. = FALSE, immediate. = TRUE)

  vsd(new_model, data, arguments, interactive, .include, ...)
}

#' @describeIn vsd Wraps \code{coxph(...)} in a survfit object (Kaplan-Meier
#' model)
#' @export
vsd.coxph <- function(fit,
                      data = NULL,
                      arguments = list(),
                      interactive = FALSE,
                      .include = c("fit", "forest", "residuals", "hazard"),
                      ...) {
  if (is.null(data)) {
    data <- eval(fit$call$data)
    if (is.null(data)) {
      stop("Original data structure couldn't be extracted, supply it to function call instead")
    }
  }

  # http://adv-r.had.co.nz/Expressions.html#capturing-call
  new_model <- survfit(fit)
  new_model$call$formula <- substitute(fit)
  vsd(new_model, data, arguments, interactive, .include, ...)
}


# @describeIn vsd Default processing, erroring out
# @export
# vsd.default <- function(fit,
#                         data = NULL,
#                         arguments = list(),
#                         interactive = FALSE,
#                         .include = c("fit", "parametric", "forest", "residuals", "hazard"),
#                         ...) {
#   types <- names(fit)
#   if (is.null(types)) {
#     types <- typeof(fit)
#   }
#
#   stop("Unknown fit type: [", types, "] ", quote(fit))
# }

#' @describeIn vsd Graphical output for survfit objects (Cox model)
#' @export
vsd.survfitcox <- function(fit,
                           data = NULL,
                           arguments = list(),
                           interactive = FALSE,
                           .include = c("fit", "forest", "residuals", "hazard"),
                           ...) {
  plots <- list()
  options <- .options(arguments, ...)
  include <- as.vector(match.arg(.include, several.ok = TRUE))

  # retrieving mid-objects
  if (is.null(data)) {
    data <- eval(eval(fit$call$formula)$call$data)

    if (is.null(data)) {
      stop("Original data structure couldn't be extracted, supply it to function call instead")
    }
  }

  cox_model <- eval(fit$call$formula, data)
  formula <- cox_model$formula
  model <- model.frame(formula, data)
  strata <- .get_strata(cox_model, model)

  # subset <- fit$subset
  surv <- as.data.frame(as.matrix(model[, 1]))
  # surv_type <- fit$type

  #### PLOT$FIT
  if (("fit" %in% include)) {
    plots$fit <-
      do.call(survminer::ggsurvplot,
              append(list(fit, data), options$fit))
  }

  #### PLOT$FOREST
  if (("forest" %in% include)) {
    forest_plots <-
      do.call(plot_forest, append(list(cox_model, data, strata), options$forest))
    plots <- append(plots, forest_plots)
  }

  #### PLOT$RESIDUALS (for coxph)
  if (("residuals" %in% include)) {
    plots$residuals <-
      do.call(survminer::ggcoxzph,
              append(list(cox.zph(cox_model)), options$residuals))
  }

  #### PLOT$HAZARD
  if (("hazard" %in% include)) {
    hazard_plots <-
      do.call(plot_hazard, append(list(surv, strata), options$hazard))
    plots <- append(plots, hazard_plots)
  }

  if (interactive && interactive()) {
    .interactive(plots)
  }
  return(plots)
}

#' @describeIn vsd Graphical output for survfit objects (Kaplan-Meier model)
#' @export
vsd.survfit <- function(fit,
                        data = NULL,
                        arguments = list(),
                        interactive = FALSE,
                        .include = c("fit", "hazard"),
                        ...) {
  plots <- list()
  options <- .options(arguments, ...)
  include <-
    as.vector(match.arg(.include, several.ok = TRUE))

  # retrieving mid-objects
  if (is.null(data)) {
    if (!is.null(fit$call$data)) {
      data <- eval(fit$call$data)
    } else if (is.call(fit$call$formula) &&
               fit$call$formula[[1]] == "coxph") {
      data <- eval(eval(fit$call$formula)$call$data)
    }

    if (is.null(data)) {
      stop("Original data structure couldn't be extracted, supply it to function call instead")
    }
  }

  formula <- fit$call$formula
  model <- model.frame(formula, data)
  strata <- .get_strata(formula, model)

  # subset <- fit$subset
  surv <- as.data.frame(as.matrix(model[, 1]))
  # surv_type <- fit$type


  #### PLOT$FIT
  if (("fit" %in% include)) {
    plots$fit <-
      do.call(survminer::ggsurvplot,
              append(list(fit, data), options$fit))
  }

  #### PLOT$HAZARD
  if (("hazard" %in% include)) {
    hazard_plots <-
      do.call(plot_hazard, append(list(surv, strata), options$hazard))
    plots <- append(plots, hazard_plots)
  }

  if (interactive && interactive()) {
    .interactive(plots)
  }
  return(plots)
}

#' @describeIn vsd Graphical output for flexsurvreg objects (parametric models)
#' @export
vsd.flexsurvreg <- function(fit,
                            data = NULL,
                            arguments = list(),
                            interactive = FALSE,
                            .include = c("fit", "parametric", "hazard"),
                            ...) {
  plots <- list()
  options <- .options(arguments, ...)
  include <-
    as.vector(match.arg(.include, several.ok = TRUE))

  if (is.null(data)) {
    if (!is.null(fit$call$data)) {
      data <- eval(fit$call$data)
    }

    if (is.null(data)) {
      stop("Original data structure couldn't be extracted, supply it to function call instead")
    }
  }

  formula <- eval(fit$call$formula, data)
  model <- model.frame(fit)
  strata <-
    .get_strata(formula, model[,!(names(model) == "(weights)")])
  surv <- as.data.frame(as.matrix(model[, 1]))

  km_fit <- survfit(formula, data)
  km_fit$call$formula <-
    eval(km_fit$call$formula, data)


  #### PLOT$FIT
  if ("fit" %in% include) {
    plot_fit <-
      do.call(survminer::ggsurvplot, append(list(km_fit, data),
                                            options$fit))
    plots$fit <- plot_fit
  }

  #### PLOT$PARAMETRIC
  if (("parametric" %in% include)) {
    parametric_plots <-
      do.call(plot_parametric, append(list(fit, km_fit, strata, data),
                                      options$parametric))
    plots <- append(plots, parametric_plots)
  }

  #### PLOT$HAZARD
  if (("hazard" %in% include)) {
    hazard_plots <-
      do.call(plot_hazard, append(list(surv, strata), options$hazard))
    plots <- append(plots, hazard_plots)
  }

  if (interactive && interactive()) {
    .interactive(plots)
  }
  return(plots)
}


.interactive <- function(plots) {
  # TODO: make choices into two lists: plots, types ?
  choices <- .unlist_plots(plots)
  whitelist <- c("fit", "parametric", "residuals", "hazard")

  repeat {
    choice <-
      utils::menu(names(choices), title = "Pick a graphic (or 0 to exit)")
    if (choice <= 0)
      break
    choice <- choices[[choice]]
    plot <- choice$plot
    type <- choice$type

    if (type %in% whitelist) {
      if (requireNamespace("plotly", quietly = TRUE)) {
        if (inherits(plot, "ggsurvplot")) {
          print(plotly::ggplotly(plot$plot))
        } else {
          print(plotly::ggplotly(plot))
        }
      } else {
        print(plot)
      }
    } else {
      print(plot)
    }
  }
}
