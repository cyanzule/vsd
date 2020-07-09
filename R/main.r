#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted survival
#' data analysis data into a graphically pleasing output, while also explaining
#' any data processing that has occured during the way. Or it will one day, for
#' now it only goes \code{'survplot()'} on most objects.
#'
#' More info will come here soon.
#'
#' @param fit The estimated survival curve for the model, defaults to a single
#' @param data Dataframe from where the model fetches its variable, if left
#'   blank will be extracted from the model, if possible
#' @param main Main title for the graphs, defaults to $fit if only one is given
#' @param xlab Label for the x axis, shared among all relevant graphs
#' @param color Color(s) used for graphs
#' @param size Line with with graphs
#' @param interactive Allows to explore the generated graphs before returning
#'   (use with \code{'plotly'} for best results)
#' @param ... Miscellaneous arguments, passed to ALL graphs
#'
#' @importFrom stats model.frame
#' @import survival
#' @import ggplot2
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
#' vsd(survfit(coxph(Surv(time, status) ~ sex + strata(rx) + adhere,
#'             data = colon)))
#'
#' # parametric models
#' vsd(flexsurv::flexsurvreg(Surv(rectime, censrec) ~ group,
#'                           data = flexsurv::bc, dist = 'gengamma'))
vsd <-
  function(fit,
           data = NULL,
           arguments = list(),
           interactive = FALSE,
           ...) {
    plots <- list()
    model <- NULL

    # constructs default-option survival models from common survival-related objects
    fit_original <- fit
    if (inherits(fit_original, "Surv")) {
      # Surv object TODO: more than just right-censored survival
      data <- as.data.frame(as.matrix(fit))

      formula <- Surv(time, status) ~ 1
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- NULL
      subset <- NULL

      fit$call$formula <- eval(fit$call$formula, data)
    } else if (is.call(fit_original)) {
      # (Assumedly) '~' call (TODO: fail first?)
      if (is.null(data)) {
        # stop('Data structure required with fit object of type call.')
      }

      formula <- eval(fit_original, data)
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- get_strata(formula, model)
      subset <- NULL

      fit$call$formula <- eval(fit$call$formula, data)
    } else if (inherits(fit_original, "coxph")) {
      # coxph object
      if (is.null(data)) {
        data <- eval(fit_original$call$data)
        if (is.null(data)) {
          stop(
            "Original data structure couldn't be extracted, supply it to function call instead"
          )
        }
      }

      fit <- survfit(formula = fit_original, data = data)
      formula <- fit_original
      model <- model.frame(formula, data)
      strata <- get_strata(formula, model)
      subset <- NULL  # TODO: is this true?
    }

    # survival fit curve
    if (inherits(fit, "survfit")) {
      # prevents a lot of 'symbol not recognized' bugs later on fit$call$formula<-eval(fit$call$formula, data)

      # retrieving mid-objects
      if (is.null(data)) {
        if (!is.null(fit$call$data)) {
          data <- eval(fit$call$data)
        } else if (is.call(fit$call$formula) &&
                   fit$call$formula[[1]] == "coxph") {
          data <- eval(eval(fit$call$formula)$call$data)
        }

        if (is.null(data)) {
          stop(
            "Original data structure couldn't be extracted, supply it to function call instead"
          )
        }
      }

      if (is.null(model)) {
        formula <- fit$call$formula

        # todo: use grepl with '^coxph(' instead
        if (inherits(eval(formula, data), "coxph")) {
          formula <- eval(formula, data)
          model <- model.frame(formula$formula, data)
        } else {
          model <- model.frame(formula, data)
        }

        strata <- get_strata(formula, model)
      }

      subset <- fit$subset
      surv <- as.data.frame(as.matrix(model[, 1]))
      surv_type <- fit$type

      #### PLOT$FIT
      plots$fit <-
        do.call(survminer::ggsurvplot,
                get_options("fit", list(fit, data), arguments, ...))


      #### PLOT$FOREST (for coxph)
      if (inherits(formula, "coxph")) {
        plots <- append(plots,
                        do.call(
                          plot_forest,
                          get_options("forest", list(formula, data, strata), arguments, ...)
                        ))
      }

      #### PLOT$HAZARD
      plots <- append(plots,
                      do.call(plot_hazard, get_options(
                        "hazard", list(surv, strata), arguments, ...
                      )))



    } else if (inherits(fit, "flexsurvreg")) {
      plots$fit.parametric <-
        do.call(
          survminer::ggflexsurvplot,
          get_options("fit.parametric", list(fit, data), arguments, ...)
        )



    } else {
      types <- names(fit_original)
      if (is.null(types)) {
        types <- typeof(fit_original)
      }

      stop("Unknown fit type: [", types, "] ", quote(fit_original))
    }

    if (interactive && interactive()) {
      choices <- unlist_plots(plots)
      whitelist <- c("fit", "hazard")

      repeat {
        choice <-
          utils::menu(names(choices), title = "Pick a graphic (or 0 to exit)")
        if (choice <= 0)
          break
        choice_name <- names(choices)[[choice]]

        plot <- choices[[choice]]

        if (choice_name %in% whitelist) {
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

    return(plots)
  }

# Generates forest plots for coxph model (more than one if strata isn't null)
plot_forest <-
  function(formula, data, strata, title, ...) {
    plots <- list()

    if (!is.null(strata)) {
      # strategy: remake formula coxph with strata removed (using call and grep, 'optional:(+\w*)?strata\(.+\)'
      # with '') then do several ggforests, using the strata to filter the *data* off

      fit_expression <- deparse(formula$call)
      fit_expression <-
        str2expression(gsub("\\+?\\s?(strata\\(.+\\)) ", "", fit_expression))

      fit_strataless <- eval(fit_expression, data)

      plots$forest <-
        survminer::ggforest(fit_strataless, data, main = title, ...)

      plots$forest.strata <- list()

      for (i in levels(strata)) {
        # does a forest for each strata, separatedly!
        subdata <- data[strata == i,]
        fit_expression[[1]]$data <- subdata

        plots$forest.strata[[i]] <- survminer::ggforest(eval(fit_expression),
                                    subdata, main = paste(title, i, sep = ", "), ...)
      }
    } else {
      plots$forest <- survminer::ggforest(formula, data, main = title, ...)
    }

    return(plots)
  }

# Generates hazard plot
plot_hazard <- function(surv, strata, size, ...) {
  plots <- list()

  if (is.null(strata)) {
    # make a simple muhaz graphic
    hazard <- muhaz::muhaz(surv$time, surv$status)
    hazard_df <-
      data.frame(
        x = hazard$est.grid,
        y = hazard$haz.est,
        strata = factor(rep("All", length(
          hazard$est.grid
        )))
      )
  } else {
    # make several separate hazard maps
    hazard_df <-
      data.frame(x = numeric(),
                 y = numeric(),
                 strata = numeric())

    hazard_count <- table(strata)

    for (i in levels(strata)) {
      # TODO: is it always ten?
      if (hazard_count[[i]] < 10) {
        warning(
          "Level ",
          i,
          " doesn't have enough datapoints to estimate the hazard function",
          call. = FALSE,
          immediate. = TRUE
        )
      } else {
        # creates a sub-table with each muhaz graphic, appends the corresponding strata
        hazard <-
          muhaz::muhaz(surv$time, surv$status, strata == i)
        hazard_df_level <-
          data.frame(
            x = hazard$est.grid,
            y = hazard$haz.est,
            strata = rep(i, length(hazard$est.grid))
          )
        hazard_df <- rbind(hazard_df, hazard_df_level)
      }
    }

    hazard_df$strata <-
      factor(hazard_df$strata, levels(strata))
  }

  plot <-
    ggplot(hazard_df, aes(x, y, color = strata)) + geom_line(size = size)

  plots$hazard <- ggpubr::ggpar(plot,
                                ...)

  return(plots)
}
