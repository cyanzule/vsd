#' @include options.R
#' @include util.R
NULL

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
#' @param arguments Collection of list of arguments, indexed by the specific type of graph they should be passed to, has priority over \dots
#' @param interactive Allows to explore the generated graphs before returning
#'   (use with \code{'plotly'} for best results)
#' @param ... Miscellaneous arguments, passed to ALL graphs
#'
#' @import survival
#' @import ggplot2
#' @importFrom stats model.frame
#' @importFrom magrittr %>%
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
vsd <-
  function(fit,
           data = NULL,
           arguments = list(),
           interactive = FALSE,
           .include = c("fit", "parametric", "forest", "residuals", "hazard"),
           ...) {
    plots <- list()
    options <- .options(arguments, ...)
    include <-
      as.vector(match.arg(.include, several.ok = TRUE))

    # constructs default-option survival models from common survival-related objects
    fit_original <- fit
    if (inherits(fit_original, "Surv")) {
      # Surv object TODO: more than just right-censored survival
      data <- as.data.frame(as.matrix(fit))

      formula <- Surv(time, status) ~ 1
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- NULL

      fit$call$formula <- eval(fit$call$formula, data)
    } else if (is.call(fit_original)) {
      # (Assumedly) '~' call (TODO: fail first?)
      if (is.null(data)) {
        # stop('Data structure required with fit object of type call.')
      }

      formula <- eval(fit_original, data)
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- .get_strata(formula, model)
      # subset <- NULL

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
      strata <- .get_strata(formula, model)
      # subset <- NULL  # TODO: is this true?
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

        strata <- .get_strata(formula, model)
      }

      # subset <- fit$subset
      surv <- as.data.frame(as.matrix(model[, 1]))
      # surv_type <- fit$type

      #### PLOT$FIT
      if (("fit" %in% include)) {
        plots$fit <-
          do.call(survminer::ggsurvplot,
                  append(list(fit, data), options$fit))
      }


      #### PLOT$FOREST (for coxph)
      if (inherits(formula, "coxph")) {
        if (("forest" %in% include)) {
          forest_plots <-
            do.call(plot_forest, append(list(formula, data, strata), options$forest))
          plots <- append(plots, forest_plots)
        }

        if (("residuals" %in% include)) {
          plots$residuals <-
            do.call(survminer::ggcoxzph,
                    append(list(cox.zph(formula)), options$residuals))
        }
      }

      #### PLOT$HAZARD
      if (("hazard" %in% include)) {
        hazard_plots <-
          do.call(plot_hazard, append(list(surv, strata), options$hazard))
        plots <- append(plots, hazard_plots)
      }



    } else if (inherits(fit, "flexsurvreg")) {
      if (is.null(data)) {
        if (!is.null(fit$call$data)) {
          data <- eval(fit$call$data)
        }

        if (is.null(data)) {
          stop(
            "Original data structure couldn't be extracted, supply it to function call instead"
          )
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
          do.call(survminer::ggsurvplot, append(list(km_fit, data), options$fit))
        plots$fit <- plot_fit
      }

      #### PLOT$PARAMETRIC
      if (("parametric" %in% include)) {
        parametric_plots <-
          do.call(plot_parametric, append(list(fit, km_fit, strata, data), options$parametric))
        plots <- append(plots, parametric_plots)
      }

      #### PLOT$HAZARD
      if (("hazard" %in% include)) {
        hazard_plots <-
          do.call(plot_hazard, append(list(surv, strata), options$hazard))
        plots <- append(plots, hazard_plots)
      }


    } else {
      types <- names(fit_original)
      if (is.null(types)) {
        types <- typeof(fit_original)
      }

      stop("Unknown fit type: [", types, "] ", quote(fit_original))
    }

    if (interactive && interactive()) {
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

    return(plots)
  }

# Generates parametric graph, with KM on the background
plot_parametric <-
  function(fit,
           km_fit,
           strata = NULL,
           data,
           alpha,
           conf.int,
           conf.int.km,
           size,
           ...) {
    if (missing(conf.int)) {
      conf.int = TRUE
    }

    plots <- list()

    summary <- summary(fit)
    if (!is.factor(strata)) {
      plots$parametric <- do.call(survminer::ggflexsurvplot,
                                  append(
                                    list(
                                      fit,
                                      data,
                                      size = size,
                                      alpha = alpha,
                                      conf.int = conf.int,
                                      conf.int.km = conf.int.km
                                    ),
                                    list(...)
                                  ))

      # summary <- summary[[1]] %>% dplyr::mutate(strata = "All")
    } else {
      for (level in levels(strata)) {
        summary[[level]] <-
          summary[[level]] %>% dplyr::mutate(strata = level)
      }

      summary <- do.call(rbind, summary)

      plot_fit <-
        do.call(survminer::ggsurvplot, append(
          list(
            km_fit,
            data,
            alpha = alpha / 2,
            size = size,
            conf.int = conf.int.km
          ),
          list(...)
        ))$plot

      plot_parametric <-
        plot_fit + geom_line(aes(time, est, color = strata),
                             data = summary,
                             size = size)

      if (conf.int) {
        plot_parametric <- plot_parametric +
          geom_line(
            aes(time, lcl, color = strata),
            data = summary,
            size = size / 2,
            linetype = "dashed"
          ) +
          geom_line(
            aes(time, ucl, color = strata),
            data = summary,
            size = size / 2,
            linetype = "dashed"
          )
      }
      plots$parametric <- plot_parametric

    }

    return(plots)
  }

# Generates forest plots for coxph model (more than one if strata isn't null)
plot_forest <-
  function(formula, data, strata = NULL, title, ...) {
    plots <- list()

    if (!is.null(strata)) {
      # strategy: remake formula coxph with strata removed (using call and grep, 'optional:(+\w*)?strata\(.+\)'
      # with '') then do several ggforests, using the strata to filter the *data* off

      fit_expression <- deparse(formula$call)
      fit_expression <-
        str2expression(gsub("\\+?\\s?(strata\\(.+\\)) ", "", fit_expression))

      fit_strataless <- eval(fit_expression, data)

      forest_plots <- list()
      class(forest_plots) <- "vsdstrata"

      forest_plots$all <-
        survminer::ggforest(fit_strataless, data, main = title, ...)

      forest_plots$strata <- list()

      for (i in levels(strata)) {
        # does a forest for each strata, separatedly!
        subdata <- data[strata == i,]
        fit_expression[[1]]$data <- subdata

        forest_plots$strata[[i]] <-
          survminer::ggforest(eval(fit_expression),
                              subdata,
                              main = paste(title, i, sep = ", "),
                              ...)
      }

      plots$forest <- forest_plots
    } else {
      plots$forest <- survminer::ggforest(formula, data, main = title)
    }

    return(plots)
  }

# Generates hazard plot
plot_hazard <- function(surv, strata = NULL, size, ...) {
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

  plots$hazard <- ggpubr::ggpar(plot, ...)

  return(plots)
}
