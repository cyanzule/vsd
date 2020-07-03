#' Visualizing Survival Data
#'
#' This function tries to guess the best format to render the inputted survival
#' data analysis data into a graphically pleasing output, while also explaining
#' any data processing that has occured during the way. Or it will one day, for
#' now it only goes \code{"survplot()"} on most objects.
#'
#' More info will come here soon.
#'
#' @param fit The estimated survival curve for the model, defaults to a single
#' @param data Dataframe from where the model fetches its variable, if left
#'   blank will be extracted from the model, if possible
#' @param main Main title for the graphs, defaults to $fit if only one is given
#' @param xlab Label for the x axis, shared among all relevant graphs
#' @param color Color(s) used for graphs
#' @param interactive Allows to explore the generated graphs before returning
#'   (use with \code{"plotly"} for best results)
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
#'                           data = flexsurv::bc, dist = "gengamma"))
vsd <-
  function(fit,
           data = NULL,
           main = NULL,
           xlab = "Time",
           size = 1,
           color = NULL,
           interactive = FALSE,
           ...) {
    plots <- list()
    model <- NULL

    # constructs default-option survival models from common survival-related objects
    fit.original <- fit
    if (inherits(fit.original, 'Surv')) {
      # Surv object
      # TODO: more than just right-censored survival
      data <- as.data.frame(as.matrix(fit))

      formula <- Surv(time, status) ~ 1
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- NULL
      subset <- NULL

      fit$call$formula <- eval(fit$call$formula, data)
    } else if (is.call(fit.original)) {
      # (Assumedly) '~' call (TODO: fail first?)
      if (is.null(data)) {
        # stop("Data structure required with fit object of type call.")
      }

      formula <- eval(fit.original, data)
      fit <- survfit(formula = formula, data = data)
      model <- model.frame(formula, data)
      strata <- getStrata(formula, model)
      subset <- NULL

      fit$call$formula <- eval(fit$call$formula, data)
    } else if (inherits(fit.original, "coxph")) {
      # coxph object
      if (is.null(data)) {
        data <- eval(fit.original$call$data)
        if (is.null(data)) {
          stop(
            "Original data structure couldn't be extracted, supply it to function call instead"
          )
        }
      }

      fit <- survfit(formula = fit.original, data = data)
      formula <- fit.original
      model <- model.frame(formula, data)
      strata <- getStrata(formula, model)
      subset <- NULL # TODO: is this true?
    }

    # survival fit curve
    if (inherits(fit, 'survfit')) {
      # prevents a lot of 'symbol not recognized' bugs later on
      # fit$call$formula<-eval(fit$call$formula, data)

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

        # todo: use grepl with "^coxph(" instead
        if (inherits(eval(formula, data), "coxph")) {
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
      # plots$fit <- survminer::ggsurvplot(fit, data, surv.median.line = "hv", main = main, xlab = xlab, ...)
      plots$fit <-
        survminer::ggsurvplot(
          fit,
          data,
          main = main,
          xlab = xlab,
          size = size,
          palette = color,
          ggtheme = ggpubr::theme_pubr(),
          ...
        )
      # plots$fit$plot <- plots$fit$plot + ggpubr::theme_pubr()


      #### PLOT$FOREST (for coxph)
      if (inherits(formula, "coxph")) {
        if (!is.null(strata)) {
          # strategy: remake formula coxph with strata removed
          # (using call and grep, 'optional:(+\w*)?strata\(.+\)' with '')
          # then do several ggforests, using the strata to filter the *data* off

          plots$forest <- list()

          fit.expression <- deparse(formula$call)
          fit.expression <-
            str2expression(gsub('\\+?\\s?(strata\\(.+\\)) ', '', fit.expression))

          fit.strataless <- eval(fit.expression)

          plots$forest <-
            survminer::ggforest(fit.strataless, data, main = "Hazard ratio (all datapoints)")

          plots$forest.strata

          for (i in levels(strata)) {
            # does a forest for each strata, separatedly!
            subdata <- data[strata == i, ]
            fit.expression[[1]]$data <- subdata
            plots$forest.strata[[i]] <-
              survminer::ggforest(eval(fit.expression),
                                  subdata,
                                  main = paste("Hazard ratio (", i, ")", sep = ""))
          }
        } else {
          plots$forest <- survminer::ggforest(formula, data)
        }
      }

      #### PLOT$HAZARD
      if (is.null(strata)) {
        # make a simple muhaz graphic
        hazard <- muhaz::muhaz(surv$time, surv$status)
        hazard.df <-
          data.frame(
            x = hazard$est.grid,
            y = hazard$haz.est,
            strata = factor(rep("All", length(hazard$est.grid)))
          )
      } else {
        # make several separate hazard maps
        hazard.df <-
          data.frame(x = numeric(),
                     y = numeric(),
                     strata = numeric())

        hazard.count <- table(strata)

        for (i in levels(strata)) {
          # TODO: is it always ten?
          if (hazard.count[[i]] < 10) {
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
            hazard.df.level <-
              data.frame(
                x = hazard$est.grid,
                y = hazard$haz.est,
                strata = rep(i, length(hazard$est.grid))
              )
            hazard.df <- rbind(hazard.df, hazard.df.level)
          }
        }

        hazard.df$strata <- factor(hazard.df$strata, levels(strata))
      }

      plot.hazard <-
        ggplot(hazard.df, aes(x = x, y = y, color = strata)) + geom_line(size = size)

      plots$hazard <-
        ggpubr::ggpar(
          plot.hazard,
          xlab = xlab,
          ylab = "Hazard rate",
          legend.title = "Strata",
          palette = color,
          ggtheme = ggpubr::theme_pubr()
        )


    } else if (inherits(fit, "flexsurvreg")) {
      plots$fit <-
        survminer::ggflexsurvplot(fit,
                                  data,
                                  xlab = xlab,
                                  size = size,
                                  palette = color,
                                  ggtheme = ggpubr::theme_pubr(),
                                  ...)


    } else {
      types <- names(fit.original)
      if (is.null(types)) {
        types <- typeof(fit.original)
      }

      stop("Unknown fit type: [", types, "] ", quote(fit.original))
    }

    if (interactive && interactive()) {
      choices <- unlistPlots(plots)
      choices.whitelist <- c("fit", "hazard")

      repeat {
        choice <-
          utils::menu(names(choices), title = "Pick a graphic (or 0 to exit)")
        if (choice <= 0)
          break
        choice.name <- names(choices)[[choice]]

        plot <- choices[[choice]]

        if(choice.name %in% choices.whitelist) {
          if (requireNamespace("plotly")) {
            if (inherits(plot, "ggsurvplot")) {
              print(plotly::ggplotly(plot$plot))
            } else {
              print(plotly::ggplotly(plot))
            }
          }
        } else {
          print(plot)
        }
      }
    }

    return(plots)
  }
