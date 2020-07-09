#' Title
#'
#' @param object String identifying the type of graph
#' @param prepend Additional arguments to put at the start of the list
#' @param arguments Graph-specific arguments (from function call)
#' @param ... General arguments (from function call)
#'
#' @return
get_options <- function(object, prepend = list(), arguments, ...)  {
  preset_default <-
    list(
      main = NULL,
      title = NULL,
      submain = NULL,
      subtitle = NULL,
      xlab = "Time",
      legend.title = "Strata",
      size = 1,
      color = NULL,
      palette = NULL,
      ggtheme = ggpubr::theme_pubr()

    )

  preset_ggsurv <-
    append(preset_default,
           list(
             censor = NULL,
             censor.shape = NULL,
             censor.size = 4.5
           ))
  preset_forest <-
    list(
      main = "Hazard ratio",
      title = NULL,
      cpositions = NULL,
      fontsize = NULL,
      refLabel = NULL,
      noDigits = NULL
    )
  preset_hazard <-
    append(preset_default,
           list(ylab = "Hazard rate"))

  preset <- switch(
    object,
    "fit" = preset_ggsurv,
    "fit.parametric" = preset_ggsurv,
    "forest" = preset_forest,
    "hazard" = preset_hazard,
    list()
  )

  ellipsis <- list(...)

  # replaces preset with ellipsis arguments, only if they're already named with presets
  ellipsis.replace <- ellipsis[names(ellipsis) %in% names(preset)]
  preset[names(ellipsis.replace)] <- ellipsis.replace

  # does the same for the sublist in arguments named after the object
  if (is.list(arguments) && object %in% names(arguments)) {
    arguments.object <- arguments[[object]]
    arguments.replace <-
      arguments.object[names(arguments.object) %in% names(preset)]
    preset[names(arguments.replace)] <- arguments.replace
  }

  preset <- append(as.list(prepend), preset)

  # cleanup: main/submain -> title/subtitle
  # ggsurv requires it to be title, ggpar allows title by default
  if (!is.null(preset$main) && "title" %in% names(preset)) {
    preset$title <- preset$main
    preset$main <- NULL
  }
  if (!is.null(preset$submain) && "subtitle" %in% names(preset)) {
    preset$subtitle <- preset$submain
    preset$submain <- NULL
  }

  # cleanup: remove NULL values
  preset[sapply(preset, is.null)] <- NULL

  # print(preset)
  return(preset)
}
