# ezGenomeTracks - Helper functions (split from helpers.R)
detect_column <- function(data, candidates, param_name, required = TRUE) {
  for (col in candidates) {
    if (col %in% colnames(data)) return(col)
  }
  if (required) {
    stop(paste0(
      "Could not find ",
      param_name,
      " column. Expected one of: ",
      paste(candidates, collapse = ", ")
    ))
  }
  return(NULL)
}

# Internal helper: Validate that required columns exist in a data frame
# Throws an error with informative message if any required column is missing.
#
# @param data A data frame to validate
# @param required_cols Character vector of required column names
# @param data_name Display name of the data (used in error message, e.g. "Data frame", "Coverage data")
validate_required_columns <- function(data, required_cols, data_name = "Data frame") {
  if (!all(required_cols %in% colnames(data))) {
    stop(
      data_name,
      " must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
}

# Internal helper: Apply consistent facet positioning to ggplot objects
# Handles both "left" and "top" positioning with appropriate theme settings.
#
# @param plot A ggplot object to modify
# @param facet_label_position Either "left" or "top"
# @param strip_placement Either "inside" or "outside" (default: "inside")
# @return Modified ggplot object with facet theming applied
apply_facet_position <- function(plot, facet_label_position, strip_placement = "inside") {
  if (facet_label_position == "left") {
    plot <- plot +
      ggplot2::facet_wrap(
        ~track,
        ncol = 1,
        scales = "free_y",
        strip.position = "left"
      ) +
      ggplot2::theme(
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
        strip.placement = strip_placement,
        panel.spacing.y = ggplot2::unit(0, "pt")
      )
  } else {
    plot <- plot +
      ggplot2::facet_wrap(~track, ncol = 1, scales = "free_y") +
      ggplot2::theme(panel.spacing.y = ggplot2::unit(0, "pt"))
  }
  return(plot)
}

# Internal helper: Apply a consistent black border to plot panel
# Used consistently across wrapper functions for uniform visual appearance.
#
# @param plot A ggplot object to modify
# @return Modified ggplot object with black border applied
apply_border_theme <- function(plot) {
  plot +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}

#' Extract signal data for a single input element
