#' Geom for nucleotide sequence tracks
#'
#' Renders a nucleotide sequence track using one of two visual styles.
#' Colors follow UCSC Genome Browser conventions by default
#' (A = green, C = blue, G = gold, T = red, N = grey).
#'
#' @inheritParams ggplot2::layer
#' @param mapping Aesthetic mappings. Required: `x` (genomic position),
#'   `label` (nucleotide character). Optional: `fill` (nucleotide color),
#'   `colour` (tile border color, only used when `style = "tile"`).
#' @param data Dataset (default: NULL).
#' @param stat Statistic to use (default: `"identity"`).
#' @param position Position adjustment (default: `"identity"`).
#' @param style Character. Visual style for the sequence track:
#'   - `"text"` (default): bold, colored nucleotide letters with no background
#'     or border. Letter color is taken from the `fill` aesthetic.
#'   - `"tile"`: colored background tiles with letters on top (controlled by
#'     `show_labels` and `label_color`).
#' @param show_labels Logical. Whether to draw nucleotide letter labels.
#'   For `style = "text"`, labels are always drawn when `show_labels = TRUE`.
#'   For `style = "tile"`, labels are drawn on top of the colored tiles.
#'   Default: `TRUE`.
#' @param label_size Numeric. Font size for nucleotide letters (in pt).
#'   Default: `3`.
#' @param label_color Character. Color of the nucleotide letter text, used
#'   only when `style = "tile"`. For `style = "text"`, letter color is taken
#'   from the `fill` aesthetic. Default: `"white"`.
#' @param tile_height Numeric (0–1). Height of each nucleotide tile as a
#'   proportion of the y-axis range. Only used when `style = "tile"`.
#'   Default: `0.8`.
#' @param tile_width Numeric. Width of each nucleotide tile in data (genomic)
#'   units. Only used when `style = "tile"`. Default: `1` (one base pair per tile).
#' @param na.rm Logical. Silently remove `NA` values. Default: `TRUE`.
#' @param ... Other arguments passed on to [ggplot2::layer()].
#'
#' @return A ggplot2 layer.
#' @export
#' @importFrom ggplot2 layer aes ggproto Geom draw_key_rect
#' @importFrom grid rectGrob textGrob grobTree gpar unit
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(
#'   position   = 1:10,
#'   nucleotide = c("A", "T", "G", "C", "A", "T", "G", "C", "N", "A"),
#'   fill       = c("#00A800", "#CC0000", "#CC9900", "#0000CC",
#'                  "#00A800", "#CC0000", "#CC9900", "#0000CC", "#999999", "#00A800")
#' )
#' ggplot(df, aes(x = position, label = nucleotide, fill = fill)) +
#'   geom_sequence() +
#'   ggplot2::scale_fill_identity()
#' }
geom_sequence <- function(
  mapping    = NULL,
  data       = NULL,
  stat       = "identity",
  position   = "identity",
  ...,
  style        = c("text", "tile"),
  show_labels  = TRUE,
  label_size   = 3,
  label_color  = "white",
  tile_height  = 0.8,
  tile_width   = 1,
  na.rm        = TRUE,
  show.legend  = NA,
  inherit.aes  = TRUE
) {
  style <- match.arg(style)

  default_mapping <- aes(
    x     = .data$position,
    label = .data$nucleotide,
    fill  = .data$fill
  )

  if (is.null(mapping)) {
    mapping <- default_mapping
  } else {
    mapping <- utils::modifyList(default_mapping, as.list(mapping))
    mapping <- do.call(aes, mapping)
  }

  ggplot2::layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomSequence,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      style       = style,
      show_labels = show_labels,
      label_size  = label_size,
      label_color = label_color,
      tile_height = tile_height,
      tile_width  = tile_width,
      na.rm       = na.rm,
      ...
    )
  )
}

#' @rdname geom_sequence
#' @format NULL
#' @usage NULL
#' @export
GeomSequence <- ggproto(
  "GeomSequence",
  Geom,
  required_aes = c("x", "label"),
  optional_aes = c("fill", "colour"),
  default_aes  = aes(
    fill      = "#999999",
    colour    = NA,
    linewidth = 0,
    linetype  = 1,
    alpha     = NA
  ),
  draw_key = draw_key_rect,

  draw_panel = function(data, panel_params, coord,
                        style       = "text",
                        show_labels = TRUE,
                        label_size  = 3,
                        label_color = "white",
                        tile_height = 0.8,
                        tile_width  = 1) {

    # ------------------------------------------------------------------
    # "text" style: bold colored letters, no background, no border
    # ------------------------------------------------------------------
    if (style == "text") {
      if (!show_labels) return(grid::nullGrob())

      text_coords <- coord$transform(
        data.frame(x = data$x, y = 0.5),
        panel_params
      )

      return(grid::textGrob(
        label = data$label,
        x     = text_coords$x,
        y     = text_coords$y,
        gp    = grid::gpar(
          col      = data$fill,
          fontsize = label_size * ggplot2::.pt,
          fontface = "bold"
        )
      ))
    }

    # ------------------------------------------------------------------
    # "tile" style: colored background tiles + optional white labels
    # ------------------------------------------------------------------
    half_w <- tile_width / 2
    coords <- coord$transform(
      data.frame(
        xmin = data$x - half_w,
        xmax = data$x + half_w,
        ymin = (1 - tile_height) / 2,
        ymax = (1 + tile_height) / 2
      ),
      panel_params
    )

    tile_grob <- grid::rectGrob(
      x      = (coords$xmin + coords$xmax) / 2,
      y      = (coords$ymin + coords$ymax) / 2,
      width  = coords$xmax - coords$xmin,
      height = coords$ymax - coords$ymin,
      gp     = grid::gpar(
        fill     = scales::alpha(data$fill, data$alpha),
        col      = data$colour,
        lwd      = data$linewidth * ggplot2::.pt,
        lty      = data$linetype
      )
    )

    if (!show_labels) {
      return(tile_grob)
    }

    text_coords <- coord$transform(
      data.frame(x = data$x, y = 0.5),
      panel_params
    )

    label_grob <- grid::textGrob(
      label = data$label,
      x     = text_coords$x,
      y     = text_coords$y,
      gp    = grid::gpar(
        col      = label_color,
        fontsize = label_size * ggplot2::.pt
      )
    )

    grid::grobTree(tile_grob, label_grob)
  }
)
