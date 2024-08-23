#' Plot Heatmap
#'
#' Creates a heatmap plot from a matrix or data frame.
#'
#' @param data A matrix or data frame to be visualized as a heatmap.
#' The row names and column names of the data will be used for labeling
#' the heatmap axes.
#' @param cell_width Numeric value indicating the width of each cell
#' in the heatmap. Default is 1.
#' @param cell_height Numeric value indicating the height of each cell
#' in the heatmap. Default is 1.
#' @param ... Additional arguments to be passed to ggplot2 functions.
#'
#' @return A ggplot2 object representing the heatmap.
#' @export
#'
#' @examples
#' # Example data
#' data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(data) <- paste0("Row", 1:10)
#' colnames(data) <- paste0("Col", 1:10)
#'
#' # Plot heatmap
#' plot <- plot_heatmap(data)
#' print(plot)
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
plot_heatmap <- function(data, cell_width = 1, cell_height = 1, ...) {
  Row <- NULL
  Column <- NULL
  Value <- NULL
  melted_data <- data.frame(
    Row = rep(rownames(data), each = ncol(data)),
    Column = rep(colnames(data), times = nrow(data)),
    Value = as.vector(t(data))
  )

  p <- ggplot(melted_data, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(width = cell_width, height = cell_height) +
    scale_fill_gradient(low = "red", high = "white") +
    labs(fill = "Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  return(p)
}
