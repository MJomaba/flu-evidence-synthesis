StatCI <- ggproto("StatCI", Stat,
                  setup_data = function(data,params) {
                    #data$ci <- as.factor(as.numeric(as.character(data$ci))                  )
                    data <- data[order(data$x,data$y),]
                    data$alpha <- rep(0.1,length(data$x))
                    data
                  },
                  compute_group = function(data, scales) {
                    lw <- seq(1,length(data$x),2)
                    hg <- rev(seq(2,length(data$x),2))
                    data[c(lw,hg),]
                  },
                  default_aes = aes(fill = ..ci.., colour=..ci..),
                  required_aes = c("x", "y","ci")
)

#' @title Draw credibility intervals
#' 
#' @description
#' Will draw credibitlity intervals with different colours for each interval
#' 
#' @inheritParams geom_polygon
#' 
stat_ci <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatCI, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}