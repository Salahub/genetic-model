## simple script to generate a plot of the toy model scheme
library(grid)

G <- c(0.06,0.5)
S <- c(0.28,0.5)
T <- c(0.50,0.5)
X <- c(0.72,0.5)
z <- c(0.94,0.5)

nodeText <- function(coords, text, type = "npc", ...) {
    grid.text(text, x = unit(coords[1], units = type),
              y = unit(coords[2], units = type), ...)
}

edge <- function(from, to, isVertical = FALSE, pad = 0.02,
                 type = "npc", a.angle = 30, a.length = 0.015, ...) {
    if (isVertical) {
        x <- c(from[1], to[1])
        y <- c(from[2] - pad, to[2] + pad)
    } else {
        x <- c(from[1] + pad, to[1] - pad)
        y <- c(from[2], to[2])
    }
    grid.lines(x = unit(x, type), y = unit(y, type),
               arrow = arrow(angle = a.angle,
                             length = unit(a.length, units = type)),
               ...)
}

edgeText <- function(text, from, to, pad = 0.1, type = "npc", ...){
    xvals <- c(from[1], to[1])
    yvals <- c(from[2], to[2])
    grid.text(label = text,
              x = unit(mean(xvals), type),
              y = unit(mean(yvals) + pad, type),
              rot = -atan(diff(yvals)/diff(xvals))*180/pi,
              ...)
}

grid.newpage()
nodeText(G, expression(bold(G)))
edge(G, S)
nodeText(S, expression(bold(S)))
edge(S, T)
nodeText(T, expression(bold(T)))
edge(T, X)
nodeText(X, expression(bold(X)))
edge(X, z)
nodeText(z, expression(bold(z)))
edgeText(expression(italic(select)), G, S, gp = gpar(cex = 0.8))
edgeText(expression(italic(annotate)), S, T, gp = gpar(cex = 0.8))
edgeText(expression(italic(encode)), T, X, gp = gpar(cex = 0.8))
edgeText(expression(italic(summarize)), X, z, gp = gpar(cex = 0.8))
# pdf is height = 1 width = 6
