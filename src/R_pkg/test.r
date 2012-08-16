structure(list(axis.line = structure(function (...) 
zeroGrob(), class = "theme", type = "any", call = theme_blank()), 
    axis.text.x = structure(function (label, x = xp, y = yp, 
        ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        colour = "grey50", size = base_size * 0.8, vjust = 1, 
        lineheight = 0.9)), axis.text.y = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        colour = "grey50", size = base_size * 0.8, hjust = 1, 
        lineheight = 0.9)), axis.ticks = structure(function (x0 = 0, 
        y0 = 0, x1 = 1, y1 = 1, ...) 
    {
        segmentsGrob(x0, y0, x1, y1, ..., default.units = "npc", 
            gp = gpar(col = colour, lty = linetype, lwd = size * 
                .pt), )
    }, class = "theme", type = "segment", call = theme_segment(colour = "grey50")), 
    axis.title.x = structure(function (label, x = xp, y = yp, 
        ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size, vjust = 0.5)), axis.title.y = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size, vjust = 0.5, angle = 90)), axis.ticks.length = structure(0.15, unit = "cm", valid.unit = 1L, class = "unit"), 
    axis.ticks.margin = structure(0.1, unit = "cm", valid.unit = 1L, class = "unit"), 
    legend.background = structure(function (x = 0.5, y = 0.5, 
        width = 1, height = 1, ...) 
    {
        rectGrob(x, y, width, height, ..., gp = gpar(lwd = size * 
            .pt, col = colour, fill = fill, lty = linetype), 
            )
    }, class = "theme", type = "box", call = theme_rect(colour = "white")), 
    legend.margin = structure(0.2, unit = "cm", valid.unit = 1L, class = "unit"), 
    legend.key = structure(function (x = 0.5, y = 0.5, width = 1, 
        height = 1, ...) 
    {
        rectGrob(x, y, width, height, ..., gp = gpar(lwd = size * 
            .pt, col = colour, fill = fill, lty = linetype), 
            )
    }, class = "theme", type = "box", call = theme_rect(fill = "grey95", 
        colour = "white")), legend.key.size = structure(1.2, unit = "lines", valid.unit = 3L, class = "unit"), 
    legend.key.height = NULL, legend.key.width = NULL, legend.text = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size * 0.8)), legend.text.align = NULL, legend.title = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        face = "bold", size = base_size * 0.8, hjust = 0)), legend.title.align = NULL, 
    legend.position = "right", legend.direction = NULL, legend.justification = "center", 
    legend.box = NULL, panel.background = structure(function (x = 0.5, 
        y = 0.5, width = 1, height = 1, ...) 
    {
        rectGrob(x, y, width, height, ..., gp = gpar(lwd = size * 
            .pt, col = colour, fill = fill, lty = linetype), 
            )
    }, class = "theme", type = "box", call = theme_rect(fill = "grey90", 
        colour = NA)), panel.border = structure(function (...) 
    zeroGrob(), class = "theme", type = "any", call = theme_blank()), 
    panel.grid.major = structure(function (x = 0:1, y = 0:1, 
        ..., default.units = "npc") 
    {
        polylineGrob(x, y, ..., default.units = default.units, 
            gp = gpar(lwd = size * .pt, col = colour, lty = linetype), 
            )
    }, class = "theme", type = "line", call = theme_line(colour = "white")), 
    panel.grid.minor = structure(function (x = 0:1, y = 0:1, 
        ..., default.units = "npc") 
    {
        polylineGrob(x, y, ..., default.units = default.units, 
            gp = gpar(lwd = size * .pt, col = colour, lty = linetype), 
            )
    }, class = "theme", type = "line", call = theme_line(colour = "grey95", 
        size = 0.25)), panel.margin = structure(0.25, unit = "lines", valid.unit = 3L, class = "unit"), 
    strip.background = structure(function (x = 0.5, y = 0.5, 
        width = 1, height = 1, ...) 
    {
        rectGrob(x, y, width, height, ..., gp = gpar(lwd = size * 
            .pt, col = colour, fill = fill, lty = linetype), 
            )
    }, class = "theme", type = "box", call = theme_rect(fill = "grey80", 
        colour = NA)), strip.text.x = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size * 0.8)), strip.text.y = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size * 0.8, angle = -90)), plot.background = structure(function (x = 0.5, 
        y = 0.5, width = 1, height = 1, ...) 
    {
        rectGrob(x, y, width, height, ..., gp = gpar(lwd = size * 
            .pt, col = colour, fill = fill, lty = linetype), 
            )
    }, class = "theme", type = "box", call = theme_rect(fill = "white", 
        colour = NA)), plot.title = structure(function (label, 
        x = xp, y = yp, ..., vjust = vj, hjust = hj, default.units = "npc") 
    {
        textGrob(label, x, y, hjust = hjust, vjust = vjust, ..., 
            default.units = default.units, gp = gpar(fontsize = size, 
                col = colour, fontfamily = family, fontface = face, 
                lineheight = lineheight), rot = angle)
    }, class = "theme", type = "text", call = theme_text(family = base_family, 
        size = base_size * 1.2)), plot.margin = structure(c(1, 
    1, 0.5, 0.5), unit = "lines", valid.unit = 3L, class = "unit")), .Names = c("axis.line", 
"axis.text.x", "axis.text.y", "axis.ticks", "axis.title.x", "axis.title.y", 
"axis.ticks.length", "axis.ticks.margin", "legend.background", 
"legend.margin", "legend.key", "legend.key.size", "legend.key.height", 
"legend.key.width", "legend.text", "legend.text.align", "legend.title", 
"legend.title.align", "legend.position", "legend.direction", 
"legend.justification", "legend.box", "panel.background", "panel.border", 
"panel.grid.major", "panel.grid.minor", "panel.margin", "strip.background", 
"strip.text.x", "strip.text.y", "plot.background", "plot.title", 
"plot.margin"), class = "options")
