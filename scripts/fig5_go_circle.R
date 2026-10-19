suppressPackageStartupMessages({ library(circlize); library(ragg) })
args <- commandArgs(TRUE)
root <- args[1]; tis <- args[2]; folder <- args[3]
n_per <- if (length(args) >= 4) as.integer(args[4]) else 6
source(file.path(root, "figure_theme.R"))

fdr_breaks <- c(0, 2, 4, 6, 8, 10, 15, 20, Inf)
fdr_labels <- c("(0,2]", "(2,4]", "(4,6]", "(6,8]", "(8,10]", "(10,15]", "(15,20]", ">20")

circle_plot <- function(dir, single = TRUE, tag = NULL) {
  x <- read.csv(file.path(root, folder, paste0(tis, "_GO_", dir, "_nonredundant.csv")))
  x <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) head(x[x$Ontology == o, ], n_per)))
  x$nlp <- -log10(x$padj)
  x$bin <- cut(x$nlp, fdr_breaks, labels = fdr_labels, right = TRUE)
  high <- if (dir == "Up") pal$up else pal$down
  fdr_cols <- setNames(colorRampPalette(c("#F4F1EE", high, if (dir == "Up") "#6B0F0A" else "#0E2F52"))(length(fdr_labels) + 1)[-1],
                       fdr_labels)
  sel_col <- pal$total
  xmax <- ceiling(max(x$SetSize) / 100) * 100
  sectors <- factor(x$GO, levels = x$GO)

  file <- file.path(root, folder, paste0("GO_circle_", tis, "_", dir, ".png"))
  if (single) {
    agg_png(file, width = 16, height = 10, units = "in", res = 400)
    layout(matrix(1:2, 1), widths = c(1, 1))
  }

  # ---- circle ----
  par(mar = c(0.5, 0.5, 2, 0.5))
  circos.clear()
  circos.par(gap.degree = 1.2, start.degree = 90, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)
  circos.initialize(sectors, xlim = c(0, xmax))

  # ring 1: GO ID + scale
  circos.track(ylim = c(0, 1), track.height = 0.09, bg.border = NA,
               panel.fun = function(xx, yy) {
                 i <- CELL_META$sector.numeric.index
                 circos.rect(0, 0, xmax, 1, col = col_ont[x$Ontology[i]], border = "white")
                 circos.text(xmax / 2, 0.5, CELL_META$sector.index, cex = 0.9,
                             facing = "bending.inside", niceFacing = TRUE)
                 circos.axis(h = "top", major.at = pretty(c(0, xmax), 4), labels.cex = 0.55,
                             major.tick.length = mm_y(0.8), labels.facing = "clockwise")
               })
  # ring 2: background genes in term, coloured by FDR
  circos.track(ylim = c(0, 1), track.height = 0.09, bg.border = NA,
               panel.fun = function(xx, yy) {
                 i <- CELL_META$sector.numeric.index
                 circos.rect(0, 0.1, x$SetSize[i], 0.9, col = fdr_cols[as.character(x$bin[i])], border = NA)
                 circos.text(x$SetSize[i] / 2, 0.5, as.character(x$SetSize[i]), cex = 0.85, col = if (x$nlp[i] > 6) "white" else "grey20",
                             facing = "bending.inside", niceFacing = TRUE, font = 2)
               })
  # ring 3: DEGs in term
  circos.track(ylim = c(0, 1), track.height = 0.07, bg.border = NA,
               panel.fun = function(xx, yy) {
                 i <- CELL_META$sector.numeric.index
                 circos.rect(0, 0.2, x$Count[i], 0.8, col = sel_col, border = NA)
                 circos.text(x$Count[i] + xmax * 0.02, 0.5, as.character(x$Count[i]), cex = 0.85, adj = c(0, 0.5),
                             facing = "bending.inside", niceFacing = TRUE)
               })
  # ring 4: rich factor
  circos.track(ylim = c(0, 1), track.height = 0.25, bg.border = NA,
               panel.fun = function(xx, yy) {
                 i <- CELL_META$sector.numeric.index
                 circos.rect(0, 0, xmax, 1, col = "grey93", border = "white")
                 for (h in c(0.25, 0.5, 0.75)) circos.lines(c(0, xmax), c(h, h), col = "white", lwd = 0.8)
                 rf <- x$Count[i] / x$SetSize[i]
                 circos.rect(xmax * 0.08, 0, xmax * 0.92, rf, col = col_ont[x$Ontology[i]], border = NA)
               })
  circos.clear()
  title(paste0(if (dir == "Up") "Upregulated" else "Downregulated", " genes (", tis, ")"), cex.main = 2)
  if (!is.null(tag)) mtext(tag, side = 3, adj = 0, line = 0.3, font = 2, cex = 2.2)
  # centre key
  legend("center", legend = c("Background genes in term", "DEGs in term", "Rich factor (0-1)"),
         pch = c(15, 15, 17), col = c(fdr_cols[5], sel_col, col_ont["BP"]), bty = "n", cex = 1.1, pt.cex = 1.6)

  # ---- key panel ----
  par(mar = c(0.5, 0, 2, 0.5))
  plot.new(); plot.window(c(0, 1), c(0, 1))
  legend(0, 1, legend = c("Biological Process", "Cellular Component", "Molecular Function"),
         fill = col_ont, border = NA, bty = "n", title = "Ontology", title.adj = 0, cex = 1.3)
  legend(0.55, 1, legend = rev(fdr_labels), fill = rev(fdr_cols), border = NA, bty = "n",
         title = expression(-log[10]~"(FDR)"), title.adj = 0, cex = 1.15, ncol = 2)
  y <- 0.76
  for (o in c("BP", "CC", "MF")) {
    xo <- x[x$Ontology == o, ]
    if (!nrow(xo)) next
    for (k in seq_len(nrow(xo))) {
      lab <- strwrap(xo$Term[k], 48)
      if (length(lab) > 2) lab <- c(lab[1], paste0(substr(lab[2], 1, 44), "..."))
      rect(0, y - 0.013, 0.028, y + 0.013, col = col_ont[o], border = NA)
      text(0.04, y, xo$GO[k], adj = 0, cex = 1.15, font = 2)
      text(0.27, y - (seq_along(lab) - 1) * 0.024, lab, adj = 0, cex = 1.15)
      y <- y - 0.037 - (length(lab) - 1) * 0.024
    }
    y <- y - 0.012
  }
  if (single) { dev.off(); cat("saved", file, "\n") }
}

circle_plot("Up")
circle_plot("Down")

# composite: A = upregulated, B = downregulated
comp <- file.path(root, folder, paste0("Composite_GO_circle_", tis, ".png"))
agg_png(comp, width = 16, height = 20, units = "in", res = 300)
layout(matrix(1:4, 2, byrow = TRUE), widths = c(1, 1))
circle_plot("Up", single = FALSE, tag = "A)")
circle_plot("Down", single = FALSE, tag = "B)")
dev.off()
cat("saved", comp, "\n")
