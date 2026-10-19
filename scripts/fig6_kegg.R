

suppressPackageStartupMessages({ library(circlize); library(ragg); library(patchwork); library(png); library(grid) })
root <- commandArgs(TRUE)[1]
source(file.path(root, "figure_theme.R"))
kd <- file.path(root, "kegg_ref")
n_bar <- 20; n_chord_path <- 7; n_genes_per_path <- 8

g2u  <- read.delim(file.path(root, "wheat_gene2uniprot.tsv"), stringsAsFactors = FALSE)
conv <- read.delim(file.path(kd, "conv_uniprot_taes.tsv"), header = FALSE, col.names = c("uniprot", "kegg"))
conv$uniprot <- sub("^up:", "", conv$uniprot)
link <- read.delim(file.path(kd, "link_pathway_taes.tsv"), header = FALSE, col.names = c("kegg", "path"))
link$path <- sub("^path:", "", link$path)
pnames <- read.delim(file.path(kd, "list_pathway_taes.tsv"), header = FALSE, col.names = c("path", "name"))
pnames$name <- sub(" - Triticum aestivum \\(bread wheat\\)$", "", pnames$name)

g2k <- merge(g2u, conv, by = "uniprot")[, c("gene", "kegg")]
g2p <- unique(merge(g2k, link, by = "kegg")[, c("gene", "path")])
overview <- grepl("^taes01[12]", g2p$path)

enrich <- function(deg, universe) {
  gp <- g2p[g2p$gene %in% universe, ]
  bg <- unique(gp$gene); hits <- intersect(deg, bg)
  p2g <- split(gp$gene, gp$path); p2g <- p2g[lengths(p2g) >= 5]
  N <- length(bg); n <- length(hits)
  res <- do.call(rbind, lapply(names(p2g), function(p) {
    k <- sum(p2g[[p]] %in% hits); M <- length(p2g[[p]])
    data.frame(Pathway = p, Count = k, SetSize = M, FoldEnrichment = (k / n) / (M / N),
               pvalue = phyper(k - 1, M, N - M, n, lower.tail = FALSE))
  }))
  res$padj <- p.adjust(res$pvalue, "BH")
  res$Name <- pnames$name[match(res$Pathway, pnames$path)]
  res$Overview <- grepl("^taes01[12]", res$Pathway)
  res$Genes <- vapply(res$Pathway, function(p) paste(intersect(p2g[[p]], hits), collapse = ";"), "")
  res <- res[order(res$padj), ]
  attr(res, "n") <- c(DEGs_in_pathways = n, background_in_pathways = N)
  res
}

run <- function(tis, folder) {
  out <- file.path(root, folder); dir.create(out, showWarnings = FALSE)
  degs <- read.csv(file.path(root, paste0(tis, "_DEGs.csv")), row.names = 1)
  universe <- readLines(file.path(root, paste0(tis, "_background.txt")))
  up <- rownames(degs)[degs$log2FC > 0]; dn <- rownames(degs)[degs$log2FC < 0]

  cols <- c("Pathway", "Name", "Overview", "Count", "SetSize", "FoldEnrichment", "pvalue", "padj", "Genes")
  res <- list(All = enrich(rownames(degs), universe), Up = enrich(up, universe), Down = enrich(dn, universe))
  for (k in names(res)) {
    write.csv(res[[k]][res[[k]]$padj < 0.05, cols], file.path(out, paste0(tis, "_KEGG_", k, "_significant.csv")), row.names = FALSE)
    cat(sprintf("%s %s: %d DEGs in pathways / %d background | %d significant pathways (%d non-overview)\n",
                tis, k, attr(res[[k]], "n")[1], attr(res[[k]], "n")[2],
                sum(res[[k]]$padj < 0.05), sum(res[[k]]$padj < 0.05 & !res[[k]]$Overview)))
  }

  # ---- A) bar plot (all DEGs, non-overview pathways) ----
  a <- res$All[res$All$padj < 0.05 & !res$All$Overview, ]
  a <- head(a, n_bar)
  a$nup <- vapply(strsplit(a$Genes, ";"), function(g) sum(g %in% up), 0)
  a$ndn <- a$Count - a$nup
  a$nlp <- -log10(a$padj)
  a$Name <- factor(a$Name, levels = rev(a$Name))
  a$lab <- sprintf("%d  (↑%d / ↓%d)", a$Count, a$nup, a$ndn)
  pA <- ggplot(a, aes(nlp, Name, fill = Count)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = lab), hjust = -0.08, size = 4.2, colour = "grey20") +
    scale_fill_gradient(low = "#D9D3E6", high = pal$total, name = "Gene count") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(title = "KEGG Pathway Enrichment", x = expression(-log[10]~"(FDR)"), y = NULL) +
    theme_pub +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 13.5, colour = "grey10"),
          axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 15),
          legend.title = element_text(size = 14), legend.text = element_text(size = 12),
          panel.grid.major.y = element_blank())
  ggsave(file.path(out, paste0(tis, "_KEGG_barplot.png")), pA, width = 11, height = 8.5, dpi = 400, device = agg_png)

  # ---- B) chord diagram ----
  sel <- head(a[order(a$padj), ], n_chord_path)
  edges <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
    g <- strsplit(sel$Genes[i], ";")[[1]]
    g <- head(g[order(degs[g, "padj"])], n_genes_per_path)
    data.frame(from = as.character(sel$Name[i]), to = g)
  }))
  genes <- unique(edges$to)
  path_cols <- setNames(c("#4C9F70", "#E0B24A", "#5B4E77", "#2A9D8F", "#D95F5F", "#8DBF7A", "#C77CC2", "#7FA7D6")[seq_len(nrow(sel))],
                        as.character(sel$Name))
  lfc <- degs[genes, "log2FC"]
  lim <- ceiling(max(abs(lfc)))
  fc_fun <- colorRamp2(c(-lim, 0, lim), c(pal$down, pal$mid, pal$up))
  grid_col <- c(path_cols, setNames(fc_fun(lfc), genes))
  edges$value <- 1

  draw_chord <- function() {
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(5, nrow(sel) - 1), 14, rep(0.6, length(genes) - 1), 14),
               canvas.xlim = c(-1.6, 1.85), canvas.ylim = c(-1.65, 1.6))
    chordDiagram(edges, order = c(as.character(sel$Name), genes), grid.col = grid_col,
                 col = adjustcolor(path_cols[edges$from], 0.55), transparency = 0,
                 annotationTrack = "grid", annotationTrackHeight = mm_h(4),
                 preAllocateTracks = list(track.height = mm_h(2)))
    lab_pos <- list()
    circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% genes) {
        # genes sit on the left half: read outward from the circle, all in the same direction
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], s, facing = "reverse.clockwise",
                    niceFacing = FALSE, adj = c(1, 0.5), cex = 1.0)
      } else {
        lab_pos[[s]] <<- circlize(CELL_META$xcenter, CELL_META$ylim[1], sector.index = s, track.index = 1)
      }
    })
    # pathway names: horizontal, just outside their colour band, aligned away from the circle;
    # labels whose boxes would overlap are pushed apart vertically (in band order)
    cex_p <- 1.35; pad <- 0.025
    L <- do.call(rbind, lapply(names(lab_pos), function(s) {
      th <- lab_pos[[s]][1, "theta"] * pi / 180; r <- lab_pos[[s]][1, "rou"] + 0.05
      lab <- paste(strwrap(s, 18), collapse = "\n")
      data.frame(s = s, lab = lab, x = r * cos(th), y = r * sin(th),
                 ax = (1 - cos(th)) / 2, ay = (1 - sin(th)) / 2,
                 w = strwidth(lab, cex = cex_p, font = 2), h = strheight(lab, cex = cex_p, font = 2))
    }))
    box <- function(i) c(L$x[i] - L$ax[i] * L$w[i], L$x[i] + (1 - L$ax[i]) * L$w[i],
                         L$y[i] - L$ay[i] * L$h[i], L$y[i] + (1 - L$ay[i]) * L$h[i])
    for (iter in 1:50) {
      moved <- FALSE
      for (i in seq_len(nrow(L))[-1]) for (j in seq_len(i - 1)) {
        a <- box(i); b <- box(j)
        if (a[1] < b[2] && b[1] < a[2] && a[3] < b[4] + pad && b[3] < a[4] + pad) {
          # i follows j clockwise (further down on the right side): push i down, j up
          ov <- (min(a[4], b[4]) - max(a[3], b[3]) + pad) / 2
          L$y[i] <- L$y[i] - ov; L$y[j] <- L$y[j] + ov; moved <- TRUE
        }
      }
      if (!moved) break
    }
    for (i in seq_len(nrow(L)))
      text(L$x[i], L$y[i], L$lab[i], adj = c(L$ax[i], L$ay[i]), cex = cex_p, font = 2,
           col = path_cols[L$s[i]], xpd = NA)
    circos.clear()
  }

  fB <- file.path(out, paste0(tis, "_KEGG_chord.png"))
  agg_png(fB, width = 13, height = 12, units = "in", res = 400)
  par(mar = c(0, 0, 3, 0))
  draw_chord()
  title("KEGG Pathway-Gene Associations", cex.main = 2.3)
  # log2FC legend
  yy <- seq(-1.55, -1.1, length.out = 101)
  rect(1.45, yy[-101], 1.53, yy[-1], col = fc_fun(seq(-lim, lim, length.out = 100)), border = NA, xpd = NA)
  text(1.57, c(yy[1], mean(range(yy)), yy[101]), c(-lim, 0, lim), adj = 0, cex = 1.4, xpd = NA)
  text(1.40, -1.0, "Gene log2FC", adj = 0, cex = 1.6, font = 2, xpd = NA)
  dev.off()

  pB <- wrap_elements(full = rasterGrob(readPNG(fB, native = TRUE), interpolate = TRUE))
  comp <- (pA / pB) + plot_layout(heights = c(1, 1.45)) +
    plot_annotation(tag_levels = "A", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 22))
  ggsave(file.path(out, paste0("Composite_KEGG_", tis, ".png")), comp, width = 13, height = 21, dpi = 300, device = agg_png)
  ggsave(file.path(out, paste0("Composite_KEGG_", tis, ".tif")), comp, width = 13, height = 21, dpi = 300,
         device = "tiff", compression = "lzw")
  write.csv(merge(edges[, c("from", "to")], data.frame(to = genes, log2FC = lfc), by = "to"),
            file.path(out, paste0(tis, "_KEGG_chord_genes.csv")), row.names = FALSE)
}

run("Root",  "Fig6_KEGG_Root_main")
run("Shoot", "FigS_KEGG_Shoot_supplementary")
