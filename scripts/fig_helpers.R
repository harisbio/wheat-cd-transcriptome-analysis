suppressPackageStartupMessages({ library(circlize); library(ragg) })

# abbreviated protein names for circular labels
abbrev <- c(
  "^(glutathione transferase|Glutathione S-transferase).*" = "GST",
  "^GST C-terminal domain-containing.*" = "GST (C-term)",
  "^GST N-terminal domain-containing.*" = "GST (N-term)",
  "^ABC transporter C family member.*" = "ABCC transporter",
  "^ABC transporter domain-containing.*" = "ABC transporter",
  "^MDR-like ABC transporter" = "ABCB (MDR)",
  "^AP2/ERF and B3 domain-containing.*" = "AP2/ERF-B3",
  "^AP2/ERF domain-containing.*" = "AP2/ERF",
  "^BZIP domain-containing.*" = "bZIP",
  "^WRKY domain-containing.*" = "WRKY",
  "^NAC domain-containing.*" = "NAC",
  "^Homeobox-leucine zipper.*" = "HD-Zip",
  "^HTH myb-type domain-containing.*" = "MYB",
  "^Transcription factor MYBS1" = "MYBS1",
  "^LOB domain-containing.*" = "LBD",
  "^DOG1 domain-containing.*" = "DOG1",
  "^BTB domain-containing.*" = "BTB protein",
  "^HSF-type DNA-binding.*" = "HSF",
  "^Bidirectional sugar transporter.*" = "SWEET",
  "^Major facilitator superfamily.*" = "MFS transporter",
  "^Protein DETOXIFICATION.*" = "MATE (DTX)",
  "^Tonoplast intrinsic protein.*" = "TIP aquaporin",
  "^Probable purine permease.*" = "Purine permease",
  "^P-type Cu\\(\\+\\) transporter.*" = "Cu+ P-ATPase",
  "^H\\(\\+\\)-exporting diphosphatase.*" = "V-PPase",
  "^H\\(\\+\\)/Pi cotransporter.*" = "Pi transporter",
  "^Uncoupling protein.*" = "UCP",
  "^Superoxide dismutase.*" = "SOD",
  "^Glutathione peroxidase.*" = "GPX",
  "^Flavone O-methyltransferase.*" = "Flavone OMT",
  "^C2H2-type domain-containing.*" = "C2H2 ZF",
  "^WRKY transcription factor.*" = "WRKY",
  "^BHLH domain-containing.*" = "bHLH",
  "^MYB-related protein.*" = "MYB-related",
  "^MADS-box domain-containing.*" = "MADS-box",
  "^Homeobox domain-containing.*" = "Homeobox",
  "^Sodium/hydrogen exchanger.*" = "NHX",
  "^Cationic amino acid transporter.*" = "CAT (amino acid)",
  "^Amino acid transporter.*" = "AA transporter",
  "^Zinc transporter.*" = "ZIP Zn transporter",
  "^Copper transporter.*" = "COPT",
  "^Delta tonoplast intrinsic protein.*" = "TIP aquaporin",
  "^STAS domain-containing.*" = "SULTR (STAS)",
  "^CRAL-TRIO domain-containing.*" = "CRAL-TRIO",
  "^Glycolipid transfer protein.*" = "GLTP",
  "^WAT1-related protein.*" = "WAT1",
  "^1-Cys peroxiredoxin.*" = "1-Cys PRX",
  "^Glutathione reductase.*" = "GR",
  "^Copper transport protein.*" = "COPT",
  "^Chloride channel.*" = "CLC",
  "^P-type Ca\\(2\\+\\) transporter.*" = "Ca2+ P-ATPase",
  "^V-type proton ATPase.*" = "V-ATPase",
  "^Peroxidase.*" = "Peroxidase",
  "^Catalase.*" = "Catalase",
  "^(Uncharacterized|Predicted protein|Protein of unknown).*" = "Uncharacterized"
)
make_short_name <- function(ann) function(g) {
  nm <- ann$protein_name[match(g, ann$gene)]
  nm[is.na(nm) | nm == ""] <- "Uncharacterized"
  for (p in names(abbrev)) nm <- sub(p, abbrev[[p]], nm, perl = TRUE)
  nm <- ifelse(nchar(nm) > 18, paste0(substr(nm, 1, 16), ".."), nm)
  paste0(nm, " | ", sub("TraesCS", "", g))
}

circ_heatmap <- function(vmat, dirn, cond, labs, file, title, legend = TRUE, lab_cex = 0.92, samp_cex = 0.74, canvas = 2.0) {
  z <- t(scale(t(vmat)))
  hc <- hclust(dist(z)); dend <- as.dendrogram(hc); ord <- order.dendrogram(dend)
  n <- nrow(z); ns <- ncol(z)
  col_fun <- colorRamp2(c(-2, 0, 2), c(pal$down, pal$mid, pal$up))
  dir_col <- c(Up = pal$up, Down = pal$down)

  agg_png(file, width = 9, height = 9, units = "in", res = 400)
  par(mar = c(0, 0, 2.5, 0))
  circos.clear()
  circos.par(gap.after = 34, start.degree = 90, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.004, 0.004), points.overflow.warning = FALSE,
             canvas.xlim = c(-canvas, canvas), canvas.ylim = c(-canvas - 0.15, canvas + 0.2))
  circos.initialize("genes", xlim = c(0, n))
  # gene labels (outermost)
  circos.track(ylim = c(0, 1), track.height = 0.01, bg.border = NA, panel.fun = function(x, y) {
    circos.text(seq_len(n) - 0.5, rep(0, n), labs[ord], facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = lab_cex)
  })
  # direction band
  circos.track(ylim = c(0, 1), track.height = 0.045, bg.border = NA, panel.fun = function(x, y) {
    circos.rect(0:(n - 1), 0, 1:n, 1, col = dir_col[dirn[ord]], border = "white", lwd = 0.5)
  })
  # heatmap: one ring per sample (outermost = first sample)
  circos.track(ylim = c(0, ns), track.height = 0.42, bg.border = NA, panel.fun = function(x, y) {
    for (j in seq_len(ns))
      circos.rect(0:(n - 1), ns - j, 1:n, ns - j + 1, col = col_fun(z[ord, j]), border = "white", lwd = 0.3)
    circos.text(rep(n, ns) + convert_x(1.2, "mm"), ns - seq_len(ns) + 0.5, colnames(z),
                col = col_cond[as.character(cond)], cex = samp_cex, font = 2, adj = c(0, 0.5),
                facing = "inside", niceFacing = TRUE)
  })
  # gene dendrogram
  h <- attr(dend, "height")
  circos.track(ylim = c(0, h), track.height = 0.11, bg.border = NA,
               panel.fun = function(x, y) circos.dendrogram(dend, max_height = h, facing = "inside"))
  circos.clear()
  title(title, cex.main = 1.75)
  if (legend) {
    yy <- seq(-1.58, -1.18, length.out = 101)
    rect(1.30, yy[-101], 1.36, yy[-1], col = col_fun(seq(-2, 2, length.out = 100)), border = NA, xpd = NA)
    text(1.39, c(yy[1], mean(range(yy)), yy[101]), c("-2", "0", "2"), adj = 0, cex = 1.05, xpd = NA)
    text(1.26, -1.09, "Row z-score", adj = 0, cex = 1.15, font = 2, xpd = NA)
    legend(-1.62, -1.12, legend = c("Control", "Cd", "Up", "Down"),
           fill = c(col_cond, dir_col), border = NA, bty = "n", cex = 1.1, ncol = 1,
           title = "Sample / Direction", title.adj = 0, xpd = NA)
  }
  dev.off()
}

# shared legend strip for circular heatmaps (sample rings, direction band, z-score); needs patchwork
circ_legend_strip <- function(title_size = 15, text_size = 14) {
  lg_df <- data.frame(x = 1:4, y = 1, z = c(-2, 0, 1, 2),
                      k = factor(c("Control", "Cd", "Up", "Down"), levels = c("Control", "Cd", "Up", "Down")))
  p <- ggplot(lg_df, aes(x, y)) +
    geom_tile(aes(fill = k)) + geom_point(aes(colour = z), alpha = 0) +
    scale_fill_manual(values = c(col_cond, Up = pal$up, Down = pal$down), name = "Sample ring / Direction band") +
    scale_colour_gradient2(low = pal$down, mid = pal$mid, high = pal$up, limits = c(-2, 2), name = "Row z-score") +
    guides(fill = guide_legend(nrow = 1, order = 1),
           colour = guide_colourbar(order = 2, barwidth = unit(6, "cm"), barheight = unit(0.5, "cm"))) +
    theme_void() + theme(legend.position = "bottom", legend.box = "horizontal",
                         legend.title = element_text(size = title_size, face = "bold"), legend.text = element_text(size = text_size),
                         legend.key.size = unit(0.7, "cm"), legend.spacing.x = unit(1.2, "cm"))
  patchwork::wrap_elements(full = cowplot::get_legend(p))
}
