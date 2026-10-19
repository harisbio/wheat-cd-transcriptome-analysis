suppressPackageStartupMessages(library(ggplot2))

pal <- list(
  control   = "#1B9E77",   # control samples
  cd        = "#D95F02",   # Cd-treated samples
  up        = "#C0392B",   # upregulated
  down      = "#2E6DA4",   # downregulated
  ns        = "grey78",    # not significant
  total     = "#5B4E77",   # totals / neutral accent
  mid       = "#F7F7F7"    # heatmap midpoint
)
col_cond <- c(Control = pal$control, Cd = pal$cd)
col_sig  <- c(Upregulated = pal$up, Downregulated = pal$down, `Not significant` = pal$ns)
heat_cols <- colorRampPalette(c(pal$down, pal$mid, pal$up))(100)  # z-score: low -> high
ctrl_shades <- colorRampPalette(c("#7FD1B9", pal$control, "#0B5A44"))
cd_shades   <- colorRampPalette(c("#FDB97D", pal$cd, "#8C2D04"))
grad_low_high <- c("#D6E4F0", pal$down)   # sequential scale (e.g. -log10 FDR, gene counts)

theme_pub <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank())

# GO ontology colours (circle plots, ontology annotations)
col_ont <- c(BP = "#2A9D8F", CC = "#E0B24A", MF = "#8DBF7A")
