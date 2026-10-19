
suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(reshape2); library(patchwork)
  library(ragg); library(circlize); library(png); library(grid)
})

args <- commandArgs(TRUE)
src <- args[1]; out_root <- args[2]
source(file.path(out_root, "figure_theme.R"))

padj_cut <- 0.05
lfc_cut  <- 1.5
n_heat   <- 1000

counts <- read.csv(file.path(src, "gene_counts_integer.csv"), row.names = 1, check.names = FALSE)
meta <- data.frame(
  row.names = colnames(counts),
  condition = c("Control","Control","Control","Control","Cd","Cd","Cd","Cd","Cd","Cd","Control","Control"),
  tissue    = c("Shoot","Root","Root","Root","Shoot","Shoot","Shoot","Root","Root","Root","Shoot","Shoot"),
  rep       = c("1","3","2","1","3","2","1","3","2","1","3","2")
)
meta$condition <- factor(meta$condition, levels = c("Control","Cd"))

circular_heatmap <- function(mat, cond, file, title, res = 600) {
  z <- t(scale(t(mat)))
  hc <- hclust(dist(z), method = "complete")
  dend <- as.dendrogram(hc)
  ord <- order.dendrogram(dend)
  n <- nrow(z); ns <- ncol(z)
  col_fun <- colorRamp2(c(-2, 0, 2), c(pal$down, pal$mid, pal$up))

  agg_png(file, width = 7.5, height = 7, units = "in", res = res)
  par(mar = c(0, 0, 2, 0))
  circos.clear()
  circos.par(gap.after = 32, start.degree = 90, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.005, 0.005), points.overflow.warning = FALSE,
             canvas.xlim = c(-1.05, 1.45), canvas.ylim = c(-1.1, 1.1))
  circos.initialize("genes", xlim = c(0, n))
  # outer band: condition of each sample ring
  circos.track(ylim = c(0, ns), track.height = 0.05, bg.border = NA,
               panel.fun = function(x, y) {
                 for (j in seq_len(ns))
                   circos.rect(0, ns - j, n, ns - j + 1, col = col_cond[as.character(cond[j])], border = "white", lwd = 0.3)
               })
  # heatmap: one ring per sample (outermost = first sample)
  circos.track(ylim = c(0, ns), track.height = 0.40, bg.border = NA,
               panel.fun = function(x, y) {
                 for (j in seq_len(ns))
                   circos.rect(0:(n - 1), ns - j, 1:n, ns - j + 1, col = col_fun(z[ord, j]), border = NA)
                 circos.text(rep(n, ns) + convert_x(1, "mm"), ns - seq_len(ns) + 0.5, colnames(z),
                             col = col_cond[as.character(cond)], cex = 0.5, adj = c(0, 0.5),
                             facing = "inside", niceFacing = TRUE)
               })
  # gene dendrogram inside
  h <- attr(dend, "height")
  circos.track(ylim = c(0, h), track.height = 0.22, bg.border = NA,
               panel.fun = function(x, y) circos.dendrogram(dend, max_height = h, facing = "inside"))
  circos.clear()
  title(title, cex.main = 1.2)
  legend(1.08, 0.95, legend = names(col_cond), fill = col_cond, border = NA,
         bty = "n", title = "Condition", cex = 0.85, title.adj = 0)
  yy <- seq(-0.45, 0.25, length.out = 101)
  rect(1.12, yy[-101], 1.18, yy[-1], col = col_fun(seq(-2, 2, length.out = 100)), border = NA)
  text(1.21, c(-0.45, -0.10, 0.25), c("-2", "0", "2"), adj = 0, cex = 0.8)
  text(1.10, 0.36, "Row z-score", adj = 0, cex = 0.85)
  dev.off()
}

run_tissue <- function(tis, folder) {
  out_dir <- file.path(out_root, folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  sp <- function(p, f, w, h) ggsave(file.path(out_dir, paste0(f, ".png")), p,
                                    width = w, height = h, dpi = 600, device = agg_png)

  m <- meta[meta$tissue == tis, ]
  m <- m[order(m$condition, m$rep), ]
  m$label <- paste0(ifelse(m$condition == "Control", "Ctrl", "Cd"), "-", m$rep, " (", rownames(m), ")")

  dds <- DESeqDataSetFromMatrix(counts[, rownames(m)], m, ~ condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "condition_Cd_vs_Control", alpha = padj_cut)
  shr <- lfcShrink(dds, coef = "condition_Cd_vs_Control", type = "apeglm", quiet = TRUE)
  vsd <- vst(dds, blind = TRUE)

  df <- data.frame(baseMean = res$baseMean, log2FC_unshrunk = res$log2FoldChange,
                   log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
                   pvalue = res$pvalue, padj = res$padj, row.names = rownames(res))
  df <- df[!is.na(df$padj), ]
  df$sig <- "Not significant"
  df$sig[df$padj < padj_cut & df$log2FC >=  lfc_cut] <- "Upregulated"
  df$sig[df$padj < padj_cut & df$log2FC <= -lfc_cut] <- "Downregulated"
  df$sig <- factor(df$sig, levels = names(col_sig))
  n_up <- sum(df$sig == "Upregulated"); n_dn <- sum(df$sig == "Downregulated")

  # DEG tables / lists for DAVID, KEGG, STRING (in the parent folder)
  deg <- df[df$sig != "Not significant", ]
  write.csv(df[order(df$padj), ], file.path(out_root, paste0(tis, "_all_genes_DESeq2.csv")))
  write.csv(deg[order(deg$padj), ], file.path(out_root, paste0(tis, "_DEGs.csv")))
  writeLines(rownames(deg)[deg$sig == "Upregulated"],   file.path(out_root, paste0(tis, "_DEGs_up.txt")))
  writeLines(rownames(deg)[deg$sig == "Downregulated"], file.path(out_root, paste0(tis, "_DEGs_down.txt")))
  writeLines(rownames(df), file.path(out_root, paste0(tis, "_background.txt")))

  # A. Density of log2 normalized counts
  nc <- log2(counts(dds, normalized = TRUE) + 1)
  colnames(nc) <- m[colnames(nc), "label"]
  long <- melt(nc, varnames = c("gene", "Sample"), value.name = "value")
  long$Sample <- factor(long$Sample, levels = m$label)
  pA <- ggplot(long, aes(value, colour = Sample)) +
    geom_density(linewidth = 0.7) +
    scale_colour_manual(values = setNames(c(ctrl_shades(3), cd_shades(3)), m$label)) +
    labs(title = "Density Plot", x = "Log2 normalized counts", y = "Density") +
    theme_pub + theme(legend.text = element_text(size = 9))

  # B. PCA (VST, top 500 variable genes)
  pc <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pv <- round(100 * attr(pc, "percentVar"))
  pB <- ggplot(pc, aes(PC1, PC2, colour = condition)) +
    geom_point(size = 4) +
    scale_colour_manual(values = col_cond, name = "Condition") +
    labs(title = "Principal Component Analysis",
         x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance")) +
    theme_pub

  # C. DEG summary
  sm <- data.frame(Category = factor(c("Downregulated", "Upregulated", "Total DEGs"),
                                     levels = c("Downregulated", "Upregulated", "Total DEGs")),
                   Count = c(n_dn, n_up, n_up + n_dn))
  pC <- ggplot(sm, aes(Category, Count, fill = Category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = format(Count, big.mark = ",")), vjust = -0.4, size = 4.5) +
    scale_fill_manual(values = c(Downregulated = pal$down, Upregulated = pal$up, `Total DEGs` = pal$total)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Differential Expression Summary", x = NULL, y = "Number of genes") +
    theme_pub + theme(legend.position = "none")

  # D. Volcano (shrunken log2FC)
  df$nlp <- -log10(df$padj)
  df$nlp[!is.finite(df$nlp)] <- max(df$nlp[is.finite(df$nlp)])
  ns <- df[df$sig == "Not significant", ]; sg <- df[df$sig != "Not significant", ]
  pD <- ggplot(df, aes(log2FC, nlp, colour = sig)) +
    geom_point(data = ns, size = 0.6, alpha = 0.5) +
    geom_point(data = sg, size = 0.8, alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    scale_colour_manual(values = col_sig, name = NULL, drop = FALSE) +
    labs(title = "Volcano Plot", x = "Shrunken log2 fold change (Cd vs Control)",
         y = "-Log10 adjusted p-value") +
    theme_pub + theme(legend.position = "top") +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  # E. MA (shrunken log2FC)
  pE <- ggplot(df, aes(baseMean, log2FC, colour = sig)) +
    geom_point(data = ns, size = 0.6, alpha = 0.5) +
    geom_point(data = sg, size = 0.8, alpha = 0.7) +
    scale_x_log10() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    scale_colour_manual(values = col_sig, name = NULL, drop = FALSE) +
    labs(title = "MA Plot", x = "Mean of normalized counts (log scale)",
         y = "Shrunken log2 fold change") +
    theme_pub + theme(legend.position = "top") +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  # F. Circular heatmap of the most variable genes (no DE selection)
  va <- assay(vsd)
  topv <- names(sort(apply(va, 1, var), decreasing = TRUE))[1:n_heat]
  mat <- va[topv, rownames(m)]
  colnames(mat) <- m$label
  fF <- file.path(out_dir, "F_circular_heatmap_top_variable.png")
  circular_heatmap(mat, m$condition, fF, paste0("Top ", n_heat, " Most Variable Genes"))

  sp(pA, "A_density", 7, 4.5); sp(pB, "B_PCA", 6, 4.5); sp(pC, "C_DEG_summary", 6, 4.5)
  sp(pD, "D_volcano", 6, 5);  sp(pE, "E_MA", 6, 5)

  fF_small <- tempfile(fileext = ".png")
  circular_heatmap(mat, m$condition, fF_small, paste0("Top ", n_heat, " Most Variable Genes"), res = 250)
  rm(dds, vsd, va, nc, long); invisible(gc())
  pF <- wrap_elements(full = rasterGrob(readPNG(fF_small), interpolate = TRUE))
  comp <- (pA + pB) / (pC + pD) / (pE + pF) +
    plot_annotation(tag_levels = "A", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 18))
  ggsave(file.path(out_dir, paste0("Composite_", tis, ".png")), comp,
         width = 15, height = 17, dpi = 300, device = agg_png)
  ggsave(file.path(out_dir, paste0("Composite_", tis, ".tif")), comp,
         width = 15, height = 17, dpi = 300, device = "tiff", compression = "lzw")

  write.csv(data.frame(sample = rownames(m), label = m$label, condition = m$condition, tissue = tis),
            file.path(out_dir, "sample_key.csv"), row.names = FALSE)
  cat(sprintf("%s: %d genes tested | %d DEGs (%d up, %d down) | PC1 %d%%, PC2 %d%%\n",
              tis, nrow(df), n_up + n_dn, n_up, n_dn, pv[1], pv[2]))
}

run_tissue("Root",  "Fig4_Root_main")
run_tissue("Shoot", "FigS_Shoot_supplementary")
