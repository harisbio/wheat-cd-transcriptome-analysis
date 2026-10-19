

suppressPackageStartupMessages({
  library(DESeq2); library(GO.db); library(circlize); library(patchwork); library(ragg)
  library(reshape2); library(png); library(grid)
})
args <- commandArgs(TRUE); src <- args[1]; root <- args[2]
source(file.path(root, "figure_theme.R"))
n_heat <- 30

groups <- list(
  `Oxidative stress` = c("GO:0016209", "GO:0006979", "GO:0004364"),
  `Transcription factors` = "GO:0003700",
  `Transporters` = "GO:0005215"
)
offspring <- c(as.list(GOMFOFFSPRING), as.list(GOBPOFFSPRING))
expand <- function(ids) unique(c(ids, unlist(offspring[ids], use.names = FALSE)))
g2g <- read.delim(file.path(root, "wheat_gene2go.tsv"), stringsAsFactors = FALSE)
ann <- read.delim(file.path(root, "wheat_gene_annotation.tsv"), stringsAsFactors = FALSE, quote = "")
group_genes <- lapply(groups, function(ids) unique(g2g$gene[g2g$go %in% expand(ids)]))

counts <- read.csv(file.path(src, "gene_counts_integer.csv"), row.names = 1, check.names = FALSE)
meta <- data.frame(
  row.names = colnames(counts),
  condition = c("Control","Control","Control","Control","Cd","Cd","Cd","Cd","Cd","Cd","Control","Control"),
  tissue    = c("Shoot","Root","Root","Root","Shoot","Shoot","Shoot","Root","Root","Root","Shoot","Shoot"),
  rep       = c("1","3","2","1","3","2","1","3","2","1","3","2"))
meta$condition <- factor(meta$condition, levels = c("Control","Cd"))

source(file.path(root, "fig_helpers.R"))
short_name <- make_short_name(ann)

run <- function(tis, folder) {
  out <- file.path(root, folder); dir.create(out, showWarnings = FALSE)
  m <- meta[meta$tissue == tis, ]; m <- m[order(m$condition, m$rep), ]
  m$label <- paste0(ifelse(m$condition == "Control", "Ctrl", "Cd"), "-", m$rep)
  dds <- DESeqDataSetFromMatrix(counts[, rownames(m)], m, ~ condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  vsd <- assay(vst(dds, blind = TRUE)); colnames(vsd) <- m$label
  degs <- read.csv(file.path(root, paste0(tis, "_DEGs.csv")), row.names = 1)

  heat <- list(); viol <- list(); summ <- list()
  for (k in names(groups)) {
    g <- intersect(group_genes[[k]], rownames(degs))
    d <- degs[g, ]; d$Direction <- ifelse(d$log2FC > 0, "Up", "Down")
    d$ProteinName <- ann$protein_name[match(g, ann$gene)]
    d <- d[order(d$padj), ]
    write.csv(data.frame(gene = rownames(d), d[, c("ProteinName", "baseMean", "log2FC", "padj", "Direction")]),
              file.path(out, paste0(tis, "_", gsub(" ", "_", k), "_DEGs.csv")), row.names = FALSE)
    nu <- sum(d$Direction == "Up"); nd <- sum(d$Direction == "Down")
    summ[[k]] <- data.frame(Tissue = tis, Group = k, Annotated_genes_tested = length(intersect(group_genes[[k]], rownames(vsd))),
                            DEGs = nrow(d), Up = nu, Down = nd, Percent_up = round(100 * nu / nrow(d), 1))

    # circular heatmap: top DEGs by padj
    top <- head(rownames(d), n_heat)
    fh <- file.path(out, paste0(tis, "_circular_heatmap_", gsub(" ", "_", k), ".png"))
    circ_heatmap(vsd[top, ], d[top, "Direction"], m$condition, short_name(top), fh,
                 sprintf("%s (top %d of %d DEGs)", k, length(top), nrow(d)), legend = FALSE)
    heat[[k]] <- wrap_elements(full = rasterGrob(readPNG(fh, native = TRUE), interpolate = TRUE))

    # violin: all group DEGs, per sample
    long <- melt(vsd[rownames(d), , drop = FALSE], varnames = c("gene", "Sample"), value.name = "VST")
    long$Sample <- factor(long$Sample, levels = m$label)
    long$Condition <- m$condition[match(long$Sample, m$label)]
    viol[[k]] <- ggplot(long, aes(Sample, VST, fill = Condition)) +
      geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
      geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "grey25") +
      scale_fill_manual(values = col_cond) +
      labs(title = sprintf("%s (n = %d; ↑%d / ↓%d)", k, nrow(d), nu, nd),
           x = NULL, y = "VST expression") +
      theme_pub + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                        axis.text = element_text(size = 12), axis.title.y = element_text(size = 13),
                        legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 12))
    ggsave(file.path(out, paste0(tis, "_violin_", gsub(" ", "_", k), ".png")), viol[[k]],
           width = 6, height = 5, dpi = 400, device = agg_png)
  }
  s <- do.call(rbind, summ); write.csv(s, file.path(out, paste0(tis, "_gene_group_summary.csv")), row.names = FALSE)
  print(s, row.names = FALSE)

  lg <- circ_legend_strip()
  comp <- (heat[[1]] | heat[[2]] | heat[[3]]) / lg / (viol[[1]] | viol[[2]] | viol[[3]]) +
    plot_layout(heights = c(1.55, 0.09, 1)) +
    plot_annotation(tag_levels = list(c("A)", "B)", "C)", "", "D)", "E)", "F)"))) &
    theme(plot.tag = element_text(face = "bold", size = 20))
  ggsave(file.path(out, paste0("Composite_GeneGroups_", tis, ".png")), comp, width = 24, height = 15.5, dpi = 300, device = agg_png)
  ggsave(file.path(out, paste0("Composite_GeneGroups_", tis, ".tif")), comp, width = 24, height = 15.5, dpi = 300,
         device = "tiff", compression = "lzw")
}

run("Root",  "Fig7_GeneGroups_Root_main")
run("Shoot", "FigS_GeneGroups_Shoot_supplementary")
