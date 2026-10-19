suppressPackageStartupMessages({ library(GO.db); library(AnnotationDbi); library(ragg); library(patchwork) })

root <- commandArgs(TRUE)[1]
source(file.path(root, "figure_theme.R"))

min_size <- 10; max_size <- 500; padj_cut <- 0.05; n_show <- 10; jac_cut <- 0.7

g2g <- read.delim(file.path(root, "wheat_gene2go.tsv"), stringsAsFactors = FALSE)
ont <- suppressMessages(AnnotationDbi::select(GO.db, unique(g2g$go), "ONTOLOGY", "GOID"))
g2g$ont <- ont$ONTOLOGY[match(g2g$go, ont$GOID)]
g2g <- g2g[!is.na(g2g$ont), ]

anc <- c(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))
names(anc) <- sub("^(BP|CC|MF)\\.", "", names(anc))
roots_go <- c("GO:0008150", "GO:0005575", "GO:0003674", "all")

# gene -> propagated GO sets, per ontology
propagate <- function(o) {
  x <- g2g[g2g$ont == o, ]
  sets <- split(x$go, x$gene)
  lapply(sets, function(t) setdiff(unique(c(t, unlist(anc[t], use.names = FALSE))), roots_go))
}
gene_sets <- list(BP = propagate("BP"), CC = propagate("CC"), MF = propagate("MF"))
term_name <- function(ids) suppressMessages(AnnotationDbi::select(GO.db, ids, "TERM", "GOID"))$TERM

enrich <- function(deg, universe, o) {
  gs <- gene_sets[[o]][intersect(names(gene_sets[[o]]), universe)]
  bg <- names(gs); hits <- intersect(deg, bg)
  long <- data.frame(gene = rep(names(gs), lengths(gs)), go = unlist(gs, use.names = FALSE))
  t2g <- split(long$gene, long$go)
  t2g <- t2g[lengths(t2g) >= min_size & lengths(t2g) <= max_size]
  N <- length(bg); n <- length(hits)
  res <- do.call(rbind, lapply(names(t2g), function(t) {
    k <- sum(t2g[[t]] %in% hits); M <- length(t2g[[t]])
    data.frame(GO = t, Count = k, SetSize = M,
               GeneRatio = k / n, FoldEnrichment = (k / n) / (M / N),
               pvalue = phyper(k - 1, M, N - M, n, lower.tail = FALSE))
  }))
  res$padj <- p.adjust(res$pvalue, "BH")
  res <- res[res$padj < padj_cut & res$Count > 0, ]
  res <- res[order(res$padj), ]
  res$Ontology <- o
  res$Term <- term_name(res$GO)
  res$Genes <- vapply(res$GO, function(t) paste(intersect(t2g[[t]], hits), collapse = ";"), "")
  attr(res, "sizes") <- c(DEGs_annotated = n, background_annotated = N)
  res
}

reduce_terms <- function(res) {
  keep <- character(0); kept_sets <- list()
  for (i in seq_len(nrow(res))) {
    g <- strsplit(res$Genes[i], ";")[[1]]
    red <- any(vapply(kept_sets, function(k) length(intersect(g, k)) / length(union(g, k)) >= jac_cut, TRUE))
    if (!red) { keep <- c(keep, res$GO[i]); kept_sets[[length(kept_sets) + 1]] <- g }
  }
  res[res$GO %in% keep, ]
}

dot_panel <- function(res, title, high_col) {
  ont_lab <- c(BP = "BP", CC = "CC", MF = "MF")
  d <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) head(res[res$Ontology == o, ], n_show)))
  d$Ontology <- factor(ont_lab[d$Ontology], levels = ont_lab)
  trm <- ifelse(nchar(d$Term) > 95, paste0(substr(d$Term, 1, 92), "..."), d$Term)
  d$Label <- vapply(strwrap(trm, 50, simplify = FALSE), paste, "", collapse = "\n")
  d$Label <- factor(paste(d$Ontology, d$Label, sep = "__"),
                    levels = rev(paste(d$Ontology, d$Label, sep = "__")))
  d$nlp <- -log10(d$padj)
  ggplot(d, aes(FoldEnrichment, Label, size = Count, colour = nlp)) +
    geom_point() +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
    scale_colour_gradient(low = "#E3E3E3", high = high_col, name = "-log10(FDR)") +
    scale_size_continuous(range = c(2.5, 8), name = "Gene count") +
    labs(title = title, x = "Fold enrichment", y = NULL) +
    theme_pub +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          strip.text.y = element_text(angle = 0, face = "bold", size = 15),
          strip.background = element_rect(fill = "grey92"),
          axis.text.y = element_text(size = 13.5, colour = "grey10"),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12.5))
}

run <- function(tis, folder) {
  out <- file.path(root, folder); dir.create(out, showWarnings = FALSE)
  universe <- readLines(file.path(root, paste0(tis, "_background.txt")))
  lists <- list(Up = readLines(file.path(root, paste0(tis, "_DEGs_up.txt"))),
                Down = readLines(file.path(root, paste0(tis, "_DEGs_down.txt"))))
  plots <- list()
  for (dir in names(lists)) {
    full <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) enrich(lists[[dir]], universe, o)))
    red <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) reduce_terms(full[full$Ontology == o, ])))
    cols <- c("Ontology", "GO", "Term", "Count", "SetSize", "GeneRatio", "FoldEnrichment", "pvalue", "padj", "Genes")
    write.csv(full[, cols], file.path(out, paste0(tis, "_GO_", dir, "_all_significant.csv")), row.names = FALSE)
    write.csv(red[, cols], file.path(out, paste0(tis, "_GO_", dir, "_nonredundant.csv")), row.names = FALSE)
    n_ann <- sum(lists[[dir]] %in% unique(g2g$gene))
    cat(sprintf("%s %s: %d DEGs (%d with GO) | significant terms BP %d, CC %d, MF %d | non-redundant %d\n",
                tis, dir, length(lists[[dir]]), n_ann, sum(full$Ontology == "BP"),
                sum(full$Ontology == "CC"), sum(full$Ontology == "MF"), nrow(red)))
    col <- if (dir == "Up") pal$up else pal$down
    plots[[dir]] <- dot_panel(red, paste0(if (dir == "Up") "Upregulated" else "Downregulated", " genes"), col)
    ggsave(file.path(out, paste0(tis, "_GO_", dir, ".png")), plots[[dir]],
           width = 10, height = 10.5, dpi = 600, device = agg_png)
  }
  comp <- (plots$Down | plots$Up) +
    plot_annotation(tag_levels = "A", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 22))
  ggsave(file.path(out, paste0("Composite_GO_", tis, ".png")), comp,
         width = 20, height = 11, dpi = 300, device = agg_png)
  ggsave(file.path(out, paste0("Composite_GO_", tis, ".tif")), comp,
         width = 20, height = 11, dpi = 300, device = "tiff", compression = "lzw")
}

run("Root",  "Fig5_GO_Root_main")
run("Shoot", "FigS_GO_Shoot_supplementary")
