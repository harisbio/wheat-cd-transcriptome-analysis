suppressPackageStartupMessages({
  library(DESeq2); library(igraph); library(ggraph); library(patchwork); library(ragg); library(png); library(grid)
})
args <- commandArgs(TRUE); src <- args[1]; root <- args[2]; tis <- args[3]; folder <- args[4]; gfolder <- args[5]
source(file.path(root, "figure_theme.R")); source(file.path(root, "fig_helpers.R"))
ann <- read.delim(file.path(root, "wheat_gene_annotation.tsv"), stringsAsFactors = FALSE, quote = "")
short_name <- make_short_name(ann)
name_only <- function(g) sub(" \\|.*$", "", short_name(g))
n_hub <- 20; n_lab <- 10; min_comp <- 3
out <- file.path(root, folder); dir.create(out, showWarnings = FALSE)
pdir <- file.path(root, "ppi", tis)

degs <- read.csv(file.path(root, paste0(tis, "_DEGs.csv")), row.names = 1)
g2s <- read.delim(file.path(pdir, "gene2string.tsv"), stringsAsFactors = FALSE)
s2g <- g2s$gene[!duplicated(g2s$string_id)]; names(s2g) <- g2s$string_id[!duplicated(g2s$string_id)]
groups <- c(`Oxidative stress` = "Oxidative_stress", `Transcription factors` = "Transcription_factors",
            Transporters = "Transporters")
group_genes <- lapply(groups, function(k) read.csv(file.path(root, gfolder, paste0(tis, "_", k, "_DEGs.csv")))$gene)

as_graph <- function(e) {
  g <- graph_from_data_frame(data.frame(from = s2g[e$stringId_A], to = s2g[e$stringId_B], score = e$score),
                             directed = FALSE)
  g <- simplify(g, edge.attr.comb = list(score = "max"))
  V(g)$log2FC <- degs[V(g)$name, "log2FC"]
  V(g)$degree <- as.numeric(igraph::degree(g))
  V(g)$label <- name_only(V(g)$name)
  g
}
read_edges <- function(name, min_score = 0.4) {
  e <- read.delim(file.path(pdir, paste0("network_", name, ".tsv")), stringsAsFactors = FALSE)
  e[e$score >= min_score, ]
}
mcc_rank <- function(G) {
  lmcc <- setNames(rep(-Inf, vcount(G)), V(G)$name)
  for (c in max_cliques(G, min = 2)) {
    v <- names(c); w <- lfactorial(length(c) - 1)
    lmcc[v] <- pmax(lmcc[v], w) + log1p(exp(-abs(lmcc[v] - w)))   # log-sum-exp accumulation
  }
  h <- data.frame(gene = V(G)$name, protein = ann$protein_name[match(V(G)$name, ann$gene)],
                  log10_MCC = lmcc / log(10), degree = V(G)$degree, log2FC = V(G)$log2FC,
                  padj = degs[V(G)$name, "padj"])
  h <- h[order(-h$log10_MCC, -h$degree), ]
  h$rank <- seq_len(nrow(h)); h$string_id <- g2s$string_id[match(h$gene, g2s$gene)]
  h
}

lfc_scale <- function(lim) scale_fill_gradient2(low = pal$down, mid = pal$mid, high = pal$up, midpoint = 0,
                                                limits = c(-lim, lim), oob = scales::squish, name = "log2FC")
net_plot <- function(g, title) {
  n_all <- sum(igraph::degree(g) > 0)
  comp <- components(g)
  g <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])   # largest connected cluster
  V(g)$degree <- as.numeric(igraph::degree(g))
  lab <- V(g)$name[order(-V(g)$degree)][seq_len(min(n_lab, vcount(g)))]
  V(g)$show <- ifelse(V(g)$name %in% lab, V(g)$label, "")
  set.seed(42)
  ggraph(g, layout = "stress") +
    geom_edge_link(aes(edge_alpha = score), colour = "grey55", edge_width = 0.5, show.legend = FALSE) +
    scale_edge_alpha(range = c(0.15, 0.6)) +
    geom_node_point(aes(fill = log2FC, size = degree), shape = 21, colour = "grey25", stroke = 0.3) +
    geom_node_text(aes(label = show), repel = TRUE, size = 5.2, fontface = "bold", colour = "grey10",
                   bg.colour = "white", bg.r = 0.12, max.overlaps = Inf, box.padding = 0.5,
                   point.padding = 0.4, min.segment.length = 0) +
    lfc_scale(6) + scale_size_continuous(range = c(2.5, 11), name = "Degree") +
    labs(title = title, subtitle = sprintf("largest cluster: %d of %d connected proteins, %d interactions", vcount(g), n_all, ecount(g))) +
    theme_void() + theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
                         plot.subtitle = element_text(size = 13, hjust = 0.5, colour = "grey30"),
                         legend.position = "none", plot.margin = margin(2, 2, 2, 2))
}

fam_cols <- c("#4C9F70", "#E0B24A", "#5B4E77", "#2A9D8F", "#D95F5F", "#8DBF7A", "#C77CC2", "#7FA7D6",
              "#B07D4F", "#E58FA8", "#9AA0A6")

transporter_family <- function(genes) {
  nm <- ann$protein_name[match(genes, ann$gene)]; nm[is.na(nm)] <- "Uncharacterized"
  rules <- c(
    "ABC transporter|MDR-like" = "ABC transporters",
    "Major facilitator|MFS|Sucrose transport|UDP-galactose|Tonoplast dicarboxylate" = "MFS & sugar/organic acid",
    "Amino acid|amino acid|Polyamine|Peptide|Oligopeptide|Urea" = "Amino acid & peptide",
    "Ammonium|Pi cotransporter|phosphate transporter|Molybdate" = "N/P nutrient transporters",
    "Copper|Cu\\(\\+\\)|Zinc|Metal transporter|HMA domain|iron transporter" = "Metal transporters",
    "Ca\\(2\\+\\)|Calcium-transporting|Phospholipid-transporting|H\\(\\+\\)-exporting|V-type proton|ATP synthase" = "P/V-type ATPases & pumps",
    "Potassium|Sodium/hydrogen|Chloride channel|anion channel|RCK|glutamate receptor|EF-hand|Ethylene-insensitive" = "Ion channels & antiporters",
    "intrinsic protein|Aquaporin" = "Aquaporins",
    "DETOXIFICATION" = "MATE",
    "^Uncharacterized" = "Uncharacterized")
  fam <- rep("Other", length(nm))
  for (p in rev(names(rules))) fam[grepl(p, nm)] <- rules[[p]]
  fam
}

circ_network <- function(g, title, file, n_fam = 10, family = NULL) {
  g <- induced_subgraph(g, V(g)[igraph::degree(g) > 0])
  fam <- if (is.null(family)) sub(" \\((C|N)-term\\)$", "", name_only(V(g)$name)) else family(V(g)$name)
  top_f <- names(sort(table(fam[!fam %in% c("Uncharacterized", "Other")]), decreasing = TRUE))[seq_len(min(n_fam, length(unique(fam))))]
  top_f <- top_f[!is.na(top_f)]
  fam[!fam %in% c(top_f, "Uncharacterized")] <- "Other"
  lev <- c(top_f, intersect(c("Other", "Uncharacterized"), fam))
  fcol <- setNames(c(fam_cols[seq_along(top_f)], "#B8B8B8", "#D9D9D9")[seq_along(lev)], lev)
  if ("Other" %in% lev) fcol["Other"] <- "#B8B8B8"; if ("Uncharacterized" %in% lev) fcol["Uncharacterized"] <- "#D9D9D9"
  d <- data.frame(name = V(g)$name, fam = factor(fam, levels = lev), deg = V(g)$degree, lfc = V(g)$log2FC)
  d <- d[order(d$fam, -d$deg), ]
  d$x <- ave(seq_len(nrow(d)), d$fam, FUN = seq_along) - 0.5
  sz <- table(d$fam)
  el <- as_edgelist(g); sc <- E(g)$score
  col_fun <- colorRamp2(c(-6, 0, 6), c(pal$down, pal$mid, pal$up))
  rad <- function(k) 0.9 + 2.6 * sqrt(k / max(d$deg))                         # node size by degree (mm)

  agg_png(file, width = 9, height = 12.4, units = "in", res = 400)
  layout(matrix(1:2, 2), heights = c(9, 3.4))
  par(mar = c(0, 0, 3.6, 0))
  circos.clear()
  circos.par(start.degree = 90, gap.after = 2.5, cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.004),
             points.overflow.warning = FALSE, canvas.xlim = c(-1.04, 1.04), canvas.ylim = c(-1.04, 1.04))
  circos.initialize(d$fam, xlim = cbind(0, as.numeric(sz)))
  # family labels
  circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = NA, panel.fun = function(x, y) {
    s <- CELL_META$sector.index; n <- sz[[s]]
    lab <- if (n / nrow(d) >= 0.07 && nchar(s) <= 18) paste0(s, " (", n, ")") else ""
    circos.text(CELL_META$xcenter, 0.05, lab, facing = "bending.inside", niceFacing = TRUE,
                adj = c(0.5, 0), cex = 1.6, font = 2,
                col = adjustcolor(fcol[[s]], red.f = 0.75, green.f = 0.75, blue.f = 0.75))
  })
  # family band
  circos.track(ylim = c(0, 1), track.height = 0.035, bg.border = NA, panel.fun = function(x, y)
    circos.rect(0, 0, CELL_META$xlim[2], 1, col = fcol[[CELL_META$sector.index]], border = NA))
  # node dots
  circos.track(ylim = c(0, 1), track.height = 0.06, bg.border = NA, panel.fun = function(x, y) {
    dd <- d[d$fam == CELL_META$sector.index, ]
    circos.points(dd$x, rep(0.5, nrow(dd)), pch = 21, bg = col_fun(dd$lfc), col = "grey25",
                  cex = rad(dd$deg) / 1.6)
  })
  # links (drawn last, inside), coloured by the family of the higher-degree partner
  pos <- d[match(c(el[, 1], el[, 2]), d$name), ]
  a <- pos[seq_len(nrow(el)), ]; b <- pos[nrow(el) + seq_len(nrow(el)), ]
  lead <- ifelse(a$deg >= b$deg, as.character(a$fam), as.character(b$fam))
  o <- order(sc)
  for (i in o) circos.link(as.character(a$fam[i]), a$x[i], as.character(b$fam[i]), b$x[i],
                           h.ratio = 0.55, lwd = 0.4 + 1.2 * (sc[i] - 0.4) / 0.6,
                           col = adjustcolor(fcol[[lead[i]]], 0.18 + 0.35 * (sc[i] - 0.4) / 0.6))
  circos.clear()
  title(title, cex.main = 2.0, line = 2.4)
  mtext(sprintf("%d connected proteins, %d interactions (STRING score \u2265 0.4)", nrow(d), nrow(el)),
        side = 3, line = 0.6, cex = 1.6, col = "grey30")
  par(mar = c(0, 1, 0, 1)); plot.new()
  legend("top", legend = paste0(names(sz), " (", as.numeric(sz), ")"), fill = fcol[names(sz)], border = NA,
         bty = "n", ncol = 2, cex = 1.55, x.intersp = 0.6)
  dev.off()
  invisible(d)
}

pg <- lapply(names(groups), function(k) {
  g <- as_graph(read_edges(groups[[k]]))
  nd <- data.frame(gene = V(g)$name, protein = name_only(V(g)$name), degree = V(g)$degree, log2FC = V(g)$log2FC)
  write.csv(nd[order(-nd$degree), ], file.path(out, paste0(tis, "_PPI_", groups[[k]], "_nodes.csv")), row.names = FALSE)
  f <- file.path(out, paste0(tis, "_PPI_", groups[[k]], ".png"))
  circ_network(g, k, f, family = if (k == "Transporters") transporter_family else NULL)
  wrap_elements(full = rasterGrob(readPNG(f, native = TRUE), interpolate = TRUE))
})

# ---- hub genes ----
e_all <- read_edges("all_DEGs", 0.4)
# supplementary: high-confidence network of all DEGs
G_all <- as_graph(e_all[e_all$score >= 0.7, ])
hub_all <- mcc_rank(G_all)
write.csv(head(hub_all, 200), file.path(out, paste0(tis, "_hub_genes_MCC_allDEGs_score0.7_top200.csv")), row.names = FALSE)
# main: stress-response network (union of the three groups)
stress <- unique(unlist(group_genes))
e_s <- e_all[s2g[e_all$stringId_A] %in% stress & s2g[e_all$stringId_B] %in% stress, ]
G <- as_graph(e_s)
cat(sprintf("%s stress-response network: %d proteins, %d interactions | all-DEG network (>=0.7): %d proteins, %d interactions\n",
            tis, vcount(G), ecount(G), vcount(G_all), ecount(G_all)))
hub <- mcc_rank(G)
hub$group <- vapply(hub$gene, function(x) paste(names(groups)[vapply(group_genes, function(s) x %in% s, TRUE)], collapse = "; "), "")
write.csv(hub, file.path(out, paste0(tis, "_hub_genes_MCC_stress_network.csv")), row.names = FALSE)
top <- head(hub, n_hub)
print(top[, c("rank", "gene", "protein", "degree", "log2FC")], row.names = FALSE)

# D: hub network, nodes on a circle in rank order, labels outside
H <- induced_subgraph(G, top$gene)
ord <- match(top$gene, V(H)$name)
ang <- pi / 2 - 2 * pi * (seq_len(n_hub) - 0.5) / n_hub   # half-step offset: no node at the very top/bottom
xy <- matrix(NA, vcount(H), 2); xy[ord, ] <- cbind(cos(ang), sin(ang))
V(H)$rank <- top$rank[match(V(H)$name, top$gene)]
V(H)$lab <- paste0(V(H)$rank, ". ", name_only(V(H)$name))
V(H)$lx <- 1.16 * xy[, 1]; V(H)$ly <- 1.16 * xy[, 2] + ifelse(abs(xy[, 1]) < 0.35, sign(xy[, 2]) * 0.05, 0)
V(H)$hj <- ifelse(xy[, 1] > 0.05, 0, ifelse(xy[, 1] < -0.05, 1, 0.5))
pD <- ggraph(H, layout = "manual", x = xy[, 1], y = xy[, 2]) +
  geom_edge_arc(strength = 0.18, colour = "grey50", edge_width = 0.45, edge_alpha = 0.55) +
  geom_node_point(aes(fill = log2FC, size = 21 - rank), shape = 21, colour = "grey20", stroke = 0.5) +
  geom_node_text(aes(x = lx, y = ly, label = lab, hjust = hj), size = 6.4, fontface = "bold", colour = "grey10") +
  lfc_scale(6) + scale_size_continuous(range = c(5, 13), guide = "none") +
  coord_fixed(clip = "off", xlim = c(-1.9, 1.9), ylim = c(-1.35, 1.35)) +
  labs(title = sprintf("Top %d hub genes of the stress-response network (MCC)", n_hub)) +
  theme_void() + theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                       legend.position = "none", plot.margin = margin(10, 10, 10, 10))

# ---- E: hub-gene expression (circular heatmap) ----
counts <- read.csv(file.path(src, "gene_counts_integer.csv"), row.names = 1, check.names = FALSE)
meta <- data.frame(row.names = colnames(counts),
  condition = factor(c("Control","Control","Control","Control","Cd","Cd","Cd","Cd","Cd","Cd","Control","Control"),
                     levels = c("Control", "Cd")),
  tissue = c("Shoot","Root","Root","Root","Shoot","Shoot","Shoot","Root","Root","Root","Shoot","Shoot"),
  rep = c("1","3","2","1","3","2","1","3","2","1","3","2"))
m <- meta[meta$tissue == tis, ]; m <- m[order(m$condition, m$rep), ]
m$label <- paste0(ifelse(m$condition == "Control", "Ctrl", "Cd"), "-", m$rep)
dds <- DESeqDataSetFromMatrix(counts[, rownames(m)], m, ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- assay(vst(dds, blind = TRUE)); colnames(vsd) <- m$label
fE <- file.path(out, paste0(tis, "_hub_genes_circular_heatmap.png"))
circ_heatmap(vsd[top$gene, ], ifelse(top$log2FC > 0, "Up", "Down"), m$condition,
             paste0(top$rank, ". ", name_only(top$gene)), fE, "Hub-gene expression", legend = FALSE,
             lab_cex = 1.45, samp_cex = 0.85, canvas = 1.75)

# ---- save panels and composite ----
ggsave(file.path(out, paste0(tis, "_hub_network.png")), pD, width = 10, height = 8, dpi = 400, device = agg_png)

lg_p <- ggplot(data.frame(x = 1:3, y = 1, l = c(-6, 0, 6), d = c(10, 50, 100)), aes(x, y)) +
  geom_point(aes(fill = l, size = d), shape = 21) + lfc_scale(6) +
  scale_size_continuous(range = c(2, 7), name = "Degree", breaks = c(10, 50, 100)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(6, "cm"))) +
  theme_void() + theme(legend.position = "bottom", legend.title = element_text(size = 20, face = "bold"),
                       legend.text = element_text(size = 18), legend.spacing.x = unit(1, "cm"))
lg1 <- wrap_elements(full = cowplot::get_legend(lg_p))
lg2 <- circ_legend_strip(20, 18)
pE <- wrap_elements(full = rasterGrob(readPNG(fE, native = TRUE), interpolate = TRUE))
comp <- (pg[[1]] | pg[[2]] | pg[[3]]) / lg1 / (pD | pE) / lg2 +
  plot_layout(heights = c(11, 0.5, 9.5, 0.5)) +
  plot_annotation(tag_levels = list(c("A)", "B)", "C)", "", "D)", "E)", ""))) &
  theme(plot.tag = element_text(face = "bold", size = 20))
ggsave(file.path(out, paste0("Composite_PPI_", tis, ".png")), comp, width = 24, height = 22.5, dpi = 300, device = agg_png)
ggsave(file.path(out, paste0("Composite_PPI_", tis, ".tif")), comp, width = 24, height = 22.5, dpi = 300,
       device = "tiff", compression = "lzw")
cat("done\n")
