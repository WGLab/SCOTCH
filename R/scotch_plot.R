#!/usr/bin/env Rscript
# SCOTCH Visualization: Coverage + Sashimi + Gene Model
#
# Standalone ggplot2 script for gene-level visualization.
# Reads intermediate TSV files produced by src/visualization.py.
#
# Usage:
#   Rscript R/scotch_plot.R \
#     --coverage  vis/GENE/coverage.tsv \
#     --junctions vis/GENE/junctions.tsv \
#     --gene_model vis/GENE/gene_model.tsv \
#     --metadata  vis/GENE/metadata.tsv \
#     --output    vis/GENE/GENE_plot.pdf \
#     [--format pdf] [--width 12] [--height 10] \
#     [--annotation_scale 0.25] [--overlay] [--separate]

# ============================================================
# Dependencies
# ============================================================
required_pkgs <- c("ggplot2", "cowplot", "data.table", "optparse")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    "\nInstall with: install.packages(c('",
    paste(missing, collapse = "', '"), "'))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(optparse)
})

# ============================================================
# Argument parsing
# ============================================================
option_list <- list(
  make_option("--coverage", type = "character", help = "Path to coverage.tsv"),
  make_option("--junctions", type = "character", help = "Path to junctions.tsv"),
  make_option("--gene_model", type = "character", help = "Path to gene_model.tsv"),
  make_option("--metadata", type = "character", help = "Path to metadata.tsv"),
  make_option("--output", type = "character", help = "Output file path"),
  make_option("--format", type = "character", default = "pdf",
              help = "Output format: pdf or png [default %default]"),
  make_option("--width", type = "double", default = 12,
              help = "Plot width in inches [default %default]"),
  make_option("--height", type = "double", default = 10,
              help = "Plot height in inches [default %default]"),
  make_option("--annotation_scale", type = "double", default = 0.25,
              help = "Relative height of gene model panel [default %default]"),
  make_option("--overlay", action = "store_true", default = FALSE,
              help = "Overlay known/novel in same track with colors"),
  make_option("--separate", action = "store_true", default = FALSE,
              help = "Separate known/novel into distinct tracks")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================
# Read data
# ============================================================
meta <- fread(opt$metadata)
gene_name  <- meta$gene_name[1]
gene_chr   <- meta$gene_chr[1]
gene_start <- meta$gene_start[1]
gene_end   <- meta$gene_end[1]
gene_strand <- meta$gene_strand[1]

coverage_dt  <- fread(opt$coverage)
junctions_dt <- fread(opt$junctions)
gm_dt        <- fread(opt$gene_model)

# ============================================================
# Color palette
# ============================================================
COL_KNOWN <- "#2166ac"
COL_NOVEL <- "#b2182b"
COL_OTHER <- "#888888"

status_colors <- c("known" = COL_KNOWN, "novel" = COL_NOVEL, "other" = COL_OTHER, "all" = COL_KNOWN)

# Auto-generate group colors when not in overlay mode
make_group_palette <- function(groups) {
  n <- length(groups)
  if (n <= 1) return(setNames(COL_KNOWN, groups))
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  setNames(hcl(h = hues, l = 55, c = 80), groups)
}

# ============================================================
# Coverage panel
# ============================================================
build_coverage_plot <- function(dt, overlay) {
  if (nrow(dt) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No coverage data") +
        theme_void()
    )
  }

  if (overlay && "status" %in% names(dt)) {
    groups <- unique(dt$group)
    p <- ggplot(dt, aes(x = pos, y = depth, fill = status)) +
      geom_area(alpha = 0.6, position = "identity") +
      scale_fill_manual(values = status_colors, name = "Mapping") +
      facet_wrap(~ group, ncol = 1, scales = "free_y",
                 strip.position = "left")
  } else {
    groups <- unique(dt$group)
    pal <- make_group_palette(groups)
    p <- ggplot(dt, aes(x = pos, y = depth, fill = group)) +
      geom_area(alpha = 0.6) +
      scale_fill_manual(values = pal, name = "Group") +
      facet_wrap(~ group, ncol = 1, scales = "free_y",
                 strip.position = "left")
  }

  p +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    coord_cartesian(xlim = c(gene_start, gene_end)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, size = 8),
      strip.placement = "outside",
      legend.position = "bottom",
      plot.margin = margin(6, 5, 2, 5)
    ) +
    labs(y = "Read depth")
}

# ============================================================
# Junction / sashimi panel
# ============================================================
build_junction_plot <- function(dt, cov_dt) {
  if (nrow(dt) == 0) {
    return(NULL)
  }

  # Compute arc midpoint and height
  dt[, mid := (start + end) / 2]
  max_count <- max(dt$count, na.rm = TRUE)
  dt[, arc_height := sqrt(count / max_count) * 0.45 + 0.05]
  dt[, linetype_val := ifelse(is_annotated, "annotated", "novel")]

  # For faceting alignment with coverage, keep same groups
  groups <- unique(dt$group)

  p <- ggplot(dt) +
    geom_curve(
      aes(x = start, xend = end,
          y = 0, yend = 0,
          linewidth = count,
          linetype = linetype_val),
      curvature = -0.4, ncp = 20, color = "grey30"
    ) +
    geom_text(
      aes(x = mid, y = arc_height * 0.6, label = count),
      size = 2.5, color = "grey20"
    ) +
    scale_linetype_manual(
      values = c("annotated" = "solid", "novel" = "dashed"),
      name = "Junction"
    ) +
    scale_linewidth_continuous(range = c(0.3, 2), guide = "none") +
    facet_wrap(~ group, ncol = 1, scales = "free_y",
               strip.position = "left") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(gene_start, gene_end), ylim = c(-0.1, 0.6)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text.y.left = element_text(angle = 0, size = 8),
      strip.placement = "outside",
      legend.position = "bottom",
      plot.margin = margin(2, 5, 2, 5)
    )

  p
}

# ============================================================
# Gene model panel
# ============================================================
build_gene_model_plot <- function(dt) {
  if (nrow(dt) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No gene model data") +
        theme_void()
    )
  }

  transcripts <- unique(dt$transcript_id)
  dt[, y_pos := match(transcript_id, transcripts)]

  exons <- dt[feature == "exon"]
  introns <- dt[feature == "intron"]

  # Isoform coloring
  dt[, color_val := ifelse(is_novel, COL_NOVEL, COL_KNOWN)]
  tx_colors <- setNames(dt[!duplicated(transcript_id)]$color_val,
                         dt[!duplicated(transcript_id)]$transcript_id)

  # Format labels
  format_label <- function(x) {
    ifelse(
      grepl("^novelIsoform_", x),
      paste0("Novel\n(", sub("^novelIsoform_", "", x), ")"),
      x
    )
  }

  y_labels <- setNames(sapply(transcripts, format_label), seq_along(transcripts))

  # Strand arrow
  strand_val <- dt$strand[1]
  arrow_dir <- if (strand_val == "+") "last" else "first"

  p <- ggplot() +
    # Intron lines with strand arrows
    geom_segment(
      data = introns,
      aes(x = start, xend = end, y = y_pos, yend = y_pos,
          color = transcript_id),
      linewidth = 0.4,
      arrow = arrow(length = unit(0.06, "inches"), ends = arrow_dir,
                    type = "open")
    ) +
    # Exon rectangles
    geom_rect(
      data = exons,
      aes(xmin = start, xmax = end,
          ymin = y_pos - 0.3, ymax = y_pos + 0.3,
          fill = transcript_id)
    ) +
    scale_fill_manual(values = tx_colors, guide = "none") +
    scale_color_manual(values = tx_colors, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq_along(transcripts),
      labels = y_labels,
      expand = expansion(mult = c(0.15, 0.15))
    ) +
    coord_cartesian(xlim = c(gene_start, gene_end)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 7, hjust = 1),
      axis.text.x = element_text(size = 8),
      plot.margin = margin(2, 5, 6, 5)
    ) +
    labs(
      x = paste0(
        gene_name, ": ", gene_chr, " ",
        format(gene_start, big.mark = ",", scientific = FALSE),
        " - ",
        format(gene_end, big.mark = ",", scientific = FALSE),
        " (", gene_strand, ")"
      )
    )

  p
}

# ============================================================
# Compose and save
# ============================================================

p_cov  <- build_coverage_plot(coverage_dt, overlay = opt$overlay)
p_junc <- build_junction_plot(junctions_dt, coverage_dt)
p_gene <- build_gene_model_plot(gm_dt)

# Calculate relative panel heights
ann_scale <- opt$annotation_scale
if (!is.null(p_junc)) {
  junc_scale <- 0.15
  cov_scale <- 1 - ann_scale - junc_scale
  panels <- plot_grid(
    p_cov  + theme(legend.position = "none"),
    p_junc + theme(legend.position = "none"),
    p_gene,
    ncol = 1,
    rel_heights = c(cov_scale, junc_scale, ann_scale),
    align = "v",
    axis = "lr"
  )
} else {
  cov_scale <- 1 - ann_scale
  panels <- plot_grid(
    p_cov + theme(legend.position = "none"),
    p_gene,
    ncol = 1,
    rel_heights = c(cov_scale, ann_scale),
    align = "v",
    axis = "lr"
  )
}

# Extract legends
cov_legend <- get_legend(p_cov)
junc_legend <- if (!is.null(p_junc)) get_legend(p_junc) else NULL

legends <- if (!is.null(junc_legend)) {
  plot_grid(cov_legend, junc_legend, nrow = 1)
} else {
  cov_legend
}

final <- plot_grid(panels, legends, ncol = 1, rel_heights = c(0.92, 0.08))

# Save
out_path <- opt$output
if (opt$format == "png") {
  ggsave(out_path, final, width = opt$width, height = opt$height, dpi = 300)
} else {
  ggsave(out_path, final, width = opt$width, height = opt$height, device = "pdf")
}

cat("Plot saved to:", out_path, "\n")
