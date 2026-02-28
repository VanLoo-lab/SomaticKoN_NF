#!/usr/bin/env Rscript
# plot_somatic_qc.R
#
# QC plots for somatic variant pipeline
# Input: outputs from make_somatic_metrics.sh
#
# Robust behavior:
# - If optional packages are missing, plots that depend on them are skipped (script still runs)
# - If an input TSV is missing/empty, the corresponding plot is skipped

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(scales)
})

# Optional: UpSetR
has_upsetr <- requireNamespace("UpSetR", quietly = TRUE)
if (has_upsetr) suppressPackageStartupMessages(library(UpSetR))

# ----------------------------
# Args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
  cat("
Usage: plot_somatic_qc.R <metrics_dir> <out_dir> [options]

Required:
  metrics_dir    Directory containing TSV outputs from make_somatic_metrics.sh
  out_dir        Output directory for plots

Options:
  --sample-id ID     Sample name for plot titles (auto-detected if not provided)
  --format FMT       Output format: png, pdf, or both (default: png)
  --width W          Plot width in inches (default: 8)
  --height H         Plot height in inches (default: 5)
  --dpi D            DPI for PNG output (default: 300)

Expected input files (in metrics_dir):
  variant_counts_by_stage.tsv
  variant_counts_wide.tsv
  upset_matrix.tsv
  caller_agreement_matrix.tsv
  caller_contribution.tsv
  k_curve.tsv
  k_curve_by_vartype.tsv

Optional (new, if generated):
  k_curve_exact.tsv
  k_curve_exact_by_vartype.tsv
  k_curve_exact_by_caller.tsv

Genomic:
  genomic_bins.tsv
  genomic_summary.tsv
")
}

if (length(args) < 2) {
  print_usage()
  quit(status = 2)
}

metrics_dir <- args[1]
out_dir <- args[2]

# Optionals
sample_id <- NULL
out_format <- "png"
plot_width <- 8
plot_height <- 5
plot_dpi <- 300

i <- 3
while (i <= length(args)) {
  if (args[i] == "--sample-id" && i < length(args)) {
    sample_id <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--format" && i < length(args)) {
    out_format <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--width" && i < length(args)) {
    plot_width <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--height" && i < length(args)) {
    plot_height <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--dpi" && i < length(args)) {
    plot_dpi <- as.integer(args[i + 1]); i <- i + 2
  } else {
    cat("Unknown option:", args[i], "\n"); i <- i + 1
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
save_plot <- function(p, basename, width = plot_width, height = plot_height) {
  if (out_format %in% c("png", "both")) {
    ggsave(file.path(out_dir, paste0(basename, ".png")), p,
           width = width, height = height, dpi = plot_dpi)
  }
  if (out_format %in% c("pdf", "both")) {
    ggsave(file.path(out_dir, paste0(basename, ".pdf")), p,
           width = width, height = height)
  }
}

file_exists_nonempty <- function(path) {
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

read_tsv_safe <- function(path) {
  if (!file_exists_nonempty(path)) return(NULL)
  tryCatch(
    read_tsv(path, show_col_types = FALSE),
    error = function(e) { message("Error reading ", path, ": ", e$message); NULL }
  )
}

# Colors (you can tweak later)
caller_colors <- c(
  "mutect2"  = "#E41A1C",
  "muse2"    = "#377EB8",
  "strelka2" = "#4DAF4A",
  "vardict"  = "#984EA3",
  "lofreq"   = "#FF7F00",
  "varscan2" = "#A65628"
)

vartype_colors <- c(
  "SNP"          = "#4DAF4A",
  "INDEL"        = "#377EB8",
  "MNV"          = "#FF7F00",
  "MULTIALLELIC" = "#984EA3"
)

vartype_levels <- c("SNP", "INDEL", "MNV", "MULTIALLELIC")
stage_levels <- c("raw", "norm", "split", "consensus", "merged", "event")

# ----------------------------
# Theme
# ----------------------------
theme_set(
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey35"),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
)

cat("Reading metrics from:", metrics_dir, "\n")

# ----------------------------
# Load TSVs
# ----------------------------
counts_long       <- read_tsv_safe(file.path(metrics_dir, "variant_counts_by_stage.tsv"))
counts_wide       <- read_tsv_safe(file.path(metrics_dir, "variant_counts_wide.tsv"))
upset_matrix      <- read_tsv_safe(file.path(metrics_dir, "upset_matrix.tsv"))
caller_agreement  <- read_tsv_safe(file.path(metrics_dir, "caller_agreement_matrix.tsv"))
caller_contrib    <- read_tsv_safe(file.path(metrics_dir, "caller_contribution.tsv"))
kcurve            <- read_tsv_safe(file.path(metrics_dir, "k_curve.tsv"))
kcurve_vartype    <- read_tsv_safe(file.path(metrics_dir, "k_curve_by_vartype.tsv"))

# New optional metrics from your updated make_somatic_metrics.sh
kcurve_exact          <- read_tsv_safe(file.path(metrics_dir, "k_curve_exact.tsv"))
kcurve_exact_vartype  <- read_tsv_safe(file.path(metrics_dir, "k_curve_exact_by_vartype.tsv"))
kcurve_exact_caller   <- read_tsv_safe(file.path(metrics_dir, "k_curve_exact_by_caller.tsv"))

genome_bins     <- read_tsv_safe(file.path(metrics_dir, "genomic_bins.tsv"))
genome_summary  <- read_tsv_safe(file.path(metrics_dir, "genomic_summary.tsv"))
event_status_by_type <- read_tsv_safe(file.path(metrics_dir, "event_status_by_type.tsv"))


# Auto sample ID
if (is.null(sample_id)) {
  if (!is.null(counts_long) && "sample" %in% names(counts_long)) {
    sample_id <- unique(counts_long$sample)[1]
  } else if (!is.null(kcurve) && "sample" %in% names(kcurve)) {
    sample_id <- unique(kcurve$sample)[1]
  } else {
    sample_id <- "Sample"
  }
}
cat("Sample ID:", sample_id, "\n")

# Available stages
available_stages <- c()
if (!is.null(counts_long) && "stage" %in% names(counts_long)) {
  s <- unique(as.character(counts_long$stage))
  available_stages <- c(intersect(stage_levels, s), setdiff(s, stage_levels))
}
cat("Available stages:", paste(available_stages, collapse = ", "), "\n")

# Available callers
available_callers <- c()
if (!is.null(counts_long) && "caller" %in% names(counts_long)) {
  available_callers <- unique(as.character(counts_long$caller))
}
cat("Available callers:", paste(available_callers, collapse = ", "), "\n")

# ----------------------------
# Fig 1A: Variant counts by stage (stacked)
# ----------------------------
if (!is.null(counts_long) && nrow(counts_long) > 0) {
  cat("Generating: Fig1A variant counts by stage\n")

  df <- counts_long %>%
    mutate(
      stage = factor(stage, levels = available_stages),
      vartype = factor(vartype, levels = vartype_levels)
    ) %>%
    group_by(sample, stage, vartype) %>%
    summarize(count = sum(count, na.rm = TRUE), .groups = "drop")

  p <- ggplot(df, aes(x = stage, y = count, fill = vartype)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = vartype_colors, drop = FALSE) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = NULL, y = "Variant count", fill = "Type",
      title = "Variant counts across pipeline stages",
      subtitle = sample_id
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_plot(p, "variant_counts_by_stage", width = 8, height = 5)
}

# ----------------------------
# Fig 1B: Variant composition by caller (for a representative stage)
# + Adds your requested plot: proportion SNP/MNV/MULTIALLELIC/INDEL per caller
# ----------------------------
if (!is.null(counts_long) && nrow(counts_long) > 0) {
  cat("Generating: Fig1B variant counts & proportions by caller\n")

  rep_stage <- intersect(c("norm", "split", "raw", "consensus"), available_stages)[1]
  if (is.na(rep_stage) || is.null(rep_stage)) rep_stage <- available_stages[1]

  df <- counts_long %>%
    filter(stage == rep_stage) %>%
    mutate(
      vartype = factor(vartype, levels = vartype_levels),
      caller = factor(caller, levels = sort(unique(caller)))
    )

  # stacked counts
  p_counts <- ggplot(df, aes(x = caller, y = count, fill = vartype)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = vartype_colors, drop = FALSE) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = NULL, y = "Variant count", fill = "Type",
      title = paste0("Variant counts by caller (stage: ", rep_stage, ")"),
      subtitle = sample_id
    )

  save_plot(p_counts, paste0("variant_counts_by_caller_", rep_stage), width = 7, height = 5)

  # stacked proportions
  df_prop <- df %>%
    group_by(caller) %>%
    mutate(prop = ifelse(sum(count) > 0, count / sum(count), NA_real_)) %>%
    ungroup()

  p_prop <- ggplot(df_prop, aes(x = caller, y = prop, fill = vartype)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = vartype_colors, drop = FALSE) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = NULL, y = "Proportion", fill = "Type",
      title = paste0("Variant composition by caller (stage: ", rep_stage, ")"),
      subtitle = sample_id
    )

  save_plot(p_prop, paste0("variant_proportion_by_caller_", rep_stage), width = 7, height = 5)
}

# ----------------------------
# Fig 3A: K-of-N curve (cumulative >=K) - Combined count and percentage
# ----------------------------
if (!is.null(kcurve) && nrow(kcurve) > 0) {
  cat("Generating: Fig2A K-of-N curve (combined)\n")
  
  # Create label with both count and percentage
  kcurve <- kcurve %>%
    mutate(
      label_text = paste0(
        format(pass_variants, big.mark = ","), "\n",
        "(", round(percent_pass, 1), "%)"
      )
    )
  
  # Calculate scaling factor for secondary axis
  max_variants <- max(kcurve$pass_variants, na.rm = TRUE)
  scale_factor <- max_variants / 100
  
  p_combined <- ggplot(kcurve, aes(x = K)) +
    # Count line (primary y-axis)
    geom_line(aes(y = pass_variants), linewidth = 1.1, color = "#377EB8") +
    geom_point(aes(y = pass_variants), size = 3, color = "#377EB8") +
    # Labels with both count and percentage
    geom_label(aes(y = pass_variants, label = label_text),
               vjust = -0.3, size = 3, color = "black", 
               fill = "white", alpha = 0.8, label.size = 0.2,
               label.padding = unit(0.15, "lines")) +
    # Scales
    scale_x_continuous(breaks = sort(unique(kcurve$K))) +
    scale_y_continuous(
      name = "Variants passing (>=K)",
      labels = label_comma(),
      limits = c(0, max_variants * 1.25),
      sec.axis = sec_axis(~ . / scale_factor, 
                          name = "Percent of union retained",
                          labels = function(x) paste0(x, "%"))
    ) +
    labs(
      x = "K threshold (min callers required)",
      title = "K-of-N consensus curve",
      subtitle = sample_id
    ) +
    theme(
      axis.title.y.left = element_text(color = "#377EB8"),
      axis.text.y.left = element_text(color = "#377EB8"),
      axis.title.y.right = element_text(color = "#E41A1C"),
      axis.text.y.right = element_text(color = "#E41A1C")
    )
  
  save_plot(p_combined, "k_curve_combined", width = 8, height = 6)
}

# ----------------------------
# Fig 3B: K-of-N curve by vartype (cumulative >=K)
# ----------------------------
if (!is.null(kcurve_vartype) && nrow(kcurve_vartype) > 0) {
  cat("Generating: Fig2B K-of-N by vartype\n")

  df <- kcurve_vartype %>%
    mutate(vartype = factor(vartype, levels = vartype_levels))

  p_line <- ggplot(df, aes(x = K, y = count, color = vartype, group = vartype)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.3) +
    scale_color_manual(values = vartype_colors, drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(df$K))) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = "K threshold",
      y = "Variants passing (>=K)",
      color = "Type",
      title = "K-of-N by variant type",
      subtitle = sample_id
    )

  save_plot(p_line, "k_curve_by_vartype", width = 8, height = 5)

  p_bar <- ggplot(df, aes(x = factor(K), y = count, fill = vartype)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = vartype_colors, drop = FALSE) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = "K threshold",
      y = "Variants passing (>=K)",
      fill = "Type",
      title = "K-of-N by type (stacked)",
      subtitle = sample_id
    )

  save_plot(p_bar, "k_curve_by_vartype_bar", width = 7, height = 5)
}

# ----------------------------
# Fig 4A: Caller agreement heatmap (counts)
# ----------------------------
if (!is.null(caller_agreement) && nrow(caller_agreement) > 1) {
  cat("Generating: Fig3A caller agreement heatmap\n")

  df <- caller_agreement %>%
    pivot_longer(-caller, names_to = "caller2", values_to = "overlap") %>%
    rename(caller1 = caller)

  p <- ggplot(df, aes(x = caller1, y = caller2, fill = overlap)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = format(overlap, big.mark = ",")), size = 3.2) +
    scale_fill_distiller(palette = "Blues", direction = 1, labels = label_comma()) +
    labs(
      x = NULL, y = NULL, fill = "Shared\nvariants",
      title = "Caller agreement matrix",
      subtitle = paste0(sample_id, " (diagonal = total per caller)")
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    coord_fixed()

  save_plot(p, "caller_agreement_heatmap", width = 6, height = 5)
}

# ----------------------------
# Fig 4B: Caller contribution (unique vs shared buckets)
# ----------------------------
if (!is.null(caller_contrib) && nrow(caller_contrib) > 0) {
  cat("Generating: Fig3B caller contribution\n")

  df <- caller_contrib %>%
    pivot_longer(
      cols = c(unique, shared_2, shared_3, shared_all),
      names_to = "category",
      values_to = "count"
    ) %>%
    mutate(
      category = factor(category,
                        levels = c("unique", "shared_2", "shared_3", "shared_all"),
                        labels = c("Unique", "Shared (2)", "Shared (3)", "Shared (all)"))
    )

  contrib_colors <- c(
    "Unique" = "#E41A1C",
    "Shared (2)" = "#FF7F00",
    "Shared (3)" = "#4DAF4A",
    "Shared (all)" = "#377EB8"
  )

  p_counts <- ggplot(df, aes(x = caller, y = count, fill = category)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = contrib_colors) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = NULL, y = "Variant count", fill = "Agreement",
      title = "Caller contribution: unique vs shared",
      subtitle = sample_id
    )

  save_plot(p_counts, "caller_contribution", width = 7, height = 5)
}

# ----------------------------
# Fig 5: UpSet plot (if UpSetR available)
# ----------------------------
if (has_upsetr && !is.null(upset_matrix) && nrow(upset_matrix) > 0) {
  cat("Generating: Fig4 UpSet plot\n")

  meta_cols <- c("sample", "CHROM", "POS", "REF", "ALT", "n_callers", "callers", "vartype")
  caller_cols <- setdiff(names(upset_matrix), meta_cols)

  if (length(caller_cols) >= 2) {
    upset_df <- upset_matrix %>%
      select(all_of(caller_cols)) %>%
      as.data.frame()

    # PDF
    pdf(file.path(out_dir, "upset_callers.pdf"), width = 8, height = 5)
    print(
      upset(
        upset_df,
        sets = caller_cols,
        order.by = "freq",
        decreasing = TRUE,
        keep.order = FALSE,
        mainbar.y.label = "Intersection size",
        sets.x.label = "Variants per caller",
        text.scale = 1.2,
        point.size = 3,
        line.size = 1
      )
    )
    dev.off()

    if (out_format %in% c("png", "both")) {
      png(file.path(out_dir, "upset_callers.png"),
          width = 8, height = 5, units = "in", res = plot_dpi)
      print(
        upset(
          upset_df,
          sets = caller_cols,
          order.by = "freq",
          decreasing = TRUE,
          keep.order = FALSE,
          mainbar.y.label = "Intersection size",
          sets.x.label = "Variants per caller",
          text.scale = 1.2,
          point.size = 3,
          line.size = 1
        )
      )
      dev.off()
    }
  }
}

# ----------------------------
# Fig 6: Variants per chromosome (grouped by caller) at a representative stage
# ----------------------------
if (!is.null(genome_summary) && nrow(genome_summary) > 0) {
  cat("Generating: Fig5A variants per chromosome\n")

  rep_stage <- intersect(c("consensus", "norm", "split", "merged"), unique(genome_summary$stage))[1]
  if (is.na(rep_stage) || is.null(rep_stage)) rep_stage <- unique(genome_summary$stage)[1]

  chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM",
                    as.character(1:22), "X", "Y", "MT")

  df <- genome_summary %>%
    filter(stage == rep_stage) %>%
    mutate(
      CHROM = factor(CHROM, levels = chrom_levels),
      caller = factor(caller, levels = sort(unique(caller)))
    )

  p <- ggplot(df, aes(x = CHROM, y = total_variants, fill = caller)) +
    geom_col(position = "dodge", width = 0.72) +
    scale_fill_manual(values = caller_colors, na.value = "grey60") +
    scale_y_continuous(labels = label_comma()) +
    labs(
      x = NULL, y = "Variant count", fill = "Caller",
      title = paste0("Variants per chromosome (stage: ", rep_stage, ")"),
      subtitle = sample_id
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  save_plot(p, "variants_per_chromosome", width = 12, height = 5)
}


# ----------------------------
# Fig 7: Venn diagram - ALWAYS generate
# Shows unique and overlap between all callers
# ----------------------------
has_venn <- requireNamespace("VennDiagram", quietly = TRUE)
cat("Generating: Fig6 Venn diagram\n")
suppressPackageStartupMessages(library(VennDiagram))
# Helper functions for Venn
create_venn_sets <- function(upset_df) {
  meta_cols <- c("sample", "CHROM", "POS", "REF", "ALT", "n_callers", "callers", "vartype")
  caller_cols <- setdiff(names(upset_df), meta_cols)
  
  if (length(caller_cols) < 2) return(NULL)
  
  upset_df <- upset_df %>% mutate(key = paste(CHROM, POS, REF, ALT, sep = ":"))
  
  sets <- lapply(caller_cols, function(cc) upset_df$key[upset_df[[cc]] == 1])
  names(sets) <- caller_cols
  
  return(sets)
}

calc_venn_regions <- function(sets) {
  n <- length(sets)
  callers <- names(sets)
  
  if (n == 2) {
    a <- sets[[1]]; b <- sets[[2]]
    only_a <- length(setdiff(a, b))
    only_b <- length(setdiff(b, a))
    both <- length(intersect(a, b))
    
    return(data.frame(
      region = c(paste0(callers[1], " only"), paste0(callers[2], " only"), 
                 paste(callers[1], callers[2], sep = " ∩ ")),
      count = c(only_a, only_b, both),
      type = c("unique", "unique", "shared_2"),
      stringsAsFactors = FALSE
    ))
  } else if (n == 3) {
    a <- sets[[1]]; b <- sets[[2]]; c <- sets[[3]]
    
    only_a <- length(setdiff(setdiff(a, b), c))
    only_b <- length(setdiff(setdiff(b, a), c))
    only_c <- length(setdiff(setdiff(c, a), b))
    ab_not_c <- length(setdiff(intersect(a, b), c))
    ac_not_b <- length(setdiff(intersect(a, c), b))
    bc_not_a <- length(setdiff(intersect(b, c), a))
    abc <- length(intersect(intersect(a, b), c))
    
    return(data.frame(
      region = c(
        paste0(callers[1], " only"), paste0(callers[2], " only"), paste0(callers[3], " only"),
        paste(callers[1], callers[2], sep = " ∩ "),
        paste(callers[1], callers[3], sep = " ∩ "),
        paste(callers[2], callers[3], sep = " ∩ "),
        paste(callers, collapse = " ∩ ")
      ),
      count = c(only_a, only_b, only_c, ab_not_c, ac_not_b, bc_not_a, abc),
      type = c("unique", "unique", "unique", "shared_2", "shared_2", "shared_2", "shared_all"),
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

if (!is.null(upset_matrix) && nrow(upset_matrix) > 0) {
  sets <- create_venn_sets(upset_matrix)
  
  if (!is.null(sets) && length(sets) >= 2) {
    n_callers <- length(sets)
    venn_regions <- calc_venn_regions(sets)
    
    # Save region counts as TSV
    if (!is.null(venn_regions)) {
      write_tsv(venn_regions, file.path(out_dir, "venn_region_counts.tsv"))
      cat("  Saved venn_region_counts.tsv\n")
    }
    
    # Use VennDiagram package if available
    if (has_venn && n_callers <= 5) {
      venn_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")[1:n_callers]
      
      venn_plot <- venn.diagram(
        x = sets,
        filename = NULL,
        main = paste0("Caller Overlap (", sample_id, ")"),
        main.cex = 1.3,
        main.fontface = "bold",
        sub = paste0("Total unique variants: ", length(unique(unlist(sets)))),
        sub.cex = 0.9,
        fill = venn_colors,
        alpha = 0.5,
        cex = 1.2,
        cat.cex = 1.1,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        margin = 0.1,
        euler.d = TRUE,
        scaled = TRUE
      )
      
      if (out_format %in% c("png", "both")) {
        png(file.path(out_dir, "venn_callers.png"), width = plot_width, height = plot_height, units = "in", res = plot_dpi)
        grid::grid.draw(venn_plot)
        dev.off()
      }
      
      if (out_format %in% c("pdf", "both")) {
        pdf(file.path(out_dir, "venn_callers.pdf"), width = plot_width, height = plot_height)
        grid::grid.draw(venn_plot)
        dev.off()
      }
      
      cat("  Generated Venn diagram with VennDiagram package\n")
    } else {
      cat("  VennDiagram package not available or >5 callers\n")
    }
    
    # Always create bar chart version of Venn regions
    if (!is.null(venn_regions)) {
      df <- venn_regions %>%
        mutate(type = factor(type, levels = c("unique", "shared_2", "shared_all")))
      
      type_colors <- c("unique" = "#E41A1C", "shared_2" = "#FF7F00", "shared_all" = "#4DAF4A")
      
      p_bar <- ggplot(df, aes(x = reorder(region, count), y = count, fill = type)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = format(count, big.mark = ",")), hjust = -0.1, size = 3.2) +
        coord_flip() +
        scale_fill_manual(values = type_colors, labels = c("Unique", "Shared (2)", "Shared (all)")) +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.18))) +
        labs(x = NULL, y = "Variant count", fill = "Type",
             title = "Venn Region Breakdown", subtitle = paste0(sample_id, " — Total: ", length(unique(unlist(sets)))))
      
      save_plot(p_bar, "venn_region_bar", width = 8, height = 5)
    }
  }
} else if (!is.null(caller_contrib) && nrow(caller_contrib) > 0) {
  # Fallback from caller_contrib
  cat("  Creating Venn-like visualization from caller_contrib\n")
  
  df <- caller_contrib %>%
    pivot_longer(cols = c(unique, shared_2, shared_3, shared_all), names_to = "category", values_to = "count") %>%
    filter(count > 0) %>%
    mutate(category = factor(category, levels = c("unique", "shared_2", "shared_3", "shared_all"),
                             labels = c("Unique", "Shared (2)", "Shared (3)", "Shared (all)")))
  
  contrib_colors <- c("Unique" = "#E41A1C", "Shared (2)" = "#FF7F00", "Shared (3)" = "#4DAF4A", "Shared (all)" = "#377EB8")
  
  p <- ggplot(df, aes(x = caller, y = count, fill = category)) +
    geom_col(width = 0.7, position = "dodge") +
    geom_text(aes(label = format(count, big.mark = ",")), position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
    scale_fill_manual(values = contrib_colors) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Variant count", fill = "Agreement",
         title = "Caller Agreement Distribution", subtitle = sample_id)
  
  save_plot(p, "caller_agreement_distribution", width = 9, height = 5)
} else {
  cat("  Skipping Venn: no data available\n")
}


# ----------------------------
# Fig8: Pie plot: proportion of SNP/MNV/MULTIALLELIC/INDEL in each caller
# ----------------------------

if (!is.null(counts_long) && nrow(counts_long) > 0) {
  message("Generating: Pie plots of variant type composition per caller")

  # choose stage to summarize (prefer norm > split > raw > consensus)
  rep_stage <- intersect(c("norm","split","raw","consensus"), unique(counts_long$stage))[1]
  if (is.na(rep_stage)) rep_stage <- unique(counts_long$stage)[1]

  df <- counts_long %>%
    filter(stage == rep_stage) %>%
    mutate(
      vartype = factor(vartype, levels = vartype_levels),
      caller = as.character(caller)
    ) %>%
    group_by(caller, vartype) %>%
    summarize(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    group_by(caller) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()

  # Pie via coord_polar; facet by caller
  p_pie <- ggplot(df, aes(x = "", y = prop, fill = vartype)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    facet_wrap(~ caller, nrow = 1) +
    scale_fill_manual(values = vartype_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = NULL, y = NULL,
      fill = "Type",
      title = paste0("Variant composition by caller (stage: ", rep_stage, ")"),
      subtitle = sample_id
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  save_plot(p_pie, paste0("pie_variant_composition_by_caller_", rep_stage),
            width = max(8, 3 * length(unique(df$caller))), height = 4.5)
}

# Fig6C: Full vs Partial agreement (simplified binary view)
if (!is.null(upset_matrix) && nrow(upset_matrix) > 0 && "n_callers" %in% names(upset_matrix)) {
  n_total_callers <- length(setdiff(names(upset_matrix), meta_cols))
  
  df <- upset_matrix %>%
    mutate(
      agreement = ifelse(n_callers == n_total_callers, "Full agreement", "Partial agreement"),
      agreement = factor(agreement, levels = c("Full agreement", "Partial agreement"))
    )
  
  df_summary <- df %>%
    count(agreement, name = "count") %>%
    mutate(
      proportion = count / sum(count),
      pct_label = sprintf("%.1f%%", proportion * 100)
    )
  
  binary_colors <- c("Full agreement" = "#2E7D32", "Partial agreement" = "#FFA726")
  
  p <- ggplot(df_summary, aes(agreement, proportion, fill = agreement)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(format(count, big.mark = ","), "\n(", pct_label, ")")), 
              vjust = -0.2, size = 4) +
    scale_fill_manual(values = binary_colors) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Proportion of consensus variants",
         title = "Full vs Partial caller agreement in consensus",
         subtitle = paste0(sample_id, " — Full = all ", n_total_callers, " callers agree")) +
    theme(legend.position = "none")
  save_plot(p, "full_vs_partial_agreement", 6, 5)
}

# Fig9: Full vs Partial agreement (simplified binary view)
plot_event_status_pies <- function(
  event_df,
  stage = NULL,
  caller = NULL,
  sample_id = "Sample",
  out_basename = "FigX_event_status_pies",
  width = 8,
  height = 4.8
) {
  if (is.null(event_df) || nrow(event_df) == 0) return(invisible(NULL))

  # Required columns
  req <- c("EVENT_TYPE", "EVENT_STATUS", "count")
  missing <- setdiff(req, names(event_df))
  if (length(missing) > 0) {
    stop("plot_event_status_pies: missing required columns: ", paste(missing, collapse = ", "))
  }

  df <- event_df

  # Optional filtering
  if (!is.null(stage) && "stage" %in% names(df)) {
    df <- df %>% dplyr::filter(stage == stage)
  }
  if (!is.null(caller) && "caller" %in% names(df)) {
    df <- df %>% dplyr::filter(caller == caller)
  }

  # Focus on your complex event types + statuses
  df <- df %>%
    dplyr::filter(EVENT_TYPE %in% c("MNV", "MULTIALLELIC")) %>%
    dplyr::filter(EVENT_STATUS %in% c("FULL", "PARTIAL"))

  if (nrow(df) == 0) return(invisible(NULL))

  # Aggregate in case table is more granular than needed
  df2 <- df %>%
    dplyr::group_by(EVENT_TYPE, EVENT_STATUS) %>%
    dplyr::summarize(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(EVENT_TYPE) %>%
    dplyr::mutate(
      total = sum(count),
      pct = ifelse(total > 0, 100 * count / total, 0),
      label = paste0(EVENT_STATUS, "\n", "N=", scales::comma(count), "\n", sprintf("%.1f%%", pct))
    ) %>%
    dplyr::ungroup()

  # Make ordering stable
  df2 <- df2 %>%
    dplyr::mutate(
      EVENT_TYPE = factor(EVENT_TYPE, levels = c("MNV", "MULTIALLELIC")),
      EVENT_STATUS = factor(EVENT_STATUS, levels = c("FULL", "PARTIAL"))
    )

  # Colors (you can change if you want)
  status_colors <- c(
    "FULL" = "#4DAF4A",
    "PARTIAL" = "#FF7F00"
  )

  subtitle <- paste(
    sample_id,
    if (!is.null(stage)) paste0("stage=", stage) else NULL,
    if (!is.null(caller)) paste0("caller=", caller) else NULL,
    sep = "  |  "
  )

  p <- ggplot2::ggplot(df2, ggplot2::aes(x = "", y = count, fill = EVENT_STATUS)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.7,
      lineheight = 0.95
    ) +
    ggplot2::facet_wrap(~ EVENT_TYPE, nrow = 1) +
    ggplot2::scale_fill_manual(values = status_colors) +
    ggplot2::labs(
      title = "Event completeness by event type",
      subtitle = subtitle,
      fill = "Event status"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey35"),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  save_plot(p, out_basename, width = width, height = height)
}

plot_mnv_multiallelic_pie <- function(
  df,
  stage = NULL,
  caller = NULL,
  title_suffix = "",
  out_basename = "FigX_mnv_multiallelic_pie"
) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Expected columns:
  # site-based: sample stage caller vartype count
  # event-based: sample stage caller EVENT_TYPE EVENT_STATUS count

  # Normalize column names
  if ("EVENT_TYPE" %in% names(df)) {
    df2 <- df %>%
      rename(vartype = EVENT_TYPE)
  } else {
    df2 <- df
  }

  # Filter to MNV / MULTIALLELIC
  df2 <- df2 %>%
    filter(vartype %in% c("MNV", "MULTIALLELIC"))

  if (!is.null(stage) && "stage" %in% names(df2)) {
    df2 <- df2 %>% filter(stage == stage)
  }

  if (!is.null(caller) && "caller" %in% names(df2)) {
    df2 <- df2 %>% filter(caller == caller)
  }

  if (nrow(df2) == 0) return(NULL)

  pie_df <- df2 %>%
    group_by(vartype) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(
      prop = count / sum(count),
      label = paste0(vartype, "\n", scales::percent(prop, accuracy = 0.1))
    )

  p <- ggplot(pie_df, aes(x = "", y = count, fill = vartype)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    scale_fill_manual(values = vartype_colors) +
    labs(
      title = "MNV vs Multiallelic variants",
      subtitle = title_suffix,
      fill = "Variant type"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )

  save_plot(p, out_basename, width = 5, height = 5)
}

plot_mnv_multiallelic_pie(
  df = event_status_by_type,
  stage = "merged",
  title_suffix = paste(sample_id, "— merged events"),
  out_basename = "mnv_vs_multiallelic_events"
)

plot_event_status_pies(
  event_df = event_status_by_type,
  stage = "merged",
  sample_id = sample_id,
  out_basename = "event_status_pies_merged"
)
# ----------------------------
# Summary
# ----------------------------
cat("\n========================================\n")
cat("Plots generated in:", out_dir, "\n")
cat("========================================\n")

plots_generated <- list.files(out_dir, pattern = "\\.(png|pdf)$", full.names = FALSE)
cat("Files created:", length(plots_generated), "\n")
for (f in sort(plots_generated)) cat("  ", f, "\n")
