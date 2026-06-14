# ============================================================
# Publication-Ready 4-Cohort Combined Enrichment Figures
# Works for GO:BP, KEGG, and Reactome
# Layout: 2x2 grid (one panel per cohort)
# ============================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(scales)

# ============================================================
# SECTION 1: FILE PATHS — edit these to your actual files
# ============================================================

# ── GO:BP files ──────────────────────────────────────────────
bp_luad <- read.csv("D:/Mcode final__LUAD/output/ORA_Results_BP.csv")
bp_hnsc <- read.csv("D:/Mcode final__HNSC/output/ORA_Results_BP.csv")
bp_kirc <- read.csv("D:/Mcode final__KIRC/output/ORA_Results_BP.csv")
bp_coad <- read.csv("D:/Mcode final__COAD/output/ORA_Results_BP.csv")

# ── KEGG files ───────────────────────────────────────────────
kegg_luad <- read.csv("D:/Mcode final__LUAD/output/ORA_Results_KEGG.csv")
kegg_hnsc <- read.csv("D:/Mcode final__HNSC/output/ORA_Results_KEGG.csv")
kegg_kirc <- read.csv("D:/Mcode final__KIRC/output/ORA_Results_KEGG.csv")
kegg_coad <- read.csv("D:/Mcode final__COAD/output/ORA_Results_KEGG.csv")

# ── Reactome files ───────────────────────────────────────────
reactome_luad <- read.csv("D:/Mcode final__LUAD/output/ORA_Results_Reactome.csv")
reactome_hnsc <- read.csv("D:/Mcode final__HNSC/output/ORA_Results_Reactome.csv")
reactome_kirc <- read.csv("D:/Mcode final__KIRC/output/ORA_Results_Reactome.csv")
reactome_coad <- read.csv("D:/Mcode final__COAD/output/ORA_Results_Reactome.csv")

# ============================================================
# SECTION 2: PARAMETERS — adjust if needed
# ============================================================

TOP_N_TERMS  <- 10      # how many top terms to show per panel
WRAP_WIDTH   <- 40      # character wrap width for term labels
POINT_RANGE  <- c(4, 12) # Adjusted dot sizes up slightly to balance with the large text sizes

# ============================================================
# SECTION 3: HELPER FUNCTIONS
# ============================================================

# Parse "6/121" → 0.0496
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"),
         function(v) as.numeric(v[1]) / as.numeric(v[2]))
}

# Prepare a single cohort dataframe for plotting
prep_data <- function(df, top_n = TOP_N_TERMS, wrap = WRAP_WIDTH) {
  df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio_num = parse_ratio(GeneRatio),
      Description   = str_wrap(Description, width = wrap),
      Description   = factor(Description, levels = rev(unique(Description)))
    )
}

# Build one dotplot panel
make_panel <- function(df, cohort_label, color_low, color_high) {
  
  # Handle case where cohort has no significant results
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0(cohort_label, "\nNo significant terms"),
                 size = 6, fontface = "bold", color = "grey50", hjust = 0.5, vjust = 0.5) +
        theme_void() +
        labs(title = cohort_label) +
        theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  }
  
  ggplot(df, aes(x = GeneRatio_num, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(
      low    = color_low,
      high   = color_high,
      name   = "Adjusted\np-value",
      labels = scientific_format(digits = 2)
    ) +
    scale_size_continuous(
      name  = "Gene\nCount",
      range = POINT_RANGE
    ) +
    labs(
      title = cohort_label,
      x     = "Gene Ratio",
      y     = NULL
    ) +
    theme_bw(base_size = 18) + # Set global base text metrics across subplots
    theme(
      # All axis and mapping text parameters converted to Size 18 and Bold
      axis.text.y        = element_text(size = 18, face = "bold", color = "black"),
      axis.text.x        = element_text(size = 18, face = "bold", color = "black"),
      axis.title.x       = element_text(size = 18, face = "bold", margin = margin(t = 12)),
      plot.title         = element_text(size = 18, face = "bold",
                                        hjust = 0.5, color = "#2C3E50", margin = margin(b = 12)),
      legend.title       = element_text(size = 18, face = "bold"),
      legend.text        = element_text(size = 18, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right"
    )
}

# Build combined 2x2 figure and save
make_combined_figure <- function(luad_df, hnsc_df, kirc_df, coad_df,
                                 db_label,
                                 color_low, color_high,
                                 out_prefix,
                                 top_n    = TOP_N_TERMS,
                                 fig_w    = 24, # Increased width from 18 to hold long bold text cleanly
                                 fig_h    = 22) # Increased height from 16 to resolve Word clipping
{
  
  # Prep all four cohorts
  d_luad <- prep_data(luad_df, top_n)
  d_hnsc <- prep_data(hnsc_df, top_n)
  d_kirc <- prep_data(kirc_df, top_n)
  d_coad <- prep_data(coad_df, top_n)
  
  # Build panels
  p_luad <- make_panel(d_luad, "LUAD", color_low, color_high)
  p_hnsc <- make_panel(d_hnsc, "HNSC", color_low, color_high)
  p_kirc <- make_panel(d_kirc, "KIRC", color_low, color_high)
  p_coad <- make_panel(d_coad, "COAD", color_low, color_high)
  
  # Combine into 2x2 grid
  combined <- (p_luad | p_hnsc) / (p_kirc | p_coad) +
    plot_annotation(
      title   = paste0(db_label, " Enrichment Analysis Across Cohorts"),
      caption = paste0(
        "Dot size: gene count; Color: BH-adjusted p-value. ",
        "Top ", top_n, " significant terms per cohort shown."
      ),
      theme = theme(
        # Global annotation labels transformed to Size 18 and Bold
        plot.title   = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        plot.caption = element_text(size = 18, face = "bold", color = "grey40", hjust = 0)
      )
    )
  
  # Save PDF (publication quality)
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot     = combined,
    width    = fig_w,
    height   = fig_h,
    units    = "in",
    device   = cairo_pdf
  )
  
  # Save PNG with high-definition 600 DPI execution specs to prevent blurriness inside Word files
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot     = combined,
    width    = fig_w,
    height   = fig_h,
    units    = "in",
    dpi      = 600
  )
  
  message("Saved: ", out_prefix, ".pdf and .png")
  invisible(combined)
}

# ============================================================
# SECTION 4: MAKE COMBINED TABLES
# ============================================================

make_combined_table <- function(luad_df, hnsc_df, kirc_df, coad_df,
                                db_label, out_prefix,
                                top_n = TOP_N_TERMS) {
  
  add_cohort <- function(df, cohort) {
    df %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n) %>%
      mutate(Cohort = cohort) %>%
      select(Cohort, ID, Description, GeneRatio, Count,
             pvalue, p.adjust, geneID)
  }
  
  combined_table <- bind_rows(
    add_cohort(luad_df, "LUAD"),
    add_cohort(hnsc_df, "HNSC"),
    add_cohort(kirc_df, "KIRC"),
    add_cohort(coad_df, "COAD")
  ) %>%
    rename(
      `Gene Ratio`   = GeneRatio,
      `Gene Count`   = Count,
      `p-value`      = pvalue,
      `Adj. p-value` = p.adjust,
      `Gene IDs`     = geneID
    )
  
  write.csv(combined_table,
            file      = paste0(out_prefix, "_combined_table.csv"),
            row.names = FALSE)
  
  message("Saved: ", out_prefix, "_combined_table.csv")
  invisible(combined_table)
}

# ============================================================
# SECTION 5: RUN — one call per database
# ============================================================

# ── GO:BP ────────────────────────────────────────────────────
make_combined_figure(
  luad_df    = bp_luad,
  hnsc_df    = bp_hnsc,
  kirc_df    = bp_kirc,
  coad_df    = bp_coad,
  db_label   = "GO Biological Process",
  color_low  = "#B22222",
  color_high = "#4682B4",
  out_prefix = "Figure_GO_BP_4cohorts"
)

make_combined_table(
  luad_df    = bp_luad,
  hnsc_df    = bp_hnsc,
  kirc_df    = bp_kirc,
  coad_df    = bp_coad,
  db_label   = "GO:BP",
  out_prefix = "Figure_GO_BP_4cohorts"
)

# ── KEGG ─────────────────────────────────────────────────────
make_combined_figure(
  luad_df    = kegg_luad,
  hnsc_df    = kegg_hnsc,
  kirc_df    = kegg_kirc,
  coad_df    = kegg_coad,
  db_label   = "KEGG Pathway",
  color_low  = "#2E8B57",
  color_high = "#9370DB",
  out_prefix = "Figure_KEGG_4cohorts"
)

make_combined_table(
  luad_df    = kegg_luad,
  hnsc_df    = kegg_hnsc,
  kirc_df    = kegg_kirc,
  coad_df    = kegg_coad,
  db_label   = "KEGG",
  out_prefix = "Figure_KEGG_4cohorts"
)

# ── Reactome ─────────────────────────────────────────────────
make_combined_figure(
  luad_df    = reactome_luad,
  hnsc_df    = reactome_hnsc,
  kirc_df    = reactome_kirc,
  coad_df    = reactome_coad,
  db_label   = "Reactome Pathway",
  color_low  = "#8B0000",
  color_high = "#1E90FF",
  out_prefix = "Figure_Reactome_4cohorts",
  top_n      = 15,   # Reactome usually has more terms, show top 15
  fig_w      = 24,   # Explicit size allocations maintained
  fig_h      = 26    # Height bumped further to handle 15 bold multi-line rows
)

make_combined_table(
  luad_df    = reactome_luad,
  hnsc_df    = reactome_hnsc,
  kirc_df    = reactome_kirc,
  coad_df    = reactome_coad,
  db_label   = "Reactome",
  out_prefix = "Figure_Reactome_4cohorts",
  top_n      = 15
)