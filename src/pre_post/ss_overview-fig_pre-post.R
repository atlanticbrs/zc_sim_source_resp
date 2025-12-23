library(gridExtra)
library(grid)
library(patchwork)
# Make OverView Figure for Southall - Scaled Source Paper
# Dot Plot --------------
source('Rscripts/ss_overview-fig_pre-post_dot-chart.R') # makes p_strip
  
# 1D: Zc93 --------------
load("~/Documents/_research/rrr/dive_eda_reproducibility/output/scaled_source_24hr/univariate_plots/ggplot_univariate_ZcTag093.rdata")
p_1d_093 <- pl_uni

# 1D: Zc 96 --------------
load("~/Documents/_research/rrr/dive_eda_reproducibility/output/scaled_source_24hr/univariate_plots/ggplot_univariate_ZcTag096.rdata")
p_1d_096 <- pl_uni

# KDE: Zc93 --------------
load("~/Documents/_research/rrr/dive_eda_reproducibility/output/scaled_source_24hr/kde_plots/ggplot_kde_ZcTag093.rdata")
p_kde_093 <- pl_kde

# KDE: Zc96 --------------
load("~/Documents/_research/rrr/dive_eda_reproducibility/output/scaled_source_24hr/kde_plots/ggplot_kde_ZcTag096.rdata")
p_kde_096 <- pl_kde


# Lay Out and Print to PDF --------------------------------------
standard_theme <- theme(
  plot.margin = margin(6, 6, 6, 6),
  axis.title.y = element_text(margin = margin(r = 6)),
  axis.title.x = element_text(margin = margin(t = 6))
)

p_kde_093  <- p_kde_093  + standard_theme
p_kde_096  <- p_kde_096  + standard_theme
p_1d_093   <- p_1d_093   + standard_theme
p_1d_096   <- p_1d_096   + standard_theme
p_strip    <- p_strip    + standard_theme

# Explicit 2x3 layout: left four panels in a 2x2 grid, strip spans both rows
layout_design <- "
ABe
CDe
"

final_plot <- (
  p_kde_093 +    # A
    p_1d_093  +    # B
    p_kde_096  +    # E
    p_1d_096 +    # D
    p_strip       # c  (spans both rows)
) +
  plot_layout(design = layout_design, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

print(final_plot)

ggsave(filename = here::here(
  'output/scaled_source_24hr',
  paste0("pre-post_Figure-03.pdf")
),
device = cairo_pdf,
dpi = 320,
width = 22, height = 18, units = 'cm') 

