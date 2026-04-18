# paper_figures.R
# Script to generate Figures for PLOS ONE submission

library(cbamm)
library(ggplot2)
library(patchwork)

# Set seed for reproducibility
set.seed(42)

# Figure 1: Multiverse Plot
message("Generating Figure 1: Multiverse Plot...")
data_mv <- simulate_cbamm_data(n_studies = 12)
config_mv <- cbamm_config(methods = list(multiverse = TRUE, parallel = FALSE, bayesian = FALSE))
results_mv <- cbamm(data_mv, config = config_mv)

# Custom plot for publication
mv_df <- results_mv$multiverse
p1 <- ggplot(mv_df, aes(x = estimate, y = estimator, color = subset, shape = hksj)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), position = position_dodge(width = 0.5), height = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Figure 1: Meta-Analysis Multiverse",
       subtitle = "Estimates across τ² methods, study subsets, and HKSJ settings",
       x = "Pooled Effect (Log Scale)",
       y = "τ² Estimator",
       color = "Data Subset",
       shape = "HKSJ Adj.") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Figure 2: Component NMA Forest Plot
message("Generating Figure 2: Component NMA Forest Plot...")
data_cnma <- data.frame(
  study_id = rep(1:8, each = 2),
  treatment = c("A", "Placebo", "B", "Placebo", "A + B", "A", "A + B", "Placebo", 
                "B", "A", "A + C", "A", "C", "Placebo", "A + B + C", "Placebo"),
  effect = c(0.45, 0, 0.32, 0, 0.78, 0.42, 0.85, 0, 0.35, 0.48, 0.65, 0.44, 0.22, 0, 1.1, 0),
  se = runif(16, 0.08, 0.15)
)

result_cnma <- network_meta_analysis(data_cnma, studies="study_id", treatments="treatment",
                                    effect_size="effect", std_error="se",
                                    reference="Placebo", model="component")

p2 <- cbamm_plot(result_cnma, type = "forest") +
  labs(title = "Figure 2: Component Network Meta-Analysis",
       subtitle = "Individual additive effects of Drug Components A, B, and C") +
  theme_minimal()

# Save figures
message("Saving figures...")
ggsave("Figure_1.png", p1, width = 10, height = 7, dpi = 300)
ggsave("Figure_2.png", p2, width = 10, height = 6, dpi = 300)

message("Figures generated successfully.")
