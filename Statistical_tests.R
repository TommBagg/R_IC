# Title: Advanced Sediment Connectivity (IC) Analysis: Temporal Shift, Lithology, & Hydrogeology
# Authors: Adapted for CAMELS-IND Peninsular India catchments

library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

###################################
# USER INPUTS & DIRECTORIES SETUP
###################################

dir_1990 <- "F:/Peninsular India/Connectivity/IC/Output_Albers/1990"
dir_2021 <- "F:/Peninsular India/Connectivity/IC/Output_Albers/2021"
csv_path <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/attributes_csv/camels_ind_geol.csv"

# Output folder for figures and summary statistics
out_dir  <- "F:/Peninsular India/Connectivity/IC/Plots/Statistical/Albers"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

col_gauge_id <- "gauge_id"
col_regime   <- "Regime"
col_geol     <- "geol_class_1st"
col_poro     <- "geol_porosity"
col_perm     <- "geol_permeability"

# Geology re-grouping dictionary
geo_grouping <- c(
  'Acid Plutonic Rocks'             = 'Igneous',    
  'Basic Volcanic Rocks'            = 'Igneous',
  'Intermediate Volcanic Rocks'     = 'Igneous',
  'Metamorphic Rocks'               = 'Metamorphic',
  'Mixed Sedimentary Rocks'         = 'Sedimentary',
  'Siliciclastic Sedimentary Rocks' = 'Sedimentary',
  'Unconsolidated Sediments'        = 'Unconsolidated'
)

# Standardized density evaluation grid
grid_n   <- 512
ic_min   <- -13
ic_max   <- 2
sample_n <- NULL # Subsample per catchment for fast density calculations


###################################
# STEP 1: LOAD METADATA & PROCESS ATTRIBUTES
###################################
cat("Loading and processing catchment attributes...\n")
meta_df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Clean IDs and apply geology grouping
meta_df$gauge_id_clean <- sprintf("%05d", as.numeric(trimws(as.character(meta_df[[col_gauge_id]]))))
meta_df$Regime         <- as.factor(trimws(as.character(meta_df[[col_regime]])))
meta_df$raw_geol       <- trimws(as.character(meta_df[[col_geol]]))
meta_df$Geology_Class  <- factor(geo_grouping[meta_df$raw_geol], 
                                 levels = c('Igneous', 'Metamorphic', 'Sedimentary', 'Unconsolidated'))

# Ensure numeric types for porosity & permeability
meta_df$Porosity     <- as.numeric(meta_df[[col_poro]])
meta_df$Permeability <- as.numeric(meta_df[[col_perm]])

files_1990 <- list.files(dir_1990, pattern = "^IC_.*\\.tif$", full.names = TRUE)
files_2021 <- list.files(dir_2021, pattern = "^IC_.*\\.tif$", full.names = TRUE)


###################################
# STEP 2: EXTRACT IC PIXELS & COMPUTE SUMMARY STATS
###################################
cat("Extracting IC distributions from 1990 and 2021 GeoTIFFs...\n")

basin_stats_list <- list()
pixel_dens_list  <- list()
matched_count    <- 0

for (i in 1:nrow(meta_df)) {
  gid <- meta_df$gauge_id_clean[i]
  
  idx_90 <- grep(paste0("IC_.*", gid, "\\.tif$"), files_1990)
  idx_21 <- grep(paste0("IC_.*", gid, "\\.tif$"), files_2021)
  
  if (length(idx_90) > 0 && length(idx_21) > 0) {
    # Load 1990 & 2021 rasters
    c_90 <- terra::values(terra::rast(files_1990[idx_90[1]]), mat = FALSE)
    c_21 <- terra::values(terra::rast(files_2021[idx_21[1]]), mat = FALSE)
    
    c_90 <- c_90[!is.na(c_90) & !is.nan(c_90) & !is.infinite(c_90)]
    c_21 <- c_21[!is.na(c_21) & !is.nan(c_21) & !is.infinite(c_21)]
    
    if (length(c_90) > 50 && length(c_21) > 50) {
      matched_count <- matched_count + 1
      
      med_90 <- median(c_90)
      med_21 <- median(c_21)
      
      # Catchment-level record
      basin_stats_list[[matched_count]] <- data.frame(
        gauge_id      = gid,
        Regime        = meta_df$Regime[i],
        Geology_Class = meta_df$Geology_Class[i],
        Porosity      = meta_df$Porosity[i],
        Permeability  = meta_df$Permeability[i],
        mean_1990     = mean(c_90),
        mean_2021     = mean(c_21),
        median_1990   = med_90,
        median_2021   = med_21,
        delta_IC      = med_21 - med_90 # Catchment-scale temporal shift
      )
      
      # Subsample for kernel density computation
      if (!is.null(sample_n)) {
        set.seed(42)
        if (length(c_90) > sample_n) c_90 <- sample(c_90, sample_n)
        if (length(c_21) > sample_n) c_21 <- sample(c_21, sample_n)
      }
      
      d_90 <- density(c_90, from = ic_min, to = ic_max, n = grid_n)
      d_21 <- density(c_21, from = ic_min, to = ic_max, n = grid_n)
      
      pixel_dens_list[[paste0(gid, "_90")]] <- data.frame(
        gauge_id = gid, Year = "1990", Regime = meta_df$Regime[i],
        Geology_Class = meta_df$Geology_Class[i], x = d_90$x, y = d_90$y
      )
      pixel_dens_list[[paste0(gid, "_21")]] <- data.frame(
        gauge_id = gid, Year = "2021", Regime = meta_df$Regime[i],
        Geology_Class = meta_df$Geology_Class[i], x = d_21$x, y = d_21$y
      )
    }
  }
}

catchment_data <- dplyr::bind_rows(basin_stats_list)
density_data   <- dplyr::bind_rows(pixel_dens_list)
cat(sprintf("Matched %d catchments with complete 1990-2021 data & geology.\n", nrow(catchment_data)))

# Save master catchment dataset
write.csv(catchment_data, file.path(out_dir, "Catchment_IC_Geology_Synthesis.csv"), row.names = FALSE)


###################################################################################################################
# ANALYSIS & PLOT 1: TEMPORAL CHANGE (DELTA IC) REGIME-WISE
###################################################################################################################
cat("\n--- Running Analysis 1: Delta IC by Regime ---\n")

# Wilcoxon test comparing Delta IC between regimes
delta_regime_test <- wilcox.test(delta_IC ~ Regime, data = catchment_data)
cat(sprintf("Wilcoxon Rank-Sum p-value for Delta IC across Regimes: %.4e\n", delta_regime_test$p.value))

regime_palette <- c("0" = "#2b5c8f", "1" = "#d95f02")

p1 <- ggplot(catchment_data, aes(x = Regime, y = delta_IC, fill = Regime)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
  # Violin plot showing distribution shape
  geom_violin(alpha = 0.45, trim = FALSE, color = "gray30") +
  # Boxplot overlaid inside violin
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.8, color = "black") +
  # Individual catchment points with jitter
  geom_jitter(width = 0.08, alpha = 0.65, size = 2, color = "gray20") +
  scale_fill_manual(values = regime_palette) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("Mann-Whitney U Test:\np = %.3e %s", 
                           delta_regime_test$p.value, 
                           ifelse(delta_regime_test$p.value < 0.05, "(Significant)", "(Not Sig.)")),
           size = 4, fontface = "bold.italic") +
  labs(
    title = expression(bold("Catchment Temporal Shift in Connectivity: ") ~ bold(Delta * "IC (2021 - 1990)")),
    subtitle = "Positive: Increased connectivity (degradation) | Negative: Decreased connectivity (greening)",
    x = "Hydrological Regime",
    y = expression(Delta * "IC (Median"[2021] * " - Median"[1990] * ") [-]")
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figure1_Delta_IC_by_Regime.png"), plot = p1, width = 7.5, height = 6.5, dpi = 300)


###################################################################################################################
# ANALYSIS & PLOT 2: IC RELATION WITH GEOLOGY CLASSES
###################################################################################################################
cat("\n--- Running Analysis 2: IC across Geological Classes ---\n")

# Filter complete geology cases
geol_data <- catchment_data %>% filter(!is.na(Geology_Class))

# 1. Kruskal-Wallis Test (Non-parametric ANOVA)
kw_test_1990 <- kruskal.test(median_1990 ~ Geology_Class, data = geol_data)
kw_test_2021 <- kruskal.test(median_2021 ~ Geology_Class, data = geol_data)
kw_test_diff <- kruskal.test(delta_IC ~ Geology_Class, data = geol_data)

cat(sprintf("Kruskal-Wallis Test across Geology Classes (2021 IC): p = %.4e\n", kw_test_2021$p.value))
cat(sprintf("Kruskal-Wallis Test across Geology Classes (Delta IC): p = %.4e\n", kw_test_diff$p.value))

# 2. Pairwise Dunn's / Wilcoxon Post-Hoc Tests
posthoc_results <- pairwise.wilcox.test(geol_data$median_2021, geol_data$Geology_Class, p.adjust.method = "holm")
print(posthoc_results)

# Geology color palette
geol_palette <- c(
  "Igneous"        = "#e41a1c", # Red/Volcanic
  "Metamorphic"    = "#984ea3", # Purple
  "Sedimentary"    = "#4daf4a", # Green
  "Unconsolidated" = "#ff7f00"  # Amber/Sediment
)

# Panel A: Boxplot with Regime interaction
p2a <- ggplot(geol_data, aes(x = Geology_Class, y = median_2021, fill = Geology_Class)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, color = "black") +
  geom_jitter(aes(shape = Regime), width = 0.15, alpha = 0.7, size = 2.2) +
  scale_fill_manual(values = geol_palette) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
           label = sprintf("Kruskal-Wallis (2021): p = %.3e", kw_test_2021$p.value),
           size = 4, fontface = "italic") +
  labs(
    title = "(a) Catchment Median IC by Geology Class (2021)",
    x = "Geological Class",
    y = expression("Median IC"[2021] ~ "[-]")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

# Panel B: Kernel Density Envelopes with 95% CIs per Geology Class (2021)
ci_geol_2021 <- density_data %>%
  filter(Year == "2021" & !is.na(Geology_Class)) %>%
  group_by(Geology_Class, x) %>%
  summarize(
    mean_y   = mean(y),
    se_y     = sd(y) / sqrt(n()),
    ci_lower = pmax(0, mean_y - 1.96 * se_y),
    ci_upper = mean_y + 1.96 * se_y,
    .groups  = "drop"
  )

p2b <- ggplot(ci_geol_2021, aes(x = x, color = Geology_Class, fill = Geology_Class)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = mean_y), size = 1.25) +
  geom_line(aes(y = ci_lower), linetype = "dashed", size = 0.45) +
  geom_line(aes(y = ci_upper), linetype = "dashed", size = 0.45) +
  scale_color_manual(values = geol_palette) +
  scale_fill_manual(values = geol_palette) +
  labs(
    title = "(b) Mean Density Profiles by Geology (with 95% CI Bands)",
    x = "Index of Connectivity (IC) [-]",
    y = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

fig2_combined <- plot_grid(p2a, p2b, nrow = 2, align = "v")
ggsave(file.path(out_dir, "Figure2_IC_vs_Geology_Classes.png"), plot = fig2_combined, width = 10, height = 9, dpi = 300)


###################################################################################################################
# ANALYSIS & PLOT 3: IC RELATION WITH POROSITY & PERMEABILITY
###################################################################################################################
cat("\n--- Running Analysis 3: Continuous Hydrogeology (Porosity & Permeability) ---\n")

hydro_data <- catchment_data %>% filter(!is.na(Porosity) & !is.na(Permeability))

# Correlations
cor_poro <- cor.test(hydro_data$Porosity, hydro_data$median_2021, method = "spearman")
cor_perm <- cor.test(hydro_data$Permeability, hydro_data$median_2021, method = "spearman")

cat(sprintf("Spearman Correlation (IC vs Porosity)    : rho = %+.3f (p = %.3e)\n", cor_poro$estimate, cor_poro$p.value))
cat(sprintf("Spearman Correlation (IC vs Permeability): rho = %+.3f (p = %.3e)\n", cor_perm$estimate, cor_perm$p.value))

# Multivariable Linear Model
mod_multivar <- lm(median_2021 ~ Porosity + log10(Permeability) + Regime + Geology_Class, data = hydro_data)
cat("\n--- Multivariable Regression Model Summary ---\n")
print(summary(mod_multivar))

# Panel 3A: IC vs Porosity
p3a <- ggplot(hydro_data, aes(x = Porosity, y = median_2021, color = Regime)) +
  geom_point(aes(shape = Geology_Class), size = 2.8, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Regime), alpha = 0.2, size = 1.1) +
  scale_color_manual(values = regime_palette) +
  scale_fill_manual(values = regime_palette) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("Spearman rho = %+.3f\np = %.3e", cor_poro$estimate, cor_poro$p.value),
           size = 4, fontface = "bold.italic") +
  labs(
    title = "(a) IC vs. Subsurface Porosity",
    x = "Geological Porosity [%]",
    y = expression("Median IC"[2021] ~ "[-]")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

# Panel 3B: IC vs Permeability (Log10 scale)
p3b <- ggplot(hydro_data, aes(x = Permeability, y = median_2021, color = Regime)) +
  geom_point(aes(shape = Geology_Class), size = 2.8, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Regime), alpha = 0.2, size = 1.1) +
  scale_x_log10() +
  scale_color_manual(values = regime_palette) +
  scale_fill_manual(values = regime_palette) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("Spearman rho = %+.3f\np = %.3e", cor_perm$estimate, cor_perm$p.value),
           size = 4, fontface = "bold.italic") +
  labs(
    title = expression(bold("(b) IC vs. Geological Permeability (Log"[10] * " Scale)")),
    x = expression("Permeability [m"^2 * " or Darcy]"),
    y = expression("Median IC"[2021] ~ "[-]")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

fig3_combined <- plot_grid(p3a, p3b, nrow = 2, align = "v")
ggsave(file.path(out_dir, "Figure3_IC_vs_Porosity_Permeability.png"), plot = fig3_combined, width = 10, height = 9, dpi = 300)

cat("\n==================================================\n")
cat("All analyses executed and figures saved successfully!\n")
cat(sprintf("Output files located in: %s\n", out_dir))



###################################################################################################################
# COMPREHENSIVE STATISTICAL TESTING SUITE (ALL TESTS INCORPORATED)
###################################################################################################################

cat("\n====================================================================\n")
cat("               COMPREHENSIVE STATISTICAL ANALYSIS REPORT             \n")
cat("====================================================================\n\n")

# ------------------------------------------------------------------------------
# TEST 1: TEMPORAL DYNAMICS (1990 vs. 2021 OVERALL)
# ------------------------------------------------------------------------------
# 1A. Paired Wilcoxon Signed-Rank Test (Do catchments shift over 31 years?)
paired_test <- wilcox.test(catchment_data$median_2021, catchment_data$median_1990, paired = TRUE)

# 1B. Two-Sample Kolmogorov-Smirnov Test (Did the distribution shape shift?)
ks_temporal <- ks.test(catchment_data$median_1990, catchment_data$median_2021)

cat("--- 1. TEMPORAL DYNAMICS (1990 vs. 2021 ACROSS ALL BASINS) ---\n")
cat(sprintf("  • 1990 Regional Median IC  : %.3f (Mean: %.3f)\n", 
            median(catchment_data$median_1990), mean(catchment_data$mean_1990)))
cat(sprintf("  • 2021 Regional Median IC  : %.3f (Mean: %.3f)\n", 
            median(catchment_data$median_2021), mean(catchment_data$mean_2021)))
cat(sprintf("  • Mean Delta IC Shift      : %+.4f\n", mean(catchment_data$delta_IC)))
cat(sprintf("  • Paired Wilcoxon p-value  : %.4e %s\n", 
            paired_test$p.value, ifelse(paired_test$p.value < 0.05, "[SIGNIFICANT]", "[NOT SIG.]")))
cat(sprintf("  • Two-Sample KS p-value    : %.4e (D = %.3f)\n\n", 
            ks_temporal$p.value, ks_temporal$statistic))


# ------------------------------------------------------------------------------
# TEST 2: REGIME-WISE COMPARISONS (REGIME 0 vs. REGIME 1)
# ------------------------------------------------------------------------------
# 2A. Mann-Whitney U test on Delta IC (Does one regime degrade/green faster?)
delta_regime_test <- wilcox.test(delta_IC ~ Regime, data = catchment_data)

# 2B. Mann-Whitney U test on static IC in 1990 and 2021
regime_90_test <- wilcox.test(median_1990 ~ Regime, data = catchment_data)
regime_21_test <- wilcox.test(median_2021 ~ Regime, data = catchment_data)

# 2C. KS-test on distribution shape differences between regimes in 2021
regimes <- levels(catchment_data$Regime)
r0_vals <- catchment_data$median_2021[catchment_data$Regime == regimes[1]]
r1_vals <- catchment_data$median_2021[catchment_data$Regime == regimes[2]]
ks_regime_2021 <- ks.test(r0_vals, r1_vals)

cat("--- 2. HYDROLOGICAL REGIME COMPARISONS (Regime 0 vs. Regime 1) ---\n")
cat(sprintf("  • Delta IC Shift p-value   : %.4e %s\n", 
            delta_regime_test$p.value, ifelse(delta_regime_test$p.value < 0.05, "[SIGNIFICANT]", "[NOT SIG.]")))
cat(sprintf("  • 1990 IC Regime Diff p-val: %.4e\n", regime_90_test$p.value))
cat(sprintf("  • 2021 IC Regime Diff p-val: %.4e\n", regime_21_test$p.value))
cat(sprintf("  • 2021 KS Distribution Diff: p = %.4e (D = %.3f)\n\n", 
            ks_regime_2021$p.value, ks_regime_2021$statistic))


# ------------------------------------------------------------------------------
# TEST 3: LITHOLOGICAL INFLUENCE (GEOLOGICAL CLASSES)
# ------------------------------------------------------------------------------
geol_data <- catchment_data %>% filter(!is.na(Geology_Class))

# 3A. Kruskal-Wallis Non-Parametric ANOVA
kw_2021 <- kruskal.test(median_2021 ~ Geology_Class, data = geol_data)
kw_diff <- kruskal.test(delta_IC ~ Geology_Class, data = geol_data)

cat("--- 3. LITHOLOGY TESTS (Across 4 Geological Classes) ---\n")
cat(sprintf("  • Kruskal-Wallis (2021 IC) : H = %.2f, p = %.4e %s\n", 
            kw_2021$statistic, kw_2021$p.value, ifelse(kw_2021$p.value < 0.05, "[SIGNIFICANT]", "[NOT SIG.]")))
cat(sprintf("  • Kruskal-Wallis (Delta IC): H = %.2f, p = %.4e %s\n", 
            kw_diff$statistic, kw_diff$p.value, ifelse(kw_diff$p.value < 0.05, "[SIGNIFICANT]", "[NOT SIG.]")))

# 3B. Post-Hoc Pairwise Comparisons (Holm adjusted)
cat("  • Pairwise Post-Hoc Comparisons (2021 IC):\n")
posthoc_geol <- pairwise.wilcox.test(geol_data$median_2021, geol_data$Geology_Class, p.adjust.method = "holm")
print(posthoc_geol$p.value)

# 3C. Factorial Interaction Test: Does Geology effect depend on Regime?
# Using two-way factorial ANOVA on ranks (Scheirer-Ray-Hare / LM interaction)
mod_interaction <- lm(rank(median_2021) ~ Geology_Class * Regime, data = geol_data)
cat("\n  • Two-Way Interaction (Geology × Regime):\n")
print(anova(mod_interaction))


# ------------------------------------------------------------------------------
# TEST 4: HYDROGEOLOGICAL CORRELATIONS & MULTIVARIABLE REGRESSION
# ------------------------------------------------------------------------------
hydro_data <- catchment_data %>% filter(!is.na(Porosity) & !is.na(Permeability))

# 4A. Overall Spearman Correlations
cor_poro_all <- cor.test(hydro_data$Porosity, hydro_data$median_2021, method = "spearman")
cor_perm_all <- cor.test(hydro_data$Permeability, hydro_data$median_2021, method = "spearman")

# 4B. Stratified Correlations by Regime
cor_poro_r0 <- cor.test(hydro_data$Porosity[hydro_data$Regime == regimes[1]], 
                        hydro_data$median_2021[hydro_data$Regime == regimes[1]], method = "spearman")
cor_poro_r1 <- cor.test(hydro_data$Porosity[hydro_data$Regime == regimes[2]], 
                        hydro_data$median_2021[hydro_data$Regime == regimes[2]], method = "spearman")

cor_perm_r0 <- cor.test(hydro_data$Permeability[hydro_data$Regime == regimes[1]], 
                        hydro_data$median_2021[hydro_data$Regime == regimes[1]], method = "spearman")
cor_perm_r1 <- cor.test(hydro_data$Permeability[hydro_data$Regime == regimes[2]], 
                        hydro_data$median_2021[hydro_data$Regime == regimes[2]], method = "spearman")

cat("\n--- 4. HYDROGEOLOGICAL PROPERTY CORRELATIONS (Spearman rho) ---\n")
cat(sprintf("  • Overall Porosity vs IC      : rho = %+.3f (p = %.3e)\n", cor_poro_all$estimate, cor_poro_all$p.value))
cat(sprintf("    - Within Regime %s           : rho = %+.3f (p = %.3e)\n", regimes[1], cor_poro_r0$estimate, cor_poro_r0$p.value))
cat(sprintf("    - Within Regime %s           : rho = %+.3f (p = %.3e)\n", regimes[2], cor_poro_r1$estimate, cor_poro_r1$p.value))
cat(sprintf("  • Overall Permeability vs IC  : rho = %+.3f (p = %.3e)\n", cor_perm_all$estimate, cor_perm_all$p.value))
cat(sprintf("    - Within Regime %s           : rho = %+.3f (p = %.3e)\n", regimes[1], cor_perm_r0$estimate, cor_perm_r0$p.value))
cat(sprintf("    - Within Regime %s           : rho = %+.3f (p = %.3e)\n", regimes[2], cor_perm_r1$estimate, cor_perm_r1$p.value))

# 4C. Multivariable Linear Model (Controlling for lithology, porosity, permeability, and regime)
mod_multivar <- lm(median_2021 ~ Porosity + log10(Permeability) + Regime + Geology_Class, data = hydro_data)
cat("\n--- 5. MULTIVARIABLE REGRESSION MODEL SUMMARY ---\n")
print(summary(mod_multivar))
cat("====================================================================\n\n")

# ------------------------------------------------------------------------------
# SAVE COMPLETE STATISTICAL SUMMARY TABLE TO CSV
# ------------------------------------------------------------------------------
stats_export <- data.frame(
  Test_Category   = c("Temporal Shift", "Temporal Shape Shift", 
                      "Regime Delta Shift", "Regime 2021 Difference", "Regime 2021 Shape Shift",
                      "Geology Influence (2021)", "Geology Delta Influence", "Geology x Regime Interaction",
                      "Porosity Correlation (All)", "Permeability Correlation (All)"),
  Test_Name       = c("Paired Wilcoxon Signed-Rank", "Two-Sample Kolmogorov-Smirnov",
                      "Mann-Whitney U on Delta IC", "Mann-Whitney U on 2021 IC", "Two-Sample KS Test",
                      "Kruskal-Wallis ANOVA", "Kruskal-Wallis ANOVA", "Two-Way Factorial ANOVA on Ranks",
                      "Spearman Rank Correlation", "Spearman Rank Correlation"),
  Statistic       = c(paste("V =", round(paired_test$statistic, 2)), paste("D =", round(ks_temporal$statistic, 3)),
                      paste("W =", round(delta_regime_test$statistic, 2)), paste("W =", round(regime_21_test$statistic, 2)),
                      paste("D =", round(ks_regime_2021$statistic, 3)), paste("H =", round(kw_2021$statistic, 2)),
                      paste("H =", round(kw_diff$statistic, 2)), paste("F =", round(anova(mod_interaction)["Geology_Class:Regime", "F value"], 2)),
                      paste("rho =", round(cor_poro_all$estimate, 3)), paste("rho =", round(cor_perm_all$estimate, 3))),
  p_value         = c(paired_test$p.value, ks_temporal$p.value,
                      delta_regime_test$p.value, regime_21_test$p.value, ks_regime_2021$p.value,
                      kw_2021$p.value, kw_diff$p.value, anova(mod_interaction)["Geology_Class:Regime", "Pr(>F)"],
                      cor_poro_all$p.value, cor_perm_all$p.value),
  Interpretation  = c(ifelse(paired_test$p.value < 0.05, "Significant regional shift", "No significant shift"),
                      ifelse(ks_temporal$p.value < 0.05, "Distribution shape shifted", "Shape unchanged"),
                      ifelse(delta_regime_test$p.value < 0.05, "Regimes change at different rates", "Regimes change similarly"),
                      ifelse(regime_21_test$p.value < 0.05, "Regimes have different median IC", "Regimes have similar median IC"),
                      ifelse(ks_regime_2021$p.value < 0.05, "Regime distribution shapes differ", "Regime shapes similar"),
                      ifelse(kw_2021$p.value < 0.05, "Lithology significantly controls IC", "Lithology has no significant effect"),
                      ifelse(kw_diff$p.value < 0.05, "Lithology controls rate of IC change", "No lithological difference in change"),
                      ifelse(anova(mod_interaction)["Geology_Class:Regime", "Pr(>F)"] < 0.05, "Significant interaction", "Independent effects"),
                      ifelse(cor_poro_all$p.value < 0.05, "Significant correlation", "No correlation"),
                      ifelse(cor_perm_all$p.value < 0.05, "Significant correlation", "No correlation"))
)

write.csv(stats_export, file.path(out_dir, "Comprehensive_Statistical_Hypothesis_Tests.csv"), row.names = FALSE)
cat(sprintf("Statistical test matrix exported to: %s\n", file.path(out_dir, "Comprehensive_Statistical_Hypothesis_Tests.csv")))


# ==============================================================================
# TEMPORAL SHIFT ANALYSIS: 1990 vs. 2021 (ALL CATCHMENTS COMBINED)
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

cat("\nGenerating Paired Temporal Shift Violin Plot (1990 vs. 2021)...\n")

# 1. Reshape catchment data from wide to long format
paired_long <- catchment_data %>%
  dplyr::select(gauge_id, median_1990, median_2021, delta_IC) %>%
  tidyr::pivot_longer(
    cols      = c(median_1990, median_2021),
    names_to  = "Year",
    values_to = "Median_IC"
  ) %>%
  dplyr::mutate(
    Year = factor(ifelse(Year == "median_1990", "1990", "2021"), levels = c("1990", "2021"))
  )

# 2. Paired Statistical Metrics
paired_test   <- wilcox.test(catchment_data$median_2021, catchment_data$median_1990, paired = TRUE)
mean_shift    <- mean(catchment_data$delta_IC)
median_shift  <- median(catchment_data$delta_IC)
n_catchments  <- nrow(catchment_data)

# Mean and median for each year
stats_90 <- c(mean = mean(catchment_data$median_1990), med = median(catchment_data$median_1990))
stats_21 <- c(mean = mean(catchment_data$median_2021), med = median(catchment_data$median_2021))

# Year color palette (Teal vs. Orange-Red)
year_palette <- c("1990" = "#1b9e77", "2021" = "#d95f02")

# 3. Create the Paired Violin Plot
p_temporal_violin <- ggplot(paired_long, aes(x = Year, y = Median_IC, fill = Year)) +
  # A. Thin gray lines connecting paired catchments from 1990 to 2021
  geom_line(aes(group = gauge_id), color = "gray65", alpha = 0.35, size = 0.35) +
  # B. Violin shape showing distribution density
  geom_violin(alpha = 0.45, trim = FALSE, color = "gray30", width = 0.7) +
  # C. Boxplot overlaid inside the violin
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "black", size = 0.6) +
  # D. Jittered catchment points
  geom_point(aes(color = Year), position = position_jitter(width = 0.05, seed = 42), 
             size = 2, alpha = 0.65) +
  scale_fill_manual(values = year_palette) +
  scale_color_manual(values = year_palette) +
  # Annotate statistical hypothesis test and key metrics
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25,
           label = sprintf("Paired Wilcoxon Test: p = %.3e %s\nN = %d Catchments\nMean Net Shift (Δ): %+.3f\nMedian Net Shift (Δ): %+.3f",
                           paired_test$p.value,
                           ifelse(paired_test$p.value < 0.05, "(Significant)", "(Not Sig.)"),
                           n_catchments, mean_shift, median_shift),
           size = 4, fontface = "bold.italic", color = "gray20") +
  # Add mean/median markers on top of the x-axis labels
  annotate("text", x = 1, y = -Inf, vjust = -1.2, 
           label = sprintf("Mean: %.2f\nMed: %.2f", stats_90["mean"], stats_90["med"]), 
           size = 3.6, fontface = "bold", color = "#1b9e77") +
  annotate("text", x = 2, y = -Inf, vjust = -1.2, 
           label = sprintf("Mean: %.2f\nMed: %.2f", stats_21["mean"], stats_21["med"]), 
           size = 3.6, fontface = "bold", color = "#d95f02") +
  labs(
    title = "Regional Shift in Sediment Connectivity (1990 vs. 2021)",
    subtitle = "Connecting lines represent individual catchment trajectories | Boxplots show median & IQR",
    x = "Year",
    y = expression("Catchment Median IC [-]")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 16),
    plot.subtitle     = element_text(size = 11, color = "gray30"),
    axis.text.x       = element_text(face = "bold", size = 13)
  )

# Display the plot
print(p_temporal_violin)

# 4. Save High-Resolution Figure
out_violin_file <- file.path(out_dir, "Figure_Temporal_Shift_1990_vs_2021_Violin.png")
ggsave(out_violin_file, plot = p_temporal_violin, width = 7.5, height = 7, dpi = 300)

cat(sprintf("Paired violin plot saved successfully to: %s\n", out_violin_file))