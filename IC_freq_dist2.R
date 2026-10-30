# Title: Multi-Temporal and Regime-Wise Sediment Connectivity (IC) Analysis: 1990 vs. 2021
# Feature: Background individual catchment curves + 95% CI ribbons + Statistical tests

library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

###################################
# USER INPUTS & DIRECTORIES
###################################

# Directories containing output IC GeoTIFFs for each year
dir_1990 <- "F:/Peninsular India/Connectivity/IC/Output_UTM/1990"
dir_2021 <- "F:/Peninsular India/Connectivity/IC/Output_UTM/2021"

# Metadata CSV containing gauge_id and Regime
csv_path <- "D:/GitHub/Peninsular_Rivers/Files/monthly_avg_discharge.csv"

col_gauge_id <- "gauge_id"
col_regime   <- "Actual_regime"

# Output folder for saving the comparison figures and stats
out_dir <- "F:/Peninsular India/Connectivity/IC/Plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Standardized evaluation grid for densities
grid_n   <- 512
ic_min   <- -15
ic_max   <- 3
x_grid   <- seq(ic_min, ic_max, length.out = grid_n)
sample_n <- NULL  # Number of random pixels per catchment (keeps runtime fast)


###################################
# STEP 1: LOAD METADATA & RESTORE 5-DIGIT IDS
###################################
cat("Loading CSV metadata...\n")
meta_df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Restore leading zeroes to ensure 5-digit string matching
meta_df$gauge_id_clean <- sprintf("%05d", as.numeric(trimws(as.character(meta_df[[col_gauge_id]]))))
meta_df$regime_clean   <- as.factor(trimws(as.character(meta_df[[col_regime]])))

# Find available TIFF files
files_1990 <- list.files(dir_1990, pattern = "^IC_.*\\.tif$", full.names = TRUE)
files_2021 <- list.files(dir_2021, pattern = "^IC_.*\\.tif$", full.names = TRUE)

cat(sprintf("Found %d TIFFs for 1990 and %d TIFFs for 2021.\n", 
            length(files_1990), length(files_2021)))


###################################
# STEP 2: EXTRACT PIXELS & EVALUATE DENSITIES
###################################
density_list  <- list()
basin_medians <- list()
counter       <- 0

cat("Extracting IC distributions for 1990 and 2021...\n")

for (i in 1:nrow(meta_df)) {
  gid    <- meta_df$gauge_id_clean[i]
  regime <- meta_df$regime_clean[i]
  
  idx_1990 <- grep(paste0("IC_.*", gid, "\\.tif$"), files_1990)
  idx_2021 <- grep(paste0("IC_.*", gid, "\\.tif$"), files_2021)
  
  # Only process catchments that have successful runs in both years
  if (length(idx_1990) > 0 && length(idx_2021) > 0) {
    
    # Process Year 1990
    r_90 <- terra::rast(files_1990[idx_1990[1]])
    v_90 <- terra::values(r_90, mat = FALSE)
    c_90 <- v_90[!is.na(v_90) & !is.nan(v_90) & !is.infinite(v_90)]
    
    # Process Year 2021
    r_21 <- terra::rast(files_2021[idx_2021[1]])
    v_21 <- terra::values(r_21, mat = FALSE)
    c_21 <- v_21[!is.na(v_21) & !is.nan(v_21) & !is.infinite(v_21)]
    
    if (length(c_90) > 50 && length(c_21) > 50) {
      counter <- counter + 1
      
    # Catchment summary statistics for statistical testing
      basin_medians[[counter]] <- data.frame(
        gauge_id    = gid,
        Regime      = regime,
        mean_1990   = mean(c_90),
        mean_2021   = mean(c_21),
        median_1990 = median(c_90),
        median_2021 = median(c_21),
        delta_IC    = mean(c_21) - mean(c_90) # or median(c_21) - median(c_90)
      )
      
      # Subsample if necessary
      if (!is.null(sample_n)) {
        set.seed(42)
        if (length(c_90) > sample_n) c_90 <- sample(c_90, sample_n)
        if (length(c_21) > sample_n) c_21 <- sample(c_21, sample_n)
      }
      
      d_90 <- density(c_90, from = ic_min, to = ic_max, n = grid_n)
      d_21 <- density(c_21, from = ic_min, to = ic_max, n = grid_n)
      
      density_list[[paste0(gid, "_1990")]] <- data.frame(
        gauge_id = gid,
        Year     = "1990",
        Regime   = regime,
        x        = d_90$x,
        y        = d_90$y
      )
      
      density_list[[paste0(gid, "_2021")]] <- data.frame(
        gauge_id = gid,
        Year     = "2021",
        Regime   = regime,
        x        = d_21$x,
        y        = d_21$y
      )
    }
  }
}

densities_df <- dplyr::bind_rows(density_list)
medians_df   <- dplyr::bind_rows(basin_medians)
cat(sprintf("Successfully matched %d paired catchments between 1990 and 2021.\n", nrow(medians_df)))


###################################
# STEP 3: STATISTICAL HYPOTHESIS TESTING
###################################
paired_test <- wilcox.test(medians_df$median_2021, medians_df$median_1990, paired = TRUE)
mean_change <- mean(medians_df$delta_IC)
median_change <- median(medians_df$delta_IC)

regimes <- levels(medians_df$Regime)
test_regime_1990 <- wilcox.test(median_1990 ~ Regime, data = medians_df)
test_regime_2021 <- wilcox.test(median_2021 ~ Regime, data = medians_df)

# Overall regional mean IC values
overall_mean_1990 <- mean(medians_df$mean_1990)
overall_mean_2021 <- mean(medians_df$mean_2021)
overall_delta_mean <- overall_mean_2021 - overall_mean_1990

# Also compute overall regional medians
overall_med_1990  <- median(medians_df$median_1990)
overall_med_2021  <- median(medians_df$median_2021)

cat("\n=================== STATISTICAL ANALYSIS RESULTS ===================\n")
cat(sprintf("1. Regional Mean IC Values:\n"))
cat(sprintf("   - 1990 Regional Mean IC      : %.4f (Median: %.4f)\n", overall_mean_1990, overall_med_1990))
cat(sprintf("   - 2021 Regional Mean IC      : %.4f (Median: %.4f)\n", overall_mean_2021, overall_med_2021))
cat(sprintf("   - Net Change (2021 - 1990)   : %+.4f\n", overall_delta_mean))
cat(sprintf("   - Paired Wilcoxon p-value    : %.4e %s\n\n", 
            paired_test$p.value, ifelse(paired_test$p.value < 0.05, "(SIGNIFICANT)", "(NOT SIGNIFICANT)")))
# 
# cat("\n=================== STATISTICAL ANALYSIS RESULTS ===================\n")
# cat(sprintf("1. Temporal Shift (1990 vs 2021 across all catchments):\n"))
# cat(sprintf("   - Mean Change in Median IC   : %+.4f\n", mean_change))
# cat(sprintf("   - Median Change in Median IC : %+.4f\n", median_change))
# cat(sprintf("   - Paired Wilcoxon p-value    : %.4e %s\n\n", 
#             paired_test$p.value, ifelse(paired_test$p.value < 0.05, "(SIGNIFICANT)", "(NOT SIGNIFICANT)")))

cat(sprintf("2. Regime Differences (Regime %s vs Regime %s):\n", regimes[1], regimes[2]))
cat(sprintf("   - 1990 Mann-Whitney U p-value: %.4e %s\n", 
            test_regime_1990$p.value, ifelse(test_regime_1990$p.value < 0.05, "(SIGNIFICANT)", "(NOT SIGNIFICANT)")))
cat(sprintf("   - 2021 Mann-Whitney U p-value: %.4e %s\n", 
            test_regime_2021$p.value, ifelse(test_regime_2021$p.value < 0.05, "(SIGNIFICANT)", "(NOT SIGNIFICANT)")))
cat("====================================================================\n\n")

write.csv(medians_df, file.path(out_dir, "Catchment_Median_Comparison_1990_2021.csv"), row.names = FALSE)


###################################
# STEP 4: PLOT 1 — 1990 vs. 2021 (WITH INDIVIDUAL GRAY CURVES)
###################################
ci_year <- densities_df %>%
  group_by(Year, x) %>%
  summarize(
    mean_y   = mean(y),
    se_y     = sd(y) / sqrt(n()),
    ci_lower = pmax(0, mean_y - 1.96 * se_y),
    ci_upper = mean_y + 1.96 * se_y,
    .groups  = "drop"
  )

year_colors <- c("1990" = "#1b9e77", "2021" = "#d95f02") # Teal vs Orange-Red

plot_task1 <- ggplot() +
  # 1. Background individual catchment lines in light-gray
  geom_line(data = densities_df, 
            aes(x = x, y = y, group = interaction(gauge_id, Year)), 
            color = "gray75", alpha = 0.22, size = 0.3) +
  # 2. Shaded 95% Confidence Interval band only
  geom_ribbon(data = ci_year, 
              aes(x = x, ymin = ci_lower, ymax = ci_upper, fill = Year), 
              alpha = 0.28, color = NA) +
  # 3. Dashed CI boundary lines
  geom_line(data = ci_year, 
            aes(x = x, y = ci_lower, color = Year), 
            linetype = "dashed", size = 0.55) +
  geom_line(data = ci_year, 
            aes(x = x, y = ci_upper, color = Year), 
            linetype = "dashed", size = 0.55) +
  # 4. Bold primary mean curves
  geom_line(data = ci_year, 
            aes(x = x, y = mean_y, color = Year), 
            size = 1.35) +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values = year_colors, name = "Year") +
  # annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
  #          label = sprintf("Catchments: %d\nPaired Wilcoxon: p = %.3e\nMean Delta IC (2021 - 1990): %+.3f", 
  #                          nrow(medians_df), paired_test$p.value, mean_change),
  #          size = 4.2, fontface = "bold", color = "gray20") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("Catchments: %d\n1990 Mean IC: %.3f\n2021 Mean IC: %.3f\nNet Shift (Δ): %+.3f\nPaired p-value: %.3e", 
                           nrow(medians_df), overall_mean_1990, overall_mean_2021, overall_delta_mean, paired_test$p.value),
           size = 4.2, fontface = "bold", color = "gray20") +
  labs(
    title = "Regional Shift in Sediment Connectivity: 1990 vs. 2021",
    subtitle = "Gray curves: Individual Catchments | Solid lines: Mean density | Shaded bands: 95% CIs",
    x = "Index of Connectivity (IC) [-]",
    y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = c(0.85, 0.8),
    legend.background = element_rect(color = "gray80"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 16)
  )

file_task1 <- file.path(out_dir, "Comparison_Task1_1990_vs_2021_AllStations.png")
ggsave(file_task1, plot = plot_task1, width = 9.5, height = 6.5, dpi = 300)
cat(sprintf("Task 1 plot saved to: %s\n", file_task1))


###################################
# STEP 5: PLOT 2 — REGIME-WISE IN A SINGLE COMBINED FIGURE
###################################
ci_year_regime <- densities_df %>%
  group_by(Year, Regime, x) %>%
  summarize(
    mean_y   = mean(y),
    se_y     = sd(y) / sqrt(n()),
    ci_lower = pmax(0, mean_y - 1.96 * se_y),
    ci_upper = mean_y + 1.96 * se_y,
    .groups  = "drop"
  )

ci_year_regime$Group <- factor(paste0(ci_year_regime$Year, " - Regime ", ci_year_regime$Regime))

# Distinct, bold palette for the 4 combinations
group_colors <- c(
  "1990 - Regime 0" = "#4575b4", # Dark Blue
  "1990 - Regime 1" = "#91bfdb", # Light Blue
  "2021 - Regime 0" = "#d73027", # Deep Red
  "2021 - Regime 1" = "#fc8d59"  # Orange
)

plot_task2 <- ggplot() +
  # 1. Background individual catchment lines in light-gray
  geom_line(data = densities_df, 
            aes(x = x, y = y, group = interaction(gauge_id, Year, Regime)), 
            color = "gray80", alpha = 0.18, size = 0.25) +
  # 2. Shaded 95% CI bands
  geom_ribbon(data = ci_year_regime, 
              aes(x = x, ymin = ci_lower, ymax = ci_upper, fill = Group), 
              alpha = 0.22, color = NA) +
  # 3. Upper & Lower dashed CI bounds
  geom_line(data = ci_year_regime, 
            aes(x = x, y = ci_lower, color = Group), 
            linetype = "dashed", size = 0.5) +
  geom_line(data = ci_year_regime, 
            aes(x = x, y = ci_upper, color = Group), 
            linetype = "dashed", size = 0.5) +
  # 4. Bold mean curves
  geom_line(data = ci_year_regime, 
            aes(x = x, y = mean_y, color = Group, linetype = Regime), 
            size = 1.3) +
  scale_color_manual(values = group_colors, name = "Year & Regime") +
  scale_fill_manual(values = group_colors, name = "Year & Regime") +
  scale_linetype_manual(values = c("solid", "longdash"), name = "Regime") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("Regime Effect (Mann-Whitney U):\n1990: p = %.3e\n2021: p = %.3e", 
                           test_regime_1990$p.value, test_regime_2021$p.value),
           size = 4, fontface = "italic", color = "gray20") +
  labs(
    title = "Spatio-Temporal Connectivity Dynamics by Hydrological Regime",
    subtitle = "Gray curves: Catchments | Solid/dashed lines: Mean densities | Shaded bands: 95% CIs",
    x = "Index of Connectivity (IC) [-]",
    y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = c(0.82, 0.72),
    legend.background = element_rect(color = "gray80"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 16)
  )

file_task2 <- file.path(out_dir, "Comparison_Task2_Regime_Dynamics_SinglePlot.png")
ggsave(file_task2, plot = plot_task2, width = 10, height = 7, dpi = 300)
cat(sprintf("Task 2 plot saved to: %s\n", file_task2))

cat("\nAll tasks completed successfully!\n")