# Title: Frequency Distribution Analysis of IC Pixel Values by Regime
# Description: Reads batch IC GeoTIFFs, restores truncated 5-digit gauge IDs, 
#              and generates distribution plots for each regime.

# 1. Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(cowplot)

###################################
# USER INPUTS & DIRECTORIES SETUP
###################################

# Directory containing your batch output IC GeoTIFFs
ic_dir   <- "F:/Peninsular India/Connectivity/IC/Output/1990"

# Path to your metadata CSV containing gauge_id and Regime
csv_path <- "D:/GitHub/Peninsular_Rivers/Files/monthly_avg_discharge.csv"

# Column names in your CSV
col_gauge_id <- "gauge_id"
col_regime   <- "Actual_regime"

# Output folder for saving the plots
plot_out_dir <- file.path(ic_dir, "Distribution_Plots")
if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)

# Pixel Sampling per catchment (for memory efficiency & speed)
# Set to e.g., 25000 to sample 25k pixels per catchment (recommended for speed)
# Set to NULL if you want to use 100% of all pixels from every catchment
sample_pixels <- NULL 


###################################
# STEP 1: LOAD CSV & FIX LEADING ZEROES
###################################
cat("Loading CSV and restoring 5-digit gauge IDs...\n")
meta_df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Check columns exist
if (!col_gauge_id %in% names(meta_df) || !col_regime %in% names(meta_df)) {
  stop("Specified column names not found in CSV! Check col_gauge_id and col_regime.")
}

# CRITICAL FIX: Format gauge_id to 5 digits with leading zeroes
meta_df$gauge_id_clean <- sprintf("%05d", as.numeric(trimws(as.character(meta_df[[col_gauge_id]]))))
meta_df$regime_clean   <- as.factor(trimws(as.character(meta_df[[col_regime]])))

cat(sprintf("Found %d records in CSV with Regimes: %s\n", 
            nrow(meta_df), paste(levels(meta_df$regime_clean), collapse = ", ")))


###################################
# STEP 2: MATCH TIFF FILES & EXTRACT PIXELS
###################################
cat("Scanning for IC GeoTIFF files...\n")
tif_files <- list.files(ic_dir, pattern = "^IC_.*\\.tif$", full.names = TRUE)

if (length(tif_files) == 0) {
  stop("No files starting with 'IC_' and ending in '.tif' found in the specified directory!")
}

# Container for all pixel values
pixel_data_list <- list()

cat("Extracting pixel values from catchments...\n")
matched_count <- 0

for (i in 1:nrow(meta_df)) {
  gid    <- meta_df$gauge_id_clean[i]
  regime <- meta_df$regime_clean[i]
  
  # Find matching TIFF file that contains this 5-digit ID
  # Matches e.g. IC_100m_5x5_03002.tif, IC_03002.tif, etc.
  match_idx <- grep(paste0("IC_.*", gid, "\\.tif$"), tif_files)
  
  if (length(match_idx) > 0) {
    matched_count <- matched_count + 1
    tif_path <- tif_files[match_idx[1]]
    
    # Read raster values using terra
    r <- terra::rast(tif_path)
    vals <- terra::values(r, mat = FALSE)
    
    # Filter out NA, NaN, and Infinite values
    vals_clean <- vals[!is.na(vals) & !is.nan(vals) & !is.infinite(vals)]
    
    if (length(vals_clean) > 0) {
      # Optional random subsampling to keep plotting fast and light
      if (!is.null(sample_pixels) && length(vals_clean) > sample_pixels) {
        set.seed(42)
        vals_clean <- sample(vals_clean, sample_pixels)
      }
      
      # Append to data frame list
      pixel_data_list[[matched_count]] <- data.frame(
        gauge_id = gid,
        Regime   = regime,
        IC_value = vals_clean,
        stringsAsFactors = FALSE
      )
    }
  }
}

cat(sprintf("Successfully matched and extracted pixels from %d / %d catchments.\n", 
            matched_count, nrow(meta_df)))

if (matched_count == 0) {
  stop("Could not match any CSV gauge IDs with GeoTIFF filenames! Check naming convention.")
}

# Combine into a single master dataframe
all_pixels <- dplyr::bind_rows(pixel_data_list)
rm(pixel_data_list); gc()


###################################
# STEP 3: SUMMARY STATISTICS PER REGIME
###################################
stats_summary <- all_pixels %>%
  group_by(Regime) %>%
  summarize(
    Count_Catchments = n_distinct(gauge_id),
    Total_Pixels     = n(),
    Mean_IC          = mean(IC_value),
    Median_IC        = median(IC_value),
    SD_IC            = sd(IC_value),
    Q25_IC           = quantile(IC_value, 0.25),
    Q75_IC           = quantile(IC_value, 0.75)
  )

print(stats_summary)
write.csv(stats_summary, file.path(plot_out_dir, "IC_Regime_Statistics.csv"), row.names = FALSE)


###################################
# STEP 4: PLOTTING FUNCTION
###################################
# Function to generate distribution plot for a specific regime
plot_regime_distribution <- function(data_subset, regime_label, color_theme) {
  
  # Calculate overall median and mean for annotation
  overall_mean   <- mean(data_subset$IC_value)
  overall_median <- median(data_subset$IC_value)
  num_basins     <- n_distinct(data_subset$gauge_id)
  
  p <- ggplot(data_subset, aes(x = IC_value)) +
    # Thin, semi-transparent density curves for individual catchments
    geom_density(aes(group = gauge_id), color = "gray60", alpha = 0.25, size = 0.3) +
    # Bold, prominent density curve representing the entire regime
    geom_density(color = color_theme, fill = color_theme, alpha = 0.35, size = 1.2) +
    # Vertical line indicating the regime median
    geom_vline(xintercept = overall_median, color = "black", linetype = "dashed", size = 0.8) +
    # Annotate summary text on the plot
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
             label = sprintf("Catchments: %d\nRegime Mean: %.2f\nRegime Median: %.2f", 
                             num_basins, overall_mean, overall_median),
             size = 4, fontface = "bold", color = "gray20") +
    labs(
      title = paste("IC Frequency Distribution — Regime", regime_label),
      subtitle = "Gray curves: Individual Catchments | Bold curve: Overall Regime",
      x = "Index of Connectivity (IC) [-]",
      y = "Kernel Density"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


###################################
# STEP 5: GENERATE & SAVE PLOTS
###################################
regimes <- levels(all_pixels$Regime)
regime_colors <- c("#2b5c8f", "#d95f02") # Blue for Regime 1, Orange for Regime 2

# Generate Plot 1
cat(sprintf("Generating plot for Regime %s...\n", regimes[1]))
data_r1 <- all_pixels %>% filter(Regime == regimes[1])
p1 <- plot_regime_distribution(data_r1, regimes[1], regime_colors[1])

# Save Plot 1
out_file1 <- file.path(plot_out_dir, paste0("IC_Distribution_Regime_", regimes[1], ".png"))
ggsave(out_file1, plot = p1, width = 8, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", out_file1))

# Generate Plot 2
cat(sprintf("Generating plot for Regime %s...\n", regimes[2]))
data_r2 <- all_pixels %>% filter(Regime == regimes[2])
p2 <- plot_regime_distribution(data_r2, regimes[2], regime_colors[2])

# Save Plot 2
out_file2 <- file.path(plot_out_dir, paste0("IC_Distribution_Regime_", regimes[2], ".png"))
ggsave(out_file2, plot = p2, width = 8, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", out_file2))

# Generate & Save Combined Comparison Plot (Side-by-Side)
cat("Generating combined side-by-side comparison figure...\n")
combined_plot <- plot_grid(p1, p2, ncol = 2, align = "h")
out_file_combined <- file.path(plot_out_dir, "IC_Distribution_Comparison_Both_Regimes.png")
ggsave(out_file_combined, plot = combined_plot, width = 16, height = 6.5, dpi = 300)
cat(sprintf("Saved Combined Plot: %s\n", out_file_combined))

cat("\n==================================================\n")
cat("All distribution plots generated and saved successfully!\n")
cat(sprintf("Output folder: %s\n", plot_out_dir))