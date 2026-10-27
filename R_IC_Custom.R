# It incorporates the exact methodology and critical fixes from the batch pipeline:
#   Replaces rgdal with terra: Uses modern spatial drivers and handles the PROJ conflict environment setup.
# Automated Reprojection to Metric Grid (target_crs at 100 m): Reprojects the DEM, NDVI, and targets into a projected coordinate system (default: UTM Zone 44N EPSG:32644).
# Upstream Pit-Filling (PitRemove): Automatically resolves all artificial depressions created during resampling before flow direction calculations, ensuring the stream network connects continuously to the outlet.
# Flow Accumulation & Stream Snapping: Runs AreaD8 and snaps the gauge/target point to the cell with the highest flow accumulation within a 
# 500
# m
# 500 m
# buffer, preventing the while-loop from terminating early on the riverbank.
# Cai et al. (2000) 
# C
# C
# -factor: Converts NDVI to vegetation cover percentage and clamps weights between 
# 0.001
# 0.001
# and 
# 1.0
# 1.0
# (avoiding negative values and division by zero).
# Modernized Longitudinal Profile Analysis: Updated without any rgdal dependencies (readOGR/writeOGR replaced by terra), compatible with flag_profile <- 1 or 0.



# Title: Custom R_IC v. 2.0 (Single Site: In-Code Reprojection, PitRemove, Snapping, & C-factor)
# Original Authors: Tommaso Baggio, Lorenzo Martini, Loris Torresani
# Updated for: terra integration, TauDEM PitRemove, stream-snapping, and Cai et al. C-factor

# 1. Resolve TauDEM PROJ environment conflict
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

# 2. Packages used
library(terra)
library(raster)
library(ggplot2)
library(cowplot)

# Allow tolerance for raster alignment checks
rasterOptions(tolerance = 0.05)

# Premise: all files will be overwritten if they already exist

###################################
# INPUTS & DIRECTORIES SETUP
###################################

# Set working / output directory (do not write "/" as last character)
wd <- "F:/Peninsular India/Connectivity/IC/Single_Output"
if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

# Input file paths
dem_path    <- "F:/Peninsular India/Connectivity/IC/Input/DEM/DEM_03002.tif"
ndvi_path   <- "F:/Peninsular India/Connectivity/IC/Input/NDVI/2021/NDVI_03002.tif"
target_path <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/shapefiles_catchment/gauge_stations.shp"

# Target Gauge ID (if using a multi-gauge shapefile; leave NA if shapefile contains only 1 feature)
gauge_id     <- "03002"
gauge_id_col <- "gauge_id"

# Coordinate System & Target Resolution
target_crs   <- "EPSG:32644" # WGS 84 / UTM Zone 44N
target_res   <- 100          # 100 meters

# Stream Snapping Radius (meters)
snap_radius  <- 500          # Search radius to snap gauge onto the main D8 stream

# Output cleanup flag: 0 = keep essential files (IC, Dup, Ddn); 1 = keep all intermediate files
flag_k <- 0

# Number of MPI cores for TauDEM
mpi_cores <- 4

###################################
# Optional: Longitudinal Analysis of IC (Torresani & Piton, 2021)
# 0 = do NOT run; 1 = RUN
flag_profile <- 0

# Input longitudinal profile points shapefile (only used if flag_profile == 1)
profile_path <- file.path(wd, "input_points.shp")

# Elevation thresholds to split longitudinal domains (descending order, or NA)
input_thresh <- c(1041)

# Automatically keep all intermediate files if longitudinal analysis is requested
if (flag_profile == 1) flag_k <- 1


###################################################################################################################
# START OF SCRIPT
###################################################################################################################
start_timeTOT <- Sys.time()

# --------------------------------------------------------------------------
# 1. Load Raw Inputs & Reproject to Metric Grid
# --------------------------------------------------------------------------
cat("1. Loading and reprojecting DEM & NDVI...\n")
dem_raw  <- terra::rast(dem_path)
ndvi_raw <- terra::rast(ndvi_path)

# Reproject DEM to target CRS at specified resolution
dem_proj_t <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)

# --------------------------------------------------------------------------
# 2. Pit-Filling (TauDEM PitRemove)
# --------------------------------------------------------------------------
cat("2. Running TauDEM PitRemove to resolve all surface depressions...\n")
system(paste("mpiexec -n", mpi_cores, "PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif"))

# Load the depressionless DEM as master grid
z_fill <- raster("dem_pitfilled.tif")
res    <- res(z_fill)[1]

# --------------------------------------------------------------------------
# 3. NDVI Reprojection & C-Factor Weighting (Cai et al., 2000)
# --------------------------------------------------------------------------
cat("3. Computing C-factor from NDVI...\n")
# Align NDVI directly to the pit-filled DEM template
wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
W_factor  <- raster(wf_proj_t)

vals <- getValues(W_factor)
vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
if (length(vals_clean) == 0) stop("Error: NDVI raster is empty (all NA) after reprojection!")

ndvi_min <- min(vals_clean)
ndvi_max <- max(vals_clean)

# Fractional vegetation cover in percent (0 to 100%)
if (ndvi_max == ndvi_min) {
  f_g <- W_factor * 0
} else {
  f_g <- ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
}

# Piecewise C-factor formulation
f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))

C_factor[f_g <= 0]     <- 1.0
C_factor[f_g > 78.3]   <- 0.001   # Minimum threshold prevents division-by-zero
C_factor[C_factor > 1] <- 1.0
C_factor[C_factor < 0.001] <- 0.001

W_factor <- C_factor
writeRaster(W_factor, "W_factor.tif", overwrite = TRUE)

# --------------------------------------------------------------------------
# 4. TauDEM D8 Flow Directions & Stream Snapping
# --------------------------------------------------------------------------
cat("4. Computing D8 Flow Directions and Stream Snapping...\n")
system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif"))
flow_dir_D8 <- raster("D8_flow_dir.tif")
slope_D8    <- raster("slope_D8.tif")

# Calculate flow accumulation for stream channel identification
system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc"))

# Load target shapefile
targets_all <- terra::vect(target_path)
if (is.na(terra::crs(targets_all)) || terra::crs(targets_all) == "") {
  terra::crs(targets_all) <- "EPSG:4326"
}

# Filter by gauge ID if provided, otherwise take the first feature
if (!is.na(gauge_id) && gauge_id_col %in% names(targets_all)) {
  shp_ids <- trimws(as.character(targets_all[[gauge_id_col]][[1]]))
  is_match <- (shp_ids == gauge_id) | (sprintf("%05d", suppressWarnings(as.numeric(shp_ids))) == gauge_id)
  target_feat <- targets_all[which(is_match)[1], ]
} else {
  target_feat <- targets_all[1, ]
}

target_feat <- terra::project(target_feat, target_crs)

# Snap target to cell with maximum flow accumulation within buffer
buf        <- terra::buffer(target_feat, snap_radius)
accD8_t    <- terra::rast("accD8.tif")
accD8_crop <- terra::crop(accD8_t, buf)
accD8_mask <- terra::mask(accD8_crop, buf)

acc_vals <- terra::values(accD8_mask)
max_idx  <- which.max(acc_vals)

if (length(max_idx) == 0 || is.na(max_idx)) {
  warning("No accumulation found within buffer; using original coordinates.")
  snapped_xy <- terra::geom(target_feat)[, c("x", "y")]
  max_acc <- 0
} else {
  snapped_xy <- terra::xyFromCell(accD8_mask, max_idx)
  max_acc    <- acc_vals[max_idx]
}

cat(sprintf("  Target snapped to X=%.1f, Y=%.1f (Max Accumulation: %.0f cells)\n", 
            snapped_xy[1, 1], snapped_xy[1, 2], max_acc))

# Release terra pointer to avoid file-locking
rm(accD8_t, accD8_crop, accD8_mask); gc()

# Create snapped SpatVector & rasterize
target_snapped <- terra::vect(snapped_xy, crs = target_crs)
target_snapped$burn <- 1

target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
target_rst   <- raster(target_rst_t)
target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))

# --------------------------------------------------------------------------
# 5. Downslope Component (Ddn) Upstream Traversal
# --------------------------------------------------------------------------
cat("5. Tracing Downslope Flow Length...\n")
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8 * target_rst

slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)

inv_CS <- 1 / (W_factor * slope_D8_r_m)
writeRaster(inv_CS, "inv_CS.tif", overwrite = TRUE)

dir8_m <- as.matrix(flow_dir_D8)
W_FL_m <- as.matrix(inv_CS)

# Boundary NA padding
coladd <- vector(mode = "numeric", length = nrow(dir8_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd, dir8_m, coladd); W_FL_m <- cbind(coladd, W_FL_m, coladd)
rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd, dir8_m, rowadd); W_FL_m <- rbind(rowadd, W_FL_m, rowadd)

Wgt_m   <- W_FL_m 
NA_data <- which(is.na(dir8_m), arr.ind = TRUE); NA_data <- as.data.frame(NA_data)
posrow  <- NA_data$row; poscol <- NA_data$col

pos1col_NA <- poscol - 1; pos1row_NA <- posrow
pos1row <- pos1row_NA[-pos1col_NA != 0]; pos1col <- pos1col_NA[-pos1col_NA != 0]
pos1dir <- dir8_m[cbind(pos1row, pos1col)]
pos1row_0 <- pos1row[which(pos1dir == 1)]; pos1col_0 <- pos1col[which(pos1dir == 1)]
W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0

pos2col_NA <- poscol - 1; pos2row_NA <- posrow + 1
pos2rowt <- pos2row_NA[-pos2col_NA != 0]; pos2colt <- pos2col_NA[-pos2col_NA != 0]
pos2row <- pos2rowt[pos2rowt <= (nrow(dir8_m))]; pos2col <- pos2colt[pos2rowt <= (nrow(dir8_m))]
pos2dir <- dir8_m[cbind(pos2row, pos2col)]
pos2row_0 <- pos2row[which(pos2dir == 2)]; pos2col_0 <- pos2col[which(pos2dir == 2)]
W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0

pos3col_NA <- poscol; pos3row_NA <- posrow + 1
pos3row <- pos3row_NA[pos3row_NA <= (nrow(dir8_m))]; pos3col <- pos3col_NA[pos3row_NA <= (nrow(dir8_m))]
pos3dir <- dir8_m[cbind(pos3row, pos3col)]
pos3row_0 <- pos3row[which(pos3dir == 3)]; pos3col_0 <- pos3col[which(pos3dir == 3)]
W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0

pos4col_NA <- poscol + 1; pos4row_NA <- posrow + 1
pos4rowt <- pos4row_NA[pos4col_NA <= (ncol(dir8_m))]; pos4colt <- pos4col_NA[pos4col_NA <= (ncol(dir8_m))]
pos4row <- pos4rowt[pos4rowt <= (nrow(dir8_m))]; pos4col <- pos4colt[pos4rowt <= (nrow(dir8_m))]
pos4dir <- dir8_m[cbind(pos4row, pos4col)]
pos4row_0 <- pos4row[which(pos4dir == 4)]; pos4col_0 <- pos4col[which(pos4dir == 4)]
W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0

pos5col_NA <- poscol + 1; pos5row_NA <- posrow
pos5row <- pos5row_NA[pos5col_NA <= (ncol(dir8_m))]; pos5col <- pos5col_NA[pos5col_NA <= (ncol(dir8_m))]
pos5dir <- dir8_m[cbind(pos5row, pos5col)]
pos5row_0 <- pos5row[which(pos5dir == 5)]; pos5col_0 <- pos5col[which(pos5dir == 5)]
W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0

pos6col_NA <- poscol + 1; pos6row_NA <- posrow - 1
pos6rowt <- pos6row_NA[pos6col_NA <= (ncol(dir8_m))]; pos6colt <- pos6col_NA[pos6col_NA <= (ncol(dir8_m))]
pos6row <- pos6rowt[-pos6rowt != 0]; pos6col <- pos6colt[-pos6rowt != 0]
pos6dir <- dir8_m[cbind(pos6row, pos6col)]
pos6row_0 <- pos6row[which(pos6dir == 6)]; pos6col_0 <- pos6col[which(pos6dir == 6)]
W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0

pos7col_NA <- poscol; pos7row_NA <- posrow - 1
pos7row <- pos7row_NA[-pos7row_NA != 0]; pos7col <- pos7col_NA[-pos7row_NA != 0]
pos7dir <- dir8_m[cbind(pos7row, pos7col)]
pos7row_0 <- pos7row[which(pos7dir == 7)]; pos7col_0 <- pos7col[which(pos7dir == 7)]
W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0 

pos8col_NA <- poscol - 1; pos8row_NA <- posrow - 1
pos8rowt <- pos8row_NA[-pos8row_NA != 0]; pos8colt <- pos8col_NA[-pos8row_NA != 0]
pos8row <- pos8rowt[-pos8colt != 0]; pos8col <- pos8colt[-pos8colt != 0]
pos8dir <- dir8_m[cbind(pos8row, pos8col)]
pos8row_0 <- pos8row[which(pos8dir == 8)]; pos8col_0 <- pos8col[which(pos8dir == 8)]
W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0

rm(pos1col, pos1col_NA, pos1row, pos1row_NA, pos1dir,
   pos3col, pos3col_NA, pos3row, pos3row_NA, pos3dir,
   pos5col, pos5col_NA, pos5row, pos5row_NA, pos5dir,
   pos7col, pos7col_NA, pos7row, pos7row_NA, pos7dir,
   pos2col, pos2col_NA, pos2row, pos2row_NA, pos2colt, pos2rowt, pos2dir,
   pos4col, pos4col_NA, pos4row, pos4row_NA, pos4colt, pos4rowt, pos4dir,
   pos6col, pos6col_NA, pos6row, pos6row_NA, pos6colt, pos6rowt, pos6dir,
   pos8col, pos8col_NA, pos8row, pos8row_NA, pos8colt, pos8rowt, pos8dir)

count <- 1
while ((length(pos1col_0) | length(pos2col_0) | length(pos3col_0) | length(pos4col_0) | 
        length(pos5col_0) | length(pos6col_0) | length(pos7col_0) | length(pos8col_0)) > 0) {
  
  pos1row_0_bis <- pos1row_0; pos1col_0_bis <- pos1col_0 + 1
  pos2row_0_bis <- pos2row_0 - 1; pos2col_0_bis <- pos2col_0 + 1
  pos3row_0_bis <- pos3row_0 - 1; pos3col_0_bis <- pos3col_0
  pos4row_0_bis <- pos4row_0 - 1; pos4col_0_bis <- pos4col_0 - 1
  pos5row_0_bis <- pos5row_0;     pos5col_0_bis <- pos5col_0 - 1
  pos6row_0_bis <- pos6row_0 + 1; pos6col_0_bis <- pos6col_0 - 1
  pos7row_0_bis <- pos7row_0 + 1; pos7col_0_bis <- pos7col_0
  pos8row_0_bis <- pos8row_0 + 1; pos8col_0_bis <- pos8col_0 + 1
  
  if (count == 1) {
    W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0; W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
    W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0; W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
    W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0; W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
    W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0; W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0
  } else {
    W_FL_m[cbind(pos1row_0, pos1col_0)] <- W_FL_m[cbind(pos1row_0_bis, pos1col_0_bis)] + (res * ((Wgt_m[cbind(pos1row_0_bis, pos1col_0_bis)] + Wgt_m[cbind(pos1row_0, pos1col_0)]) / 2))
    W_FL_m[cbind(pos2row_0, pos2col_0)] <- W_FL_m[cbind(pos2row_0_bis, pos2col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos2row_0_bis, pos2col_0_bis)] + Wgt_m[cbind(pos2row_0, pos2col_0)]) / 2))
    W_FL_m[cbind(pos3row_0, pos3col_0)] <- W_FL_m[cbind(pos3row_0_bis, pos3col_0_bis)] + (res * ((Wgt_m[cbind(pos3row_0_bis, pos3col_0_bis)] + Wgt_m[cbind(pos3row_0, pos3col_0)]) / 2))
    W_FL_m[cbind(pos4row_0, pos4col_0)] <- W_FL_m[cbind(pos4row_0_bis, pos4col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos4row_0_bis, pos4col_0_bis)] + Wgt_m[cbind(pos4row_0, pos4col_0)]) / 2))
    W_FL_m[cbind(pos5row_0, pos5col_0)] <- W_FL_m[cbind(pos5row_0_bis, pos5col_0_bis)] + (res * ((Wgt_m[cbind(pos5row_0_bis, pos5col_0_bis)] + Wgt_m[cbind(pos5row_0, pos5col_0)]) / 2))
    W_FL_m[cbind(pos6row_0, pos6col_0)] <- W_FL_m[cbind(pos6row_0_bis, pos6col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos6row_0_bis, pos6col_0_bis)] + Wgt_m[cbind(pos6row_0, pos6col_0)]) / 2))
    W_FL_m[cbind(pos7row_0, pos7col_0)] <- W_FL_m[cbind(pos7row_0_bis, pos7col_0_bis)] + (res * ((Wgt_m[cbind(pos7row_0_bis, pos7col_0_bis)] + Wgt_m[cbind(pos7row_0, pos7col_0)]) / 2))
    W_FL_m[cbind(pos8row_0, pos8col_0)] <- W_FL_m[cbind(pos8row_0_bis, pos8col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos8row_0_bis, pos8col_0_bis)] + Wgt_m[cbind(pos8row_0, pos8col_0)]) / 2))
  }
  
  pos0row_bis <- c(pos1row_0, pos2row_0, pos3row_0, pos4row_0, pos5row_0, pos6row_0, pos7row_0, pos8row_0)
  pos0col_bis <- c(pos1col_0, pos2col_0, pos3col_0, pos4col_0, pos5col_0, pos6col_0, pos7col_0, pos8col_0)
  
  poscolt <- pos0col_bis - 1; posrowt <- pos0row_bis
  posrow <- posrowt[-poscolt != 0]; poscol <- poscolt[-poscolt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos1row_0 <- posrow[which(posdir == 1)]; pos1col_0 <- poscol[which(posdir == 1)]
  
  poscolt <- pos0col_bis - 1; posrowt <- pos0row_bis + 1
  posrow <- posrowt[-poscolt != 0]; poscol <- poscolt[-poscolt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos2row_0 <- posrow[which(posdir == 2)]; pos2col_0 <- poscol[which(posdir == 2)]
  
  poscol <- pos0col_bis; posrow <- pos0row_bis + 1
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos3row_0 <- posrow[which(posdir == 3)]; pos3col_0 <- poscol[which(posdir == 3)]
  
  poscol <- pos0col_bis + 1; posrow <- pos0row_bis + 1
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos4row_0 <- posrow[which(posdir == 4)]; pos4col_0 <- poscol[which(posdir == 4)]
  
  poscol <- pos0col_bis + 1; posrow <- pos0row_bis
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos5row_0 <- posrow[which(posdir == 5)]; pos5col_0 <- poscol[which(posdir == 5)]
  
  poscolt <- pos0col_bis + 1; posrowt <- pos0row_bis - 1
  posrow <- posrowt[-posrowt != 0]; poscol <- poscolt[-posrowt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos6row_0 <- posrow[which(posdir == 6)]; pos6col_0 <- poscol[which(posdir == 6)]
  
  poscolt <- pos0col_bis; posrowt <- pos0row_bis - 1
  posrow <- posrowt[-posrowt != 0]; poscol <- poscolt[-posrowt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos7row_0 <- posrow[which(posdir == 7)]; pos7col_0 <- poscol[which(posdir == 7)]
  
  poscolt <- pos0col_bis - 1; posrowt <- pos0row_bis - 1
  posrowt1 <- posrowt[-posrowt != 0]; poscolt1 <- poscolt[-posrowt != 0]
  posrow <- posrowt1[-poscolt1 != 0]; poscol <- poscolt1[-poscolt1 != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos8row_0 <- posrow[which(posdir == 8)]; pos8col_0 <- poscol[which(posdir == 8)]
  
  count <- count + 1
}

cat(sprintf("  Flow length loop completed in %d iterations.\n", count))

# Remove matrix boundary padding
W_FL_m <- W_FL_m[-1, ]; W_FL_m <- W_FL_m[-nrow(W_FL_m), ]
dir8_m <- dir8_m[-1, ]; dir8_m <- dir8_m[-nrow(dir8_m), ]
W_FL_m <- W_FL_m[, -1]; W_FL_m <- W_FL_m[, -ncol(W_FL_m)] 
dir8_m <- dir8_m[, -1]; dir8_m <- dir8_m[, -ncol(dir8_m)]

crs_out <- CRS("+init=epsg:32644")
ext <- as.matrix(extent(inv_CS))
W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs_out)
writeRaster(W_FL_fin, "Weighted_Flow_length.tif", overwrite = TRUE)

dif_WFL <- W_FL_fin - inv_CS
dif_WFL[dif_WFL == 0] <- NA
mask_cells <- dif_WFL * 0 + 1

Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
crs(Ddn) <- crs_out
writeRaster(Ddn, "Ddn.tif", overwrite = TRUE)

rm(Wgt_m, W_FL_m, dir8_m, dif_WFL, W_FL, dir8, inv_CS, flow_dir_D8); gc()

# --------------------------------------------------------------------------
# 6. Upslope Component (Dup)
# --------------------------------------------------------------------------
cat("6. Computing Upslope Component (Dup)...\n")
system(paste("mpiexec -n", mpi_cores, "DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_pitfilled.tif"))
ang <- raster("ang.tif") * target_rst
writeRaster(ang, "ang.tif", overwrite = TRUE)
gc()

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc"))
accS <- raster("accS.tif")

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca sca.tif -nc"))
sca <- raster("sca.tif")

constant_res <- z_fill * 0 + res
Acc_fin      <- sca / constant_res
smean        <- (accS + slope_D8_r_m) / Acc_fin

rm(slope_D8_r_m, accS, target_rst); gc()

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc"))
accW <- raster("accW.tif") * res
cmean <- (accW + W_factor) / Acc_fin

Dup <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
crs(Dup) <- crs_out
writeRaster(Dup, "Dup.tif", overwrite = TRUE)

rm(constant_res, ang, Acc_fin, smean, accW, cmean); gc()

# --------------------------------------------------------------------------
# 7. Final Index of Connectivity (IC)
# --------------------------------------------------------------------------
cat("7. Calculating Final IC...\n")
IC <- log10(Dup / Ddn)
IC <- IC * mask_cells
crs(IC) <- crs_out

pal <- colorRampPalette(c("blue", "lightyellow1", "red"))

png("IC.png", width = 800, height = 800)
plot(IC, col = pal(100), main = "Index of Connectivity (IC)")
dev.off()

plot(IC, col = pal(100), main = "Index of Connectivity (IC)")
writeRaster(IC, paste0("IC_", res, "m.tif"), overwrite = TRUE)

# Cleanup intermediate files if flag_k == 0
if (flag_k == 0) {
  tmp_files <- c("accS.tif", "accW.tif", "ang.tif", "D8_flow_dir.tif", 
                 "sca.tif", "slope_D8.tif", "slope_D8_rec_m.tif", 
                 "slp.tif", "Weighted_Flow_length.tif", "W_factor.tif", 
                 "dem_raw_proj.tif", "dem_pitfilled.tif", "accD8.tif")
  suppressWarnings(file.remove(tmp_files[file.exists(tmp_files)]))
}

end_timeTOT <- Sys.time()
cat(sprintf("\nConnectivity run completed in %.2f minutes.\n", 
            as.numeric(difftime(end_timeTOT, start_timeTOT, units = "mins"))))


###################################################################################################################
# 8. Optional Longitudinal Profile Analysis (Torresani & Piton, 2021)
###################################################################################################################

if (flag_profile == 1) {
  cat("\nRunning Longitudinal Profile Analysis...\n")
  
  pointprofile_vect <- terra::vect(profile_path)
  pointprofile_vect <- terra::project(pointprofile_vect, target_crs)
  
  # Remove extraneous attribute columns
  if (ncol(pointprofile_vect) > 0) {
    pointprofile_vect <- pointprofile_vect[, 0]
  }
  
  system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accW8D.tif -nc"))
  Acc_fin <- raster("accW8D.tif")
  
  # Stack all metric raster layers
  rasStack <- stack(W_factor, Dup, Ddn, IC, z_fill, slope_D8, Acc_fin)
  
  # Extract components along points
  pointprofile_sp <- as(pointprofile_vect, "Spatial")
  pointprofile$POINT_X <- pointprofile_sp@coords[, 1]
  pointprofile$POINT_Y <- pointprofile_sp@coords[, 2]
  
  pointprofile2 <- raster::extract(rasStack, pointprofile_sp, method = "simple", df = TRUE, sp = TRUE)
  names(pointprofile2) <- c("POINT_X", "POINT_Y", "W", "Dup", "Ddn", "IC", "Z", "S8", "A")
  
  # Compute Abscissa (Distance along river channel)
  pointprofile2@data <- pointprofile2@data[order(pointprofile2@data$Z), ]
  pointprofile2@data$Absc[1] <- 0
  for (i in 2:length(pointprofile2@data$Z)) {
    pointprofile2@data$Absc[i] <- pointprofile2@data$Absc[i-1] + 
      ((pointprofile2@data$POINT_X[i] - pointprofile2@data$POINT_X[i-1])^2 + 
         (pointprofile2@data$POINT_Y[i] - pointprofile2@data$POINT_Y[i-1])^2)^0.5
  }
  
  # Split profile into domains
  pointprofile2@data$Domain <- NA
  if (length(input_thresh) == 1 && !is.na(input_thresh[1])) {
    pointprofile2@data$Domain[pointprofile2@data$Z > input_thresh[1]] <- "A"
    pointprofile2@data$Domain[pointprofile2@data$Z <= input_thresh[1]] <- "B"
  } else if (length(input_thresh) > 1) {
    pointprofile2@data$Domain[pointprofile2@data$Z > input_thresh[1]] <- LETTERS[1]
    for (i in 1:length(input_thresh)) {
      pointprofile2@data$Domain[pointprofile2@data$Z <= input_thresh[i]] <- LETTERS[i + 1]
    }
  } else {
    pointprofile2@data$Domain <- "A"
  }
  pointprofile2@data$Domain <- as.factor(pointprofile2@data$Domain)
  
  DATA <- as.data.frame(pointprofile2)
  write.csv(DATA, file = "DATA.csv", row.names = FALSE)
  terra::writeVector(terra::vect(pointprofile2), "point_profile_IC.shp", overwrite = TRUE)
  
  # Fit cubic polynomial
  lm_pol <- lm(DATA$IC ~ poly(DATA$Absc, 3, raw = TRUE))
  rsqr   <- round(summary(lm_pol)$adj.r.squared, 3)
  
  # Summary Longitudinal Plots
  number_ticks <- function(n) { function(limits) pretty(limits, n) }
  
  Profile_plot <- ggplot(DATA, aes(x = Absc, y = Z, colour = Domain)) + 
    geom_line() + scale_colour_brewer(palette = "Set1") + theme_bw() +
    theme(legend.position = c(.25, .8), legend.direction = "horizontal", text = element_text(size = 14)) +
    scale_x_continuous(breaks = number_ticks(10)) + labs(y = "Z [m.a.s.l]", x = "", tag = "(a)")
  
  Slope_plot <- ggplot(DATA, aes(x = Absc, y = S8, color = Domain)) + 
    geom_line() + scale_colour_brewer(palette = "Set1") + theme_bw() +
    theme(legend.position = "none", text = element_text(size = 14)) +
    scale_x_continuous(breaks = number_ticks(10)) + labs(x = "", y = "Slope [m/m]", tag = "(b)")
  
  Area_plot <- ggplot(DATA, aes(x = Absc, y = A / 10^6, color = Domain)) + 
    geom_line() + scale_colour_brewer(palette = "Set1") + theme_bw() +
    theme(legend.position = "none", text = element_text(size = 14)) +
    scale_x_continuous(breaks = number_ticks(10)) + labs(x = "", y = expression("A [km"^2~"]"), tag = "(c)")
  
  IC_plot <- ggplot(DATA, aes(x = Absc, y = IC, color = Domain)) + 
    geom_line() + scale_colour_brewer(palette = "Set1") + theme_bw() +
    theme(legend.position = "none", text = element_text(size = 14)) +
    scale_x_continuous(breaks = number_ticks(10)) + labs(x = "Distance [m]", y = "IC [-]", tag = "(d)") +
    geom_line(aes(x = DATA$Absc, y = predict(lm_pol, data.frame(x = DATA$Absc))), col = "black", lwd = 0.7) +
    annotate("text", x = mean(DATA$Absc), y = max(DATA$IC, na.rm = TRUE) * 0.9, 
             label = bquote("R"^2 == .(rsqr)), size = 3.5)
  
  png("Summary_plot.png", width = 18, height = 24, units = "cm", res = 300)
  print(cowplot::plot_grid(Profile_plot, Slope_plot, Area_plot, IC_plot, nrow = 4, align = "hv"))
  dev.off()
}

cat("\nCongratulations!! Script execution fully completed.\n")