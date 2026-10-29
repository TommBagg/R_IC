# Title: Custom R_IC v. 2.3 (Single Site: Dual Outlet Method, PitRemove, & Corrected 8-Dir Traversal)
# Author details: Tommaso Baggio, Lorenzo Martini, Loris Torresani

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

rasterOptions(tolerance = 0.05)

# Premise: all files will be overwritten if they already exist

###################################
# INPUTS & DIRECTORIES SETUP
###################################

# Set working / output directory
wd <- "F:/Peninsular India/Connectivity/IC/Single_Output"
if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

gauge_id <- "09012"

# Input file paths
dem_path    <- paste0("F:/Peninsular India/Connectivity/IC/Input/DEM/DEM_", gauge_id, ".tif")
ndvi_path   <- paste0("F:/Peninsular India/Connectivity/IC/Input/NDVI/2021/NDVI_", gauge_id, ".tif")
target_path <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/shapefiles_catchment/gauge_stations.shp"
gauge_id_col <- "gauge_id"

# Coordinate System & Target Resolution
target_crs   <- "EPSG:32644" # WGS 84 / UTM Zone 44N
target_res   <- 100          # 100 meters

# -------------------------------------------------------------------------
# OUTLET METHOD SELECTOR:
# "catchment_outlet" = Automatically finds the true river outlet (max accumulation)
# "gauge_snapping"   = Uses the gauge station shapefile point snapped within snap_radius
# -------------------------------------------------------------------------
outlet_method <- "catchment_outlet"
snap_radius   <- 2500 # Used only if outlet_method == "gauge_snapping"

# Output cleanup flag: 0 = keep essential files (IC, Dup, Ddn); 1 = keep all intermediate files
flag_k <- 0
mpi_cores <- 4

# Optional: Longitudinal Analysis of IC (Torresani & Piton, 2021)
flag_profile <- 0
profile_path <- file.path(wd, "input_points.shp")
input_thresh <- c(1041)

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

dem_proj_t <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)

# --------------------------------------------------------------------------
# 2. Pit-Filling Upfront (TauDEM PitRemove)
# --------------------------------------------------------------------------
cat("2. Running TauDEM PitRemove to resolve all surface depressions...\n")
system(paste("mpiexec -n", mpi_cores, "PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif"))

z_fill <- raster("dem_pitfilled.tif")
res    <- res(z_fill)[1]

# --------------------------------------------------------------------------
# 3. NDVI Reprojection & C-Factor Weighting (Cai et al., 2000)
# --------------------------------------------------------------------------
cat("3. Computing C-factor from NDVI...\n")
wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
W_factor  <- raster(wf_proj_t)

vals <- getValues(W_factor)
vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
if (length(vals_clean) == 0) stop("Error: NDVI raster is empty (all NA) after reprojection!")

ndvi_min <- min(vals_clean)
ndvi_max <- max(vals_clean)

if (ndvi_max == ndvi_min) {
  f_g <- W_factor * 0
} else {
  f_g <- ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
}

f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))

C_factor[f_g <= 0]     <- 1.0
C_factor[f_g > 78.3]   <- 0.001   # Minimum threshold prevents division-by-zero
C_factor[C_factor > 1] <- 1.0
C_factor[C_factor < 0.001] <- 0.001

W_factor <- C_factor
writeRaster(W_factor, "W_factor.tif", overwrite = TRUE)

# --------------------------------------------------------------------------
# 4. TauDEM D8 Flow Directions & Outlet Selection
# --------------------------------------------------------------------------
cat("4. Computing D8 Flow Directions & Setting Target...\n")
system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif"))
flow_dir_D8 <- raster("D8_flow_dir.tif")
slope_D8    <- raster("slope_D8.tif")

system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc"))
accD8 <- raster("accD8.tif")

if (outlet_method == "catchment_outlet") {
  # METHOD A: True catchment outlet via maximum D8 accumulation
  cat("  Using Method: Catchment Outlet (Max Flow Accumulation)...\n")
  max_cell   <- which.max(getValues(accD8))
  snapped_xy <- xyFromCell(accD8, max_cell)
  max_accum  <- max(getValues(accD8), na.rm = TRUE)
  
  cat(sprintf("  Outlet located at X=%.1f, Y=%.1f (Total Drainage: %d cells)\n", 
              snapped_xy[1], snapped_xy[2], as.integer(max_accum)))
  
} else if (outlet_method == "gauge_snapping") {
  # METHOD B: Gauge station shapefile with buffer snapping
  cat(sprintf("  Using Method: Gauge Station Snapping (Radius: %dm)...\n", snap_radius))
  all_targets <- terra::vect(target_path)
  if (is.na(terra::crs(all_targets)) || terra::crs(all_targets) == "") terra::crs(all_targets) <- "EPSG:4326"
  
  shp_ids <- trimws(as.character(all_targets[[gauge_id_col]][[1]]))
  is_match <- (shp_ids == gauge_id) | (sprintf("%05d", suppressWarnings(as.numeric(shp_ids))) == gauge_id)
  target_feat <- all_targets[which(is_match)[1], ]
  target_feat <- terra::project(target_feat, target_crs)
  
  buf        <- terra::buffer(target_feat, snap_radius)
  accD8_t    <- terra::rast("accD8.tif")
  accD8_crop <- terra::crop(accD8_t, buf)
  accD8_mask <- terra::mask(accD8_crop, buf)
  
  acc_vals <- terra::values(accD8_mask)
  max_idx  <- which.max(acc_vals)
  
  if (length(max_idx) == 0 || is.na(max_idx)) {
    warning("No accumulation found within buffer; using original coordinates.")
    snapped_xy <- terra::geom(target_feat)[, c("x", "y")]
    max_accum  <- 0
  } else {
    snapped_xy <- terra::xyFromCell(accD8_mask, max_idx)
    max_accum  <- acc_vals[max_idx]
  }
  
  cat(sprintf("  Gauge %s snapped to X=%.1f, Y=%.1f (Drainage: %d cells)\n", 
              gauge_id, snapped_xy[1], snapped_xy[2], as.integer(max_accum)))
  
  rm(accD8_t, accD8_crop, accD8_mask); gc()
}

# Create target point SpatVector & rasterize
target_snapped      <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs)
target_snapped$burn <- 1

target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
target_rst   <- raster(target_rst_t)
target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))

# --------------------------------------------------------------------------
# 5. Downslope Component (Ddn) Upstream Traversal
# --------------------------------------------------------------------------
cat("5. Tracing Downslope Flow Length from Outlet...\n")
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8 * target_rst

slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)

inv_CS <- 1 / (W_factor * slope_D8_r_m)
writeRaster(inv_CS, paste0("inv_CS_100m_5x5_", gauge_id, ".tif"), overwrite = TRUE)

# FIX: Find target index BEFORE adding matrix padding (prevents boundary artifact)
dir8_raw_m <- as.matrix(flow_dir_D8)
target_idx <- which(is.na(dir8_raw_m), arr.ind = TRUE)

coladd <- vector(mode = "numeric", length = nrow(dir8_raw_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd, dir8_raw_m, coladd)
W_FL_m <- cbind(coladd, as.matrix(inv_CS), coladd)

rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd, dir8_m, rowadd)
W_FL_m <- rbind(rowadd, W_FL_m, rowadd)

Wgt_m <- W_FL_m 

# Sinks to trace from: ONLY the true target cell, NOT the outer border!
posrow <- target_idx[, "row"] + 1
poscol <- target_idx[, "col"] + 1

# Initial search around the outlet across all 8 directions
r1 <- posrow;     c1 <- poscol - 1; idx1 <- which(dir8_m[cbind(r1, c1)] == 1); pos1row_0 <- r1[idx1]; pos1col_0 <- c1[idx1]; W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0
r2 <- posrow + 1; c2 <- poscol - 1; idx2 <- which(dir8_m[cbind(r2, c2)] == 2); pos2row_0 <- r2[idx2]; pos2col_0 <- c2[idx2]; W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
r3 <- posrow + 1; c3 <- poscol;     idx3 <- which(dir8_m[cbind(r3, c3)] == 3); pos3row_0 <- r3[idx3]; pos3col_0 <- c3[idx3]; W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0
r4 <- posrow + 1; c4 <- poscol + 1; idx4 <- which(dir8_m[cbind(r4, c4)] == 4); pos4row_0 <- r4[idx4]; pos4col_0 <- c4[idx4]; W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
r5 <- posrow;     c5 <- poscol + 1; idx5 <- which(dir8_m[cbind(r5, c5)] == 5); pos5row_0 <- r5[idx5]; pos5col_0 <- c5[idx5]; W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0
r6 <- posrow - 1; c6 <- poscol + 1; idx6 <- which(dir8_m[cbind(r6, c6)] == 6); pos6row_0 <- r6[idx6]; pos6col_0 <- c6[idx6]; W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
r7 <- posrow - 1; c7 <- poscol;     idx7 <- which(dir8_m[cbind(r7, c7)] == 7); pos7row_0 <- r7[idx7]; pos7col_0 <- c7[idx7]; W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0
r8 <- posrow - 1; c8 <- poscol - 1; idx8 <- which(dir8_m[cbind(r8, c8)] == 8); pos8row_0 <- r8[idx8]; pos8col_0 <- c8[idx8]; W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0

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
  
  # Search all 8 directions with correct respective offsets (fixes Direction 3, 4, 5, 7 bug)
  r1 <- pos0row_bis;     c1 <- pos0col_bis - 1; idx1 <- which(dir8_m[cbind(r1, c1)] == 1); pos1row_0 <- r1[idx1]; pos1col_0 <- c1[idx1]
  r2 <- pos0row_bis + 1; c2 <- pos0col_bis - 1; idx2 <- which(dir8_m[cbind(r2, c2)] == 2); pos2row_0 <- r2[idx2]; pos2col_0 <- c2[idx2]
  r3 <- pos0row_bis + 1; c3 <- pos0col_bis;     idx3 <- which(dir8_m[cbind(r3, c3)] == 3); pos3row_0 <- r3[idx3]; pos3col_0 <- c3[idx3]
  r4 <- pos0row_bis + 1; c4 <- pos0col_bis + 1; idx4 <- which(dir8_m[cbind(r4, c4)] == 4); pos4row_0 <- r4[idx4]; pos4col_0 <- c4[idx4]
  r5 <- pos0row_bis;     c5 <- pos0col_bis + 1; idx5 <- which(dir8_m[cbind(r5, c5)] == 5); pos5row_0 <- r5[idx5]; pos5col_0 <- c5[idx5]
  r6 <- pos0row_bis - 1; c6 <- pos0col_bis + 1; idx6 <- which(dir8_m[cbind(r6, c6)] == 6); pos6row_0 <- r6[idx6]; pos6col_0 <- c6[idx6]
  r7 <- pos0row_bis - 1; c7 <- pos0col_bis;     idx7 <- which(dir8_m[cbind(r7, c7)] == 7); pos7row_0 <- r7[idx7]; pos7col_0 <- c7[idx7]
  r8 <- pos0row_bis - 1; c8 <- pos0col_bis - 1; idx8 <- which(dir8_m[cbind(r8, c8)] == 8); pos8row_0 <- r8[idx8]; pos8col_0 <- c8[idx8]
  
  count <- count + 1
}

cat(sprintf("  Flow length loop finished in %d iterations!\n", count))

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
writeRaster(Ddn, paste0("Ddn_100m_5x5_", gauge_id, ".tif"), overwrite = TRUE)

# Cleanup Ddn memory
rm(Wgt_m, W_FL_m, dir8_m, dif_WFL, inv_CS, flow_dir_D8); gc()

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

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc"))
accW <- raster("accW.tif") * res

constant_res <- z_fill * 0 + res
Acc_fin      <- sca / constant_res
smean        <- (accS + slope_D8_r_m) / Acc_fin
cmean        <- (accW + W_factor) / Acc_fin

Dup <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
crs(Dup) <- crs_out
writeRaster(Dup, paste0("Dup_100m_5x5_", gauge_id, ".tif"), overwrite = TRUE)

rm(constant_res, ang, Acc_fin, smean, accW, cmean); gc()

# --------------------------------------------------------------------------
# 7. Final Index of Connectivity (IC)
# --------------------------------------------------------------------------
cat("7. Calculating Final Index of Connectivity...\n")
ratio <- Dup / Ddn
IC    <- log10(ratio)
IC    <- IC * mask_cells
crs(IC) <- crs_out

pal <- colorRampPalette(c("blue", "lightyellow1", "red"))

png(paste0("IC_100m_5x5_", gauge_id, ".png"), width = 800, height = 800)
plot(IC, col = pal(100), main = paste("Index of Connectivity (IC) - Gauge", gauge_id))
dev.off()

plot(IC, col = pal(100), main = paste("Index of Connectivity (IC) - Gauge", gauge_id))
writeRaster(IC, paste0("IC_", res, "m_5x5_", gauge_id, ".tif"), overwrite = TRUE)

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
  
  if (ncol(pointprofile_vect) > 0) {
    pointprofile_vect <- pointprofile_vect[, 0]
  }
  
  system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accW8D.tif -nc"))
  Acc_fin <- raster("accW8D.tif")
  
  rasStack <- stack(W_factor, Dup, Ddn, IC, z_fill, slope_D8, Acc_fin)
  
  pointprofile_sp <- as(pointprofile_vect, "Spatial")
  pointprofile_sp$POINT_X <- pointprofile_sp@coords[, 1]
  pointprofile_sp$POINT_Y <- pointprofile_sp@coords[, 2]
  
  pointprofile2 <- raster::extract(rasStack, pointprofile_sp, method = "simple", df = TRUE, sp = TRUE)
  names(pointprofile2) <- c("POINT_X", "POINT_Y", "W", "Dup", "Ddn", "IC", "Z", "S8", "A")
  
  pointprofile2@data <- pointprofile2@data[order(pointprofile2@data$Z), ]
  pointprofile2@data$Absc[1] <- 0
  for (i in 2:length(pointprofile2@data$Z)) {
    pointprofile2@data$Absc[i] <- pointprofile2@data$Absc[i - 1] + 
      ((pointprofile2@data$POINT_X[i] - pointprofile2@data$POINT_X[i - 1])^2 + 
         (pointprofile2@data$POINT_Y[i] - pointprofile2@data$POINT_Y[i - 1])^2)^0.5
  }
  
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
  
  lm_pol <- lm(DATA$IC ~ poly(DATA$Absc, 3, raw = TRUE))
  rsqr   <- round(summary(lm_pol)$adj.r.squared, 3)
  
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