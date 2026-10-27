# ==============================================================================
# STANDALONE DEBUG SCRIPT FOR GAUGE 03002 (WITH STREAM SNAPPING)
# ==============================================================================

if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

library(terra)
library(raster)

rasterOptions(tolerance = 0.05)

debug_dir <- "F:/Peninsular India/Connectivity/IC/Debug"
if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)
setwd(debug_dir)

dem_path    <- "F:/Peninsular India/Connectivity/IC/Input/DEM/DEM_03002.tif"
ndvi_path   <- "F:/Peninsular India/Connectivity/IC/Input/NDVI/2021/NDVI_03002.tif"
target_path <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/shapefiles_catchment/gauge_stations.shp"
gauge_id    <- "03002"
target_crs  <- "EPSG:32644" # UTM Zone 44N
target_res  <- 100

check_raster <- function(r, name) {
  if (inherits(r, "SpatRaster")) {
    v <- terra::values(r)
  } else {
    v <- raster::getValues(r)
  }
  v_clean <- v[!is.na(v) & !is.nan(v) & !is.infinite(v)]
  cat(sprintf("--------------------------------------------------\n"))
  cat(sprintf("STATUS: %s\n", name))
  cat(sprintf("  Total Cells : %d\n", length(v)))
  cat(sprintf("  Valid Cells : %d (%.1f%%)\n", length(v_clean), (length(v_clean)/length(v))*100))
  cat(sprintf("  NA/NaN Cells: %d\n", sum(is.na(v) | is.nan(v))))
  if (length(v_clean) > 0) {
    cat(sprintf("  Min / Max   : %.4f / %.4f\n", min(v_clean), max(v_clean)))
    cat(sprintf("  Mean        : %.4f\n", mean(v_clean)))
  } else {
    cat("  >>> CRITICAL ERROR: RASTER IS COMPLETELY BLANK (100% NA)! <<<\n")
  }
}

cat("\n=================== STARTING DIAGNOSTICS ===================\n\n")

# 1. Reproject DEM & NDVI
cat("1. Reprojecting DEM & NDVI to UTM 44N...\n")
dem_raw   <- terra::rast(dem_path)
ndvi_raw  <- terra::rast(ndvi_path)
dem_proj  <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
ndvi_proj <- terra::project(ndvi_raw, dem_proj, method = "bilinear")

z_fill   <- raster(dem_proj)
W_factor <- raster(ndvi_proj)
writeRaster(z_fill, "dem_filled.tif", overwrite = TRUE)

# 2. C-factor Calculation
cat("2. Computing C-factor...\n")
vals <- getValues(W_factor)
vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
ndvi_min <- min(vals_clean)
ndvi_max <- max(vals_clean)

f_g <- ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))
C_factor[f_g <= 0]     <- 1.0
C_factor[f_g > 78.3]   <- 0.001
C_factor[C_factor > 1] <- 1.0
C_factor[C_factor < 0.001] <- 0.001

W_factor <- C_factor
writeRaster(W_factor, "W_factor.tif", overwrite = TRUE)

# 3. TauDEM D8 Flow Direction & Flow Accumulation (for Snapping)
cat("3. Running TauDEM D8 FlowDir & Accumulation...\n")
system("mpiexec -n 4 D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_filled.tif")
flow_dir_D8 <- raster("D8_flow_dir.tif")
slope_D8    <- raster("slope_D8.tif")

# Compute D8 area to find where the stream channel actually runs
system("mpiexec -n 4 AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc")
accD8 <- raster("accD8.tif")

# 4. Target Point Extraction & Stream Snapping
cat("4. Snapping Gauge Point to Stream Channel...\n")
all_targets <- terra::vect(target_path)
shp_ids <- trimws(as.character(all_targets$gauge_id))
target_feat <- all_targets[shp_ids == gauge_id | sprintf("%05d", suppressWarnings(as.numeric(shp_ids))) == gauge_id, ]
target_feat <- terra::project(target_feat, target_crs)

# Snap gauge point to cell with highest flow accumulation within 500m (5 cells)
pt_coords <- terra::geom(target_feat)[, c("x", "y")]
cat(sprintf("  Original Gauge Coords: X=%.1f, Y=%.1f\n", pt_coords[1], pt_coords[2]))

# # Search buffer around the gauge
# snap_radius <- 500 # 500 meters search radius
# buf <- terra::buffer(target_feat, snap_radius)
# accD8_crop <- crop(accD8, extent(buf))
# accD8_mask <- mask(accD8_crop, raster(buf))
# 
# max_cell <- which.max(getValues(accD8_mask))
# snapped_xy <- xyFromCell(accD8_mask, max_cell)
# cat(sprintf("  Snapped Stream Coords: X=%.1f, Y=%.1f (Max Accumulation: %.0f cells)\n", 
#             snapped_xy[1], snapped_xy[2], max(getValues(accD8_mask), na.rm = TRUE)))
# 
# # Create snapped target point
# target_snapped <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs)
# target_snapped$burn <- 1
# 
# # Rasterize snapped target
# target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
# target_rst   <- raster(target_rst_t)
# target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))


# Search buffer around the gauge using terra natively
snap_radius <- 500 # 500 meters search radius
buf <- terra::buffer(target_feat, snap_radius)

accD8_t    <- terra::rast("accD8.tif")
accD8_crop <- terra::crop(accD8_t, buf)
accD8_mask <- terra::mask(accD8_crop, buf)

# Find coordinates of maximum flow accumulation cell
acc_vals <- terra::values(accD8_mask)
max_idx  <- which.max(acc_vals)
snapped_xy <- terra::xyFromCell(accD8_mask, max_idx)

max_acc <- acc_vals[max_idx]
cat(sprintf("  Snapped Stream Coords: X=%.1f, Y=%.1f (Max Accumulation: %.0f cells)\n", 
            snapped_xy[1, 1], snapped_xy[1, 2], max_acc))

# Create snapped target point as SpatVector
target_snapped <- terra::vect(snapped_xy, crs = target_crs)
target_snapped$burn <- 1

# Rasterize snapped target
target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
target_rst   <- raster(target_rst_t)
target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))


# 5. Downslope Component (Ddn)
cat("5. Tracing Downslope Flow Length...\n")
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8 * target_rst

slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE) # FIX: Written to disk!

inv_CS <- 1 / (W_factor * slope_D8_r_m)
writeRaster(inv_CS, "inv_CS.tif", overwrite = TRUE)

res <- res(z_fill)[1]
dir8_m <- as.matrix(flow_dir_D8)
W_FL_m <- as.matrix(inv_CS)

coladd <- vector(mode = "numeric", length = nrow(dir8_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd, dir8_m, coladd)
W_FL_m <- cbind(coladd, W_FL_m, coladd)
rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd, dir8_m, rowadd)
W_FL_m <- rbind(rowadd, W_FL_m, rowadd)

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
    W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0
    W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
    W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0
    W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
    W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0
    W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
    W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0
    W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0
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

cat(sprintf("  Loop completed in %d iterations.\n", count))

W_FL_m <- W_FL_m[-1, ]; W_FL_m <- W_FL_m[-nrow(W_FL_m), ]
dir8_m <- dir8_m[-1, ]; dir8_m <- dir8_m[-nrow(dir8_m), ]
W_FL_m <- W_FL_m[, -1]; W_FL_m <- W_FL_m[, -ncol(W_FL_m)] 
dir8_m <- dir8_m[, -1]; dir8_m <- dir8_m[, -ncol(dir8_m)]

ext <- as.matrix(extent(inv_CS))
W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = CRS("+init=epsg:32644"))

dif_WFL <- W_FL_fin - inv_CS
dif_WFL[dif_WFL == 0] <- NA
mask_cells <- dif_WFL * 0 + 1
check_raster(mask_cells, "Watershed Cell Mask (mask_cells)")

Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
writeRaster(Ddn, "Ddn.tif", overwrite = TRUE)
check_raster(Ddn, "Downslope Factor (Ddn)")

# 6. Upslope Component (Dup)
cat("6. Running TauDEM D-inf Area...\n")
system("mpiexec -n 4 DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_filled.tif")
ang <- raster("ang.tif") * target_rst
writeRaster(ang, "ang.tif", overwrite = TRUE)

system("mpiexec -n 4 AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc")
accS <- raster("accS.tif")

system("mpiexec -n 4 AreaDinf -ang ang.tif -sca sca.tif -nc")
sca <- raster("sca.tif")

system("mpiexec -n 4 AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc")
accW <- raster("accW.tif") * res

constant_res <- z_fill * 0 + res
Acc_fin      <- sca / constant_res
smean        <- (accS + slope_D8_r_m) / Acc_fin
cmean        <- (accW + W_factor) / Acc_fin
Dup          <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
writeRaster(Dup, "Dup.tif", overwrite = TRUE)
check_raster(Dup, "Upslope Factor (Dup)")

# 7. Final IC
cat("7. Calculating Final IC...\n")
ratio <- Dup / Ddn
check_raster(ratio, "Ratio (Dup / Ddn)")

IC <- log10(ratio)
IC <- IC * mask_cells
writeRaster(IC, "IC.tif", overwrite = TRUE)
check_raster(IC, "Final IC Map")

cat("\n=================== DIAGNOSTICS COMPLETE ===================\n")