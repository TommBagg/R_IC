# ==============================================================================
# STANDALONE DEBUG SCRIPT (WITH PITREMOVE, TRUE OUTLET, & CLEAN MASK_CELLS)
# ==============================================================================

# 1. Resolve TauDEM PROJ environment conflict
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

library(terra)
library(raster)

rasterOptions(tolerance = 0.05)

# Set working / output directory
debug_dir <- "F:/Peninsular India/Connectivity/IC/Debug"
if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)
setwd(debug_dir)

# File Paths & Parameters (Change gauge_id as needed)
gauge_id    <- "09012"
dem_path    <- paste0("F:/Peninsular India/Connectivity/IC/Input/DEM/DEM_", gauge_id, ".tif")
ndvi_path   <- paste0("F:/Peninsular India/Connectivity/IC/Input/NDVI/2021/NDVI_", gauge_id, ".tif")
target_crs  <- "EPSG:32644" # WGS 84 / UTM Zone 44N
target_res  <- 100          # 100 meters
mpi_cores   <- 4

cat(sprintf("\n=================== RUNNING DEBUG WITH MASK_CELLS: GAUGE %s ===================\n\n", gauge_id))

# ------------------------------------------------------------------------------
# STEP 1: Load Raw Inputs & Reproject DEM
# ------------------------------------------------------------------------------
cat("1. Loading raw DEM & NDVI, reprojecting DEM to 100m metric grid...\n")
dem_raw  <- terra::rast(dem_path)
ndvi_raw <- terra::rast(ndvi_path)

dem_proj_t <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 2: Pit-Filling Upfront (TauDEM PitRemove)
# ------------------------------------------------------------------------------
cat("2. Running TauDEM PitRemove to resolve all surface depressions...\n")
system(paste("mpiexec -n", mpi_cores, "PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif"))

# Master elevation raster used for all further calculations
z_fill <- raster("dem_pitfilled.tif")
res    <- res(z_fill)[1]

# ------------------------------------------------------------------------------
# STEP 3: Align NDVI & Compute C-factor (Cai et al., 2000)
# ------------------------------------------------------------------------------
cat("3. Computing C-factor from NDVI...\n")
wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
W_factor  <- raster(wf_proj_t)

vals <- getValues(W_factor)
vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
ndvi_min <- min(vals_clean)
ndvi_max <- max(vals_clean)

# Fractional vegetation cover in percent (0 to 100%)
if (ndvi_max == ndvi_min) {
  f_g <- W_factor * 0
} else {
  f_g <- ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
}

f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))

C_factor[f_g <= 0]     <- 1.0
C_factor[f_g > 78.3]   <- 0.001   # Minimum threshold avoids division-by-zero
C_factor[C_factor > 1] <- 1.0
C_factor[C_factor < 0.001] <- 0.001

W_factor <- C_factor
writeRaster(W_factor, "W_factor.tif", overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 4: Flow Directions & True Catchment Outlet Identification
# ------------------------------------------------------------------------------
cat("4. Computing D8 Flow Directions & Locating Catchment Outlet...\n")
system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif"))
flow_dir_D8 <- raster("D8_flow_dir.tif")
slope_D8    <- raster("slope_D8.tif")

# Calculate flow accumulation
system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc"))
accD8 <- raster("accD8.tif")

# Find the true catchment outlet cell (maximum D8 accumulation)
max_cell   <- which.max(getValues(accD8))
snapped_xy <- xyFromCell(accD8, max_cell)
max_accum  <- max(getValues(accD8), na.rm = TRUE)

cat(sprintf("  Catchment Outlet located at X=%.1f, Y=%.1f (Total Drainage: %d cells)\n", 
            snapped_xy[1], snapped_xy[2], as.integer(max_accum)))

target_snapped      <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs)
target_snapped$burn <- 1

target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
target_rst   <- raster(target_rst_t)
target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))

# ------------------------------------------------------------------------------
# STEP 5: Impedance & Downslope Flow Length Algorithm (Boundary Artifact Fixed)
# ------------------------------------------------------------------------------
cat("5. Tracing Downslope Flow Length from Outlet...\n")
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8 * target_rst

slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)

inv_CS <- 1 / (W_factor * slope_D8_r_m)
writeRaster(inv_CS, "inv_CS.tif", overwrite = TRUE)

# FIX FOR BOUNDARY VALUES: Find target index BEFORE padding matrix
dir8_raw_m <- as.matrix(flow_dir_D8)
target_idx <- which(is.na(dir8_raw_m), arr.ind = TRUE)

# Add 1-cell boundary padding to prevent index-out-of-bounds
coladd <- vector(mode = "numeric", length = nrow(dir8_raw_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd, dir8_raw_m, coladd)
W_FL_m <- cbind(coladd, as.matrix(inv_CS), coladd)

rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd, dir8_m, rowadd)
W_FL_m <- rbind(rowadd, W_FL_m, rowadd)

Wgt_m <- W_FL_m 

# Sinks to trace from: ONLY the true outlet (shifted +1 for padding), NOT the outer border!
posrow <- target_idx[, "row"] + 1
poscol <- target_idx[, "col"] + 1

pos1col_NA <- poscol - 1; pos1row_NA <- posrow
pos1row <- pos1row_NA[pos1col_NA != 0]; pos1col <- pos1col_NA[pos1col_NA != 0]
pos1dir <- dir8_m[cbind(pos1row, pos1col)]
pos1row_0 <- pos1row[which(pos1dir == 1)]; pos1col_0 <- pos1col[which(pos1dir == 1)]
W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0

pos2col_NA <- poscol - 1; pos2row_NA <- posrow + 1
pos2rowt <- pos2row_NA[pos2col_NA != 0]; pos2colt <- pos2col_NA[pos2col_NA != 0]
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
pos6row <- pos6rowt[pos6rowt != 0]; pos6col <- pos6colt[pos6rowt != 0]
pos6dir <- dir8_m[cbind(pos6row, pos6col)]
pos6row_0 <- pos6row[which(pos6dir == 6)]; pos6col_0 <- pos6col[which(pos6dir == 6)]
W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0

pos7col_NA <- poscol; pos7row_NA <- posrow - 1
pos7row <- pos7row_NA[pos7row_NA != 0]; pos7col <- pos7col_NA[pos7row_NA != 0]
pos7dir <- dir8_m[cbind(pos7row, pos7col)]
pos7row_0 <- pos7row[which(pos7dir == 7)]; pos7col_0 <- pos7col[which(pos7dir == 7)]
W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0 

pos8col_NA <- poscol - 1; pos8row_NA <- posrow - 1
pos8rowt <- pos8row_NA[pos8row_NA != 0]; pos8colt <- pos8col_NA[pos8row_NA != 0]
pos8row <- pos8rowt[pos8rowt != 0]; pos8col <- pos8colt[pos8rowt != 0]
pos8dir <- dir8_m[cbind(pos8row, pos8col)]
pos8row_0 <- pos8row[which(pos8dir == 8)]; pos8col_0 <- poscol[which(pos8dir == 8)]
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
  posrow <- posrowt[poscolt != 0]; poscol <- poscolt[poscolt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos1row_0 <- posrow[which(posdir == 1)]; pos1col_0 <- poscol[which(posdir == 1)]
  
  poscolt <- pos0col_bis - 1; posrowt <- pos0row_bis + 1
  posrow <- posrowt[poscolt != 0]; poscol <- poscolt[poscolt != 0]
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
  posrow <- posrowt[posrowt != 0]; poscol <- poscolt[posrowt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos6row_0 <- posrow[which(posdir == 6)]; pos6col_0 <- poscol[which(posdir == 6)]
  
  poscolt <- pos0col_bis; posrowt <- pos0row_bis - 1
  posrow <- posrowt[posrowt != 0]; poscol <- poscolt[posrowt != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos7row_0 <- posrow[which(posdir == 7)]; pos7col_0 <- poscol[which(posdir == 7)]
  
  poscolt <- pos0col_bis - 1; posrowt <- pos0row_bis - 1
  posrowt1 <- posrowt[posrowt != 0]; poscolt1 <- poscolt[posrowt != 0]
  posrow <- posrowt1[poscolt1 != 0]; poscol <- poscolt1[poscolt1 != 0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos8row_0 <- posrow[which(posdir == 8)]; pos8col_0 <- poscol[which(posdir == 8)]
  
  count <- count + 1
}

cat(sprintf("  Flow length loop finished in %d iterations!\n", count))

# Remove boundary padding
W_FL_m <- W_FL_m[-1, ]; W_FL_m <- W_FL_m[-nrow(W_FL_m), ]
dir8_m <- dir8_m[-1, ]; dir8_m <- dir8_m[-nrow(dir8_m), ]
W_FL_m <- W_FL_m[, -1]; W_FL_m <- W_FL_m[, -ncol(W_FL_m)] 
dir8_m <- dir8_m[, -1]; dir8_m <- dir8_m[, -ncol(dir8_m)]

crs_out <- CRS("+init=epsg:32644")
ext <- as.matrix(extent(inv_CS))
W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs_out)
writeRaster(W_FL_fin, "Weighted_Flow_length.tif", overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 6: Compute Watershed Mask (mask_cells)
# ------------------------------------------------------------------------------
dif_WFL <- W_FL_fin - inv_CS
dif_WFL[dif_WFL == 0] <- NA
mask_cells <- dif_WFL * 0 + 1
writeRaster(mask_cells, "mask_cells.tif", overwrite = TRUE)

Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
crs(Ddn) <- crs_out
writeRaster(Ddn, "Ddn.tif", overwrite = TRUE)

# Cleanup Ddn memory
rm(Wgt_m, W_FL_m, dir8_m, dif_WFL, inv_CS, flow_dir_D8); gc()

# ------------------------------------------------------------------------------
# STEP 7: Upslope Component (Dup)
# ------------------------------------------------------------------------------
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
writeRaster(Dup, "Dup.tif", overwrite = TRUE)

rm(constant_res, ang, Acc_fin, smean, accW, cmean); gc()

# ------------------------------------------------------------------------------
# STEP 8: Final Index of Connectivity (IC) WITH mask_cells
# ------------------------------------------------------------------------------
cat("7. Calculating Final Index of Connectivity with mask_cells...\n")
ratio <- Dup / Ddn

# Apply mask_cells to neatly clip to the watershed
IC <- log10(ratio)
IC <- IC * mask_cells
crs(IC) <- crs_out

# Preview Plot
pal <- colorRampPalette(c("blue", "lightyellow1", "red"))
png(paste0("IC_", gauge_id, ".png"), width = 800, height = 800)
plot(IC, col = pal(100), main = paste("Index of Connectivity (IC) - Catchment", gauge_id))
dev.off()

plot(IC, col = pal(100), main = paste("Index of Connectivity (IC) - Catchment", gauge_id))
writeRaster(IC, paste0("IC_", gauge_id, ".tif"), overwrite = TRUE)

cat("\n=================== FINAL STATUS ===================\n")
valid_mask <- sum(!is.na(getValues(mask_cells)))
valid_ic   <- sum(!is.na(getValues(IC)))
total_cells <- length(getValues(IC))

cat(sprintf("Valid Cells in mask_cells : %d (%.1f%% of DEM)\n", valid_mask, (valid_mask / total_cells) * 100))
cat(sprintf("Valid Cells in Final IC   : %d (%.1f%% of DEM)\n", valid_ic, (valid_ic / total_cells) * 100))
cat("IC Range                  : [", min(getValues(IC), na.rm = TRUE), "to", max(getValues(IC), na.rm = TRUE), "]\n")
cat("Run successfully completed!\n")