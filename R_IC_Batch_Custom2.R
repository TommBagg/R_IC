# Title: Batch R_IC v. 2.2 (Single Target Shapefile, DEM/NDVI 5-Digit ID Matching)
# Author details: Tommaso Baggio, Lorenzo Martini, Loris Torresani

'''
# 1. Point R back to its own internal PROJ database
if (dir.exists(system.file("proj", package = "rgdal"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "rgdal"))
} else if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else if (dir.exists(system.file("proj", package = "sf"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "sf"))
} else {
  # If none of the above, unset it so R doesnt look at TauDEMs folder
  Sys.unsetenv("PROJ_LIB")
}

# Packages used
require("raster")
require("rgdal")
require("sp")
'''

# 1. Point R to terra's PROJ database (resolves the TauDEM PROJ warning)
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

# 2. Load modern packages (no rgdal needed)
require("raster")
require("terra")

# Premise: all files will be overwritten if they already exist

###################################
# INPUTS & DIRECTORIES SETUP
###################################

# Output working directory
wd <- "F:/Peninsular India/Connectivity/IC/Output"
if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

# Directories for DEM and NDVI weight rasters
dem_dir      <- "F:/Peninsular India/Connectivity/IC/Input/DEM"       # Folder containing DEM_xxxxx.tif files
w_factor_dir <- "F:/Peninsular India/Connectivity/IC/Input/NDVI/1990"      # Folder containing NDVI_xxxxx.tif files

# Path to the SINGLE master target shapefile and the ID column name
target_shp_path <- "F:/Peninsular India/Connectivity/IC/Input/Catchment/Catchments.shp"
gauge_id_col    <- "gauge_id"                           # Name of the attribute column

# Moving window size (retained for nomenclature and function metadata)
mw <- 5

# Flag: save all outputs? 
# 0 = keep essential files (IC, Dup, Ddn)
# 1 = keep all intermediate files (accS, accW, ang, D8 flow dir, etc.)
flag_k <- 0

# Number of MPI cores to use for TauDEM functions
mpi_cores <- 4


###################################
# LOAD MASTER TARGET SHAPEFILE & MATCH FILES
###################################

cat("Loading master target shapefile...\n")
# all_targets <- shapefile(target_shp_path)
all_targets <- terra::vect(target_shp_path)

# List all DEM files starting with DEM_ and ending with 5 digits
dem_file_list <- list.files(dem_dir, pattern = "^DEM_[0-9]{5}\\.(tif|TIF)$", full.names = TRUE)

# Fallback: if filenames have variations, find all .tif with 'DEM_'
if (length(dem_file_list) == 0) {
  dem_file_list <- list.files(dem_dir, pattern = "\\.(tif|TIF)$", full.names = TRUE)
  dem_file_list <- dem_file_list[grepl("DEM_", basename(dem_file_list), ignore.case = TRUE)]
}

batch_list <- data.frame(
  gauge_id      = character(),
  dem_path      = character(),
  w_factor_path = character(),
  stringsAsFactors = FALSE
)

# Available gauge IDs in the master shapefile (trimmed and formatted)
# shp_gauge_ids <- trimws(as.character(all_targets[[gauge_id_col]]))
shp_gauge_ids <- trimws(as.character(all_targets[[gauge_id_col]][[1]]))

# Match DEM, NDVI, and shapefile feature
for (df in dem_file_list) {
  fname <- basename(df)
  
  # Extract the 5-digit numeric ID (xxxxx)
  gid <- sub(".*DEM_([0-9]{5})\\..*", "\\1", fname, ignore.case = TRUE)
  
  # Expected corresponding NDVI weight file
  wf_path <- file.path(w_factor_dir, paste0("NDVI_", gid, ".tif"))
  
  # Verify if NDVI file and matching gauge_id in shapefile exist
  wf_exists  <- file.exists(wf_path)
  shp_exists <- any(shp_gauge_ids == gid | 
                      sprintf("%05d", suppressWarnings(as.numeric(shp_gauge_ids))) == gid)
  
  if (wf_exists && shp_exists) {
    batch_list <- rbind(batch_list, data.frame(
      gauge_id      = gid,
      dem_path      = df,
      w_factor_path = wf_path,
      stringsAsFactors = FALSE
    ))
  } else {
    if (!wf_exists)  warning(paste("Missing NDVI file for ID:", gid, "- skipping!"))
    if (!shp_exists) warning(paste("Gauge ID", gid, "not found in master shapefile - skipping!"))
  }
}

cat(sprintf("Found %d matched sets ready for processing.\n", nrow(batch_list)))
print(head(batch_list))


###################################################################################################################
# START OF THE FUNCTION
###################################################################################################################

IC_function <- function(z_fill, target_feat, W_factor, mw, flag_k, site_id) {
  
  start_timeTOT <- Sys.time()
  
  #### Check if x and y cellsize are equal
  resNS <- res(z_fill)[1]; resWE <- res(z_fill)[2]
  ifelse(resNS == resWE, print("Resolution NS/WE: OK"), print("Error!! NS and WE resolution are not equal"))
  res <- res(z_fill)[1]
  
  #### Check if resolution dem = resolution weighting factor
  res_WF <- res(W_factor)[1]
  ifelse(res_WF == res, print("WF resolution OK"), print("Error!! Weighting factor and DTM resolution do not match"))
  
  ################### Weighting factor calculation
  # Exact calculation (negative values strictly retained)
  W_factor <- 0.431 - (0.805 * W_factor)
  writeRaster(W_factor, filename = "W_factor.tif", overwrite = TRUE)
  
  ######################################################################
  # CREATE MASK FOR TARGETS
  # Ensure target feature CRS matches DEM CRS
  #if (!compareCRS(target_feat, z_fill)) {
  #  target_feat <- spTransform(target_feat, crs(z_fill))
  #}
  
  #target_rst <- rasterize(target_feat, z_fill, background = 2, field = 1)
  #target_rst <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))
  
  if (crs(target_feat) != crs(z_fill)) {
    target_feat <- terra::project(target_feat, crs(z_fill))
  }
  
  # raster::rasterize natively accepts terra vector objects
  target_rst <- rasterize(target_feat, z_fill, background = 2, field = 1)
  target_rst <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))
  
  ######################################################################
  # COMPUTING Ddown
  writeRaster(z_fill, "dem_filled.tif", overwrite = TRUE)
  
  # Flow directions 8 directions
  system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_filled.tif"))
  flow_dir_D8 <- raster("D8_flow_dir.tif")
  slope_D8    <- raster("slope_D8.tif")
  
  # Integrate targets in flow directions
  flow_dir_D8[is.na(flow_dir_D8)] <- -9999
  flow_dir_D8 <- flow_dir_D8 * target_rst
  
  # Reclassified values < 0.005 m/m and > 1 m/m
  slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
  writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)
  
  inv_CS <- 1 / (W_factor * slope_D8_r_m)
  writeRaster(inv_CS, paste0("inv_CS_", res, "m_", mw, "x", mw, "_", site_id, ".tif"), overwrite = TRUE)
  
  gc()
  
  ######################### FLOW WEIGHTED LENGTH
  dir8 <- flow_dir_D8
  W_FL <- inv_CS
  
  W_FL_m <- as.matrix(W_FL)
  dir8_m <- as.matrix(dir8)
  
  ## Add boundary padding
  coladd <- vector(mode = "numeric", length = nrow(dir8_m)); coladd[1:length(coladd)] <- NA
  dir8_m <- cbind(coladd, dir8_m, coladd)
  W_FL_m <- cbind(coladd, W_FL_m, coladd)
  
  rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
  dir8_m <- rbind(rowadd, dir8_m, rowadd)
  W_FL_m <- rbind(rowadd, W_FL_m, rowadd)
  
  start_time <- Sys.time()
  
  Wgt_m <- W_FL_m 
  NA_data <- which(is.na(dir8_m), arr.ind = TRUE); NA_data <- as.data.frame(NA_data)
  posrow <- NA_data$row; poscol <- NA_data$col
  gc()
  
  ########## Flow direction 1
  pos1col_NA <- poscol - 1; pos1row_NA <- posrow
  pos1row <- pos1row_NA[-pos1col_NA != 0]; pos1col <- pos1col_NA[-pos1col_NA != 0]
  pos1dir <- dir8_m[cbind(pos1row, pos1col)]
  pos1row_0 <- pos1row[which(pos1dir == 1)]; pos1col_0 <- pos1col[which(pos1dir == 1)]
  W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0
  
  ########## Flow direction 2
  pos2col_NA <- poscol - 1; pos2row_NA <- posrow + 1
  pos2rowt <- pos2row_NA[-pos2col_NA != 0]; pos2colt <- pos2col_NA[-pos2col_NA != 0]
  pos2row <- pos2rowt[pos2rowt <= (nrow(dir8_m))]; pos2col <- pos2colt[pos2rowt <= (nrow(dir8_m))]
  pos2dir <- dir8_m[cbind(pos2row, pos2col)]
  pos2row_0 <- pos2row[which(pos2dir == 2)]; pos2col_0 <- pos2col[which(pos2dir == 2)]
  W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
  
  ########## Flow direction 3
  pos3col_NA <- poscol; pos3row_NA <- posrow + 1
  pos3row <- pos3row_NA[pos3row_NA <= (nrow(dir8_m))]; pos3col <- pos3col_NA[pos3row_NA <= (nrow(dir8_m))]
  pos3dir <- dir8_m[cbind(pos3row, pos3col)]
  pos3row_0 <- pos3row[which(pos3dir == 3)]; pos3col_0 <- pos3col[which(pos3dir == 3)]
  W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0
  
  ########## Flow direction 4
  pos4col_NA <- poscol + 1; pos4row_NA <- posrow + 1
  pos4rowt <- pos4row_NA[pos4col_NA <= (ncol(dir8_m))]; pos4colt <- pos4col_NA[pos4col_NA <= (ncol(dir8_m))]
  pos4row <- pos4rowt[pos4rowt <= (nrow(dir8_m))]; pos4col <- pos4colt[pos4rowt <= (nrow(dir8_m))]
  pos4dir <- dir8_m[cbind(pos4row, pos4col)]
  pos4row_0 <- pos4row[which(pos4dir == 4)]; pos4col_0 <- pos4col[which(pos4dir == 4)]
  W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
  
  ########## Flow direction 5
  pos5col_NA <- poscol + 1; pos5row_NA <- posrow
  pos5row <- pos5row_NA[pos5col_NA <= (ncol(dir8_m))]; pos5col <- pos5col_NA[pos5col_NA <= (ncol(dir8_m))]
  pos5dir <- dir8_m[cbind(pos5row, pos5col)]
  pos5row_0 <- pos5row[which(pos5dir == 5)]; pos5col_0 <- pos5col[which(pos5dir == 5)]
  W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0
  
  ########## Flow direction 6
  pos6col_NA <- poscol + 1; pos6row_NA <- posrow - 1
  pos6rowt <- pos6row_NA[pos6col_NA <= (ncol(dir8_m))]; pos6colt <- pos6col_NA[pos6col_NA <= (ncol(dir8_m))]
  pos6row <- pos6rowt[-pos6rowt != 0]; pos6col <- pos6colt[-pos6rowt != 0]
  pos6dir <- dir8_m[cbind(pos6row, pos6col)]
  pos6row_0 <- pos6row[which(pos6dir == 6)]; pos6col_0 <- pos6col[which(pos6dir == 6)]
  W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
  
  ########## Flow direction 7
  pos7col_NA <- poscol; pos7row_NA <- posrow - 1
  pos7row <- pos7row_NA[-pos7row_NA != 0]; pos7col <- pos7col_NA[-pos7row_NA != 0]
  pos7dir <- dir8_m[cbind(pos7row, pos7col)]
  pos7row_0 <- pos7row[which(pos7dir == 7)]; pos7col_0 <- pos7col[which(pos7dir == 7)]
  W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0 
  
  ########## Flow direction 8
  pos8col_NA <- poscol - 1; pos8row_NA <- posrow - 1
  pos8rowt <- pos8row_NA[-pos8row_NA != 0]; pos8colt <- pos8col_NA[-pos8row_NA != 0]
  pos8row <- pos8rowt[-pos8colt != 0]; pos8col <- pos8colt[-pos8colt != 0]
  pos8dir <- dir8_m[cbind(pos8row, pos8col)]
  pos8row_0 <- pos8row[which(pos8dir == 8)]; pos8col_0 <- pos8col[which(pos8dir == 8)]
  W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0
  
  rm(pos1col, pos1col_NA, pos1row, pos1row_NA, pos1dir)
  rm(pos3col, pos3col_NA, pos3row, pos3row_NA, pos3dir)
  rm(pos5col, pos5col_NA, pos5row, pos5row_NA, pos5dir)
  rm(pos7col, pos7col_NA, pos7row, pos7row_NA, pos7dir)
  rm(pos2col, pos2col_NA, pos2row, pos2row_NA, pos2colt, pos2rowt, pos2dir)
  rm(pos4col, pos4col_NA, pos4row, pos4row_NA, pos4colt, pos4rowt, pos4dir)
  rm(pos6col, pos6col_NA, pos6row, pos6row_NA, pos6colt, pos6rowt, pos6dir)
  rm(pos8col, pos8col_NA, pos8row, pos8row_NA, pos8colt, pos8rowt, pos8dir)
  
  ## Compute weighted sum following direction
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
  
  gc()
  
  # Remove padding
  W_FL_m <- W_FL_m[-1, ]; W_FL_m <- W_FL_m[-nrow(W_FL_m), ]
  dir8_m <- dir8_m[-1, ]; dir8_m <- dir8_m[-nrow(dir8_m), ]
  W_FL_m <- W_FL_m[, -1]; W_FL_m <- W_FL_m[, -ncol(W_FL_m)] 
  dir8_m <- dir8_m[, -1]; dir8_m <- dir8_m[, -ncol(dir8_m)]
  
  crs <- projection(W_FL)
  ext <- as.matrix(extent(W_FL))
  
  W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs)
  
  # Remove cells
  dif_WFL <- W_FL_fin - inv_CS
  dif_WFL[dif_WFL == 0] <- NA
  mask_cells <- dif_WFL * 0 + 1
  
  ##### Computing final Downslope component
  Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
  writeRaster(Ddn, paste0("Ddn_", res, "m_", mw, "x", mw, "_", site_id, ".tif"), overwrite = TRUE)
  
  rm(Wgt_m, W_FL_m, dir8_m, dif_WFL, W_FL, dir8, inv_CS, flow_dir_D8); gc()
  
  ######################################################################
  # COMPUTING Dup
  system(paste("mpiexec -n", mpi_cores, "DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_filled.tif"))
  ang <- raster("ang.tif")
  ang <- ang * target_rst
  writeRaster(ang, "ang.tif", overwrite = TRUE)
  gc()
  
  # Computing D-inf area weighted for slope
  system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc"))
  accS <- raster("accS.tif")
  
  # Computing D-inf area
  system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca sca.tif -nc"))
  sca <- raster("sca.tif")
  
  constant_res <- z_fill * 0 + res
  Acc_fin      <- sca / constant_res
  smean        <- (accS + slope_D8_r_m) / Acc_fin
  
  rm(slope_D8_r_m, accS, target_rst); gc()
  
  # Computing D-inf Area weighted for W_factor (acc W)
  system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc"))
  accW <- raster("accW.tif") * res
  cmean <- (accW + W_factor) / Acc_fin
  
  ### Computing final Upslope component
  Dup <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
  writeRaster(Dup, paste0("Dup_", res, "m_", mw, "x", mw, "_", site_id, ".tif"), overwrite = TRUE)
  
  rm(constant_res, ang, Acc_fin, smean, accW, cmean); gc()
  
  ################################## Computing IC ############################################
  IC <- log10(Dup / Ddn)
  IC <- IC * mask_cells
  
  pal <- colorRampPalette(c("blue", "lightyellow1", "red"))
  
  # Save preview PNG
  png(paste0("IC_", res, "m_", mw, "x", mw, "_", site_id, ".png"))
  plot(IC, col = pal(50), main = paste("IC -", site_id))
  dev.off()
  
  writeRaster(IC, paste0("IC_", res, "m_", mw, "x", mw, "_", site_id, ".tif"), overwrite = TRUE)
  
  # Clean temporary files if flag_k == 0
  if (flag_k == 0) {
    tmp_files <- c("accS.tif", "accW.tif", "ang.tif", "D8_flow_dir.tif", 
                   "sca.tif", "slope_D8.tif", "slope_D8_rec_m.tif", 
                   "slp.tif", "Weighted_Flow_length.tif", "W_factor.tif", "dem_filled.tif")
    suppressWarnings(file.remove(tmp_files[file.exists(tmp_files)]))
  }
  
  end_timeTOT <- Sys.time()
  cat(sprintf("Done with site %s in %.2f minutes.\n", site_id, as.numeric(difftime(end_timeTOT, start_timeTOT, units = "mins"))))
  
  return(IC)
}


###################################################################################################################
# RUN BATCH LOOP OVER ALL SITES
###################################################################################################################

data_IC <- batch_list
data_IC$mw        <- mw
data_IC$max_IC    <- NA
data_IC$min_IC    <- NA
data_IC$median_IC <- NA
data_IC$SD_IC     <- NA

for (i in 1:nrow(batch_list)) {
  
  gid <- batch_list$gauge_id[i]
  cat(sprintf("\n==================================================\nProcessing gauge [%d/%d]: %s\n==================================================\n", 
              i, nrow(batch_list), gid))
  
  tryCatch({
    # 1. Load DEM and NDVI rasters
    z_fill   <- raster(batch_list$dem_path[i])
    W_factor <- raster(batch_list$w_factor_path[i])
    
    # 2. Extract the single corresponding feature from the master shapefile
    #is_match <- (shp_gauge_ids == gid) | 
    #  (sprintf("%05d", suppressWarnings(as.numeric(shp_gauge_ids))) == gid)
    #target_feat <- all_targets[which(is_match)[1], ]
    
    is_match <- (shp_gauge_ids == gid) | 
      (sprintf("%05d", suppressWarnings(as.numeric(shp_gauge_ids))) == gid)
    target_feat <- all_targets[which(is_match)[1], ]
    
    # 3. Run IC calculation
    IC <- IC_function(z_fill = z_fill, target_feat = target_feat, W_factor = W_factor, mw = mw, flag_k = flag_k, site_id = gid)
    
    # 4. Extract summary statistics
    IC_val <- getValues(IC)
    data_IC[i, "max_IC"]    <- max(IC_val, na.rm = TRUE)
    data_IC[i, "min_IC"]    <- min(IC_val, na.rm = TRUE)
    data_IC[i, "median_IC"] <- median(IC_val, na.rm = TRUE)
    data_IC[i, "SD_IC"]     <- sd(IC_val, na.rm = TRUE)
    
  }, error = function(e) {
    cat(sprintf("ERROR encountered in gauge %s: %s\nSkipping to next gauge...\n", gid, e$message))
  })
}

# Export summary statistics table
write.csv(data_IC, file = file.path(wd, "stat_IC_batch.csv"), row.names = FALSE)

print("Congratulations!! Batch processing successfully completed.")