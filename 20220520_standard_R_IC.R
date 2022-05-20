# Title: Standard R_IC v. 1.0
# Author details: Tommaso Baggio, Lorenzo Martini, Loris Torresani

# Script and data info: This script performs the Index of Connectivity (Cavalli et al.,2013)

# Copyright statement: This script is the product of the work of Tommaso Baggio, Lorenzo Martini and Loris Torresani

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# Packages used
require("raster")
require("shapefiles")
require("rgdal")

# Premise: all files will be overwritten if they already exist

###################################
# INPUTS
# Set working directory to your location
# !!! WARNING !!!   do not write "/" as last character
wd <- ("D:/scripts/sed_connectivity/final_script/github_20220126/training_data/standard")
setwd(wd)

# input dem -> pit removed - user defined
z_fill=raster("vergini_5m_fill.tif"); plot(z_fill)

# input moving window size for roughness calculation
mw <- 3

# Input targets -> polygon shapefile - user defined
targets <- readOGR(dsn = wd, "target_vergini", stringsAsFactors = F); plot(targets, add=TRUE)

# Flag: save all outputs? 
# 0=keep essential files (IC, Ddup, Ddown)
# 1=keep all the files used
flag_k <- 0

###################################
# CHECK INPUT DATA
# !!! WARNING !!!   run until start of the script (line 54) to double check one by one the inputs

####Check if x and y cellsize are equal
resNS <- res(z_fill)[1]; resWE <- res(z_fill)[2]
ifelse(resNS==resWE, dtm_chk <- "DTM resolution ok", dtm_chk <- "Error!! NS and WE resolution are not equal")
remove(resNS,resWE)
#get raster and check if NS and EW resolution are equal
res <- res(z_fill)[1]

print(dtm_chk)

###################################################################################################################
# START OF THE SCRIPT
###################################################################################################################

start_timeTOT <- Sys.time()
################### Weighting factor

### smooth the dtm ###
dem_smooth <- focal(x = z_fill, fun=mean, w=matrix(1,mw,mw), pad=TRUE, na.rm=TRUE)
diff <- dem_smooth - z_fill

### Moving window algorithm to compute RI ###
crs_diff <- projection(diff)
ext_diff <- as.matrix(extent(diff))
diff_mat <- as.matrix(diff)

#Add row and column with NA
for (i in 1:mw) {
  coladd <- vector(mode = "numeric", length=nrow(diff_mat))
  coladd [1:length(coladd)] <- NA
  diff_mat <- cbind(coladd,diff_mat,coladd)
}

for (i in 1:mw) {
  rowadd <- vector(mode = "numeric", length=ncol(diff_mat))
  rowadd [1:length(rowadd)] <- NA
  diff_mat <- rbind(rowadd,diff_mat,rowadd)
}

#get position excluding NA values
XY_diff <- which(!is.na(diff_mat),arr.ind=T)
X <- XY_diff[,2]; Y <- XY_diff[,1]
RI_mat <- diff_mat

# Release memory
rm(XY_diff); gc()

start_time <- Sys.time()
for (i in 1:length(X)) { 
  pos_X <- X[i]
  pos_Y <- Y[i]
  val_MW <- (diff_mat[(((-mw+1)/2)+pos_Y):(((mw-1)/2)+pos_Y),(((-mw+1)/2)+pos_X):(((mw-1)/2)+pos_X)])
  val_MW <- val_MW[!is.na(val_MW)]
  mean_val_MW <- mean(val_MW)
  RI_mat[pos_Y,pos_X] <- sqrt(mean((val_MW-mean_val_MW)^2))
}
end_time <- Sys.time()
end_time - start_time

# remove the added rows and last columns
RI_mat <- RI_mat[-c(1:mw),]
for (i in 1:mw) {
  RI_mat <- RI_mat[-nrow(RI_mat),]
}
RI_mat <- RI_mat[,-c(1:mw)]
for (i in 1:mw) {
  RI_mat <- RI_mat[,-ncol(RI_mat)]
}

RI <- raster(RI_mat ,xmn=ext_diff[1,1],xmx=ext_diff[1,2],ymn=ext_diff[2,1],ymx=ext_diff[2,2], crs=crs_diff)
writeRaster(RI, filename = "RI.tif", overwrite=TRUE)

W_factor <- 1-(RI/(maxValue(RI))); plot(W_factor)
W_factor <- reclassify(W_factor, c(-1,0.001,0.001))
writeRaster(W_factor, filename = "W_factor.tif", overwrite=TRUE)

# Release memory
rm(diff_mat,RI_mat,diff,X,Y); gc()

######################################################################
# CREATE MASK FOR TARGETS

#polygon to raster
target_rst <- rasterize(targets, z_fill, background=2, field=1)
target_rst <- reclassify(target_rst, c(0.9,1.1,NA, 1.9,2.1,1))

######################################################################
# COMPUTING Ddown
writeRaster(z_fill, "dem_filled.tif", overwrite=TRUE)
#flow directions 8 directions
system("mpiexec -n 8 D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_filled.tif")
flow_dir_D8=raster("D8_flow_dir.tif")
slope_D8=raster("slope_D8.tif")

# integrate targets in flow directions
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8*target_rst


# reclassified values < 0.005 m/m and > 1 m/m
slope_D8_r_m <- reclassify(slope_D8, c(-2,0.005,0.005, 1,10^50,1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite=TRUE)

W_factor_r <- reclassify(W_factor, c(-2,0.001,0.001)) 

inv_CS <- 1/(W_factor_r * slope_D8_r_m)
writeRaster(inv_CS, "inv_CS.tif", overwrite=TRUE)

# Release memory
rm(W_factor_r);gc()
######################### FLOW WEIGHTED LENGTH
# Inputs

dir8 <- flow_dir_D8
W_FL <- inv_CS

# Start of the algorithm
W_FL_m <- as.matrix(W_FL)
dir8_m <- as.matrix(dir8)

############################################
## Add a row up and a column on the right
coladd <- vector(mode = "numeric", length=nrow(dir8_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd,dir8_m,coladd)
W_FL_m <- cbind(coladd,W_FL_m,coladd)

rowadd <- vector(mode = "numeric", length=ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd,dir8_m,rowadd)
W_FL_m <- rbind(rowadd,W_FL_m,rowadd)

start_time <- Sys.time()
# initialize flow length for cells draining to No DATA assigning 0
####################

# make a copy of the weights
Wgt_m <- W_FL_m 

NA_data <- which(is.na(dir8_m),arr.ind=T); NA_data <- as.data.frame(NA_data)
posrow <- NA_data$row; poscol <- NA_data$col

# Release memory
gc()

########## Flow direction 1
pos1col_NA <- poscol - 1
pos1row_NA <- posrow
pos1row <- pos1row_NA[-pos1col_NA!=0]; pos1col <- pos1col_NA[-pos1col_NA!=0]
pos1dir <- dir8_m[cbind(pos1row, pos1col)]
pos1row_0 <- pos1row[which(pos1dir==1)]; pos1col_0 <- pos1col[which(pos1dir==1)]
#pos1row_0 <- pos1row_0; pos1col_0 <- pos1col_0+1
W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0

########## Flow direction 2
pos2col_NA <- poscol - 1
pos2row_NA <- posrow + 1
pos2rowt <- pos2row_NA[-pos2col_NA!=0]; pos2colt <- pos2col_NA[-pos2col_NA!=0]
pos2row <- pos2rowt[pos2rowt<=(nrow(dir8_m))]; pos2col <- pos2colt[pos2rowt<=(nrow(dir8_m))]
pos2dir <- dir8_m[cbind(pos2row, pos2col)]
pos2row_0 <- pos2row[which(pos2dir==2)]; pos2col_0 <- pos2col[which(pos2dir==2)]
#pos2row_0 <- pos2row_0-1; pos3col_0 <- pos2col_0+1
W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0

########## Flow direction 3
pos3col_NA <- poscol
pos3row_NA <- posrow + 1
pos3row <- pos3row_NA[pos3row_NA<=(nrow(dir8_m))]; pos3col <- pos3col_NA[pos3row_NA<=(nrow(dir8_m))]
pos3dir <- dir8_m[cbind(pos3row, pos3col)]
pos3row_0 <- pos3row[which(pos3dir==3)]; pos3col_0 <- pos3col[which(pos3dir==3)]
#pos3row_0 <- pos3row_0-1; pos3col_0 <- pos3col_0
W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0

########## Flow direction 4
pos4col_NA <- poscol + 1
pos4row_NA <- posrow + 1
pos4rowt <- pos4row_NA[pos4col_NA<=(ncol(dir8_m))]; pos4colt <- pos4col_NA[pos4col_NA<=(ncol(dir8_m))]
pos4row <- pos4rowt[pos4rowt<=(nrow(dir8_m))]; pos4col <- pos4colt[pos4rowt<=(nrow(dir8_m))]
pos4dir <- dir8_m[cbind(pos4row, pos4col)]
pos4row_0 <- pos4row[which(pos4dir==4)]; pos4col_0 <- pos4col[which(pos4dir==4)]
#pos4row_0 <- pos4row_0-1; pos4col_0 <- pos4col_0-1
W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0

########## Flow direction 5
pos5col_NA <- poscol + 1
pos5row_NA <- posrow
pos5row <- pos5row_NA[pos5col_NA<=(ncol(dir8_m))]; pos5col <- pos5col_NA[pos5col_NA<=(ncol(dir8_m))]
pos5dir <- dir8_m[cbind(pos5row, pos5col)]
pos5row_0 <- pos5row[which(pos5dir==5)]; pos5col_0 <- pos5col[which(pos5dir==5)]
#pos5row_0 <- pos5row_0; pos5col_0 <- pos5col_0-1
W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0

########## Flow direction 6
pos6col_NA <- poscol + 1
pos6row_NA <- posrow - 1
pos6rowt <- pos6row_NA[pos6col_NA<=(ncol(dir8_m))]; pos6colt <- pos6col_NA[pos6col_NA<=(ncol(dir8_m))]
pos6row <- pos6rowt[-pos6rowt!=0]; pos6col <- pos6colt[-pos6rowt!=0]
pos6dir <- dir8_m[cbind(pos6row, pos6col)]
pos6row_0 <- pos6row[which(pos6dir==6)]; pos6col_0 <- pos6col[which(pos6dir==6)]
#pos6row_0 <- pos6row_0+1; pos6col_0 <- pos6col_0-1
W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0

########## Flow direction 7
pos7col_NA <- poscol
pos7row_NA <- posrow - 1
pos7row <- pos7row_NA[-pos7row_NA!=0]; pos7col <- pos7col_NA[-pos7row_NA!=0]
pos7dir <- dir8_m[cbind(pos7row, pos7col)]
pos7row_0 <- pos7row[which(pos7dir==7)]; pos7col_0 <- pos7col[which(pos7dir==7)]
#pos7row_0 <- pos7row_0+1; pos7col_0 <- pos7col_0
W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0 

########## Flow direction 8
pos8col_NA <- poscol - 1
pos8row_NA <- posrow - 1
pos8rowt <- pos8row_NA[-pos8row_NA!=0]; pos8colt <- pos8col_NA[-pos8row_NA!=0]
pos8row <- pos8rowt[-pos8colt!=0]; pos8col <- pos8colt[-pos8colt!=0]
pos8dir <- dir8_m[cbind(pos8row, pos8col)]
pos8row_0 <- pos8row[which(pos8dir==8)]; pos8col_0 <- pos8col[which(pos8dir==8)]
#pos8row_0 <- pos8row_0+1; pos8col_0 <- pos8col_0+1
W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0


## remove non necessary vectors
rm(pos1col, pos1col_NA, pos1row, pos1row_NA, pos1dir)
rm(pos3col, pos3col_NA, pos3row, pos3row_NA, pos3dir)
rm(pos5col, pos5col_NA, pos5row, pos5row_NA, pos5dir)
rm(pos7col, pos7col_NA, pos7row, pos7row_NA, pos7dir)
rm(pos2col, pos2col_NA, pos2row, pos2row_NA, pos2colt, pos2rowt, pos2dir)
rm(pos4col, pos4col_NA, pos4row, pos4row_NA, pos4colt, pos4rowt, pos4dir)
rm(pos6col, pos6col_NA, pos6row, pos6row_NA, pos6colt, pos6rowt, pos6dir)
rm(pos8col, pos8col_NA, pos8row, pos8row_NA, pos8colt, pos8rowt, pos8dir)


## Compute weighted sum following the direction 
######################
count=1
while ( (length(pos1col_0) | length(pos2col_0) | length(pos3col_0) | length(pos4col_0) | length(pos5col_0) | length(pos6col_0) | length(pos7col_0) | length(pos8col_0)) > 0) {
  
  # check the moves
  pos1row_0_bis <- pos1row_0; pos1col_0_bis <- pos1col_0+1
  pos2row_0_bis <- pos2row_0-1; pos2col_0_bis <- pos2col_0+1
  pos3row_0_bis <- pos3row_0-1; pos3col_0_bis <- pos3col_0
  pos4row_0_bis <- pos4row_0-1; pos4col_0_bis <- pos4col_0-1
  pos5row_0_bis <- pos5row_0; pos5col_0_bis <- pos5col_0-1
  pos6row_0_bis <- pos6row_0+1; pos6col_0_bis <- pos6col_0-1
  pos7row_0_bis <- pos7row_0+1; pos7col_0_bis <- pos7col_0
  pos8row_0_bis <- pos8row_0+1; pos8col_0_bis <- pos8col_0+1
  
  
  if (count==1) {
    W_FL_m[cbind(pos1row_0,pos1col_0)] <- 0
    W_FL_m[cbind(pos2row_0,pos2col_0)] <- 0
    W_FL_m[cbind(pos3row_0,pos3col_0)] <- 0
    W_FL_m[cbind(pos4row_0,pos4col_0)] <- 0
    W_FL_m[cbind(pos5row_0,pos5col_0)] <- 0
    W_FL_m[cbind(pos6row_0,pos6col_0)] <- 0
    W_FL_m[cbind(pos7row_0,pos7col_0)] <- 0
    W_FL_m[cbind(pos8row_0,pos8col_0)] <- 0
  } 
  
  else {
    W_FL_m[cbind(pos1row_0,pos1col_0)] <- W_FL_m[cbind(pos1row_0_bis,pos1col_0_bis)] + (res*((Wgt_m[cbind(pos1row_0_bis,pos1col_0_bis)] + Wgt_m[cbind(pos1row_0,pos1col_0)])/2))
    W_FL_m[cbind(pos2row_0,pos2col_0)] <- W_FL_m[cbind(pos2row_0_bis,pos2col_0_bis)] + (res*sqrt(2)*((Wgt_m[cbind(pos2row_0_bis,pos2col_0_bis)] + Wgt_m[cbind(pos2row_0,pos2col_0)])/2))
    W_FL_m[cbind(pos3row_0,pos3col_0)] <- W_FL_m[cbind(pos3row_0_bis,pos3col_0_bis)] + (res*((Wgt_m[cbind(pos3row_0_bis,pos3col_0_bis)] + Wgt_m[cbind(pos3row_0,pos3col_0)])/2))
    W_FL_m[cbind(pos4row_0,pos4col_0)] <- W_FL_m[cbind(pos4row_0_bis,pos4col_0_bis)] + (res*sqrt(2)*((Wgt_m[cbind(pos4row_0_bis,pos4col_0_bis)] + Wgt_m[cbind(pos4row_0,pos4col_0)])/2))
    W_FL_m[cbind(pos5row_0,pos5col_0)] <- W_FL_m[cbind(pos5row_0_bis,pos5col_0_bis)] + (res*((Wgt_m[cbind(pos5row_0_bis,pos5col_0_bis)] + Wgt_m[cbind(pos5row_0,pos5col_0)])/2))
    W_FL_m[cbind(pos6row_0,pos6col_0)] <- W_FL_m[cbind(pos6row_0_bis,pos6col_0_bis)] + (res*sqrt(2)*((Wgt_m[cbind(pos6row_0_bis,pos6col_0_bis)] + Wgt_m[cbind(pos6row_0,pos6col_0)])/2))
    W_FL_m[cbind(pos7row_0,pos7col_0)] <- W_FL_m[cbind(pos7row_0_bis,pos7col_0_bis)] + (res*((Wgt_m[cbind(pos7row_0_bis,pos7col_0_bis)] + Wgt_m[cbind(pos7row_0,pos7col_0)])/2))
    W_FL_m[cbind(pos8row_0,pos8col_0)] <- W_FL_m[cbind(pos8row_0_bis,pos8col_0_bis)] + (res*sqrt(2)*((Wgt_m[cbind(pos8row_0_bis,pos8col_0_bis)] + Wgt_m[cbind(pos8row_0,pos8col_0)])/2))
  }
  
  #reconstruting the position and moving up
  
  pos0row_bis <- c(pos1row_0, pos2row_0, pos3row_0, pos4row_0, pos5row_0, pos6row_0, pos7row_0, pos8row_0)
  pos0col_bis <- c(pos1col_0, pos2col_0, pos3col_0, pos4col_0, pos5col_0, pos6col_0, pos7col_0, pos8col_0)
  
  ########## Flow direction 1
  poscolt <- pos0col_bis - 1
  posrowt <- pos0row_bis
  posrow <- posrowt[-poscolt!=0]; poscol <- poscolt[-poscolt!=0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos1row_0 <- posrow[which(posdir==1)]; pos1col_0 <- poscol[which(posdir==1)]
  
  ########## Flow direction 2
  poscolt <- pos0col_bis - 1
  posrowt <- pos0row_bis + 1
  posrow <- posrowt[-poscolt!=0]; poscol <- poscolt[-poscolt!=0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos2row_0 <- posrow[which(posdir==2)]; pos2col_0 <- poscol[which(posdir==2)]
  
  ########## Flow direction 3
  poscol <- pos0col_bis
  posrow <- pos0row_bis + 1
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos3row_0 <- posrow[which(posdir==3)]; pos3col_0 <- poscol[which(posdir==3)]
  
  ########## Flow direction 4
  poscol <- pos0col_bis + 1
  posrow <- pos0row_bis + 1
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos4row_0 <- posrow[which(posdir==4)]; pos4col_0 <- poscol[which(posdir==4)]
  
  ########## Flow direction 5
  poscol <- pos0col_bis + 1
  posrow <- pos0row_bis
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos5row_0 <- posrow[which(posdir==5)]; pos5col_0 <- poscol[which(posdir==5)]
  
  ########## Flow direction 6
  poscolt <- pos0col_bis + 1
  posrowt <- pos0row_bis - 1
  posrow <- posrowt[-posrowt!=0]; poscol <- poscolt[-posrowt!=0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos6row_0 <- posrow[which(posdir==6)]; pos6col_0 <- poscol[which(posdir==6)]
  
  ########## Flow direction 7
  poscolt <- pos0col_bis
  posrowt <- pos0row_bis - 1
  posrow <- posrowt[-posrowt!=0]; poscol <- poscolt[-posrowt!=0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos7row_0 <- posrow[which(posdir==7)]; pos7col_0 <- poscol[which(posdir==7)]
  
  ########## Flow direction 8
  poscolt <- pos0col_bis - 1
  posrowt <- pos0row_bis - 1
  posrowt1 <- posrowt[-posrowt!=0]; poscolt1 <- poscolt[-posrowt!=0]
  posrow <- posrowt1[-poscolt1!=0]; poscol <- poscolt1[-poscolt1!=0]
  posdir <- dir8_m[cbind(posrow, poscol)]
  pos8row_0 <- posrow[which(posdir==8)]; pos8col_0 <- poscol[which(posdir==8)]
  
  count = count + 1
}

# Release memory
gc()

end_time <- Sys.time()
end_time - start_time
# remove the added first row and last column
W_FL_m <- W_FL_m[-1,]
W_FL_m <- W_FL_m[-nrow(W_FL_m),]
dir8_m <- dir8_m[-1,]
dir8_m <- dir8_m[-nrow(dir8_m),]

W_FL_m <- W_FL_m[,-1]
W_FL_m <- W_FL_m[,-ncol(W_FL_m)] 
dir8_m <- dir8_m[,-1]
dir8_m <- dir8_m[,-ncol(dir8_m)]


crs <- projection(W_FL)
ext <- as.matrix(extent(W_FL))

W_FL_fin <- raster(W_FL_m ,xmn=ext[1,1],xmx=ext[1,2],ymn=ext[2,1],ymx=ext[2,2], crs=crs); plot(W_FL_fin)
writeRaster(W_FL_fin, "Weighted_Flow_length.tiff", overwrite=TRUE)


#remove cells
dif_WFL <- W_FL_fin - inv_CS
dif_WFL[dif_WFL==0] <- NA
mask_cells <- dif_WFL*0 + 1

##### Computing final Downslope component
Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
writeRaster(Ddn, "Ddn.tif", overwrite=TRUE)

# Release memory
rm(Wgt_m,W_FL_m,dir8_m,dif_WFL,W_FL,dir8,inv_CS,flow_dir_D8);gc()

######################################################################
# COMPUTING Dup

#computing d-inf directions
system("mpiexec -n 8 DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_filled.tif")
ang=raster("ang.tif")

##
# integrating targets
ang <- ang *target_rst
writeRaster(ang, "ang.tif", overwrite=TRUE)

# Release memory
gc()

##
slp_dinf=raster("slp.tif")


#computing d-inf-area weighted (acc S) 
system("mpiexec -n 8 AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc")
accS=raster("accS.tif")

#computing d-inf-area
system("mpiexec -n 8 AreaDinf -ang ang.tif -sca sca.tif -nc")
sca=raster("sca.tif")

constant_res <- z_fill*0 + res

Acc_fin <- sca / constant_res
writeRaster(Acc_fin, "Acc_fin.tif", overwrite=TRUE)

smean <- (accS + slope_D8_r_m) / Acc_fin

# Release memory
rm(slope_D8_r_m,accS,target_rst); gc()


##computing d-inf-area weighted for W_factor(acc W) 
system("mpiexec -n 8 AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc")
accW=raster("accW.tif") * res

cmean <- (accW + W_factor) / Acc_fin


### computing final Upslope component
Dup <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
writeRaster(Dup, "Dup.tif", overwrite=TRUE)

# Release memory
rm(constant_res,ang,Acc_fin,smean,accW,cmean); gc()

################################## Computing IC ############################################

IC <- log10(Dup/Ddn)
IC <- IC * mask_cells

pal <- colorRampPalette(c("blue","lightyellow1","red"))

# Save IC plot as png image
png("IC.png")
plot(IC,col = pal(50))
dev.off()

plot(IC,col = pal(50))

writeRaster(IC, paste("IC", res, "m", mw,"x",mw, ".tif", sep = "_"), overwrite=TRUE)

# remove files

if (flag_k==0) {
  file.remove("accS.tif","accW.tif","ang.tif","D8_flow_dir.tif","sca.tif","slope_D8.tif","slope_D8_rec_m.tif","slp.tif","Weighted_Flow_length.tif", "W_factor.tif","dem_filled.tif","Acc_fin.tif","IC.tif")
}

end_timeTOT <- Sys.time()
end_timeTOT - start_timeTOT

print("Congratulations!! Script completed")
