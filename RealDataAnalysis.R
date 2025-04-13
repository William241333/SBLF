# Clean
rm(list=ls(all=TRUE))


# Library
library(ggplot2)
library(reshape)
library(dplyr)
library(lattice)
library(gridExtra)

# Pathways
mainpath <- "E:/c++code/"
datpath <- paste0(mainpath, "Data/", sep="")
resultpath <- paste0(mainpath, "Result20-8/", sep="")

# Functions
source('E:/c++code/biom13420-sup-0002-code/SouceCode/RCode/Functions.R', echo=TRUE)

# Input
Xtrain <- as.matrix(read.table(paste0(datpath, "ROI_dat_mat.txt"), quote="\"", comment.char=""))
Xtest <- as.matrix(read.table(paste0(datpath, "ROI_dat_mat_test.txt"), quote="\"", comment.char=""))
Ztrain <- as.matrix(read.table(paste0(datpath, "ROI_task.txt"), quote="\"", comment.char=""))
Ztest <- as.matrix(read.table(paste0(datpath, "ROI_task_test.txt"), quote="\"", comment.char=""))
voxel_loc <- as.matrix(read.table(paste0(datpath, "ROI_axes.txt"), quote="\"", comment.char=""))

image_len <- ncol(Ztrain) # number of voxels per image
n_predictors <- ncol(Xtrain) / image_len # number of imaging predictors
n_train <- nrow(Xtrain)
n_test <- nrow(Xtest)

# Model fit
res_linear <- model_linear(Xtrain, Ztrain, Xtest, Ztest)
res_voxel <- model_voxel(Xtrain, Ztrain, Xtest, Ztest)
res_bayes <- model_bayes(Xtrain, Ztrain, Xtest, Ztest, scenario=NULL, resultpath)

# MSE, MSPE
model_names <- c("Linear", "Voxelwise", "SpatialBayesLatent")
mse_mspe <- rbind(res_linear$err, res_voxel$err, res_bayes$err)
rownames(mse_mspe) = model_names
mse_mspe


model_compare <- function(Ztrue, Zpred_m1, Zpred_m2){
  mse1 <- apply(Zpred_m1 - Ztrue, 1, function(x) mean(x^2))
  mse2 <- apply(Zpred_m2 - Ztrue, 1, function(x) mean(x^2))
  return (round((sum(mse1 >= mse2 ) / length(mse2)), 4))
}


# Model comparison
proportions <- c(
  model_compare(Ztrain, res_linear$pred_tr, res_bayes$pred_tr),
  model_compare(Ztest, res_linear$pred_ts, res_bayes$pred_ts),
  model_compare(Ztrain, res_voxel$pred_tr, res_bayes$pred_tr),
  model_compare(Ztest, res_voxel$pred_ts, res_bayes$pred_ts)
)


proportions

## Figures
slice_z <- sort(unique(voxel_loc[, 3])) # z-axis values
predictor_idx <- sample(1:n_predictors, 1) # predictor id

# Axial view, subject from training set
obs_tr_idx <- sample(n_train, 1) # subject id from training set
fig_dat <- Ztrain[obs_tr_idx, ]
fig_dat2 <- res_bayes$pred_tr[obs_tr_idx, ]
title <- paste0("Observed Left Amygdala Image (Axial View)\nPredictor ID=", predictor_idx)
show.axial(fig_dat,  all_coords=voxel_loc, 
           slice_list=slice_z, 
           col_lim = range(fig_dat),
           layout=c(length(slice_z), 1), 
           x_lim = range(voxel_loc[, 1]), 
           y_lim = range(voxel_loc[, 2]), 
           cexaxis=1.0, cexcol=1.1, cexmain=1.0, cexlab=1.3,
           x_lab="", y_lab="", main=title)

title <- paste0("Estimated Left Amygdala Image (Axial View)\nPredictor ID=", predictor_idx)
show.axial(fig_dat2,  all_coords=voxel_loc, 
           slice_list=slice_z, col_lim = range(fig_dat2),
           layout=c(length(slice_z), 1), 
           x_lim = range(voxel_loc[, 1]), 
           y_lim = range(voxel_loc[, 2]), 
           cexaxis=1.0, cexcol=1.1, cexmain=1.0, cexlab=1.3,
           x_lab="", y_lab="", main=title)

# Axial view, subject from test set
obs_ts_idx <- sample(n_test, 1) # subject id from training set
fig_dat <- Ztest[obs_ts_idx, ]
fig_dat2 <- res_bayes$pred_ts[obs_ts_idx, ]
title <- paste0("Observed Left Amygdala Image (Axial View)\nPredictor ID=", predictor_idx)
par(mfrow=c(2, 1))
show.axial(fig_dat,  all_coords=voxel_loc, 
           slice_list=slice_z, 
           col_lim = range(fig_dat),
           layout=c(length(slice_z), 1), 
           x_lim = range(voxel_loc[, 1]), 
           y_lim = range(voxel_loc[, 2]), 
           cexaxis=1.0, cexcol=1.1, cexmain=1.0, cexlab=1.3,
           x_lab="", y_lab="", main=title)

title <- paste0("Predicted Left Amygdala Image (Axial View)\nPredictor ID=", predictor_idx)
show.axial(fig_dat2,  all_coords=voxel_loc, 
           slice_list=slice_z, 
           col_lim = range(fig_dat2),
           layout=c(length(slice_z), 1), 
           x_lim = range(voxel_loc[, 1]), 
           y_lim = range(voxel_loc[, 2]), 
           cexaxis=1.0, cexcol=1.1, cexmain=1.0, cexlab=1.3,
           x_lab="", y_lab="", main=title)
