## Read results from spatial Bayesian latent factor model
read_sim_result <- function(scenario="Linear", resultpath){
  # result pathway
  resultpath_temp <- NULL
  if(scenario == "Linear"){
    resultpath_temp <- paste0(resultpath, "Linear/", sep="")
  }else if(scenario == "Voxelwise"){
    resultpath_temp <- paste0(resultpath, "Voxelwise/", sep="")
  }else if(scenario == "SpatialBayesLatent"){
    resultpath_temp <- paste0(resultpath, "SpatialBayesLatent/", sep="")
  }
  
  # read results
  pred_train <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Out_train.txt", sep=""), 
                                     quote="\"", comment.char="")) 
  pred_test <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Out_test.txt", sep=""), 
                                    quote="\"", comment.char="")) 
  if(scenario == "SpatialBayesLatent"){
    latent <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Latent.txt"), quote="\"", comment.char=""))
    out <- list(pred_train=pred_train, pred_test=pred_test, latent=latent)
  }else{
    out <- list(pred_train=pred_train, pred_test=pred_test)
  }

  return(out)
}

read_real_result <- function(resultpath){
  pred_train <- as.matrix(read.table(paste0(resultpath, "PostMean_Out_train.txt", sep=""), 
                                     quote="\"", comment.char="")) 
  pred_test <- as.matrix(read.table(paste0(resultpath, "PostMean_Out_test.txt", sep=""), quote="\"", comment.char="")) 
  
  latent <- as.matrix(read.table(paste0(resultpath, "PostMean_Latent.txt"), quote="\"", comment.char=""))
  
  return(list(pred_train=pred_train, pred_test=pred_test, latent=latent))
}

## Fit linear model
model_linear <- function(Xtrain, Ztrain, Xtest, Ztest){
  n_train = nrow(Xtrain)
  n_test = nrow(Xtest)
  n_predictors = ncol(Xtrain) / ncol(Ztrain)
  image_len = ncol(Ztrain)
  ## training
  beta_linear <- matrix(NA, nrow=n_train, ncol=n_predictors+1)
  pred_tr_linear <- matrix(NA, nrow=n_train, ncol=image_len)
  for(i in 1:n_train){ # only one parcel
    xi <- matrix(Xtrain[i, ], ncol = n_predictors, byrow=F)
    zi <- matrix(as.numeric(Ztrain[i, ]), ncol = 1)
    beta_linear[i, ] <- as.numeric(lm(zi~xi)$coef)
    pred_tr_linear[i, ] <- cbind(1, xi) %*% beta_linear[i, ]
  }
  beta_avg_linear <- matrix(colMeans(beta_linear), ncol=1)
  err_tr_linear <- pred_tr_linear - Ztrain
  mse_tr_linear <- mean(err_tr_linear^2) 
  
  ## test
  pred_ts_linear <- matrix(NA, nrow=n_test, ncol=image_len)
  for(i in 1:n_test){
    xi_new <- matrix(Xtest[i, ], ncol = n_predictors, byrow=F)
    pred_ts_linear[i, ] <- cbind(1, xi_new) %*% beta_avg_linear
  }
  err_ts_linear <- pred_ts_linear - Ztest
  mspe_ts_linear <- mean(err_ts_linear^2)
  
  err <- matrix(c(mse_tr_linear, mspe_ts_linear), ncol=2)
  colnames(err) <- c("MSE(training)", "MSPE(test)")
  return(list(err=err, pred_tr=pred_tr_linear, pred_ts=pred_ts_linear))
}

## Fit voxelwise regression model
model_voxel <- function(Xtrain, Ztrain, Xtest, Ztest){
  n_voxels <- ncol(Ztrain)
  n_predictors <- ncol(Xtrain) / n_voxels
  n_train <- nrow(Xtrain)
  n_test <- nrow(Xtest)
  pred_tr_voxel <- matrix(NA, ncol=n_voxels, nrow=n_train)
  pred_ts_voxel <- matrix(NA, ncol=n_voxels, nrow=n_test)
  for(i in 1:n_voxels){
    ## training
    zi <- Ztrain[, i]
    xi <- Xtrain[, ((0:(n_predictors-1))*n_voxels)+i]
    model <- lm(zi ~ xi)
    pred_tr_voxel[, i] <- predict(model, as.data.frame(xi))
    ## testing
    xi_new <- Xtest[, ((0:(n_predictors-1))*n_voxels)+i]
    pred_ts_voxel[, i] <- cbind(1, xi_new) %*% matrix(model$coef, ncol=1)
  }
  err_tr_voxel <- pred_tr_voxel - Ztrain
  mse_tr_voxel <- mean(err_tr_voxel^2) 
  err_ts_voxel <- pred_ts_voxel - Ztest
  mspe_ts_voxel <- mean(err_ts_voxel^2)
  
  err <- matrix(c(mse_tr_voxel, mspe_ts_voxel), ncol=2)
  colnames(err) <- c("MSE(training)", "MSPE(test)")
  return(list(err=err, pred_tr=pred_tr_voxel, pred_ts=pred_ts_voxel))
}

## Fit spatial Bayes latent factor model
model_bayes <- function(Xtrain, Ztrain, Xtest, Ztest,
                        scenario="Linear", resultpath){
  if(is.null(scenario)) {
    postmean <- read_real_result(resultpath)
  }else{
    postmean <- read_sim_result(scenario, resultpath)
  }
  
  err_tr_bayes <- postmean$pred_train - Ztrain
  mse_tr_bayes <- mean(err_tr_bayes^2)  
  
  err_ts_bayes <- postmean$pred_test - Ztest
  mspe_ts_bayes <- mean(err_ts_bayes^2)
  
  err <- matrix(c(mse_tr_bayes, mspe_ts_bayes), ncol=2)
  colnames(err) <- c("MSE(training)", "MSPE(test)")
  return(list(err=err, pred_tr=postmean$pred_train, pred_ts=postmean$pred_test, 
              latent=postmean$latent))
  
}

## Model comparison
# proportions of observations having less estimation/predition errors 
# fitted by model 2 to model 1
model_compare <- function(Ztrue, Zpred_m1, Zpred_m2){
  mse1 <- apply(Zpred_m1 - Ztrue, 1, function(x) mean(x^2))
  mse2 <- apply(Zpred_m2 - Ztrue, 1, function(x) mean(x^2))
  return (round((sum(mse1 >= mse2 ) / length(mse2)), 4))
}


## Distributions of Latent
elbow_latent <- function(latent, par=list(mfrow=c(2, 3))){
  # range
  ranges <- apply(latent, 2, function(x) max(x)-min(x))
  ranges_sort <- as.numeric(sort(ranges, decreasing = T))
  # max absolute value
  max_abs <- apply(latent, 2, function(x) max(abs(x)))
  max_abs_sort <- as.numeric(sort(max_abs, decreasing = T))
  # standard deviation
  sds <- apply(latent, 2, sd)
  sds_sort <- as.numeric(sort(sds, decreasing = T))
  # 95% CI
  cuts <- as.numeric(quantile(as.numeric(latent), probs = c(0.025, 0.975)))
  n_cuts <- apply(latent, 2, function(x) sum(x<cuts[1]|x>cuts[2]))
  n_cuts_sort <- as.numeric(sort(n_cuts, decreasing = T))
  # 90% CI
  cuts <- as.numeric(quantile(as.numeric(latent), probs = c(0.05, 0.95)))
  n_cuts <- apply(latent, 2, function(x) sum(x<cuts[1]|x>cuts[2]))
  n_cuts_sort <- as.numeric(sort(n_cuts, decreasing = T))
  # 68% CI
  cuts <- as.numeric(quantile(as.numeric(latent), probs = c(0.16, 0.84)))
  n_cuts <- apply(latent, 2, function(x) sum(x<cuts[1]|x>cuts[2]))
  n_cuts_sort <- as.numeric(sort(n_cuts, decreasing = T))
  
  if (!is.null(par)) { 
    on.exit(par(opar)) 
    opar <- par(par) 
  } 

  plot(ranges_sort, type="b", pch=16, 
       xlab="Latent Factors", ylab="", main="Range (max-min)")
  plot(max_abs_sort, type="b", pch=16, 
       xlab="Latent Factors", ylab="", main="Max absolute value")
  plot(sds_sort, type="b", pch=16, 
       xlab="Latent Factors", ylab="", main="Standard Deviation")
  plot(n_cuts_sort, type="b", pch=16,
       xlab="Latent Factor", ylab="", main="NO. of values outside 95% CI")
  plot(n_cuts_sort, type="b", pch=16,
       xlab="Latent Factor", ylab="", main="NO. of values outside 90% CI")
  plot(n_cuts_sort, type="b", pch=16,
       xlab="Latent Factor", ylab="", main="NO. of values outside 68% CI")

}



##### Figures #####
# colors
blue2red_cols = c("#000080", "#000083", "#000087", "#00008B", "#00008F", "#000093", "#000097", "#00009B",
                  "#00009F", "#0000A3", "#0000A7", "#0000AB", "#0000AF", "#0000B3", "#0000B7", "#0000BB",
                  "#0000BF", "#0000C3", "#0000C7", "#0000CB", "#0000CF", "#0000D3", "#0000D7", "#0000DB",
                  "#0000DF", "#0000E3", "#0000E7", "#0000EB", "#0000EF", "#0000F3", "#0000F7", "#0000FB",
                  "#0004FF", "#0008FF", "#000CFF", "#0010FF", "#0014FF", "#0018FF", "#001CFF", "#0020FF",
                  "#0024FF", "#0028FF", "#002CFF", "#0030FF", "#0034FF", "#0038FF", "#003CFF", "#0040FF",
                  "#0044FF", "#0048FF", "#004CFF", "#0050FF", "#0054FF", "#0058FF", "#005CFF", "#0060FF",
                  "#0064FF", "#0068FF", "#006CFF", "#0070FF", "#0074FF", "#0078FF", "#007CFF", "#0080FF",
                  "#0083FF", "#0087FF", "#008BFF", "#008FFF", "#0093FF", "#0097FF", "#009BFF", "#009FFF",
                  "#00A3FF", "#00A7FF", "#00ABFF", "#00AFFF", "#00B3FF", "#00B7FF", "#00BBFF", "#00BFFF",
                  "#00C3FF", "#00C7FF", "#00CBFF", "#00CFFF", "#00D3FF", "#00D7FF", "#00DBFF", "#00DFFF",
                  "#00E3FF", "#00E7FF", "#00EBFF", "#00EFFF", "#00F3FF", "#00F7FF", "#00FBFF", "#00FFFF",
                  "#04FFFB", "#08FFF7", "#0CFFF3", "#10FFEF", "#14FFEB", "#18FFE7", "#1CFFE3", "#20FFDF",
                  "#24FFDB", "#28FFD7", "#2CFFD3", "#30FFCF", "#34FFCB", "#38FFC7", "#3CFFC3", "#40FFBF",
                  "#44FFBB", "#48FFB7", "#4CFFB3", "#50FFAF", "#54FFAB", "#58FFA7", "#5CFFA3", "#60FF9F",
                  "#64FF9B", "#68FF97", "#6CFF93", "#70FF8F", "#74FF8B", "#78FF87", "#7CFF83", "#80FF80",
                  "#83FF7C", "#87FF78", "#8BFF74", "#8FFF70", "#93FF6C", "#97FF68", "#9BFF64", "#9FFF60",
                  "#A3FF5C", "#A7FF58", "#ABFF54", "#AFFF50", "#B3FF4C", "#B7FF48", "#BBFF44", "#BFFF40",
                  "#C3FF3C", "#C7FF38", "#CBFF34", "#CFFF30", "#D3FF2C", "#D7FF28", "#DBFF24", "#DFFF20",
                  "#E3FF1C", "#E7FF18", "#EBFF14", "#EFFF10", "#F3FF0C", "#F7FF08", "#FBFF04", "#FFFF00",
                  "#FFFB00", "#FFF700", "#FFF300", "#FFEF00", "#FFEB00", "#FFE700", "#FFE300", "#FFDF00",
                  "#FFDB00", "#FFD700", "#FFD300", "#FFCF00", "#FFCB00", "#FFC700", "#FFC300", "#FFBF00",
                  "#FFBB00", "#FFB700", "#FFB300", "#FFAF00", "#FFAB00", "#FFA700", "#FFA300", "#FF9F00",
                  "#FF9B00", "#FF9700", "#FF9300", "#FF8F00", "#FF8B00", "#FF8700", "#FF8300", "#FF8000",
                  "#FF7C00", "#FF7800", "#FF7400", "#FF7000", "#FF6C00", "#FF6800", "#FF6400", "#FF6000",
                  "#FF5C00", "#FF5800", "#FF5400", "#FF5000", "#FF4C00", "#FF4800", "#FF4400", "#FF4000",
                  "#FF3C00", "#FF3800", "#FF3400", "#FF3000", "#FF2C00", "#FF2800", "#FF2400", "#FF2000",
                  "#FF1C00", "#FF1800", "#FF1400", "#FF1000", "#FF0C00", "#FF0800", "#FF0400", "#FF0000",
                  "#FB0000", "#F70000", "#F30000", "#EF0000", "#EB0000", "#E70000", "#E30000", "#DF0000",
                  "#DB0000", "#D70000", "#D30000", "#CF0000", "#CB0000", "#C70000", "#C30000", "#BF0000",
                  "#BB0000", "#B70000", "#B30000", "#AF0000", "#AB0000", "#A70000", "#A30000", "#9F0000",
                  "#9B0000", "#970000", "#930000", "#8F0000", "#8B0000", "#870000", "#830000", "#800000")


# image functions
show.axial <- function(img, all_coords, slice_list,
                       col_lim=range(img, na.rm=T), cols=blue2red_cols, 
                       bgcol="white", main="", layout=NULL, x_lim=NULL, y_lim=NULL, 
                       cexaxis=1.5, cexcol=1.5, cexmain=1.5, cexlab=1.5, 
                       x_lab="", y_lab=""){
  all_coords <- as.data.frame(all_coords)
  names(all_coords) = c("x","y","z")
  idx <- which(is.element(all_coords$z, slice_list))
  z <- factor(paste0("z=",all_coords$z[idx], sep=""),
              levels=paste0("z=",sort(unique(all_coords$z[idx])), sep=""))
  
  if(is.null(x_lim)){
    x_lim = range(all_coords$x)
  }
  
  if(is.null(x_lim)){
    y_lim = range(all_coords$y)
  }
  
  fig <- levelplot(img~all_coords$x[idx]+all_coords$y[idx] | z, 
                   at=seq(col_lim[1],col_lim[2],length=length(cols)+1), 
                   col.regions = cols,cuts = length(cols), 
                   par.settings=list(panel.background=list(col=bgcol)), 
                   main=list(main, cex=cexmain), 
                   xlab=list(label=x_lab, cex=cexlab),
                   ylab=list(label=y_lab, cex=cexlab), 
                   aspect=1.0, layout=layout,
                   xlim=x_lim,ylim=y_lim, colorkey=list(labels=list(cex=cexcol)),
                   scales=list(x=list(cex=cexaxis),y=list(cex=cexaxis)))
  return(fig)
}

show.sagittal <- function(img, all_coords, slice_list,
                          col_lim=range(img, na.rm=T), cols=blue2red_cols, 
                          bgcol="white", main="", layout=NULL){
  all_coords <- as.data.frame(all_coords)
  names(all_coords) = c("x","y","z")
  idx <- which(is.element(all_coords$x, slice_list))
  x <- factor(paste0("x=",all_coords$x[idx], sep=""),
              levels=paste0("x=", sort(unique(all_coords$x[idx])), sep=""))
  
  fig <- levelplot(img~all_coords$y[idx]+all_coords$z[idx] | x, 
                   at=seq(col_lim[1],col_lim[2],length=length(cols)+1), 
                   col.regions = cols,cuts = length(cols), 
                   par.settings=list(panel.background=list(col=bgcol)), 
                   xlab="y", ylab="z", main=main, aspect=1.0, layout=layout, 
                   xlim=range(all_coords$y),ylim=range(all_coords$z))
  return(fig)
}

show.coronal <- function(img, all_coords, slice_list,
                         col_lim=range(img, na.rm=T), cols=blue2red_cols, 
                         bgcol="white", main="", layout=NULL){
  all_coords <- as.data.frame(all_coords)
  names(all_coords) = c("x","y","z")
  idx <- which(is.element(all_coords$y, slice_list))
  y <- factor(paste0("y=",all_coords$y[idx], sep=""),
              levels=paste0("y=",sort(unique(all_coords$y[idx])), sep=""))
  
  fig <- levelplot(img~all_coords$x[idx]+all_coords$z[idx] | y, 
                   at=seq(col_lim[1],col_lim[2],length=length(cols)+1), 
                   col.regions = cols,cuts = length(cols), 
                   par.settings=list(panel.background=list(col=bgcol)), 
                   xlab="x", ylab="z", main=main, aspect=1.0, layout=layout, 
                   xlim=range(all_coords$x),ylim=range(all_coords$z))
  return(fig)
}



