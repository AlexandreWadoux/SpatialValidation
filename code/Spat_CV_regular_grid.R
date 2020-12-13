############
############  Load packages (install from CRAN if necessary)
############

require(geoR)
require(caret)
require(raster)
require(sp)
require(rgdal)
require(pracma)
require(rgeos)
require(data.table)
require(ggplot2)
require(ranger)
require(spcosa)
require(readr)
require(gstat)
require(reshape2)

############
############  Initialize script
############

root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# make a raster stack 
files <- list.files(path = paste0(root_dir, '/data/'), recursive = FALSE, pattern = "\\.tif$")
s <- stack(paste0(root_dir, '/data/', files))

# remove water bodies 
s[[1]][s[[1]] ==0] <- NA

# convert to data frame
s_df <- as.data.frame(s, xy=T, na.rm=T)

############
############  Initialize script
############

# sample size
SampSize <- 500

# load functions
source('./Functions_spat_CV.R')

# load polygon of the area
load('./polygon.Rdata')

############
############  Generate semivariogram on AGB to define val.dist
############

# # take a large sample
# sp_s <- s_df[sample(nrow(s_df), 10000), ]
# 
# # generate sample variogram
# coordinates(sp_s) <- ~x+y
# lzn.vgm =variogram(ABG1~1, data = sp_s)
# plot(lzn.vgm$dist, lzn.vgm$gamma)

# Value chosen for distance
val.dist <- 350

############
############  Prepare for modelling
############

# set model predictors
predList_modelfull = c("AI_glob","CC_am","Clay","Elev","ETP_Glob","G_mean","NIR_mean","OCS","Prec_am","Prec_Dm","Prec_seaso",
                       "Prec_Wm","R_mean","Sand","Sha_EVI","Slope","Soc","solRad_m","SolRad_sd","SWIR1_mean","SWIR2_mean","T_am","T_mdq",
                       "T_mwarmq","T_seaso","Terra_PP","Vapor_m","Vapor_sd")

#### Set model response variable
response.name="ABG1"

# How many times a sample is selected?
nIter2 <- 500

## prepare output table
outtab.val.ME <- data.frame(Population = NA,
                            DesignBased = NA, 
                            RandonKFold = NA, 
                            SpatialKFold = NA, 
                            BLOOCV = NA, 
                            Stat = rep('ME', 1)) 
outtab.val.RMSE <- data.frame(Population = NA,
                              DesignBased = NA, 
                              RandonKFold = NA, 
                              SpatialKFold = NA, 
                              BLOOCV = NA, 
                              Stat = rep('RMSE', 1)) 
outtab.val.r2 <- data.frame(Population = NA,
                            DesignBased = NA, 
                            RandonKFold = NA, 
                            SpatialKFold = NA, 
                            BLOOCV = NA, 
                            Stat = rep('r2', 1)) 
outtab.val.MEC <- data.frame(Population = NA,
                             DesignBased = NA, 
                             RandonKFold = NA, 
                             SpatialKFold = NA, 
                             BLOOCV = NA, 
                             Stat = rep('MEC', 1)) 

for (sampling in 1:nIter2){
  
  ############
  ############  Generate the systematic (regualar) grid sample
  ############
  
  pt.reg <- spsample(pp, n=SampSize, type = 'regular')
  
  # extract covariate values and make valuetable
  valuetable.ptClus <- raster::extract(s, pt.reg, sp=T, df=T, na.rm = T)
  valuetable <- na.omit(as.data.frame(valuetable.ptClus))
  names(valuetable)[names(valuetable) == 'x1'] <- 'x'
  names(valuetable)[names(valuetable) == 'x2'] <- 'y'
  
  ############
  ############  Generate population validation statistics
  ############
  
  # formula for the RF model
  form_RF = as.formula(paste(response.name, "~", paste(predList_modelfull, collapse='+')))
  
  # build the RF model
  RF <- ranger(formula = form_RF, data = valuetable)
  
  # predict at all locations using the RF model
  pred_obj = predict(RF, data = s_df, type="response", )$predictions
  
  # compute population validation statistics
  eval(obs =s_df$ABG1 , pred = pred_obj)
  
  # store the population validation statistics
  outtab.val.ME$Population <- eval(obs =s_df$ABG1 , pred = pred_obj)$ME
  outtab.val.RMSE$Population <- eval(obs =s_df$ABG1 , pred = pred_obj)$RMSE
  outtab.val.r2$Population <- eval(obs =s_df$ABG1 , pred = pred_obj)$r2
  outtab.val.MEC$Population <- eval(obs =s_df$ABG1 , pred = pred_obj)$MEC
  
  
  ############
  ############  take probability sample for validation
  ############
  
  # make a simple random sample of size 500
  val_SRS <- as.data.frame(sampleRandom(s, size = 500))
  
  # build a RF model using the random sampling 
  RF <- ranger(formula = form_RF, data = valuetable)
  
  # predict at validation locations
  pred_obj = predict(RF, data = val_SRS, type="response")$predictions
  
  # store the design-based estimates of the population statistics
  outtab.val.ME$DesignBased <-  eval(obs =val_SRS$ABG1 , pred = pred_obj)$ME
  outtab.val.RMSE$DesignBased <-  eval(obs =val_SRS$ABG1 , pred = pred_obj)$RMSE
  outtab.val.r2$DesignBased <-  eval(obs =val_SRS$ABG1 , pred = pred_obj)$r2
  outtab.val.MEC$DesignBased <-  eval(obs =val_SRS$ABG1 , pred = pred_obj)$MEC
  
  ############
  ############ Prepare for spatial K-fold Cross-Validation (Spatial K-fold CV), following the implementation of Ploton et al., (2020)
  ############
  
  # make a distance matrix
  #set.seed(1)
  mdist <- dist(valuetable[c("x","y")])
  
  # hierachical clustering 
  hc <- hclust(mdist, method="complete")
  #plot(hc)
  d = val.dist * 1000  # the maximum distance between pixels within clusters (val.dist m * 1000 m = val.dist km)                   
  valuetable$Clust_val.distkm = cutree(hc, h=d) 
  #plot(valuetable$x, valuetable$y, col=as.factor(valuetable$Clust_val.distkm), xlab="Longitude", ylab="Latitude")
  
  
  ################################################################################################
  
  ############
  ############ Train and validate the AGB prediction model with Random and Spatial K-fold CVs
  ############ 
  
  #### Classical Random K-fold CV
  # create random folds
  flds <- createFolds(valuetable$ABG1, k = 10, list = TRUE, returnTrain = FALSE)  # with k = 10 the number of folds
  # perform CV
  for (j in 1:length(flds))
  {
    print(paste("PROCESSING CLUSTER n�", j, "out of", length(flds)))
    id = flds[[j]]
    training_data= valuetable[-id,]
    validation_data = valuetable[id,]  
    
    # formula
    form_RF = as.formula(paste(response.name, "~", paste(predList_modelfull, collapse='+')))
    
    # rf model
    RF <- ranger(formula = form_RF, data = training_data)
    
    # prediction
    pred_obj = predict(RF, data = validation_data, type="response")$predictions
    PRED.TEST = data.frame(AGB = validation_data$ABG1,
                           predRF = pred_obj)
    
    # store results
    if(j == 1) { PRED.TEST_tot = PRED.TEST}
    if(j > 1)  { PRED.TEST_tot = rbind(PRED.TEST_tot, PRED.TEST)}
    
  } # end of random K-fold CV loop
  # store results
  Res_random_CV = PRED.TEST_tot  
  # # plot
  # plot(Res_random_CV$AGB, Res_random_CV$predRF, pch=19, xlab="Observed AGB (in Mg/ha)", ylab="Predicted AGB (in Mg/ha)", main="Random K-fold CV")
  # eval(obs = Res_random_CV$AGB, pred = Res_random_CV$predRF)
  # abline(0,1)
  
  outtab.val.ME$RandonKFold <-  eval(obs = Res_random_CV$AGB, pred = Res_random_CV$predRF)$ME
  outtab.val.RMSE$RandonKFold <- eval(obs = Res_random_CV$AGB, pred = Res_random_CV$predRF)$RMSE
  outtab.val.r2$RandonKFold <-   eval(obs = Res_random_CV$AGB, pred = Res_random_CV$predRF)$r2
  outtab.val.MEC$RandonKFold <-   eval(obs = Res_random_CV$AGB, pred = Res_random_CV$predRF)$MEC
  
  
  #### Classical Spatial K-fold CV
  # create list of spatial folds
  flds = unique(valuetable$Clust_val.distkm)
  # perform CV
  for (j in 1:length(flds))
  {
    print(paste("PROCESSING CLUSTER n�", j, "out of", length(flds)))
    id = which(valuetable$Clust_val.distkm == levels(as.factor(flds))[j])
    training_data= valuetable[-id,]
    validation_data = valuetable[id,]  
    
    # formula
    form_RF = as.formula(paste(response.name, "~", paste(predList_modelfull, collapse='+')))
    
    # rf model
    RF <- ranger(formula = form_RF, data = training_data)
    
    # prediction
    pred_obj = predict(RF, data = validation_data, type="response")$predictions
    PRED.TEST = data.frame(AGB = validation_data$ABG1,
                           predRF = pred_obj)
    
    # store results
    if(j == 1) { PRED.TEST_tot = PRED.TEST}
    if(j > 1)  { PRED.TEST_tot = rbind(PRED.TEST_tot, PRED.TEST)}
    
  } # end of spatial K-fold CV loop 
  # store results
  Res_spatial_CV = PRED.TEST_tot  
  # # plot
  # plot(Res_spatial_CV$AGB, Res_spatial_CV$predRF, pch=19, xlab="Observed AGB (in Mg/ha)", ylab="Predicted AGB (in Mg/ha)", main="Spatial K-fold CV")
  # eval(obs = Res_spatial_CV$AGB, pred = Res_spatial_CV$predRF)
  # abline(0,1)
  
  outtab.val.ME$SpatialKFold <-  eval(obs = Res_spatial_CV$AGB, pred = Res_spatial_CV$predRF)$ME
  outtab.val.RMSE$SpatialKFold <-  eval(obs = Res_spatial_CV$AGB, pred = Res_spatial_CV$predRF)$RMSE
  outtab.val.r2$SpatialKFold <-  eval(obs = Res_spatial_CV$AGB, pred = Res_spatial_CV$predRF)$r2
  outtab.val.MEC$SpatialKFold <-  eval(obs = Res_spatial_CV$AGB, pred = Res_spatial_CV$predRF)$MEC
  
  
  ################################################################################################
  
  ############
  ############  Buffered Leave-One-Out CV (B-LOO CV), same as in Ploton et al., (2020)
  ############
  ####
  #### SECTION 1 : Generate Buffered-LOO CV output table
  ####
  
  #### Parameters of the analysis
  R_list = c(val.dist) # the list of r values (in kilometers). In the manuscript we used the following list : seq(0, val.dist, 10)
  NB_groups = 10       # For each r value, the number of groups of test pixels on which validation statistics will be computed (manuscript value : 10)
  NB_test_pixels = 100 # For each r value and group, the number of test pixels (manuscript value : 100)
  
  #### Prepare (empty) output table
  data_res = data.frame(R = NA,                              # Exclusion radius
                        N_training = NA,                     # "nb_training"
                        N_cell_lost = NA,                    # Number of pixels excluded
                        N_cell_training = NA,                # Number of training pixels conserved
                        AGB = NA,                            # "Observed" (i.e. field-derived) AGB of the test pixel
                        Pred_RF_FULL = NA,                   # AGB prediction of RF-RSE for the test pixel
                        AGB_calibrationDATA = NA             # Mean AGB of the training set (for null model)
  )
  a = 0  # "a" is just an incremental index used to store the results at the proper line of the output table "data_res"
  
  #### Models formulas
  response.name = "ABG1"
  form_RF_full = as.formula(paste(response.name, "~", paste(predList_modelfull, collapse='+')))  # Formula of the RF-RS model
  
  #### start of LOOP 1
  NB_iteration = NB_test_pixels*NB_groups
  for (j in 1:NB_iteration)     # from 1 to the number of test pixels in all groups
  {
    tic()
    
    Point_in_range = FALSE
    while(Point_in_range == FALSE)
    {
      ### STEP 1
      id_focal = sample(seq(1,nrow(valuetable),1),1)
      
      ### STEP 2
      # Generate the training set (i.e. remove test pixel)
      training_tmp = valuetable[-id_focal,]
      # Generate the test set
      focal_points = valuetable[id_focal,]
      # Set the longest exclusion radius size (in km)
      Ri = R_list[length(R_list)] * 1000
      # Generate the exclusion buffer with radius 'Ri'
      d = data.frame(lat = focal_points$y, lon = focal_points$x)
      coordinates(d) = ~lon + lat
      proj4string(d) = raster::crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
      buf = rgeos::gBuffer(d[1,], width = Ri)
      # Optional plot
      plot(buf)
      points(training_tmp$x, training_tmp$y, pch=19, col="grey", cex=1.2)
      plot(d, col="red", pch=19, cex=2, add=T)
      # Remove training pixels that fall within the exclusion buffer from the training set
      tmp = point.in.polygon(training_tmp$x, training_tmp$y, buf@polygons[[1]]@Polygons[[1]]@coords[,1], buf@polygons[[1]]@Polygons[[1]]@coords[,2])
      if(length(which(tmp == 1)) != 0) training_tmp = training_tmp[which(tmp != 1),]
      # Optional plot
      points(training_tmp$x, training_tmp$y, pch=19, col="orange", cex=1.2)
      # Check if the test pixelfall within the range of all predictors in the training set
      tmp_train = training_tmp[,names(training_tmp) %in% predList_modelfull]
      tmp_val = focal_points[,names(focal_points) %in% predList_modelfull]
      range_col = apply(tmp_train, 2, range)
      tmp_val2 = cbind(t(tmp_val), t(range_col))
      tmp_val2 = as.data.table(tmp_val2) ; names(tmp_val2) = c("focal", "lim_lw", "lim_up")
      if((nrow(tmp_val2) == nrow(tmp_val2[focal %between% list(lim_lw,lim_up)])) == FALSE) { print("RETRY test pixel selection : out of predictors' range")}
      if((nrow(tmp_val2) == nrow(tmp_val2[focal %between% list(lim_lw,lim_up)])) == TRUE)
      { # Record the number of pixels in the training dataset (at the longest r value)
        nb_training = nrow(training_tmp)
        # break the while loop
        Point_in_range = TRUE }
    }
    print(paste("<!> Test pixel number ", j, " selected <!>"))
    
    #### start of LOOP 2
    # Re-Generate training and test sets
    focal_points = valuetable[id_focal,]
    training_tmp_j = valuetable[-id_focal,]
    for (i in 1:length(R_list))
    {
      print(paste("PROCESSING exclusion radius of", R_list[i], "km"))
      
      ### STEPs 3 & 4
      # Set exclusion radius "i" from the exclusion radius list "R_list"
      Ri = R_list[i] * 1000
      # if Ri == 0, the exclusion buffer does not remove any training pixels at the vicinity of the test pixel
      # We just implement step 4
      if(Ri == 0) { training_tmp = training_tmp_j[sample(seq(1,nrow(training_tmp_j),1), nb_training),] }
      # if Ri > 0, generate the exclusion buffer with radius Ri & remove training pixels from the training set that fall within the exclusion buffer
      if(Ri > 0)
      {
        d = data.frame(lat = focal_points$y, lon = focal_points$x)
        coordinates(d) = ~lon + lat
        proj4string(d) = raster::crs("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
        buf = rgeos::gBuffer(d[1,], width = Ri)
        tmp = point.in.polygon(training_tmp_j$x, training_tmp_j$y, buf@polygons[[1]]@Polygons[[1]]@coords[,1], buf@polygons[[1]]@Polygons[[1]]@coords[,2])
        
        # if some training pixels fall in the exclusion buffer
        if(length(which(tmp == 1)) != 0)
        { # ... remove training pixels in the neighborhood of test pixel
          training_tmp = training_tmp_j[which(tmp != 1),]
          # ... and implement STEP 4
          nb_training_dif = nrow(training_tmp) - nb_training
          if(nb_training_dif > 0) { training_tmp = training_tmp[-sample(seq(1,nrow(training_tmp),1), nb_training_dif),] }
        }
        
        # if no training pixels fall in the exclusion buffer
        if(length(which(tmp == 1)) == 0)
        {  #... implement STEP 4
          training_tmp = training_tmp_j
          nb_training_dif = nrow(training_tmp) - nb_training
          if(nb_training_dif > 0) { training_tmp = training_tmp[-sample(seq(1,nrow(training_tmp),1), nb_training_dif),] }
        }
      }
      
      # Record the number of excluded training pixels (not used in the results)
      loss_cells = nrow(training_tmp_j) - nrow(training_tmp)
      
      ### STEP 4 & 5: RF
      
      # RF
      RF <- ranger(formula = form_RF_full, data = training_tmp)
      
      pred_RF_full = predict(RF, data = focal_points, type="response")$predictions
      
      ### STEP 6 : store results
      a = a+1
      data_res[a,] = NA
      data_res$R[a]    = R_list[i]
      data_res$N_training[a] = nb_training
      data_res$N_cell_lost[a] = loss_cells
      data_res$N_cell_training[a] = nrow(training_tmp)
      data_res$AGB[a] = focal_points$ABG1
      data_res$Pred_RF_FULL[a] = pred_RF_full
      data_res$AGB_calibrationDATA[a] = mean(training_tmp$ABG1)
      print(paste("Exclusion radius ", i, " DONE"))
    }      # end of LOOP 2
    toc()  # gives computation time per test pixel
    print(paste0("<!><!><!> R list on test pixel ", j, " is DONE <!><!><!>"))
  }  # end of LOOP 1
  
  ####
  #### SECTION 2 : Re-format the output table to compute validation statistics (RMSPE)
  ####
  
  ### Re-format the output table
  
  # List of exclusion buffer radius values (X axis)
  x = unique(data_res$R)
  
  # Generate empty vectors for mean and sd of R� by groups of test pixels (with n = 100 by groups)
  R2_RF_FULL_mean_list = c()
  # Generate empty vectors for mean and sd of RMSPE by groups of test pixels (with n = 100 by groups)
  RMSE_null_mean_list = c()
  RMSE_RF_FULL_mean_list = c()
  # Generate empty vectors for mean and sd of R� by groups of test pixels (with n = 100 by groups)
  ME_RF_FULL_mean_list = c()
  # Generate empty vectors for mean and sd of R� by groups of test pixels (with n = 100 by groups)
  MEC_RF_FULL_mean_list = c()
  
  for (i in 1:length(x))   # compute validation statistics by exclusion buffer radius r
  {
    # Select data corresponding to exlusion radius r at iteration i
    tmp = data_res[which(data_res$R == x[i]),]
    # Allocate processed test pixels to groups
    flds <- createFolds(row.names(tmp), k = NB_groups, list = F, returnTrain = FALSE)
    tmp = as.data.table(cbind(tmp, flds))
    # Generate empty vectors for mean and sd of R� by groups of test pixels FOR EXCLUSION RADIUS AT ITERATION i
    R2_RF_FULL_list = c()
    # Generate empty vectors for mean and sd of RMSPE by groups of test pixels FOR EXCLUSION RADIUS AT ITERATION i
    RMSE_null_list = c()
    RMSE_RF_FULL_list = c()
    ME_RF_FULL_list = c()
    MEC_RF_FULL_list = c()
    # Compute the R� and RMSPE by groups of test pixels
    for (j in 1:max(flds))  # from 1 to the number of groups
    {
      # RMSE Null model
      RMSE = eval(obs = tmp$AGB_calibrationDATA[which(tmp$flds == j)], pred = tmp$AGB[which(tmp$flds == j)])$RMSE
      RMSE_null_list = c(RMSE_null_list, RMSE)
      
      # RF-FULL
      map_quality_all = eval(obs = tmp$AGB[which(tmp$flds == j)], pred = tmp$Pred_RF_FULL[which(tmp$flds == j)])
      R2_RF_FULL_list = c(R2_RF_FULL_list, map_quality_all$r2)
      ME_RF_FULL_list = c(ME_RF_FULL_list,  map_quality_all$ME)
      RMSE_RF_FULL_list = c(RMSE_RF_FULL_list,  map_quality_all$RMSE)
      MEC_RF_FULL_list = c(MEC_RF_FULL_list,  map_quality_all$MEC)
      
    }
    # Generate mean and sd of R� over the groups
    R2_RF_FULL_mean_list = c(R2_RF_FULL_mean_list, mean(R2_RF_FULL_list))
    ME_RF_FULL_mean_list = c(ME_RF_FULL_mean_list, mean(ME_RF_FULL_list))
    MEC_RF_FULL_mean_list = c(MEC_RF_FULL_mean_list, mean(MEC_RF_FULL_list))
    
    # Generate mean and sd of RMSPE over the groups
    RMSE_null_mean_list = c(RMSE_null_mean_list, mean(RMSE_null_list))
    RMSE_RF_FULL_mean_list = c(RMSE_RF_FULL_mean_list, mean(RMSE_RF_FULL_list))
  }
  
  outtab.val.ME$BLOOCV <-  ME_RF_FULL_mean_list
  outtab.val.RMSE$BLOOCV <-  RMSE_RF_FULL_mean_list
  outtab.val.r2$BLOOCV <-  R2_RF_FULL_mean_list
  outtab.val.MEC$BLOOCV <- MEC_RF_FULL_mean_list
  
  outtab.val.ME.sub <- outtab.val.ME
  outtab.val.ME.sub[,c(2:5)] <- outtab.val.ME.sub[,c(2:5)] - outtab.val.ME.sub$Population
  outtab.val.ME.sub <- outtab.val.ME.sub[,-1]
  
  outtab.val.RMSE.sub <- outtab.val.RMSE
  outtab.val.RMSE.sub[,c(2:5)] <- outtab.val.RMSE.sub[,c(2:5)] - outtab.val.RMSE.sub$Population
  outtab.val.RMSE.sub <- outtab.val.RMSE.sub[,-1]
  
  outtab.val.r2.sub <- outtab.val.r2
  outtab.val.r2.sub[,c(2:5)] <- outtab.val.r2.sub[,c(2:5)] - outtab.val.r2.sub$Population
  outtab.val.r2.sub <- outtab.val.r2.sub[,-1]
  
  outtab.val.MEC.sub <- outtab.val.MEC
  outtab.val.MEC.sub[,c(2:5)] <- outtab.val.MEC.sub[,c(2:5)] - outtab.val.MEC.sub$Population
  outtab.val.MEC.sub <- outtab.val.MEC.sub[,-1]
  
  if (sampling==1){
    res.ME <- data.frame(DesignBased = NA, 
                         RandonKFold = NA, 
                         SpatialKFold = NA, 
                         BLOOCV = NA, 
                         Stat = rep('ME', nIter2))
    
    
    res.RMSE <- data.frame(DesignBased = NA, 
                           RandonKFold = NA, 
                           SpatialKFold = NA, 
                           BLOOCV = NA, 
                           Stat = rep('RMSE', nIter2))
    
    res.r2 <- data.frame(DesignBased = NA, 
                         RandonKFold = NA, 
                         SpatialKFold = NA, 
                         BLOOCV = NA, 
                         Stat = rep('r2', nIter2))
    
    res.MEC <- data.frame(DesignBased = NA, 
                          RandonKFold = NA, 
                          SpatialKFold = NA, 
                          BLOOCV = NA, 
                          Stat = rep('MEC', nIter2))
    
  }
  
  res.ME[sampling,][1:4] <- colMeans(outtab.val.ME.sub[,c(1:4)], na.rm = T)
  res.RMSE[sampling,][1:4] <- colMeans(outtab.val.RMSE.sub[,c(1:4)], na.rm = T)
  res.r2[sampling,][1:4] <- colMeans(outtab.val.r2.sub[,c(1:4)], na.rm = T)
  res.MEC[sampling,][1:4] <- colMeans(outtab.val.MEC.sub[,c(1:4)], na.rm = T)
  
  print(paste0(sampling, ' DONE, OUT OF ', nIter2))
  
  # save a the (temp, for each sampling) final results in a Rdata file
  save(res.ME, res.RMSE, res.r2, res.MEC, file = paste0('res_regular_testnew', sampling, '.Rdata'))
  
}

save(res.ME, res.RMSE, res.r2, res.MEC, file = 'res_regular_500.Rdata')


# plot the results of the validation statistics for the random sampling
outtab.val.ME2 <- melt(res.ME, 'Stat')
p1 <- ggplot(outtab.val.ME2) + geom_boxplot(aes(fill=variable, y = value, x = Stat))

outtab.val.RMSE2 <- melt(res.RMSE, 'Stat')
p2 <- ggplot(outtab.val.RMSE2) + geom_boxplot(aes(fill=variable, y = value, x = Stat))

outtab.val.r22 <- melt(res.r2, 'Stat')
p3 <- ggplot(outtab.val.r22) + geom_boxplot(aes(fill=variable, y = value, x = Stat))

outtab.val.MEC2 <- melt(res.MEC, 'Stat')
p4 <- ggplot(outtab.val.MEC2) + geom_boxplot(aes(fill=variable, y = value, x = Stat))

library(ggpubr)
ggarrange(p1, p2, p3, p4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE)


