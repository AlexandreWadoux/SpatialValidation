
## functions for the validation statistics
eval <- function(obs, pred, ...){
  
  # mean error
  ME <- round(mean(pred - obs, na.rm = TRUE), digits = 2)
  
  # root mean square error
  RMSE <-   round(sqrt(mean((pred - obs)^2, na.rm = TRUE)), digits = 2)
  
  # Pearson's correlation squared
  r2 <-  round((cor(pred, obs, method = 'spearman', use = 'pairwise.complete.obs')^2), digits = 2)
  
  # coefficient of determination - 'MEC' - 'AVE'
  SSE <- sum((pred - obs) ^ 2, na.rm = T)
  SST <- sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
  MEC <- round((1 - SSE/SST), digits = 2)
  
  return(data.frame(ME = ME, RMSE = RMSE, r2 = r2, MEC = MEC))
}

## function for the two-stage clustered sampling
twostage <- function(sframe,psu,n,m) {
  units <- sample.int(nrow(sframe),size=n,replace=TRUE)
  mypsusample <- sframe[units,psu]
  ssunits <- NULL
  for (psunit in mypsusample) {
    ssunit <- sample(x = which(sframe[,psu]==psunit),size=m,replace=FALSE)
    ssunits <- c(ssunits, ssunit)
  }
  psudraw <- rep(c(1:n),each=m)
  mysample <- data.frame(sframe[ssunits,],psudraw)
  return(mysample)
}