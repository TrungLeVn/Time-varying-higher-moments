rm(list = ls())
# Load libraries
source("Cryto/Code/utility.R")
#-----
theta = c(0.01,0.025,0.05)
noCores = detectCores()
print(noCores)
cluster <- makePSOCKcluster(noCores-1,outfile = "")
registerDoParallel(cluster)

horizonFor <- function(idxData,horizon = 1,As = TRUE,q = 0.01,windowLength = 1260){
  n = nrow(idxData)
  start = windowLength
  endpoint = seq(start,n,by = 10)
  if(As){
    estPars = matrix(NA,ncol = length(endpoint),nrow = 7)
  } else{
    estPars = matrix(NA,ncol = length(endpoint),nrow = 6)
  }
  forecast <- foreach(w = 1:length(endpoint),.combine = "rbind",
                      .packages = c("MidasQuantR","dplyr","foreach"),
                      .export = c("As","horizon","idxData","q","windowLength","idxData","theta"))%dopar%{
    d <- try(VarEs_jointAL_cav(y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]],
                               horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                               numInitialsRand = 30000, numInitials = 10,  As = As),silent = TRUE)
    if(inherits(d,"try-error")){
      ans = data_frame(date = idxData$date[seq(endpoint[w]-horizon+1,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon))
    } else if(d$conv == 1){
      # First try is to increase the numInitialsRand and numInitials
      d <- try(VarEs_jointAL_cav(y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]],
                                 horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                                 numInitialsRand = 30000, numInitials = 10,  As = As,
                                 MainSolver = "bobyqa",multiSol = FALSE),silent = TRUE)
      # Second try is to change the solvers
      if(inherits(d,"try-error")){
        ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon))
      } else if(d$conv == 1){
        d <- try(VarEs_jointAL_cav(y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]],
                                   horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                                   numInitialsRand = 30000, numInitials = 10,  As = As,
                                   MainSolver = "neldermead",multiSol = FALSE),silent = TRUE)
        if(inherits(d,"try-error")){
          ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon))
        } else if(d$conv == 1){
          ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon))
        } else{
          dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
          estPars[,w] = c(unname(d$meanFit$coef),d$estPars)
          ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)))
        }
      } else{
        dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
        estPars[,w] = c(unname(d$meanFit$coef),d$estPars)
        ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)))
      }
    } else{
      dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
      estPars[,w] = c(unname(d$meanFit$coef),d$estPars)
      ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)))
    }
    ans = ans
  }
  return(list(estPars = estPars, forecast = forecast))
}

al_as <- foreach(i = 1:4,.multicombine = T, .combine = "list")%do%{
  quantEst = foreach(ii = 1:3,.combine = "list",.multicombine = TRUE)%do%{
    idxData <- RawRet[,i]
    idxData <- as_tibble(idxData)
    idxData$date = as.Date(index(RawRet))
    colnames(idxData) = c("ret","date")
    d1for <- horizonFor(idxData = idxData, horizon = 1, As = TRUE, q = theta[ii])
  }
  names(quantEst) = as.character(theta)
  quantEst
}

saveRDS(al_as,"Cryto/Estimate/as.rds")

al_sav <- foreach(i = 1:4,.multicombine = T, .combine = "list")%do%{
  quantEst = foreach(ii = 1:3,.combine = "list",.multicombine = TRUE)%do%{
    idxData <- RawRet[,i]
    idxData <- as_tibble(idxData)
    idxData$date = as.Date(index(RawRet))
    colnames(idxData) = c("ret","date")
    d1for <- horizonFor(idxData = idxData, horizon = 1, As = FALSE, q = theta[ii])
  }
  names(quantEst) = as.character(theta)
  quantEst
}

saveRDS(al_sav,"Cryto/Estimate/sav.rds")