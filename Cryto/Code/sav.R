rm(list = ls())
# Load libraries
source("Cryto/Code/utility.R")
#-----
theta = c(0.01,0.025,0.05)
noCores = detectCores()
print(noCores)
cluster <- makePSOCKcluster(noCores-1,outfile = "")
registerDoParallel(cluster)

horizonFor <- function(idxData,horizon = 1,As = TRUE,q = 0.01,windowLength = 1260,B = 1000){
  n = nrow(idxData)
  start = windowLength
  endpoint = seq(start,n,by = 10)
  forecast <- foreach(w = 1:length(endpoint),.combine = "list",.multicombine = TRUE,
                      .packages = c("MidasQuantR","dplyr","foreach"),
                      .export = c("As","horizon","idxData","q","windowLength","idxData"))%dopar%{
                        y = idxData$ret[1:endpoint[w]]; yDate = idxData$date[1:endpoint[w]]
                        d <- try(VarEs_jointAL_cav(y = y,yDate = yDate,horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                                                   numInitialsRand = 30000, numInitials = 10,  As = As),silent = TRUE)
                        if(inherits(d,"try-error")){
                          ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon),
                                           uEst = rep(NA,10/horizon))
                        } else if(d$conv == 1){
                          # First try is to increase the numInitialsRand and numInitials
                          d <- try(VarEs_jointAL_cav(y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]],
                                                     horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                                                     numInitialsRand = 30000, numInitials = 10,  As = As,
                                                     MainSolver = "bobyqa",multiSol = FALSE),silent = TRUE)
                          # Second try is to change the solvers
                          if(inherits(d,"try-error")){
                            ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon),
                                             uEst = rep(NA,10/horizon))
                          } else if(d$conv == 1){
                            d <- try(VarEs_jointAL_cav(y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]],
                                                       horizon = horizon, armaOrder = c(1,0,0), q = q, GetSe = FALSE, forecastLength = 10,
                                                       numInitialsRand = 30000, numInitials = 10,  As = As,
                                                       MainSolver = "neldermead",multiSol = FALSE),silent = TRUE)
                            if(inherits(d,"try-error")){
                              ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon),
                                               uEst = rep(NA,10/horizon))
                            } else if(d$conv == 1){
                              ans = data_frame(date = idxData$date[seq(endpoint[w]-9,endpoint[w],by = horizon)],VaR = rep(NA,10/horizon),ES = rep(NA,10/horizon),
                                               uEst = rep(NA,10/horizon))
                            } else{
                              estPars = c(d$estPars)
                              dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
                              VaREst = dforFull$OutSample$condVaR
                              EsEst = dforFull$OutSample$condES
                              VaR_in = d$condVaR
                              y_in = d$yLowFreq
                              VaR_out = tail(dforFull$OutSample$condVaR,10/horizon)
                              ES_out = tail(dforFull$OutSample$condES,10/horizon)
                              y_out = tail(dforFull$OutSample$yLowFreq,10/horizon)
                              stdResid = (y_in - VaR_in)/abs(VaR_in)
                              meanExceed = mean(stdResid[stdResid<0])
                              y_booted = matrix(NA,nrow = B,ncol = 10/horizon)
                              for(i in 1:B){
                                stdResid_booted = sample(stdResid,size = 10/horizon,replace = T)
                                y_B = stdResid_booted*abs(VaR_out) + VaR_out
                                if(length(which(y_B<VaR_out)) > 0){
                                  y_B[y_B < VaR_out] = VaR_out[y_B < VaR_out] + stdResid_booted[y_B < VaR_out]*(ES_out[y_B < VaR_out] - VaR_out[y_B < VaR_out]) / meanExceed
                                }
                                y_booted[i,] = y_B
                              }
                              u = vector()
                              for(i in 1:ncol(y_booted)){
                                u[i] = mean(y_booted[,i] < y_out[i])
                              }
                              ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), 
                                               VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),
                                               ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)),
                                               uEst = u)        
                              }
                          } else{
                            estPars = c(d$estPars)
                            dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
                            VaREst = dforFull$OutSample$condVaR
                            EsEst = dforFull$OutSample$condES
                            VaR_in = d$condVaR
                            y_in = d$yLowFreq
                            VaR_out = tail(dforFull$OutSample$condVaR,10/horizon)
                            ES_out = tail(dforFull$OutSample$condES,10/horizon)
                            y_out = tail(dforFull$OutSample$yLowFreq,10/horizon)
                            stdResid = (y_in - VaR_in)/abs(VaR_in)
                            meanExceed = mean(stdResid[stdResid<0])
                            y_booted = matrix(NA,nrow = B,ncol = 10/horizon)
                            for(i in 1:B){
                              stdResid_booted = sample(stdResid,size = 10/horizon,replace = T)
                              y_B = stdResid_booted*abs(VaR_out) + VaR_out
                              if(length(which(y_B<VaR_out)) > 0){
                                y_B[y_B < VaR_out] = VaR_out[y_B < VaR_out] + stdResid_booted[y_B < VaR_out]*(ES_out[y_B < VaR_out] - VaR_out[y_B < VaR_out]) / meanExceed
                              }
                              y_booted[i,] = y_B
                            }
                            u = vector()
                            for(i in 1:ncol(y_booted)){
                              u[i] = mean(y_booted[,i] < y_out[i])
                            }
                            ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), 
                                             VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),
                                             ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)),
                                             uEst = u)      
                            }
                        } else{
                        estPars = c(d$estPars)
                        dforFull =VarEs_jointAL_cav_for(object = d,y = idxData$ret[1:endpoint[w]], yDate = idxData$date[1:endpoint[w]])
                        VaREst = dforFull$OutSample$condVaR
                        EsEst = dforFull$OutSample$condES
                        VaR_in = d$condVaR
                        y_in = d$yLowFreq
                        VaR_out = tail(dforFull$OutSample$condVaR,10/horizon)
                        ES_out = tail(dforFull$OutSample$condES,10/horizon)
                        y_out = tail(dforFull$OutSample$yLowFreq,10/horizon)
                        stdResid = (y_in - VaR_in)/abs(VaR_in)
                        meanExceed = mean(stdResid[stdResid<0])
                        y_booted = matrix(NA,nrow = B,ncol = 10/horizon)
                        for(i in 1:B){
                          stdResid_booted = sample(stdResid,size = 10/horizon,replace = T)
                          y_B = stdResid_booted*abs(VaR_out) + VaR_out
                          if(length(which(y_B<VaR_out)) > 0){
                            y_B[y_B < VaR_out] = VaR_out[y_B < VaR_out] + stdResid_booted[y_B < VaR_out]*(ES_out[y_B < VaR_out] - VaR_out[y_B < VaR_out]) / meanExceed
                          }
                          y_booted[i,] = y_B
                        }
                        u = vector()
                        for(i in 1:ncol(y_booted)){
                          u[i] = mean(y_booted[,i] < y_out[i])
                        }
                        ans = data_frame(date = tail(dforFull$OutSample$yDate,10/horizon), 
                                              VaR = as.numeric(tail(dforFull$OutSample$condVaR,10/horizon)),
                                              ES = as.numeric(tail(dforFull$OutSample$condES,10/horizon)),
                                              uEst = u)
                        }
                      }
  return(forecast)
}

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
stopCluster(cluster)