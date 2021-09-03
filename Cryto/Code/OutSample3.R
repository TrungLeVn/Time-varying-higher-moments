rm(list = ls())
# Load libraries
source("Cryto/Code/utility.R")
noCores = detectCores()
print(noCores)
cluster <- makePSOCKcluster(noCores-1,outfile = "")
registerDoParallel(cluster)

#===================
# Out sample estimation
# Compute out-of-sample forecasts using rolling windows of 1250 (5 years) 
# Parameters are re-estimated for each 10 days 
# In each ten days, get 10 1-day (re-iterative), 2 5-days and 1 10-day VaR and Expected Shortfall
#===================
#out_garch_Sged <- acdrollsim(spec = garch,data = RawRet[,3],horizon = 1,m.sim = 10000,n.start = 1250,
#                          refit.every = 10,refit.window = "recursive",window.size = 1250,fixARMA = FALSE,
#                          fixGARCH = FALSE,cluster = cluster,solver.control = list(trace = TRUE),
#                          VaR.alpha = c(0.001,0.005,0.01,0.025,0.05,0.1))
#saveRDS(out_garch_Sged,"Cryto/Estimate/coin3-outgarchSged.rds")

#rm(list = "out_garch_Sged")
#out_garch_Acd<- acdrollsim(spec = garch_acd,data = RawRet[,3],horizon = 1,m.sim = 10000,n.start = 1250,
#                           refit.every = 10,refit.window = "recursive",window.size = 1250,fixARMA = FALSE,
#                           fixGARCH = FALSE,cluster = cluster,solver.control = list(trace = TRUE),
#                           VaR.alpha = c(0.001,0.005,0.01,0.025,0.05,0.1))
#saveRDS(out_garch_Acd,"Cryto/Estimate/coin3-outgarchAcd.rds")
#rm(list = "out_garch_Acd")


#out_gjr_Sged <- acdrollsim(spec = gjr,data = RawRet[,3],horizon = 1,m.sim = 10000,n.start = 1250,
#                             refit.every = 10,refit.window = "recursive",window.size = 1250,fixARMA = FALSE,
#                             fixGARCH = FALSE,cluster = cluster,solver.control = list(trace = TRUE),
#                           VaR.alpha = c(0.001,0.005,0.01,0.025,0.05,0.1))
#saveRDS(out_gjr_Sged,"Cryto/Estimate/coin3-outgjrSged.rds")

#rm(list = "out_gjr_Sged")

#out_gjr_Acd<- acdrollsim(spec = gjr_acd,data = RawRet[,3],horizon = 1,m.sim = 10000,n.start = 1250,
#                           refit.every = 10,refit.window = "recursive",window.size = 1250,fixARMA = FALSE,
#                           fixGARCH = FALSE,cluster = cluster,solver.control = list(trace = TRUE),
#                         VaR.alpha = c(0.001,0.005,0.01,0.025,0.05,0.1))
#saveRDS(out_gjr_Acd,"Cryto/Estimate/coin3-outgjrAcd.rds")
#rm(list = "out_gjr_Acd")
#

############################
out_fit = getFit(spec = garch_norm,ret = RawRet[,3])
out_garch_norm <- foreach(i = 1:length(out_fit),.combine = "rbind",.packages = c("zoo","dplyr"))%dopar%{
  ans = getVaR(fit = out_fit[[i]])
}
saveRDS(out_garch_norm,"Cryto/Estimate/coin3-outgarchNorm.rds")
rm(list = "out_fit")
out_fit = getFit(spec = garch_std,ret = RawRet[,3])
out_garch_std <- foreach(i = 1:length(out_fit),.combine = "rbind",
                         .packages = c("zoo","dplyr"))%dopar%{
                           ans = getVaR(fit = out_fit[[i]])
                         }
saveRDS(out_garch_std,"Cryto/Estimate/coin3-outgarchStd.rds")
rm(list = "out_fit")
out_fit = getFit(spec = gjr_norm,ret = RawRet[,3])
out_gjr_norm <- foreach(i = 1:length(out_fit),.combine = "rbind")%dopar%{
  ans = getVaR(fit = out_fit[[i]])
}
saveRDS(out_gjr_norm,"Cryto/Estimate/coin3-outgjrNorm.rds")
rm(list = "out_fit")

out_fit = getFit(spec = gjr_std,ret = RawRet[,3])
out_gjr_std <- foreach(i = 1:length(out_fit),.combine = "rbind")%dopar%{
  ans = getVaR(fit = out_fit[[i]])
}
saveRDS(out_gjr_std,"Cryto/Estimate/coin3-outgjrStd.rds")
stopCluster(cluster)
