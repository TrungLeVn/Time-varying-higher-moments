setwd("~/2021/Higher-moments")
rm(list = ls())
# Load data
source(file = "Cryto/Code/utility.R")
noCores = detectCores()
cluster <- makePSOCKcluster(noCores-1,outfile = "")
registerDoParallel(cluster)

#===================
# Out sample estimation
# Compute out-of-sample forecasts using rolling windows of 1250 (5 years) 
# Parameters are re-estimated for each 10 days 
# In each ten days, get 10 1-day (re-iterative), 2 5-days and 1 10-day VaR and Expected Shortfall
#===================

out_garch <- acdrollsim(spec = garch,data = RawRet[,3],horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                        window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.025,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_garch,"Cryto/Estimate/Garch-out-coint3.rds")

out_gjr <- acdrollsim(spec = gjr,data = RawRet[,3],horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                      window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.025,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_gjr,"Cryto/Estimate/Gjr-out-coint3.rds")

out_garchACD <- acdrollsim(spec = garch_acd,data = RawRet[,3],horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                           window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.025,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_garchACD,"Cryto/Estimate/GarchACD-out-coint3.rds")

out_gjrACD <- acdrollsim(spec = gjr_acd,data = RawRet[,3],horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                         window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.025,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_gjrACD,"Cryto/Estimate/GjrACD-out-coint3.rds")
stopCluster(cluster)
