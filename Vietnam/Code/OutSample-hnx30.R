setwd("~/2021/Higher-moments")
rm(list = ls())
# Load data
source(file = "Vietnam/Code/utility.R")
noCores = detectCores()
cluster <- makePSOCKcluster(noCores-1,outfile = "")
registerDoParallel(cluster)

#===================
# Out sample estimation
# Compute out-of-sample forecasts using rolling windows of 1250 (5 years) 
# Parameters are re-estimated for each 10 days 
# In each ten days, get 10 1-day (re-iterative), 2 5-days and 1 10-day VaR and Expected Shortfall
#===================

out_garch <- acdrollsim(spec = garch,data = RawRet$hnx30$ret,horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                        window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_garch,"Vietnam/Estimate/hnx30-out-garch.rds")

out_gjr <- acdrollsim(spec = gjr,data = RawRet$hnx30$ret,horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                      window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_gjr,"Vietnam/Estimate/hnx30-out-gjr.rds")

out_garchACD <- acdrollsim(spec = garch_acd,data = RawRet$hnx30$ret,horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                           window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_garchACD,"Vietnam/Estimate/hnx30-out-garchACD.rds")

out_gjrACD <- acdrollsim(spec = gjr_acd,data = RawRet$hnx30$ret,horizon = c(1,5,10),m.sim = 10000,n.start = 1250,refit.every = 10,burn = 0,
                         window.size = 1250,calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05),cluster = cluster,refit.window = "moving")
saveRDS(out_gjrACD,"Vietnam/Estimate/hnx30-out-gjrACD.rds")
