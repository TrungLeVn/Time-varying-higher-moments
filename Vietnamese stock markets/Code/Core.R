# Load data
library(readxl)
library(dplyr)
library(SgtAcd)
library(foreach)
library(doParallel)
library(rugarch)
cluster <- makePSOCKcluster(5)
registerDoParallel(cluster)

RawDat <- read_excel("Cryto/5coins.xlsx", 
                     col_types = c("date", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric"))
#====================
# Transform to return
#====================
CoinNames <- colnames(RawDat)[2:6]
RawRet <- foreach(i = 1:5,.combine = "cbind")%do%{
  dat <- as.numeric(as.matrix(RawDat[,CoinNames[i]]))
  ret <- diff(log(dat))
  ans <- xts::xts(x = ret,order.by = RawDat$Date[2:nrow(RawDat)])
}

#===================
# Racd test 
#===================
garch_acd <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                     distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = c(1,0,1),shape1model = "quad",skewOrder = c(1,0,1),skewmodel = "quad"))
gjr_acd <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
                   distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = c(1,0,1),shape1model = "quad",skewOrder = c(1,0,1),skewmodel = "quad"))
gjr <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
               distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = NULL,skewOrder = NULL))
garch <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                 distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = NULL,skewOrder = NULL))

garch_est <- foreach(i = 1:5,.combine = "list",.multicombine = TRUE,.packages = c("SgtAcd","rugarch"))%dopar%{
  ans = acdfit(garch,RawRet[,i],solver = "msoptim")
}

gjr_est <- foreach(i = 1:5,.combine = "list",.multicombine = TRUE,.packages = c("SgtAcd","rugarch"))%dopar%{
  ans = acdfit(gjr,RawRet[,i],solver = "msoptim")
}

garch_acd_est <- foreach(i = 1:5,.combine = "list",.multicombine = TRUE,.packages = c("SgtAcd","rugarch"))%dopar%{
  ans = acdfit(garch,RawRet[,i],solver = "msoptim")
}

gjr_acd_est <- foreach(i = 1:5,.combine = "list",.multicombine = TRUE,.packages = c("SgtAcd","rugarch"))%dopar%{
  ans = acdfit(garch,RawRet[,i],solver = "msoptim")
}

InSample_est <- list(garch = garch_est, gjr = gjr_est, garch_acd = garch_acd_est, gjr_acd = gjr_acd_est)
saveRDS(InSample_est,"Cryto/InSample.rds")