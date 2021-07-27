#==============================
# Load Libraries
#==============================
libraries<- c("readxl", "dplyr", "xts", "lubridate","fExtremes",
              "zoo","SgtAcd","foreach","doParallel")
#install packages if not already installed 
packages_to_install <- libraries[!(libraries %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) install.packages(packages_to_install)
#load libraries
lapply(libraries, library, character.only = TRUE)

#===================
# Model Specification
#===================
garch_acd <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                     distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = c(1,0,1),shape1model = "quad",skewOrder = c(1,0,1),skewmodel = "quad"))
gjr_acd <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
                   distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = c(1,0,1),shape1model = "quad",skewOrder = c(1,0,1),skewmodel = "quad"))
gjr <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
               distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = NULL,skewOrder = NULL))
garch <- acdspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                 distribution.model = list(model = "sged",shape2Order = NULL,shape1Order = NULL,skewOrder = NULL))

#====================
# Load Data

#====================

RawDat <- read_excel("Cryto/Dat/5coins.xlsx", 
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
colnames(RawRet) = CoinNames