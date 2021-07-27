#==============================
# Load Libraries
#==============================
libraries<- c("readr", "dplyr", "xts", "lubridate","fExtremes",
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

vni <- read_csv("Vietnam/Data/vni.csv", col_types = cols(Vol. = col_skip(), 
                                                         `Change %` = col_skip()))
vni30 <- read_csv("Vietnam/Data/vni30.csv", col_types = cols(Vol. = col_skip(), 
                                                             `Change %` = col_skip()))
vniSmall <- read_csv("Vietnam/Data/vni-small.csv", col_types = cols(Vol. = col_skip(), 
                                                                    `Change %` = col_skip()))
vniMid <- read_csv("Vietnam/Data/vni-mid.csv", col_types = cols(Vol. = col_skip(), 
                                                                `Change %` = col_skip()))
hnx <- read_csv("Vietnam/Data/hnx.csv", col_types = cols(Vol. = col_skip(), 
                                                         `Change %` = col_skip()))
hnx30 <- read_csv("Vietnam/Data/hnx30.csv", col_types = cols(Vol. = col_skip(), 
                                                             `Change %` = col_skip()))
rawDat <- list(vni = vni, vni30 = vni30, vniSmall = vniSmall,vniMid = vniMid, hnx = hnx, hnx30 = hnx30)
#====================
# Transform to return series
#====================

RawRet <- foreach(i = 1:length(rawDat),.combine = "list",.multicombine = TRUE)%do%{
  dat = rawDat[[i]]
  dat$Date <- as.Date.character(dat$Date,format = "%b %d,%Y") # Transform to Date series
  dat <- dat %>% arrange(Date) %>% mutate(ret = c(NA,diff(log(dat$Price))))
  dat = dat[2:nrow(dat),]
  ans <- xts::xts(x = dat[,2:ncol(dat)],order.by = dat$Date)
}
names(RawRet) = names(rawDat)

