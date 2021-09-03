rm(list = ls())
source("Vietnam/Code/utility.R")
#===============#
# Gather out-of-sample estimation
#===============#
vni.garchNorm.Est = readRDS("Vietnam/Estimate/vni-outgarch_norm.rds")
vni.garchStd.Est = readRDS("Vietnam/Estimate/vni-outgarch_std.rds")
vni.garchSged.Est = readRDS("Vietnam/Estimate/vni-out-garch.rds") # Had to resume sim for this object
vni.garchSged.Est <- acdresumeSim(vni.garchSged.Est)
saveRDS(vni.garchSged.Est,"Vietnam/Estimate/vni-out-garch.rds")
vni.garchAcd.Est = readRDS("Vietnam/Estimate/vni-out-garchACD.rds") # Had to resume sim for this object
vni.garchAcd.Est <- acdresumeSim(vni.garchAcd.Est)
saveRDS(vni.garchAcd.Est,"Vietnam/Estimate/vni-out-garchACD.rds")
vni.gjrNorm.Est = readRDS("Vietnam/Estimate/vni-outgjr_norm.rds")
vni.gjrStd.Est = readRDS("Vietnam/Estimate/vni-outgjr_std.rds")
vni.gjrSged.Est = readRDS("Vietnam/Estimate/vni-out-gjr.rds")
vni.gjrAcd.Est = readRDS("Vietnam/Estimate/vni-out-gjrACD.rds")

MixedFreqDat <- function(DataY,DataYdate,period,ovlap = FALSE,simpleRet = FALSE){
  nobs <- length(DataY)
  nobsShort = nobs-period+1
  DataYlowfreq = matrix(NA,ncol = 1, nrow = nobsShort)
  DataYDateLow = matrix(NA,ncol = 1, nrow = nobsShort)
  for(t in 1:(nobs - period + 1)){
    DataYlowfreq[t,1] = sum(DataY[t:(t+period-1)]);
    DataYDateLow[t,1] = DataYdate[t];
  }
  if(!ovlap){
    DataYlowfreq = DataYlowfreq[seq(1,nobsShort,period),1]
    DataYDateLow = DataYDateLow[seq(1,nobsShort,period),1]
  }
  output = tibble(date = as.Date.numeric(DataYDateLow,origin = "1970-01-01"), ret = DataYlowfreq)
  return(output)
}

horizon = c(1,5,10)

AllFor <- foreach(i = 1:3,.combine = "list",.multicombine = TRUE)%do%{
  # Garch-Norm 
  vni.garchNorm <- foreach(w = 1:length(vni.garchNorm.Est),.combine = "rbind")%do%{
    ans = vni.garchNorm.Est[[w]][[i]]
  }
  dateEst <- row.names(vni.garchNorm)
  vni.garchNorm <- as_tibble(vni.garchNorm)
  colnames(vni.garchNorm) = c("mu",paste(c("uEst","var001","var005","es001","es005"),"garchNorm",sep = "_"))
  vni.garchNorm <- vni.garchNorm %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu)
  
  # Garch-std
  vni.garchStd <- foreach(w = 1:length(vni.garchStd.Est),.combine = "rbind")%do%{
    ans = vni.garchStd.Est[[w]][[i]]
  }
  
  dateEst <- row.names(vni.garchStd)
  vni.garchStd <- as_tibble(vni.garchStd)
  colnames(vni.garchStd) = c("mu",paste(c("uEst","var001","var005","es001","es005"),"garchStd",sep = "_"))
  vni.garchStd <- vni.garchStd %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu)
  
  # Garch-Sged
  vni.garchSged <- foreach(w = 1:length(vni.garchSged.Est@forecast),.combine = "rbind")%do%{
    ans = vni.garchSged.Est@forecast[[w]][[i]]
  }
  
  dateEst <- row.names(vni.garchSged)
  vni.garchSged <- as_tibble(vni.garchSged)
  colnames(vni.garchSged) = c("mu","sigma","tm","skewness","kurtosis",paste(c("uEst","var001","var005","es001","es005"),"garchSged",sep = "_"))
  vni.garchSged <- vni.garchSged %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu,-sigma,-tm,-skewness,-kurtosis)
  
  # Garch-Acd
  vni.garchAcd <- foreach(w = 1:length(vni.garchAcd.Est@forecast),.combine = "rbind")%do%{
    ans = vni.garchAcd.Est@forecast[[w]][[i]]
  }
  
  dateEst <- row.names(vni.garchAcd)
  vni.garchAcd <- as_tibble(vni.garchAcd)
  colnames(vni.garchAcd) = c("mu","sigma","tm","skewness","kurtosis",paste(c("uEst","var001","var005","es001","es005"),"garchAcd",sep = "_"))
  vni.garchAcd <- vni.garchAcd %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu,-sigma,-tm,-skewness,-kurtosis)
  
  # GJR-Norm 
  vni.gjrNorm <- foreach(w = 1:length(vni.gjrNorm.Est),.combine = "rbind")%do%{
    ans = vni.gjrNorm.Est[[w]][[i]]
  }
  dateEst <- row.names(vni.gjrNorm)
  vni.gjrNorm <- as_tibble(vni.gjrNorm)
  colnames(vni.gjrNorm) = c("mu",paste(c("uEst","var001","var005","es001","es005"),"gjrNorm",sep = "_"))
  vni.gjrNorm <- vni.gjrNorm %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu)
  
  # Gjr-std
  vni.gjrStd <- foreach(w = 1:length(vni.gjrStd.Est),.combine = "rbind")%do%{
    ans = vni.gjrStd.Est[[w]][[i]]
  }
  
  dateEst <- row.names(vni.gjrStd)
  vni.gjrStd <- as_tibble(vni.gjrStd)
  colnames(vni.gjrStd) = c("mu",paste(c("uEst","var001","var005","es001","es005"),"gjrStd",sep = "_"))
  vni.gjrStd <- vni.gjrStd %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu)
  
  # Gjr-Sged
  vni.gjrSged <- foreach(w = 1:length(vni.gjrSged.Est@forecast),.combine = "rbind")%do%{
    ans = vni.gjrSged.Est@forecast[[w]][[i]]
  }
  
  dateEst <- row.names(vni.gjrSged)
  vni.gjrSged <- as_tibble(vni.gjrSged)
  colnames(vni.gjrSged) = c("mu","sigma","tm","skewness","kurtosis",paste(c("uEst","var001","var005","es001","es005"),"gjrSged",sep = "_"))
  vni.gjrSged <- vni.gjrSged %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu,-sigma,-tm,-skewness,-kurtosis)
  
  # Gjr-Acd
  vni.gjrAcd <- foreach(w = 1:length(vni.gjrAcd.Est@forecast),.combine = "rbind")%do%{
    ans = vni.gjrAcd.Est@forecast[[w]][[i]]
  }
  
  dateEst <- row.names(vni.gjrAcd)
  vni.gjrAcd <- as_tibble(vni.gjrAcd)
  colnames(vni.gjrAcd) = c("mu","sigma","tm","skewness","kurtosis",paste(c("uEst","var001","var005","es001","es005"),"gjrAcd",sep = "_"))
  vni.gjrAcd <- vni.gjrAcd %>% mutate(date = as.Date.character(dateEst)) %>% select(-mu,-sigma,-tm,-skewness,-kurtosis)
  
  # Ret
  ret = RawRet$vni$ret
  n = length(ret)
  
  retFreq <- MixedFreqDat(DataY = ret[1251:n], DataYdate = index(ret)[1251:n], period = horizon[i], ovlap = FALSE, simpleRet = FALSE)
    
  ans <- left_join(vni.garchNorm,vni.garchStd,by = "date") %>% left_join(.,vni.garchSged,by = "date") %>% left_join(.,vni.garchAcd,by = "date") %>%
    left_join(.,vni.gjrNorm,by = "date") %>% left_join(.,vni.gjrStd,by = "date") %>% left_join(.,vni.gjrSged,by = "date") %>% 
    left_join(.,vni.gjrAcd,by = "date") %>% left_join(.,retFreq,by = "date")
}
names(AllFor) = c("1-day","5-day","10-day")