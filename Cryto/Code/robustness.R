#---------------------
# Robustness 
#---------------------
# Compare with quantile regression based
as <- readRDS("D:/OneDrive - hvnh.edu.vn/Research/Articles/2021/Higher-moments/Cryto/Estimate/as.rds")
sav <- readRDS("D:/OneDrive - hvnh.edu.vn/Research/Articles/2021/Higher-moments/Cryto/Estimate/sav.rds")

thetaMain = c(0.01,0.025,0.05)
allEst_quantile <- foreach(i = 1:4,.combine = "list",.multicombine = TRUE)%do%{
  coinEst <- foreach(ii = 1:3,.combine = "cbind")%do%{
    quantEst_as = as[[i]][[ii]]$forecast[,2:3] 
    quantEst_sav = sav[[i]][[ii]]$forecast[,2:3] 
    names(quantEst_as) = paste(c(paste("VaR",thetaMain[ii],sep = ""),paste("ES",thetaMain[ii],sep = "")),"as",sep = "_")
    names(quantEst_sav) = paste(c(paste("VaR",thetaMain[ii],sep = ""),paste("ES",thetaMain[ii],sep = "")),"sav",sep = "_")
    quantEst = as_tibble(cbind(quantEst_as,quantEst_sav))
  }
  coinEst$date = as[[i]][[1]]$forecast$date
  coinEst
}
names(allEst_quantile) = CoinNames

MegaEst = list()
for(i in 1:4){
  MegaEst[[i]] = left_join(allEst[[i]],allEst_quantile[[i]],by = "date")
}

modelName = c(modelName,"as","sav")
