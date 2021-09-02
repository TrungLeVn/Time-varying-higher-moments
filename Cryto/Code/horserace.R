rm(list = ls())
source("Cryto/Code/utility.R")
#===============#
# Gather out-of-sample estimation
#===============#
modelName = c("garchNorm","garchStd","garchSged","garchAcd","gjrNorm","gjrStd","gjrSged","gjrAcd")
theta_test = c(0.001,0.005,0.01,0.025,0.05,0.1)

allEst <- foreach(i = 1:4,.combine = "list",.multicombine = TRUE)%do%{
  coinEst <- foreach(ii = 1:length(modelName),.combine = "list",.multicombine = TRUE)%do%{
    ans = readRDS(paste("Cryto/Estimate/coin",i,"-out",modelName[ii],".rds",sep = ""))
    if(ii %in% c(1,2,5,6)){
      ans <- as_tibble(ans)
      colnames(ans) = c("date",paste(c("uEst",paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = "")),modelName[ii],sep = "_"))
    } else{
      if(i == 4 && ii == 8){
        extract_est <- foreach(iii = 1:length(ans@forecast),.combine = "rbind")%do%{
          extract = ans@forecast[[iii]]$forecast$`1`[,6:18]
        }
      }else{
        extract_est <- foreach(iii = 1:length(ans@forecast),.combine = "rbind")%do%{
          extract = ans@forecast[[iii]]$`1`[,6:18]
        }
      }
      dateEst <- as.Date.character(row.names(extract_est))
      ans = as_tibble(extract_est)
      ans$date = dateEst
      colnames(ans) = c(paste(c("uEst",paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = "")),modelName[ii],sep = "_"),"date")
    }
    ans
  }
  names(coinEst) = modelName
  ret = window(RawRet[,i],start = coinEst$garchNorm$date[1],end = tail(coinEst$garchNorm$date,1))
  ret <- as_tibble(ret)
  ret$date = coinEst$garchNorm$date
  colnames(ret) = c("ret","date")
  
  ans <- left_join(coinEst$garchNorm,coinEst$garchStd,by = "date") %>% left_join(.,coinEst$garchSged,by = "date") %>% left_join(.,coinEst$garchAcd,by = "date") %>%
    left_join(.,coinEst$gjrNorm,by = "date") %>% left_join(.,coinEst$gjrStd,by = "date") %>% left_join(.,coinEst$gjrSged,by = "date") %>% 
    left_join(.,coinEst$gjrAcd,by = "date") %>% left_join(.,ret,by = "date")
}
names(allEst) = CoinNames
saveRDS(allEst,"Cryto/Estimate/allGarchEst")
#---------------------
# Backtest results
#---------------------
GetNumberFail <- function(testResult){
  ans = vector()
  for(i in 1:nrow(testResult)){
    ans[i] = length(which(testResult[i,] < 0.05))
  }
  return(ans)
}

# Absolute Test
Absolute_Test <- foreach(i = 1:length(theta_test),.combine = "list",.multicombine = T)%do%{
  absolute <- foreach(ii = 1:4,.combine = "list",.multicombine = TRUE)%do%{
    realized = allEst[[ii]][,"ret"]$ret
    testResult <- foreach(m = 1:length(modelName),.combine = "cbind")%do%{
      VaREs = as.matrix(allEst[[ii]][,paste(c(paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = ""),"uEst"),modelName[m],sep = "_")])
      # in case has NA, intepolate
      if(length(which(is.na(VaREs[,1]))) > 0){
        for(v in 1:ncol(VaREs)){
          VaREs[,v] <- as.numeric(forecast::na.interp(as.ts(VaREs[,v])))
        }
      }
      var_results = var_backtest(realized = realized,VaR = VaREs[,i],alpha = theta_test[i],lags = 4)
      var_results = c(var_results$AE,var_results$Kupiec[2],var_results$DQVaR[1,1])
      es_Mc = es_backtest(realized = realized,VaR = VaREs[,i],ES = VaREs[,length(theta_test)+i],alpha =theta_test[i] ,B = 1000,lags = 5)
      es_Mc = c(es_Mc$UnC,es_Mc$Acerbi[2])
      es_Du = EStest_DuEs(utj = VaREs[,ncol(VaREs)],nout = nrow(VaREs),jalp = theta_test[i],ddls = 5)
      es_Du = c(es_Du$UES,es_Du$CES[4])
      ans = c(var_results,es_Mc,es_Du)
    }
  }
}
names(Absolute_Test) = as.character(theta_test)
# Loss function
Loss_rank <- foreach(i = 1:length(theta_test),.combine = "list",.multicombine = T)%do%{
  loss <- foreach(ii = 1:4,.combine = "cbind")%do%{
    realized = allEst[[ii]][,"ret"]$ret
    testResult <- foreach(m = 1:length(modelName),.combine = "rbind")%do%{
      VaREs = as.matrix(allEst[[ii]][,paste(c(paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = ""),"uEst"),modelName[m],sep = "_")])
      # in case has NA, intepolate
      if(length(which(is.na(VaREs[,1]))) > 0){
        for(v in 1:ncol(VaREs)){
          VaREs[,v] <- as.numeric(forecast::na.interp(as.ts(VaREs[,v])))
        }
      }
      loss1 = Loss(y = realized,VaR = VaREs[,i],ES = VaREs[,length(theta_test)+i],alpha = theta_test[i],choice = 1)$score
      loss2 = Loss(y = realized,VaR = VaREs[,i],ES = VaREs[,length(theta_test)+i],alpha = theta_test[i],choice = 2)$score
      ans = c(loss1,loss2)
    }
  }
}
names(Loss_rank) = as.character(theta_test)

# Mcc test 
# Seem not working ok for only 4 index
#----------------------------
# Load required packages
#----------------------------
library(dplyr)
library(doParallel)
library(rugarch)
library(foreach)
cl<-makePSOCKcluster(4)
registerDoParallel(cl)

MCS <- foreach(i = 1:length(theta_test),.combine = "list",.multicombine = T)%do%{
  L1 <- foreach(ii = 1:4,.combine = "cbind",.packages = c("rugarch","foreach","forecast"),
                .export = c("allEst","modelName"))%dopar%{
                  print(ii)
                  realized = allEst[[ii]][,"ret"]$ret
                  AE <- foreach(m = 1:length(modelName),.combine = "cbind")%do%{
                    VaREs = as.matrix(allEst[[ii]][,paste(c(paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = ""),"uEst"),modelName[m],sep = "_")])
                    # in case has NA, intepolate
                    if(length(which(is.na(VaREs[,1]))) > 0){
                      for(v in 1:ncol(VaREs)){
                        VaREs[,v] <- as.numeric(forecast::na.interp(as.ts(VaREs[,v])))
                      }
                    }
                    loss = Loss(y = realized,VaR = VaREs[,i],ES = VaREs[,length(theta_test)+i],alpha = theta_test[i],choice = 1)$loss
                  }
                  ans <- matrix(0,ncol = 1, nrow = length(modelName))
                  mcs <- rugarch::mcsTest(losses = AE,alpha = 0.05,nboot = 10000,nblock = 4,boot="block")
                  ans[mcs$includedR] = 1
                  ans
                }
  
  L2 <- foreach(ii = 1:4,.combine = "cbind",.packages = c("rugarch","foreach"),
                .export = c("allEst","modelName"))%dopar%{
                  realized = allEst[[ii]][,"ret"]$ret
                  AE <- foreach(m = 1:length(modelName),.combine = "cbind")%do%{
                    VaREs = as.matrix(allEst[[ii]][,paste(c(paste("VaR",theta_test,sep = ""),paste("ES",theta_test,sep = ""),"uEst"),modelName[m],sep = "_")])
                    # in case has NA, intepolate
                    if(length(which(is.na(VaREs[,1]))) > 0){
                      for(v in 1:ncol(VaREs)){
                        VaREs[,v] <- as.numeric(forecast::na.interp(as.ts(VaREs[,v])))
                      }
                    }
                    loss = Loss(y = realized,VaR = VaREs[,i],ES = VaREs[,length(theta_test)+i],alpha = theta_test[i],choice = 2)$loss
                  }
                  ans <- matrix(0,ncol = 1, nrow = 8)
                  mcs <- rugarch::mcsTest(losses = AE,alpha = 0.05,nboot = 10000,nblock = 5,boot="block")
                  ans[mcs$includedR] = 1
                  ans
                }
  ans = list(L1 = L1, L2 = L2)
}