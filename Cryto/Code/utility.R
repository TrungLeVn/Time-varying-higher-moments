#==============================
# Load Libraries
#==============================
libraries<- c("readxl", "dplyr", "xts", "lubridate",
              "zoo","SgtAcd","foreach","doParallel")
lapply(libraries, library, character.only = TRUE)

# Load data
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
garch_norm <- ugarchspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                         distribution.model = "norm")
gjr_norm <- ugarchspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
                       distribution.model = "norm")
garch_std <- ugarchspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "sGARCH"),
                        distribution.model = "std")
gjr_std <- ugarchspec(mean.model = list(armaOrder = c(1,0)),variance.model = list(model = "gjrGARCH"),
                      distribution.model = "std")

getFit <- function(spec,ret,windowLength = 1250,refit.every = 10){
   n = length(ret)
   lotEnd = seq(windowLength + refit.every,n,by = refit.every)
   ans = foreach(i = 1:length(lotEnd),.combine = "list",.multicombine = TRUE,.maxcombine = length(lotEnd),
                          .packages = "rugarch",.export = c("spec","ret","windowLength","lotEnd"))%dopar%{
                            fit = try(rugarch::ugarchfit(spec = spec, data = ret[1:lotEnd[i]], out.sample = refit.every),silent = TRUE)
                            if(inherits(fit,"try-error")){
                              ans = NULL
                            } else if(fit@fit$convergence == 1){
                              ans = NULL
                            } else{
                              ans = fit
                            }
                          }
  return(ans)
   }

getVaR <- function(fit, m.sim = 10000, VaR.alpha = c(0.001,0.005,0.01,0.025,0.05,0.1), refit.every = 10){
  if(is.null(fit)){
    forecast = NA
  } else{
    mx = fit@model$maxOrder
    fspec <- getspec(fit)
    fspec <- `setfixed<-`(fspec,fit@fit$coef)
    n.old = fit@model$modeldata$T
    dataflt = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
    flt <- ugarchfilter(spec = fspec,data = dataflt,n.old = n.old)
    sigmafilter 	= flt@filter$sigma
    resfilter 		= flt@filter$residuals
    uEst = matrix(NA,nrow = refit.every)
    VaRest = matrix(NA,nrow = refit.every,ncol = length(VaR.alpha))
    ESest = matrix(NA,nrow = refit.every,ncol = length(VaR.alpha))
    for(ii in seq(0,refit.every-1,by = 1)){
      presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
      prereturns = as.numeric(tail(dataflt[1:(n.old+ii)], mx))
      preresiduals = tail(resfilter[1:(n.old+ii)],mx)
      sim <- rugarch::ugarchsim(fit = fit,n.sim =1,m.sim = 10000,presigma = presig,
                                prereturns = prereturns,preresiduals = preresiduals,
                                startMethod = "sample")
      return = as.numeric(sim@simulation$seriesSim[1,])
      uEst[ii+1,] = mean(return < as.numeric(dataflt[n.old+ii+1]))
      rm(list = c("sim"))
      VaR <- quantile(return,VaR.alpha)
      ES <- VaR
      for(a in 1:length(VaR)){
        ES[a] <- mean(return[return < VaR[a]])
      }
      VaRest[ii+1,] = VaR
      ESest[ii+1,] = ES
    }
    dateEst = as.Date(tail(index(dataflt),refit.every))
    forecast = as_tibble(cbind(dateEst,uEst,VaRest,ESest))
    colnames(forecast) = c("date","uEst",paste("VaR",VaR.alpha,sep = ""),paste("ES",VaR.alpha,sep = ""))
    forecast$date <- as.Date.numeric(forecast$date,origin = "1970-01-01")
  }
  return(forecast)
}
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
CoinNames <- colnames(RawDat)[c(2,4:6)]
RawRet <- foreach(i = 1:length(CoinNames),.combine = "cbind")%do%{
  dat <- as.numeric(as.matrix(RawDat[,CoinNames[i]]))
  ret <- diff(log(dat))
  ret[ret == 0] = mean(ret)
  ret <- xts::xts(x = ret,order.by = RawDat$Date[2:nrow(RawDat)])
}
colnames(RawRet) = CoinNames
#=======================
# Back test functions
#=======================
Loss = function(y, VaR, ES, alpha, choice){
  # Check input
  if(length(y) != length(VaR)){
    stop('Length of y and VaR must be the same')
  }
  if(length(y) != length(ES)){
    stop('If ES is provided, it should be at the same length as VaR')
  }
  if(sum(is.infinite(VaR)) > 0){ 
    stop('VaR can not be Inf - Recheck the VaR estimates')
  }
  if(sum(is.infinite(ES)) > 0){ 
    stop('ES can not be Inf - Recheck the ES estimates')
  }
  if(alpha < 0|alpha > 1){
    stop('Quantile level must be between 0 and 1')
  }
  if(!is.element(choice,c(1,2,3))){
    stop('Choice of loss function need to be either 1, 2 or 3')
  }
  if(choice == 1){
    loss = (y - VaR)*(alpha - (y <= VaR))
    score = mean(loss)
  }else if(choice == 2){
    loss1 = ((y <= VaR) - alpha)*VaR - (y <= VaR)*y
    loss2 = (exp(ES)/(1 + exp(ES)))*(ES - VaR + (y <= VaR)*((VaR-y)/alpha))
    loss3 = log(2/(1 + exp(ES)))
    loss = loss1 + loss2 + loss3
    score = mean(loss1 + loss2 + loss3)
  }else{
    loss = -log((alpha - 1)/ES) - (((y - VaR)*(alpha - (y <= VaR)))/(alpha*ES)) + y/ES;
    score = mean(loss);
  }
  return(list(loss = loss, score = score))
}

KupiecTest = function(Hit,alpha){
  N = length(Hit)
  x = sum(Hit)
  rate = x/N
  test = -2 * log(((1 - alpha)^(N - x) * alpha^x)/((1 - rate)^(N - x) * rate^x))
  if (is.nan(test)) 
    test = -2 * ((N - x) * log(1 - alpha) + x * log(alpha) - 
                   (N - x) * log(1 - rate) - x * log(rate))
  pvalue = 1 - stats::pchisq(test, df = 1)
  LRpof = c(test, pvalue)
  names(LRpof) = c("Test", "Pvalue")
  return(LRpof)  
}

DQtest_VaR = function(y,VaR,alpha,lags = 4){
  cT = length(y)
  vHit = numeric(cT)
  vHit[y <= VaR] = 1 - alpha
  vHit[y > VaR] = -alpha
  vConstant = rep(1, (cT - lags))
  vHIT = vHit[(lags + 1):cT]
  vVaRforecast = VaR[(lags + 1):cT]
  mZ = matrix(0, cT - lags, lags)
  vY2_lag = y[lags:(cT - 1)]^2
  for (st in 1:lags) {
    mZ[, st] = vHit[st:(cT - (lags + 1L - st))]
  }
  mX = cbind(vConstant, vVaRforecast, mZ)
  dDQstatOut = (t(vHIT) %*% mX %*% MASS::ginv(t(mX) %*% mX) %*% 
                  t(mX) %*% (vHIT))/(alpha * (1 - alpha))
  dDQpvalueOut = 1 - stats::pchisq(dDQstatOut, ncol(mX))
  return(dDQpvalueOut)
}

var_backtest = function(realized,VaR,alpha,lags = 5){
  Hit = numeric(length(realized))
  Hit[which(realized <= VaR)] = 1L
  N = length(Hit)
  x = sum(Hit)
  rate = x/N
  AE = rate/alpha
  Kupiec = KupiecTest(Hit,alpha)
  DQVaR = DQtest_VaR(realized,VaR,alpha,lags)
  return(list(AE = AE,Kupiec = Kupiec,DQVaR = DQVaR))
}

es_McNeil = function(realized,VaR,ES,B=1000){
  # Extract standardized exceedances
  stdExceed = ((realized - ES)/VaR)[realized <= VaR]
  # Build the modified stdExceed that have zero mean to build the bootstrap population
  # Refer to Efron and Tibshirani (1993)
  newstdExceed = stdExceed - mean(stdExceed) # The newstdExceed has the mean = 1
  f = function(x)  mean(x) / (stats::sd(x) / sqrt(length(x)))
  t0 = f(stdExceed)
  t = c()
  for (i in 1:1000){ 
    newsample <- sample(newstdExceed, length(newstdExceed), replace=T)
    t <- c(t,f(newsample))
  }
  return(mean(abs(t) > abs(t0)))
}

es_Patton = function(y,VaR,ES,alpha,lags = 4){
  cT = length(y)
  eHit = numeric(cT)
  eHit[y <= VaR] = (1/alpha) * (y/ES)[y <= VaR] - 1
  eHit[y > VaR] = -1
  eConstant = rep(1, (cT - lags))
  eHIT = eHit[(lags + 1):cT]
  eESforecast = ES[(lags + 1):cT]
  mZ = matrix(0, cT - lags, lags)
  for (st in 1:lags) {
    mZ[, st] = eHit[st:(cT - (lags + 1L - st))]
  }
  reg = lm(eHIT ~ 0 + eConstant + eESforecast + mZ)
  Hnull = names(coefficients(reg))
  for(i in 1:length(Hnull)) Hnull[i] = paste(Hnull[i],"=0",sep = "")
  fTest = car::linearHypothesis(reg,Hnull)
  dDQpvalueOut = fTest$`Pr(>F)`[2]
  return(dDQpvalueOut)
}

# Define the ES backtest of Acerbi (2014) with regards to the Z2 test statistic
es_Acerbi = function(y,VaR,ES,alpha,B = 10000){
  z2 = -((1/(length(y)*alpha)) * sum((y*(y <= VaR))/ES)) + 1
  stdResid <- (y - VaR)/abs(VaR)
  meanExceed = mean(stdResid[stdResid<0])
  z2_booted = vector()
  for(i in 1:B){
    stdResid_booted <- sample(stdResid,replace = T)
    y_booted = stdResid_booted*abs(VaR) + VaR
    exceeded_lot = which(y_booted < VaR)
    y_booted[exceeded_lot] = VaR[exceeded_lot] + stdResid_booted[exceeded_lot]*(ES[exceeded_lot] - VaR[exceeded_lot]) / meanExceed
    z2_booted[i] = -((1/(length(y_booted)*alpha)) * sum((y_booted*(y_booted < VaR))/ES)) + 1
  }
  return(c(z2,mean(z2_booted > z2)))
}

es_backtest= function(realized,VaR,ES,alpha,lags = 4, B = 10000){
  UnC_backtest = es_McNeil(realized,VaR,ES,B)
  #Patton_backtest = es_Patton(realized,VaR,ES,alpha,lags)
  Acerbi_backtest = es_Acerbi(realized,VaR,ES,alpha,B)
  return(list(UnC = UnC_backtest, Acerbi = Acerbi_backtest))#, Patton = Patton_backtest))
}

#------------------------
# ES-test by Du and Escanciano (MS, 2017)
#-----------------------
library(fGarch)
EStest_DuEs <- function(utj,nout,jalp,ddls = 15){
  sqdd <- (nout-1):(nout-ddls)
  indalp<- utj<=jalp           
  utalp<-rep(0, nout)
  utalp[indalp]<-utj[indalp]-jalp
  mutalp<-utalp+jalp^2/2
  invar <- 1/var(mutalp)
  UES<-sum(mutalp)*sqrt(invar)/sqrt(nout)  #Basic Unconditional Backtest for ES
  pUES<-2*pnorm(-abs(UES))                 #p-value of Basic Unconditional Backtest for ES
  x<- c(0, mutalp[1:nout-1])
  mutldd <- Triang(toeplitz(x))[1:ddls,]     
  rhojalp<-invar*mutldd%*%mutalp/sqrt(sqdd)
  rhoj2alp <-rhojalp^2
  CES<-cumsum(rhoj2alp)                    #Basic Conditional Backtest for ES at lag 1 to 15
  pCES<- 1-pchisq(CES,(1:ddls))  
  return(list(UES = pUES,CES = pCES))
}