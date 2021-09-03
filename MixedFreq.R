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
