# SSAGA v1.1 portfolio test 1 - constructing the efficient frontier 

################### FITNESS FUNCTION ################### 

fitness_function_Markowitz_EF = function(weights, returns, risk, lambda) {
  return(lambda * t(weights) %*% returns - (1 - lambda) * t(weights) %*% risk %*% weights)
}

################### DATA IMPORT ################### 

NASDAQdata = read.csv('StockDataV2.csv')

NASDAQdata = NASDAQdata[,-95]
nRows = as.numeric(dim(NASDAQdata)[1])
nCols = as.numeric(dim(NASDAQdata)[2])

NASDAQ_No_Date = as.matrix(NASDAQdata[1:nRows, 2:nCols])
NASDAQ_No_Date[is.na(NASDAQ_No_Date)] = 0

#Generate returns data
ReturnData = matrix(1, nrow = (nRows-1), ncol = (nCols-1))

for(i in 2:(nRows-1)) {
  ReturnData[(i-1),] = try(log(NASDAQ_No_Date[i,]) - log(NASDAQ_No_Date[(i-1),]))
}

colnames(ReturnData) = colnames(NASDAQ_No_Date)
ReturnData[is.na(ReturnData)] = 0

################### RISK AND RETURNS MATRICES ################### 

HistoricalReturns = function(Returns, TimeIndex) {
  alpha = numeric(nCols-1)
  for (i in 1:(nCols-1)){
    alpha[i] = mean(Returns[1:TimeIndex,i])
  }
  alpha[is.infinite(alpha)] = 0
  alpha[is.na(alpha)] = 0
  return(alpha)
}

CovarianceMatrix = function(Returns, TimeIndex) {
  omega = matrix(rep(0, (nCols-1)^2), nCols-1, nCols-1) 
  for (j in 1:(nCols-1)){
    for (k in 1:(nCols-1)){
      omega[j,k] = cov(Returns[1:TimeIndex,j], Returns[1:TimeIndex,k])    
    }
  }
  omega[is.infinite(omega)] = 0
  omega[is.na(omega)] = 0
  return(omega)
}

#Compute returns from time series matrix
alpha = HistoricalReturns(ReturnData, as.numeric(dim(ReturnData)[1]))
omega = CovarianceMatrix(ReturnData, as.numeric(dim(ReturnData)[1]))

################### THE EFFICIENT FRONTIER ################### 
library(compiler)
enableJIT(2)
fitness_function_M_EF = cmpfun(fitness_function_Markowitz_EF) 
SSAGAcomp = cmpfun(SSAGA12_PE)

EfficientFrontier = matrix(0, nrow = 101, ncol = 2)

for(i in 0:100) {
  Lambda = i/100
  OptimizationOutput = SSAGAcomp(fitness_function_M_EF, 2000, 10, 100, dim(ReturnData)[2], alpha, omega, Lambda)
  EfficientFrontier[i,1] = OptimizationOutput$Returns
  EfficientFrontier[i,2] = OptimizationOutput$Risk
}
