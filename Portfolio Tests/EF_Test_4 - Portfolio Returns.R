# SSAGA v1.3 portfolio test 4 - comparing Markowitz to Extended

################### Packages ################### 

install.packages('beepr')
install.packages('compiler')
install.packages('ggplot2')
install.packages('PerformanceAnalytics')

library(beepr)
library(compiler)
library(ggplot2)
library(PerformanceAnalytics)
enableJIT(2)

################### Set WD ################### 

setwd("C:\\Users\\DIEHLAT\\Desktop")

################### FITNESS FUNCTION ################### 

fitness_function_Markowitz_EF = function(weights, returns, risk, lambda) {
  return(lambda * t(weights) %*% returns - (1 - lambda) * t(weights) %*% risk %*% weights)
}

fitness_function_Markowitz_EFHOM = function(weights, returns, risk, skew, kurt, lambda1, lambda2, lambda3, lambda4) {
  m1 = t(weights) %*% returns
  sq = t(weights) %*% risk %*% weights
  m2 = sq + m1^2
  m3 = t(weights) %*% skew %*% (weights %x% weights)
  m4 = t(weights) %*% kurt %*% (weights %x% weights %x% weights)
  sp = 1/(sq * sqrt(sq)) * (m3 - 3*m2*m1 + 2*m1^3)
  kp = 1/(sq^2) * (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4) - 3
  return(lambda1 * m1 + lambda2 * sp - lambda3 * sq - lambda4 * kp)
}

################### DATA IMPORT ################### 

Dowdata = read.csv("DowData.csv")

nRows = as.numeric(dim(Dowdata)[1])
nCols = as.numeric(dim(Dowdata)[2])

Dow_No_Date = as.matrix(Dowdata[1:nRows-1, 2:nCols])

#Generate returns data
ReturnData = matrix(1, nrow = (nRows-1), ncol = (nCols-1))

for(i in 2:(nRows-1)) {
  ReturnData[(i-1),] = try(log(Dow_No_Date[i,]) - log(Dow_No_Date[(i-1),]))
}

colnames(ReturnData) = colnames(Dow_No_Date)
ReturnData[is.na(ReturnData)] = 0
ReturnData = ReturnData[-(nRows-1),]

ReturnDimRow = as.numeric(dim(ReturnData)[1])
ReturnDimCol = as.numeric(dim(ReturnData)[2])

################### RISK AND RETURNS MATRICES ################### 

HistoricalReturns = function(Returns, TimeIndex) {
  alpha = numeric(ReturnDimCol)
  for (i in 1:(ReturnDimCol)){
    alpha[i] = mean(Returns[1:TimeIndex,i])
  }
  alpha[is.infinite(alpha)] = 0
  alpha[is.na(alpha)] = 0
  return(alpha)
}

CovarianceMatrix = function(Returns, TimeIndex) {
  omega = matrix(rep(0, (ReturnDimCol)^2), ReturnDimCol, ReturnDimCol) 
  for (j in 1:(ReturnDimCol)){
    for (k in 1:(ReturnDimCol)){
      omega[j,k] = cov(Returns[1:TimeIndex,j], Returns[1:TimeIndex,k])    
    }
  }
  omega[is.infinite(omega)] = 0
  omega[is.na(omega)] = 0
  return(omega)
}

################### PORTFOLIO TRACKING ################### 

#Compile functions
fitness_function_M_EF = cmpfun(fitness_function_Markowitz_EF) 
fitness_function_M_EFHOM = cmpfun(fitness_function_Markowitz_EFHOM) 
SSAGA_v12_comp = cmpfun(SSAGA12_PE)
SSAGA_v13_comp = cmpfun(SSAGA13_PE)

#Create vectors to hold returns
Markowitz_Returns = numeric(length(seq(252, ReturnDimRow, 20)))
Expanded_Returns = numeric(length(seq(252, ReturnDimRow, 20)))
PortfolioTracker_Mark = numeric(length(seq(252, ReturnDimRow, 20)) + 1)
PortfolioTracker_Exp = numeric(length(seq(252, ReturnDimRow, 20)) + 1)
PortfolioValueM = 1000000
PortfolioValueE = 1000000
PortfolioTracker_Mark[1] = PortfolioValueM
PortfolioTracker_Exp[1] = PortfolioValueE

Lambda = 0.5
lambda1 = 0.45
lambda2 = 0.05
lambda3 = 0.495
lambda4 = 0.005

#Create initial portfolios from 1 year of data
j = 252
CurrentPriceIndex = j + 1
ReferencePriceIndex = j - 19
alpha = HistoricalReturns(ReturnData, j)
omega = CovarianceMatrix(ReturnData, j)
CoskewTensor = M3.MM(ReturnData[1:j,], unbiased = T)
CokurtTensor = M4.MM(ReturnData[1:j,], unbiased = T)

#Run the algos
x = SSAGA_v12_comp(fitness_function_M_EF, 1500, 30, 100, dim(ReturnData)[2], alpha, omega, Lambda)
y = SSAGA_v13_comp(fitness_function_M_EFHOM, 1500, 30, 100, dim(ReturnData)[2], alpha, omega, CoskewTensor, CokurtTensor, lambda1, lambda2, lambda3, lambda4)

SharesInvested_M = floor((PortfolioValueM * x$Solution)/Dow_No_Date[CurrentPriceIndex,])
SharesInvested_HOM = floor((PortfolioValueE * y$Solution)/Dow_No_Date[CurrentPriceIndex,])

index = 1

#Iterate through, starting at one year and in increments of one month
begin = Sys.time()
for(i in seq(272, ReturnDimRow, 20)) {
  
  #Calculate P&Ls 
  CurrentPriceIndex = i + 1
  ReferencePriceIndex = i - 19
  
  Markowitz_Returns[index] = sum(SharesInvested_M * Dow_No_Date[CurrentPriceIndex,]) - sum(SharesInvested_M * Dow_No_Date[ReferencePriceIndex,])
  Expanded_Returns[index] = sum(SharesInvested_HOM * Dow_No_Date[CurrentPriceIndex,]) - sum(SharesInvested_HOM * Dow_No_Date[ReferencePriceIndex,])
  
  PortfolioValueM = PortfolioValueM + Markowitz_Returns[index]
  PortfolioValueE = PortfolioValueE + Expanded_Returns[index]
  PortfolioTracker_Mark[index] = PortfolioValueM
  PortfolioTracker_Exp[index] = PortfolioValueE
  
  index = index + 1
  
  #Set new data
  print(paste('Percent complete:', as.character(i/ReturnDimRow)))
  alpha = HistoricalReturns(ReturnData, i)
  omega = CovarianceMatrix(ReturnData, i)
  CoskewTensor = M3.MM(ReturnData[1:i,], unbiased = T)
  CokurtTensor = M4.MM(ReturnData[1:i,], unbiased = T)
  
  #Run the algos and update portfolio weights
  x = SSAGA_v12_comp(fitness_function_M_EF, 1500, 30, 100, dim(ReturnData)[2], alpha, omega, Lambda)
  y = SSAGA_v13_comp(fitness_function_M_EFHOM, 1500, 30, 100, dim(ReturnData)[2], alpha, omega, CoskewTensor, CokurtTensor, lambda1, lambda2, lambda3, lambda4)
  
  SharesInvested_M = floor((PortfolioValueM * x$Solution)/Dow_No_Date[CurrentPriceIndex,])
  SharesInvested_HOM = floor((PortfolioValueE * y$Solution)/Dow_No_Date[CurrentPriceIndex,])
  
}
end = Sys.time()

write.csv(Markowitz_Returns, "MarkowitzReturns.csv")
write.csv(Expanded_Returns, "ExpandedReturns.csv")
write.csv(PortfolioTracker_Mark, "MarkowitzPortfolio.csv")
write.csv(PortfolioTracker_Exp, "ExpandedPortfolio.csv")
write.csv(end-begin, "time.csv")

beep(); Sys.sleep(0.5); beep(); Sys.sleep(0.5); beep(); Sys.sleep(0.5); beep(); Sys.sleep(0.5); beep()
