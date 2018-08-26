rm(list=ls())
set.seed(21)

library("MASS")   ############################### Package to sample from the multivariate Normal

################################################# Function to generate fractionally integrated processes using the fast Fourier transform
################################################# See Jensen and Nielsen (2014)
#################################################
fracdiff = function(x, dd){
  iT = length(x)
  n = nextn(2*iT - 1, 2)
  k = 1:(iT-1)
  b = c(1, cumprod((k - dd - 1)/k))
  dx = fft(fft(c(x, rep(0, n - iT))) * fft(c(b, rep(0, n - iT))), inverse = T) / n;
  return(Re(dx[1:iT]))
}
################################################# Overall parameters
#################################################
#################################################

iR = 10^4 ############################### Number of replications
dy = 0.25 
d1 = 0.20
d2 = 0.15
d3 = 0.10

################################################# Uncorrelated case
#################################################
#################################################
#################################################
################################################# 

## T = 100

iT = 10^2 ############################### Sample size

vmus = rep(0,3) ############################### Mean of multivariate noise processes
Msigma = matrix(c(1,0,0, 0,0.75, 0, 0, 0, 0.4),3,3) ##### Variance of multivariate noise processes

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)

for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)  ######## Multivariate noise processes
  Y = fracdiff(rnorm(iT[1])*2,-dy) ######################### Regressand 
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1] #################### Regressors
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05 ##### Rejection F-test
  r2s[i]= summary(reg)$r.squared                  ################# R-squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05  ################ Rejection t-test
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals) ##Durbin-Watson
  
}


table_100 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 1000

iT = 10^3

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_1000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 10000

iT = 10^4

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

##

table_uncorr = cbind(table_100,table_1000,table_10000)

################################################# 2 variables correlated
#################################################
#################################################
#################################################
################################################# 


## T = 100

iT = 10^2

vmus = rep(0,3)
Msigma = matrix(c(1,0.4,0, 0.4,0.75, 0,0, 0, 0.4),3,3)

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}


table_100 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 1000

iT = 10^3

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_1000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 10000

iT = 10^4

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

##

table_2corr = cbind(table_100,table_1000,table_10000)



################################################# 3 variables correlated
#################################################
#################################################
#################################################
################################################# 


## T = 100

iT = 10^2

vmus = rep(0,3)
Msigma = matrix(c(1,0.4,0.6, 0.4,0.75, 0.3,0.6, 0.3, 0.4),3,3)

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}


table_100 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 1000

iT = 10^3

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_1000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

## T = 10000

iT = 10^4

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  u = mvrnorm(n=iT, vmus, Msigma)
  Y = fracdiff(rnorm(iT[1])*2,-dy)
  x1 = fracdiff(rnorm(iT[1]),-d1)+u[,1]
  x2 = fracdiff(rnorm(iT[1])*0.75,-d2)+u[,2]
  x3 = fracdiff(rnorm(iT[1])*0.75,-d3)+u[,3]
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

##

table_3corr = cbind(table_100,table_1000,table_10000)


################################################# Overall results
#################################################
#################################################
#################################################
################################################# 

table_final = cbind(table_uncorr,table_2corr,table_3corr)
table_final