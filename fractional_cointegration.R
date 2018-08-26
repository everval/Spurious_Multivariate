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

iR = 10^4
du = 0.2    ########################## Memory of linear combination


#################################################
#################################################
#################################################
#################################################
################################################# dx=0.25
dx = 0.25

## T = 100

iT = 10^2

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))


table_lowD = cbind(table_100,table_1000,table_10000)

#################################################
#################################################
#################################################
#################################################
################################################# dx=0.35
dx = 0.35

## T = 100

iT = 10^2

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

table_midD = cbind(table_100,table_1000,table_10000)

#################################################
#################################################
#################################################
#################################################
################################################# 3 var dx=0.45
dx = 0.45

## T = 100

iT = 10^2

fts = rep(0,iR)
r2s = rep(0,iR)
dws = rep(0,iR)
tabcof = matrix(0,nrow=iR,ncol=4)


for(i in 1:iR){
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
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
  
  x1 = fracdiff(rnorm(iT[1]),-dx)
  x2 = fracdiff(rnorm(iT[1])*0.75,-dx)
  x3 = fracdiff(rnorm(iT[1])*0.75,-dx)
  Y = 0.70+0.2*x1+0.3*x2+0.4*x3+fracdiff(rnorm(iT[1])*0.5,-du)
  
  reg = lm(Y~x1+x2+x3)
  
  fts[i] = pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)<0.05
  r2s[i]= summary(reg)$r.squared
  tabcof[i,] = summary(reg)$coefficients[,4]<0.05
  dws[i] =(reg$residuals[2:iT]-reg$residuals[1:(iT-1)])%*%(reg$residuals[2:iT]-reg$residuals[1:(iT-1)]) / (reg$residuals%*%reg$residuals)
  
}

table_10000 = c(colMeans(tabcof),mean(r2s),mean(fts),mean(dws))

table_highD = cbind(table_100,table_1000,table_10000)

################################################# Overall results
#################################################
#################################################
#################################################
################################################# 

table_final = cbind(table_lowD,table_midD,table_highD)
table_final