## This file contains an example of implementing ICATemporalNet.R and ICATemporalNetBoot.R


################################################
#######       Data  Generation          ########
################################################

N=100
Ntime = 200

## Set temporal network A
node=5
A = matrix(0,node,node)
diag(A)=0.8

A[1,2]=0.2
A[3,2]=0.2
A[4,3]=0.2
A[4,5]=0.2

sum(A!=0)  ## 30
diag(A)=0.8

trueA = A

## Set Omega and contemporaneous network Gamma
node=5
pcorS=matrix(0,node,node)
pcorS[2,1]=0.3
pcorS[2,4]=0.3
pcorS[3,5]=0.3
pcorS=pcorS+t(pcorS)
diag(pcorS)=0

sigmas = rep(0.1,5)

Gamma = matrix(0,node,node)
for (s in 1:node){
  for (r in 1:node){
    if (s == r){
      Gamma[s,s] = sigmas[s]
    } else{
      Gamma[s,r]=-pcorS[s,r]*sqrt(sigmas[s])*sqrt(sigmas[r])
    }
    
  }
}

Omega = solve(Gamma)
tOmega = Omega
tGamma = Gamma

## G0 covariance (use identity matrix)
cov = matrix(0,node,node)
diag(cov)=1


source("genData.R")

## Generate Gaussian components, including G(t), e(t)
GaussianData = genGaussian(N, Ntime, trueA, tOmega, node, cov)

## Set w(t)
nIC=3

set.seed(100)
w0 = matrix(rnorm(node*nIC, mean =5), nrow = node, ncol = nIC)


wts = list()

for (t in 1:Ntime){
  scalet = t/Ntime
  
  wt = w0
  
  wt[1,1] = wt[1,1] +  5*scalet
  wt[2,2] = wt[2,2] + 5*scalet^2
  wt[3,3] = wt[3,3] + 5*sin(2*scalet)
  wt[4,3] = wt[4,3] + 5*cos(3*scalet)

  wts[[t]] = wt
  
}

## Add non-Gaussian components U(t)=w(t)S to obtain final data Y(t)
min=1
max=3
set.seed(100)
Data = GenICA(GaussianData$Rts, wts, min, max, N, nIC)




################################################
#######      Estimation                 ########
################################################


source("estICA.R") # estimation for non-Gaussian components
source('temporalNet.R') # estimation for temporal network A
source('ICATemporalNet.R') # entire estimation procedures for non-Gaussian components, A, and contemporanoues network Gamma
source('ICATemporalNetBoot.R')

library(mgcv)

ncomp = 5 # maximum number of nIC
node = 5

result <- ICATemporalNet(Data$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)))
estA = result$A
estGamma = result$Gamma
nIC = result$nIC

AMSE = sum((estA - trueA)^2)
GammaMSE = sum((estGamma - tGamma)^2)


#### estimation  with bootstraps
boot.result <- ICATemporalNet.boot(Data$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)), nb = 100)

estA = boot.result$A
estGamma = boot.result$Gamma
AMSE = sum((estA - trueA)^2)
GammaMSE = sum((estGamma - tGamma)^2)

seA = boot.result$seA
seGamma = boot.result$Gamma

A_upperb = boot.result$A_upperb
A_lowerb = boot.result$A_lowerb
length_A = boot.result$length_A


Gamma_upperb = boot.result$Gamma_upperb
Gamma_lowerb = boot.result$Gamma_lowerb
length_Gamma = boot.result$length_Gamma







