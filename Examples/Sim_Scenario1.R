## This contains simulations of Scenario 1

########################################################
#######       Simulation Scenario 1 node = 5    ########
########################################################

################################################
#######       Data  Generation          ########
################################################

N = 100
Ntime = 200

## Set temporal network A
node = 5
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

library(doParallel)
library(foreach)


detectCores()

registerDoParallel(3)
getDoParWorkers()

## Generate Gaussian components, including G(t), e(t)

set.seed(123)
nSim = 100
GaussianDatas = foreach(iSim=1:nSim) %dopar% {genGaussian(N, Ntime, trueA, tOmega, node, cov)}



## Set w(t)
nIC=3

### most of coefficients are timewiseant
set.seed(100)
w0 = matrix(rnorm(node*nIC, mean = 5), nrow = node, ncol = nIC)


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
min = 1
max = 3

set.seed(100)
Data = foreach(iSim=1:nSim) %dopar% {GenICA(GaussianDatas[[iSim]]$Rts, wts, min, max, N, nIC)}


#################
##  Estimation
#################

source("estICA.R") # estimation for non-Gaussian components
source('temporalNet.R') # estimation for temporal network A
source('ICATemporalNet.R') # entire estimation procedures for non-Gaussian components, A, and contemporanoues network Gamma
source('ICATemporalNetBoot.R')

library(mgcv)

ncomp = 5 # maximum number of nIC
node = 5

result <- ICATemporalNet(Data[[1]]$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)))
estA = result$A
estGamma = result$Gamma
nIC = result$nIC

AMSE = sum((estA - trueA)^2)
GammaMSE = sum((estGamma - tGamma)^2)


#### estimation  with bootstraps

set.seed(123)
boot.results <-foreach(iSim=1:nSim) %dopar% {ICATemporalNet.boot(Data[[iSim]]$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)), nb = 100, trueA = trueA, tGamma = tGamma)}



#################
##  Evaluation
#################

library(pracma)

Tedge = (trueA!=0)*1
diag(Tedge) = 0
Fedge = (trueA==0)*1


TPedge = (tGamma!=0)*1
diag(TPedge) = 0
FPedge = (tGamma==0)*1
diag(FPedge) = 0

AUCs = NULL
AMSEs = NULL
GammaMSEs = NULL
length_As = NULL
coverA_bs = NULL
length_Gammas = NULL
coverGamma_bs = NULL
GammaAUCs = NULL

for (iSim in 1:nSim){
  cat(iSim)
  
  boot.result = boot.results[[iSim]]
  estA = boot.result$A
  estGamma = boot.result$Gamma
  AMSE = sum((estA - trueA)^2)
  GammaMSE = sum((estGamma - tGamma)^2)
  
  AMSEs = c(AMSEs, AMSE)
  GammaMSEs = c(GammaMSEs, GammaMSE)
  
  length_A = boot.result$length_A
  length_As = rbind(length_As, length_A)
  
  coverA_b = boot.result$coverA_b
  coverA_bs = rbind(coverA_bs, c(coverA_b))
  
  Gamma_upperb = boot.result$Gamma_upperb
  Gamma_lowerb = boot.result$Gamma_lowerb
  length_Gamma = boot.result$length_Gamma
  coverGamma_b = boot.result$coverGamma_b
  
  length_Gammas = rbind(length_Gammas, length_Gamma)
  coverGamma_bs = rbind(coverGamma_bs, c(coverGamma_b))
  
  ## AUC of A
  cutoff = seq(0.15,0.25,length.out = 5)
  Aave = estA
  diag(Aave)=0
  
  TPs=NULL
  FPs=NULL
  TNs=NULL
  FNs=NULL

  for (cut in cutoff){
    Acut = Aave
    Acut[abs(Acut)<cut]=0
    estAnet = (Acut!=0)*1
    estATP = estAnet*Tedge
    estAFP = estAnet*Fedge
    
    TPs=c(TPs,sum(estATP))
    
    FPs=c(FPs,sum(estAFP))
    
    TNs=c(TNs,sum((estAnet==0)*Fedge))
    
    FNs=c(FNs,sum((estAnet==0)*Tedge))

  }
  TPRs=TPs/(TPs+FNs)
  
  
  FPRs=FPs/(TNs+FPs)
  
  fpr=FPRs
  fpr=sort(fpr)
  fpr=c(0,fpr,1)
  
  tpr=sort(TPRs)
  tpr=c(0,tpr,1)
  
  auc=trapz(fpr,tpr)
  AUCs = c(AUCs,auc)
  
  ### Calculate AUC of Gamma
  cutoff = seq(0.01, 0.05,length.out = 5)
  
  GammaTPs=NULL
  GammaFPs=NULL
  GammaTNs=NULL
  GammaFNs=NULL
  
  for (cut in cutoff){
    
    Gammacut=c(matrix(Gamma, node, node))
    
    Gammacut[abs(Gammacut)<cut]=0
    estGammanet = (Gammacut!=0)*1
    estGammaTP = estGammanet*TPedge
    estGammaFP = estGammanet*FPedge
    
    GammaTPs = c(GammaTPs,sum(estGammaTP))
    
    GammaFPs = c(GammaFPs,sum(estGammaFP))
    
    GammaTNs = c(GammaTNs,sum((estGammanet==0)*FPedge))
    
    GammaFNs = c(GammaFNs,sum((estGammanet==0)*TPedge))
    
  }

  
  GammaTPRs = GammaTPs/(GammaTPs+GammaFNs)
  GammaFPRs = GammaFPs/(GammaTNs+GammaFPs)
  
  Gammafpr = GammaFPRs
  Gammafpr = sort(Gammafpr)
  Gammafpr = c(0,Gammafpr,1)
  
  Gammatpr = sort(GammaTPRs)
  Gammatpr = c(0,Gammatpr,1)
  
  Gammaauc = trapz(Gammafpr,Gammatpr)
  
  GammaAUCs = c(GammaAUCs, Gammaauc)
}


## coverage probablity of A
cpAb = matrix(apply(coverA_bs,2,sum),node,node)/nSim
cpAb_True = cpAb*Tedge
cpAb_False = cpAb*Fedge
cpAb = sum(cpAb_True+cpAb_False)/(sum(Tedge)+sum(Fedge))  
cpAb

## coverage length of A
Avelength = apply(length_As,2,mean)
Avelength95 = (sum(Avelength*c(Tedge))+sum(Avelength*c(Fedge)))/sum(Tedge+Fedge)

## A's results
AResults =c(mean(AMSEs), mean(AUCs), cpAb, Avelength95)
names(AResults)=c("MSE","AUC", "basic coverage probability",  "basic coverage length")

AResults

## coverage probablity of Gamma
cpGammab = matrix(apply(coverGamma_bs,2,sum),node,node)/nSim

cpGamma95_b_True = cpGammab*TPedge
cpGamma95_b_False = cpGammab*FPedge

cpGammab = sum(cpGamma95_b_True + cpGamma95_b_False)/(sum(TPedge)+sum(FPedge))  ## 0.9610263
cpGammab

## coverage length of Gamma
Avelength_Gamma = apply(length_Gammas, 2, mean)
Avelength_Gamma95 = (sum(Avelength_Gamma*c(TPedge))+sum(Avelength_Gamma*c(FPedge)))/sum(TPedge+FPedge)

## Gamma's results
GammaResults =c(mean(GammaMSEs), mean(GammaAUCs), cpGammab, Avelength_Gamma95)
names(GammaResults)=c("MSE","AUC", "basic coverage probability",  "basic coverage length")

GammaResults


################################################################
#######       Simulation Scenario 1 node = 10   Setting ########
################################################################

node=10
A = matrix(0,node,node)
diag(A) = 0.8

## sparse
A[1,2]=0.2
A[3,2]=0.2
A[4,5]=0.2
A[6,8]=0.2
A[6,9]=0.2

sum(A!=0)  ## 30
diag(A)=0.8

trueA=A

N = 100
Ntime = 200

node = 10
pcorS = matrix(0,node,node)
pcorS[2,1] = 0.3
pcorS[2,4] = 0.3
pcorS[3,5] = 0.3
pcorS = pcorS+t(pcorS)
diag(pcorS) = 0

sigma = 0.1
sigmas = c(rep(0.1,5), rep(0.5, 5))

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

eigen(Gamma)$value

tGamma = Gamma
Omega = solve(Gamma)
tOmega = Omega


## G0 covariance (use identity matrix)
cov=matrix(0,node,node)
diag(cov)=1


################################################################
#######       Simulation Scenario 1 node = 20   Setting ########
################################################################

node=20
A=matrix(0,node,node)
diag(A)=0.8

A[1,5]=0.2
A[3,2]=0.2
A[3,8]=0.2
A[4,5]=0.2
A[4,6]=0.2
A[5,10]=0.2
A[6,8]=0.2
A[6,9]=0.2
A[7,10]=0.2
A[7,1]=0.2
A[8,9]=0.2
A[8,20]=0.2 
A[3,11]=0.2
A[11,12]=0.2
A[11,15]=0.2
A[13,12]=0.2
A[13,18]=0.2
A[14,15]=0.2
A[14,16]=0.2
A[15,20]=0.2
A[16,18]=0.2
A[16,19]=0.2
A[17,20]=0.2
A[18,19]=0.2

sum(A!=0) 
diag(A)=0.8


trueA=A

N=100
Ntime = 200

node=20
pcorS=matrix(0,node,node)
pcorS[2,1]=0.3
pcorS[2,4]=0.3
pcorS[3,5]=0.3

pcorS[2+5,1+5]=0.3
pcorS[2+5,4+5]=0.3
pcorS[3+5,5+5]=0.3


pcorS[2+10,1+10]=0.3
pcorS[2+10,4+10]=0.3
pcorS[3+10,5+10]=0.3
pcorS=pcorS+t(pcorS)
diag(pcorS)=0

sigma = 0.1
sigmas = c(rep(0.1, 5), rep(0.4,5), rep(0.4, 5), rep(0.5, 5))


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


eigen(Gamma)$value

tGamma = Gamma
Sigma = solve(Gamma)
tSigma = Sigma


