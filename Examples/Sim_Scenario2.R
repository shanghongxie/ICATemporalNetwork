## This file contains simulations of Scenario 2

########################################################
#######       Simulation Scenario 2 node = 5    ########
########################################################

#################
## Generate Data 
#################

## get the directory from Julia app
JULIA_HOME = "/Applications/Julia-1.4.app/Contents/Resources/julia/bin"

library(JuliaCall)
library(diffeqr)

julia_setup(JULIA_HOME = "/Applications/Julia-1.4.app/Contents/Resources/julia/bin")

de = diffeqr::diffeq_setup()

node=5
B=matrix(0,node,node)
diag(B)=-0.2

B[1,2]=0.3
B[3,2]=0.3
B[4,3]=0.3
B[4,5]=0.3

trueB = B
library(expm)
trueA = expm(B)



### Systems of SDEs with Diagonal Noise
f <- JuliaCall::julia_eval("
function f(dr,r,p,t)
  dr[1] = 0.3*r[2]-0.2*r[1]
  dr[2] = -0.2*r[2]
  dr[3] = 0.3*r[2] - 0.2*r[3]
  dr[4] = 0.3*r[3] - 0.2*r[4] + 0.3*r[5]
  dr[5] =  - 0.2*r[5]
end")


g <- JuliaCall::julia_eval("
function g(dr,r,p,t)
  dr[1] = 0.1
  dr[2] = 0.1
  dr[3] = 0.1
  dr[4] = 0.1
  dr[5] = 0.1
end")


N=100

Ntime=200
tspan <- list(0.0,200.0)
u0 = c(1,rep(0,node-1))

## Generate Gaussian components from SDE: including G(t), e(t)
solveSDE <- function(f,g, N, Ntime, tspan, u0){
  
  # f: specify in SDE systems, which represents dG = BG
  # g: white noise in SDE, which represents e(t) 
  # N: sample size
  # Ntime: total number of time points
  # tspan: time span for SDE
  # u0: initial for SDE

  library(expm)
  
  Ri=list()
  for (i in 1:N){
    
    JuliaCall::julia_assign("u0", u0)
    JuliaCall::julia_assign("tspan", tspan)
    prob <- JuliaCall::julia_eval("SDEProblem(f,g,u0,tspan)")
    soli =de$solve(prob,saveat=1)
    Ri[[i]] = as.data.frame(t(sapply(soli$u,identity)))
    
  }
  Rts = list()
  
  for (t in 1:(Ntime)){
    Rt = NULL
    for (i in 1:N){
      Rt = rbind(Rt,Ri[[i]][t+1,])  ## exclude 1st time point
    }
    Rts[[t]] = Rt
  }
  
  Rts = lapply(Rts,function(x){x = as.matrix(x)})
  
  return("Rts" = Rts)
  
}

Rts = solveSDE(f,g, N, Ntime, tspan, u0)

### Add non-Gaussian components
nIC=3
########################
set.seed(100)
w0 = matrix(rnorm(node*nIC, mean =5), nrow = node, ncol = nIC)

wts=list()

for (t in 1:Ntime){
  scalet = t/Ntime
  
  wt = w0
  
  wt[1,1] = wt[1,1] +  5*scalet
  
  wt[2,2] = wt[2,2] + 5*scalet^2
  
  wt[3,3] = wt[3,3] + 5*sin(2*scalet)
  
  wt[4,3] = wt[4,3] + 5*cos(3*scalet)
  
  wts[[t]] = wt
  
}


source("genData.R")
set.seed(100)
min = 1
max = 3
Data = GenICA(Rts, wts, min = min, max = max, N, nIC)



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

result <- ICATemporalNet(Data$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)))
estA = result$A
estB = logm(estA)

nIC = result$nIC

AMSE = sum((estA - trueA)^2)
BMSE = sum((estB - trueB)^2)



#### estimation  with bootstraps
set.seed(123)
boot.result <- ICATemporalNet.boot(Data$Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)), nb = 100)

estA = boot.result$A
estB = logm(estA)

AMSE = sum((estA - trueA)^2)
BMSE = sum((estB - trueB)^2)


Aave = estA
Bave = estB
Tedge = (trueB!=0)*1
diag(Tedge) = 0
Fedge = (trueB==0)*1


ATedge = (trueA!=0)*1
diag(ATedge) = 0
AFedge = (trueA==0)*1

cutoff = seq(0.15,0.25,length.out = 5)

Acutoff = seq(0.02,0.3,length.out = 5)
library(pracma)
diag(Bave)=0

### AUC of B
TPs=NULL
FPs=NULL
TNs=NULL
FNs=NULL
for (cut in cutoff){
  Bcut = Bave
  Bcut[abs(Bcut)<cut]=0
  estBnet = (Bcut!=0)*1
  estBTP = estBnet*Tedge
  estBFP = estBnet*Fedge
  
  TPs=c(TPs,sum(estBTP))
  
  FPs=c(FPs,sum(estBFP))
  
  TNs=c(TNs,sum((estBnet==0)*Fedge))
  
  FNs=c(FNs,sum((estBnet==0)*Tedge))
  
}
TPRs=TPs/(TPs+FNs)
FPRs=FPs/(TNs+FPs)

fpr=FPRs
fpr=sort(fpr)
fpr=c(0,fpr,1)

tpr = sort(TPRs)
tpr = c(0,tpr,1)

BAUC = trapz(fpr,tpr)

### AUC of A
diag(Aave) = 0

TPs = NULL
FPs = NULL
TNs = NULL
FNs = NULL

for (cut in Acutoff){
  Acut = Aave
  Acut[abs(Acut)<cut] = 0
  estAnet = (Acut!=0)*1
  estATP = estAnet*ATedge
  estAFP = estAnet*AFedge
  
  TPs=c(TPs,sum(estATP))
  
  FPs=c(FPs,sum(estAFP))
  
  TNs=c(TNs,sum((estAnet==0)*AFedge))
  
  FNs=c(FNs,sum((estAnet==0)*ATedge))
  
}
TPRs=TPs/(TPs+FNs)
FPRs=FPs/(TNs+FPs)

fpr=FPRs
fpr=sort(fpr)
fpr=c(0,fpr,1)

tpr=sort(TPRs)
tpr=c(0,tpr,1)

AAUC = trapz(fpr,tpr)



########################################################
#######       Simulation Scenario 2 node = 10   ########
########################################################

node=10
B=matrix(0,node,node)
diag(B)=-0.2

B[1,2]=0.3
B[1,5]=0.3
B[3,2]=0.3
B[3,8]=0.3
B[4,5]=0.3
B[4,6]=0.3
B[5,10]=0.3
B[6,8]=0.3
B[6,9]=0.3
B[7,1]=0.3
B[8,9]=0.3

########################################################
#######       Simulation Scenario 2 node = 20   ########
########################################################

node = 20
B = matrix(0,node,node)
diag(B) = -0.2

B[1,5]=0.3
B[3,2]=0.3
B[3,11]=0.3
B[4,5]=0.3
B[4,6]=0.3
B[5,10]=0.3
B[6,8]=0.3
B[6,9]=0.3
B[7,1]=0.3
B[7,10]=0.3
B[8,9]=0.3
B[8,20]=0.3
B[11,12]=0.3
B[11,15]=0.3
B[13,12]=0.3
B[13,18]=0.3
B[14,15]=0.3
B[14,16]=0.3
B[15,20]=0.3
B[16,18]=0.3
B[17,20]=0.3
B[18,19]=0.3
