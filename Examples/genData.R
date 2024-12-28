
library(MASS)

# Generate Gaussian components, including G(t) and Gaussian noise e(t)
genGaussian <- function(N, Ntime, A, Omega, node, cov){
  
  # N: sample size
  # Ntime: total number of time points
  # A: temporal network A
  # Omega: cov(et)
  # node: number of nodes
  # cov: G0 covariance
  
  
  library(MASS)
  G0 = mvrnorm(n = N,mu = rep(0,node), Sigma = cov)
  Gt = G0
  
  ets=list()
  Rts=list()
  Gts=list()
  
  for (t in 1:Ntime){
    et = mvrnorm(n=N,mu=rep(0,node),Sigma = Omega)
    ets[[t]] = et
    Rt = Gt + et/sqrt(Ntime)  
    Rts[[t]] = Rt
    Gts[[t]] = Gt
    
    Gt = Gt%*%t(A)
  }
  
  return(list("Rts"=Rts, "Gts"=Gts, "G0" = G0, "ets" = ets))
}


# Add non-Gaussian components
GenICA <-function(Rts, wt, min, max, N, nIC){
  
  # Rts: generated Gaussian data from genGaussian() ,  Gt + et/sqrt(Ntime) 
  # wt: w(t)
  # min: generate S from a uniform distribution. minimal value
  # max: maximum value
  # nIC: number of ICs (M in the paper)


  ## Add two independent uniform components S (not depend on t)
  
  S = matrix(runif(nIC*N, min = min, max = max), N, nIC)  ## N *nIC
  
  Yts=list()

  Uts = list()
  for (t in 1:Ntime){
    Ut = S%*%t(wt[[t]]) # non-Gaussian components
    Yt = Rts[[t]] + S%*%t(wt[[t]]) # final data
    
    Uts[[t]] = Ut
    Yts[[t]] = Yt
    
  }

  return(list("S" = S,"wt" = wt,"Yts" = Yts, "Uts" = Uts))
}


