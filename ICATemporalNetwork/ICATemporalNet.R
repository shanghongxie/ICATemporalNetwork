ICATemporalNet <- function(Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime))){
  
  # Yts: raw data
  # N: sample size
  # Ntime: total number of time points
  # ncomp: maximum number of ICs to be chosen
  # Ta: use t<=Ta time points to estimate  temporal network A
  # Tc: use t>Tc time points to estimate contemporaneous network Gamma
  
  library(fastICA)
  library(gdata)
  
  ## Remove non-Gaussian components
  estIC = estICA(Yts, N, node, Ntime, ncomp)
  estRts = estIC$estRts # residuals after removing non-Gaussian components
  estS = estIC$estS # independent components S
  estUts = estIC$estUts # non-Gaussian signal U(t)
  estWt = estIC$coff # w(t)
  nIC = estIC$nIC

  ## Use residuals to estimate temporal network A  by moment estimations
  Results = temporalNet(estRts, N, Ntime, node, Ta)
  
  A = Results$A
  
  
  ## Estimate contemporaneous network Gamma by moment estimations
  tt = A%x%A-diag(node*node)
  
  IA=apply(tt,1,function(x){
    m=matrix(x,node,node)
    m=m+t(m)
    diag(m)=diag(m)/2
    return(mt=lowerTriangle(m,diag=T))
  })
  IA=t(IA)
  
  indext = which(lower.tri(matrix(0,node,node),diag=T)==TRUE)
  
  Alongt=IA[indext,]
  
  rm(IA)
  
  
  count=1
  Qs=NULL
  Q = NULL
  
  nT2 = Ntime-Tc
  
  covtt = matrix(0,node,node)
  covtt1 = matrix(0,node,node) 
  
  for (count in Tc:(Ntime-1)){
    covtt = covtt + cov(estRts[[count]])
    covtt1 = covtt1 + cov(estRts[[count+1]])
  }
  
  covtt = covtt/nT2
  covtt1 = covtt1/nT2
  
  Q = A%*%covtt%*%t(A)-covtt1
  q = lowerTriangle(Q, diag = T)
  AveOmega =  Ntime*solve(Alongt)%*%q
  
  Omega = matrix(0,node,node)
  
  Omega[lower.tri(Omega,diag=T)]=AveOmega
  Omega = Omega + t(Omega)
  diag(Omega) = diag(Omega)/2
  
  # Gamma: contemponranoues network
  Gamma = ginv(Omega)
  
  
  return (list("estIC" = estIC, "estRts" = estRts, "estS"=estS, "estUts" = estUts, "estWt" = estWt, "nIC" = nIC, "A" = A, "Omega" = Omega, "Gamma" = Gamma))
  
}