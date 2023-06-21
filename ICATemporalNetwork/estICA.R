# Estimate non-Gaussian components and remove from Y(t)
estICA <- function(Yts, N, node, Ntime, ncomp){

  # Yts: raw data
  # N: sample size
  # node: number of nodes
  # Ntime: total time points
  # ncomp: maximum number of ICs to be chosen

  
  ### ICA estimation
  Yave=matrix(0,N,node)

  Y=NULL
  for (t in 1:Ntime){
    Yt=Yts[[t]]

    ### normalize at each time point for each region separately
    
    ## first take average
    meanYt=apply(Yt,2,mean)
    
    centerYt=apply(Yt,1,function(x){x-meanYt})
    centerYt=t(centerYt)
    
    ## standardize centered Yt
    sdYt=apply(centerYt,2,sd)
    
    stdYt=apply(centerYt,1,function(x){x/sdYt})
    
    stdYt=t(stdYt)
    Yave=Yave+stdYt
  }
  Yave=Yave/Ntime
  
  
  Y = NULL
  for (i in 1:N){
    Y = rbind(Y, t(sapply(Yts,function(x){x[i,]}, simplify = "array")))  ## Ntime * node, i = 1: t= 1, ..., T
  }
  
  
  
  ## choose number of components by PCA based on MDL
  pca.fit <- prcomp(Yave, center=T, rank = 10)
  # summary(pca.fit)
  
  eigen = pca.fit$sdev^2
  num = length(eigen)
  
  
  MDLs=NULL
  for (k in 1:min(ncomp,node-1)){
   
    likelihood = (node-k)*N*log(prod((eigen[(k+1):node])^(1/(node-k)))/(sum(eigen[(k+1):node])/(node-k)))
    MDL = -likelihood+0.5*k*(2*node-k)*log(N)
    MDLs=c(MDLs,MDL)

  }
  
  k = 3
  repeat{

    if (abs(MDLs[k]-MDLs[k-1])<0.1*abs(MDLs[k-1]-MDLs[k-2])) break
    
    if (k == length(MDLs)) break
    k = k +1
    
  }
  
  k1 = k
  
  k2 = which.min(MDLs)
  
  if (k1 <= k2){
    nIC = k1-1
  } else{
    nIC = k2
  }
  
 
  
  ### perform ICA on the average Y
  library(fastICA)
  res <- fastICA(Yave, nIC)
  
  ### mixing matrix
  B=res$A
  ### unmixing matrix
  W=res$W
  ### independent component
  S=res$S   ### subject*component
  
  
  ### linear regression
  
  ## for each region separately
  estRts=list()  ## include G(t)+1/sqrt(T)\epsilon(t)
  cofft=list()
  estUts=list()
  for (t in 1:Ntime){
    Yt=Yts[[t]]
    
    lmfit.t <- lm(as.matrix(Yt) ~ S-1)
    estUts[[t]] = lmfit.t$fitted.values
    cofft[[t]]=lmfit.t$coefficients
    estRt=lmfit.t$residuals
    estRts[[t]]=estRt
  }
  
  # residuals estRts will be used to estimate temporal and contemporaneous networks
  
  return(list("estS"=S,"estB"=B,"estW"=W,"estRts"=estRts, "cofft"=cofft,"estUts"=estUts, "nIC" = nIC, "MDLs" = MDLs, "eigen" = eigen,  "pca.fit" = pca.fit))
  
}


