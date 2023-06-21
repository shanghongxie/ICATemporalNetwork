### only calculate temporal network A
temporalNet <- function(Yts, N, Ntime, node, Ta){
  
  # Yts: input ICA-denoised data, the user must supply a list of Yts, where each elements is a N*K data matrix at time t. N is sample size, K is the number of nodes.
  # N: sample size
  # Ntime: number of time points
  # node: number of nodes
  # Ta: first Ta time points used for estimating temporal network A
  
  ### lag 1
  Ats = list()
  covs = list()
  covsts = NULL
  covsls = NULL
  
  interval = 1 #time interval 
  
  timeseq = seq(1,Ntime,by=interval)
  
  covts =  matrix(0, node, node)
  covls = matrix(0, node, node)
  
  # use the first Ta time points to estimate A
  for (count in 1:Ta){
    
    # cov(Y(t))
    Yt = Yts[[timeseq[count]]]
    cYt = apply(Yt,2,function(x){x-mean(x)})
    
    # cov(Y(t+1))
    Yt1 = Yts[[timeseq[count+1]]]
    cYt1 = apply(Yt1,2,function(x){x-mean(x)})
    
    # cov(Y(t+2))
    Yt2 = Yts[[timeseq[count+2]]]
    cYt2 = apply(Yt2,2,function(x){x-mean(x)})
    
    covts = covts + (t(cYt2)%*%cYt)/(N-1)
    covls = covls + (t(cYt1)%*%cYt)/(N-1)
    
    covtt=cov(Yt)
  
    covs[[count]]=covtt
  
  }
  
  covts = covts/Ta
  covls = covls/Ta
  
  A = covts%*%solve(covls) # average of A from Ta time points
  
  results <- list("A" = A, "covts" = covts, "covls" = covls)
  
  return(results) 
}