ICATemporalNet.boot <- function(Yts, N, Ntime, node, ncomp, Ta = 5, Tc = ceiling(sqrt(Ntime)), nb = 100, trueA = NA, tGamma = NA){
  
  # Yts: raw data
  # N: sample size
  # Ntime: total number of time points
  # ncomp: maximum number of ICs to be chosen
  # Ta: use t<=Ta time points to estimate temporal network A
  # Tc: use t>Tc time points to estimate contemporaneous network Gamma
  # nb: number of bootstraps
  # trueA: the true value of temporal network A
  # tGamma: the true value of contemporaneous network Gamma
  
  
  library(MASS)
  library(fastICA)
  library(gdata)


  estNet <- ICATemporalNet(Yts, N, Ntime, node, ncomp, Ta, Tc)
  
  estRts = estNet$estRts
  estUts = estNet$estUts
  estWt = estNet$estWt
  estS = estNet$estS
  
  
  A = estNet$A
  Omega = estNet$Omega
  Gamma = c(estNet$Gamma)

  
  ### bootstrap residuals estRts instead of Yts
  Ab=NULL
  indexbs=NULL
  Omegab= NULL
  Gammab = NULL

  for (b in 1:nb){
    
    cat(b)
    
    ## generate bootstrap data
    indexb = sample(1:N,size = N,replace = T)
    indexbs = rbind(indexbs, indexb)
    
    ## Si: permutate S
    newS = estS[indexb,]
    
    ## regenerate e from a multivarate normal distribution with covariance matrix being estimated Omega
    ei = mvrnorm(n = N, mu = rep(0,node), Sigma = matrix(Omega, node, node))
    
    ## obtain G0
    newGi0 = Yts[[1]][indexb,] - newS%*%estWt[[1]] - ei/sqrt(Ntime)
    
    Utb = list()
    Gtb = list()
    Ytb = list()
    etb = list()
    Rtb = list()
    Gb = newGi0
    
    ## permutate Yt
    Ytb[[1]] = Yts[[1]][indexb,]
    
    ## obtain new Rt
    Rtb[[1]] = Yts[[1]][indexb,] - newS%*%estWt[[1]]
    
    for (t in 2:Ntime){
      
      Gb = Gb%*%t(matrix(A, node, node))
      Utb[[t]] = newS%*%estWt[[t]]
      
      Gtb[[t]] = Gb
      
      eb = mvrnorm(n = N, mu = rep(0,node), Sigma = matrix(Omega, node, node))
      etb[[t]] = eb/sqrt(Ntime)
      
      Rtb[[t]] = Gtb[[t]] + etb[[t]]
      
      Ytb[[t]] = Utb[[t]] + Gtb[[t]] + etb[[t]]
      
    }
    
    
    estNet.b <- ICATemporalNet(Ytb, N, Ntime, node, ncomp, Ta, Tc)
    
    Ab = rbind(Ab, c(estNet.b$A))
    
    Omegab= rbind(Omegab, c(estNet.b$Omega))

    Gammab = rbind(Gammab, c(estNet.b$Gamma))
  }
  
  # standard error of bootstrapped A, Gamma
  seA = apply(Ab,2,sd)
  seOmega = apply(Omegab, 2, sd)
  seGamma = apply(Gammab, 2, sd)
  
  # empirical coverage probability and coverage length of A 
  A_uqi = apply(Ab,2,function(x){quantile(x,prob=0.975)})
  A_lqi = apply(Ab,2,function(x){quantile(x,prob=0.025)})
  
  ### basic bootstrap (reverse percentile interval)
  A_upperb = 2*A - A_lqi
  A_lowerb = 2*A - A_uqi
  
  length_A = A_uqi- A_lqi
  
  if (!is.null(trueA)){
    coverA_b = (A_upperb>=trueA)&(A_lowerb<=trueA)
  }else{
    coverA_b = NA
  }
  
  Gamma_uqi = apply(Gammab,2,function(x){quantile(x,prob=0.975)})
  Gamma_lqi = apply(Gammab,2,function(x){quantile(x,prob=0.025)})
  
  length_Gamma = Gamma_uqi - Gamma_lqi
  
  
  ### basic bootstrap (reverse percentile interval)
  Gamma_upperb = 2*c(Gamma) - Gamma_lqi
  Gamma_lowerb = 2*c(Gamma) - Gamma_uqi

  if (!is.null(tGamma)){
    coverGamma_b = (Gamma_upperb>= tGamma)&(Gamma_lowerb<= tGamma)
  }else{
    coverGamma_b = NA
  }
  
  
  return (list("estNet" = estNet, "estRts" = estRts, "estS"=estS, "estUts" = estUts, "estWt" = estWt,  "A" = A, "Omega" = Omega, "Gamma" = Gamma, "seA" = seA, "seOmega" = seOmega, "seGamma" = seGamma,
               "A_upperb" = A_upperb, "A_lowerb" = A_lowerb, "length_A" = length_A, "coverA_b" = coverA_b,
               "Gamma_upperb" = Gamma_upperb, "Gamma_lowerb" = Gamma_lowerb, "length_Gamma" = length_Gamma, "coverGamma_b" = coverGamma_b))
  
  
}
