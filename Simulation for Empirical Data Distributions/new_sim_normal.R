generateX2 = function(Sigma_e,W,rhox){
  n = dim(W)[1]
  p = dim(W)[2]
  #Step 1: Sample a lot of WW from the empirical distribution of W and sample a lot of EE from the estimated distribution of measurement error independently, 
  NN=10000
  WW=EE=XX=matrix(NA,nrow=NN,ncol=p)
  for (j in 1:p){
  WW[,j] = sample(1:NN,n,replace = T)
  }
  
  EE = mvrnorm(n=NN,mu=rep(0,p),Sigma=Sigma_e)
  
  for (j in 1:p){
    XX[,j] = WW[,j] - EE[,j]
  }
  
  #Step 2:Sample n individual Z_i from multivariate normal distribution N(0,AR(\rho))
  Z<-mvrnorm(n=n,mu=rep(0,p),Sigma=rhox^abs(outer(1:p,1:p,"-")))
  #Step 3: For each i, j, compute X_{ij}= F_j^{-1}(\Phi(Z_{ij}))
  
  Xtilde=matrix(NA,nrow=n,ncol=p)
  for (i in 1:n) {
    for (j in 1:p) {
      Xtilde[i, j] <- quantile(XX[, j], pnorm(Z[i, j]))
    }
  }
  #F_XX[[j]] = ecdf(XX[,j])
  return(Xtilde)
}

