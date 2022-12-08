## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE----------------------------------------------------------------
l2NormSq = function(b){
  return( crossprod(b,b)[1,1] )
}


## -----------------------------------------------------------------------------
library(Rcpp)

## ----eval=TRUE----------------------------------------------------------------
Splicing = function(beta, d, setA, setI, kmax, TAWs,    n,p,x,y){
  sizeA=length(setA)
  sizeI=length(setI)
  x=as.matrix(x)
  beta = as.matrix(beta)
  
  L0=L= l2NormSq( y - matmultC(x,beta))/2/n

  SacBack=SacFort=c()
  for (j in 1:p) {
    xj=x[,j]
    temp = (t(xj)%*%xj)[1,1]
    SacBack=c(SacBack,temp/2/n*beta[j]^2)
    SacFort=c(SacFort,temp/2/n*(d[j]/temp*n)^2)
  }


  for (k in 1:kmax) {
    setAk=c()
    for (j in setA) {
      if( sum( SacBack[j]>=SacBack[setA] )<=k ){
        setAk=c(setAk,j)
      }
    }
    setIk=c()
    for (j in setI) {
      if( sum( SacFort[j]<=SacFort[setI] )<=k ){
        setIk=c(setIk,j)
      }
    }
    
    setAk_wav= union(setdiff(setA,setAk),setIk)
    setIk_wav= union(setdiff(setI,setIk),setAk)

    beta_wav=as.matrix(rep(0,p))
    xak_wav=as.matrix(x[,setAk_wav])
    beta_wav[setAk_wav]= matmultC(matmultC(solve( matmultC(t(xak_wav),xak_wav) ),t(xak_wav)),y)
    beta_wav[setIk_wav]=0
    d_wav= matmultC(t(x), ( y - matmultC(x, beta_wav)  ) )/n
    Lnbeta_wav = l2NormSq( y - matmultC(x, beta_wav)  )/2/n

    if( L>Lnbeta_wav ){
      beta_hat=beta_wav
      d_hat=d_wav
      setA_hat=setAk_wav
      setI_hat=setIk_wav
      
      L=Lnbeta_wav
    }
  }
  
  if( L0-L<TAWs ){
    beta_hat=beta
    d_hat=d
    setA_hat=setA
    setI_hat=setI    
  }
  
  return( list( "bet"=beta_hat[,1],"d"=d_hat,"A"=setA_hat,"I"=setI_hat  ,"L"=L ) )
}


## ----eval=TRUE----------------------------------------------------------------
BessFix = function(s,x,y,kmax,TAWs,   n,p,corrsY){
  setA = c()
  for (j in 1:p) {
    if( sum( corrsY[j]<=corrsY ) <= s   ){
      setA = c(setA,j)
    }
  }
  setI = setdiff(c(1:p),setA)

  x = as.matrix(x)
  beta=d=as.matrix(rep(0,p))
  xA = as.matrix(x[,setA])
  xI = as.matrix(x[,setI])
  
  beta[setA] = (matmultC(matmultC((solve(  matmultC(t(xA),xA))  ) , t(xA)) , y ) )[1,1]
  d[setI] = ( matmultC( t(xI),(y-matmultC(x,beta)) ) )[1,1] /n
  
  while (1) {
    news = Splicing(beta=beta,d=d,setA = setA,setI = setI,kmax=kmax,TAWs = TAWs,   n,p,x,y)

    beta1 = news[[1]]
    d1=news[[2]]
    setA1=news[[3]]
    setI1=news[[4]]
    L1=news[[5]]


    if( setequal(setA1,setA) & setequal(setI1,setI) ){
      break
    }
    else{
      beta = beta1
      d = d1
      setA = setA1
      setI = setI1
    }
  }
  
  return(  list("bet"=beta1,"d"=d1,"A"=setA1,"I"=setI1,  "L"=L1) )
}

myABESS = function(x,y,n,p, 
                   smax=min(p,floor( n/log(p)/log(log(n)) ) )
                   ){
  corrsY=vector(length = p)
  for (i in 1:p) {
    xi=x[,i]
    corrsY[i] = abs( crossprod(xi,y)[1,1] )/sqrt( crossprod(xi,xi)[1,1] )
  }

  
  data = vector(mode="list",length = smax)
  SIC = Ls = vector(length = smax)
  for (s in 1:smax) {
    taws=1e-2*s*log(p)*log(log(n))/n
    kmax = s
    data[[s]] = BessFix(s = s, x = x, y = y, kmax = kmax, TAWs = taws, n = n, p = p,corrsY = corrsY)
    Ls[s] = (data[[s]])[[5]]
    SIC[s] = n*log(Ls[s]) + length( (data[[s]])[[3]] )*log(p)*log(log(n))
  }
  
  best_s = which.min(SIC)
  return( list( "result"=data[[ best_s ]], "SIC"=SIC, "BestSize"=best_s  ) )
  
}


## ----eval=TRUE----------------------------------------------------------------
library(mvtnorm)

seed=0
set.seed(seed)
n=60
eps_sig=1   # for epsilon

p = 20
beta= c(16,0,0,0,0,
        18,0,0,14,0,
        rep(0,p-10)
        )

sig=matrix(0, nrow=p,ncol=p)  # matrix, for X the matrix
sig= 0.5^( abs( row(sig)-col(sig) ) )

x=rmvnorm(n,mean=rep(0,p),sigma=sig)
x=scale(x,center = T,scale = F)
y = rnorm(n, mean=x %*% beta, eps_sig)

x=as.matrix(x)
y=as.matrix(y)
syn_sim_data = list( "x"=x,"y"=y,"beta"=beta )
dat = cbind.data.frame("y" = syn_sim_data[["y"]],
                        syn_sim_data[["x"]])



## ----eval=TRUE----------------------------------------------------------------
res=myABESS(x,y,n,p,smax = p-1)
res$result[-2]  # hide d;
print(paste("Best size: ",res$BestSize))


## -----------------------------------------------------------------------------
library(abess)
library(microbenchmark)


ts <- microbenchmark('with my package'=myABESS(x,y,n,p,smax = p-1),
                     'with Official package'=abess(x,y,family = "gaussian"))
summary(ts)[,c(1,3,5,6)]


## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

