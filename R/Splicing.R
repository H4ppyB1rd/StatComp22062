#' @title Splicing
#' @description splicing method for adaptive best subset selection.
#' @param beta the regression coefficient vector that is imported into the splicing method
#' @param d the vector value of t(x)(y-x beta)/n that is imported into the splicing method
#' @param setA a vector representing the set of positions of non zero values of regression coefficient that is imported into the splicing method
#' @param setI a vector representing the set of positions of zero values of regression coefficient that is imported into the splicing method
#' @param kmax an integer representing the maximum splicing size
#' @param TAWs a numeric value representing the threshold chosen by user into the splicing method
#' @param n an integer representing number of observations
#' @param p an integer representing dimensions of regression coefficient vector
#' @param x the predictors of observations
#' @param y the responses of observations
#' @import Rcpp microbenchmark knitr rmarkdown kableExtra boot bootstrap DAAG mediation mvtnorm abess
#' @export
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
