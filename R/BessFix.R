#' @title BessFix
#' @description this function is for best subset selection with support size s given and fixed.
#' @param s an integer, the maximum support size of coefficient we allow, given and fixed
#' @param x the predictors of observations
#' @param y the responses of observations
#' @param kmax an integer representing the maximum splicing size
#' @param TAWs a numeric value representing the threshold chosen by user into the splicing method
#' @param n an integer representing number of observations
#' @param p an integer representing dimensions of regression coefficient vector
#' @param corrsY correlations imported into BessFix in order to help determine the initial value for setA_0
#' @import Rcpp
#' @export
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
