#' @title abess
#' @description the function for adaptive best-subset selection.
#' @param x the input x, predictors of observations
#' @param y the input y, responses of observations
#' @param n an integer representing number of observations
#' @param p an integer representing dimensions of regression coefficient vector
#' @param smax the maximum support size that can be allowed, default set to floor( n/log(p)/log(log(n)) )
#' @examples
#' \dontrun{
#' library(mvtnorm)
#' set.seed(123)
#' n=60
#' eps_sig=1
#' p = 20
#' beta= c(0,0,10, rep(0,p-3))
#' sig=matrix(0, nrow=p,ncol=p)
#' sig= 0.8^( abs( row(sig)-col(sig) ) )
#' 
#' x=rmvnorm(n,mean=rep(0,p),sigma=sig)
#' x=scale(x,center = T,scale = F)
#' y = rnorm(n, mean=x %*% beta, eps_sig)
#' 
#' x=as.matrix(x)
#' y=as.matrix(y)
#' syn_sim_data = list( "x"=x,"y"=y,"beta"=beta )
#' dat = cbind.data.frame("y" = syn_sim_data[["y"]],syn_sim_data[["x"]])
#' myABESS(x,y,n,p,smax = p-1)
#' }
#' @export
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