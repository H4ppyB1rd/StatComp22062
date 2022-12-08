\name{StatComp22062-package}
\alias{StatComp22062-package}
\alias{StatComp22062}
\docType{package}
\title{
  StatComp homework and an implement on abess algorithm
}
\description{
  A package containing all my solutions to the homework of this term's Statistical Computation and an implement of  the methods of the Adaptive Best-Subset Selection (ABESS) algorithm according to its thesis, A polynomial algorithm for best-subset selection problem. 
}
\details{
  The myABESS function in this package can be used to implement an Adaptive Best-Subset Selection to find the best subset in certain cases given the proper predictors and responses. BessFix is a step contained in its implement of abess algorithm, while Splicing is a step contained in BessFix. And l2NormSq is to get square of l2 norm of a vector while matmultC helps the calculating run faster. 
}
\author{
Qi Yi.

Maintainer: Qi Yi <lizardplanet@mail.ustc.edu.cn>
}
\references{
  A polynomial algorithm for best-subset selection problem.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
      library(mvtnorm)
      set.seed(0)
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
      res=myABESS(x,y,n,p,smax = p-1)
  }
}
