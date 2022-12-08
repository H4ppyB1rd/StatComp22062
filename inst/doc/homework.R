## -----------------------------------------------------------------------------
mytext = "hello R!"
print(mytext)


## -----------------------------------------------------------------------------
x <- rnorm(100)
y <- rnorm(100)
plot(x,y,col = "darkgreen")


## -----------------------------------------------------------------------------
library(kableExtra)
dt <- state.x77[1:10,]
kable(dt)


## -----------------------------------------------------------------------------
# clear all
rm(list = ls())

## -----------------------------------------------------------------------------
invtr_pareto = function(a,b,n){      #inverse transformation for pareto(a,b)
  number = n
  u = runif(number)     #uniform(0,1)
  x = (1/(1-u))^(1/a)*b      #inversed function for F()
  hist(x,prob = TRUE,main = expression( "pareto distribution" ),col = "green" )
  y = seq(0, 200, 0.1)
  lines(y,b^a*a*y^(-a-1),col="blue")   #density corresponding to F()
}

invtr_pareto(a=2,b=2,n=5000)


## -----------------------------------------------------------------------------
beta_pdf = function(x,a,b){   #pdf for beta distribution
  x^(a-1)*(1-x)^(b-1)/beta(a,b)
}

ar_beta = function(a,b,n){
  k = 0   #counter for acception
  j = 0   #iterations
  x = numeric(0)   #initiation
  
  max_p = (a-1)/(a+b-2)   # max point, if it exists
  max_v = beta_pdf(max_p,a,b)
  
  while(k<n){
    u = runif(1)
    j=j+1
    y = runif(1)    # random variate from g()
    if( beta_pdf(y,a,b)/max_v>u ){   # we accept it
      k = k+1
      x[k] = y
    }
  }
  
  #make a plot:
  plot(density(x), main = bquote("AR method vs True"))
  this.range <- seq(0, 1, 0.01)
  f <- function(x) beta_pdf(x,a,b)
  lines(this.range,f(this.range), col="red")
  legend("topleft", legend = c("A-R", "True"), col=c("black", "red"), lty=1)

}

ar_beta(3,2,1000)

## -----------------------------------------------------------------------------
exp_gam_mix = function(n,r,beta){   # my generator
  lambda = rgamma(n,r,beta)
  x = rexp(n,lambda)
  x
}

x = exp_gam_mix(n=1e3,r=4,beta=2)
# observations generated.

## -----------------------------------------------------------------------------
bet=2;r=4
x = exp_gam_mix(n=1e3,r=4,bet=2)
plot(ecdf(x),xlim=c(-0.1,9),col="orange",lty=1,main="comparison")
y = seq(0,8,.01)
lines(y,1-(bet/(bet+y))^r,col="green",lty=2)
legend("bottomright", legend = c("empirical", "pareto"), col=c("orange", "green"), lty=c(1,2))


## -----------------------------------------------------------------------------
# clear all
rm(list = ls())

## -----------------------------------------------------------------------------
my_qsort = function(x){   # function for recursively quick sort.
  num = length(x)
  if(num==0||num==1){
    return(x)
  }
  else{
    a = x[1]
    y = x[-1]
    left = y[y<a]
    right = y[y>=a]
    return(c(my_qsort(left),a,my_qsort(right)))}
}

numbers = c(1e4,2*1e4,4*1e4,6*1e4,8*1e4)    # how many numbers
an = numeric(5)       # averged runtime for length n
i = 1
print("Here is the result of runtimes:")
for(number in numbers){
  runtime = numeric(100)
  for(j in 1:100){
    mysample = sample(1:number)     # get randomly permuted numbers
    runtime[j] = system.time(my_qsort(mysample))[1]
  }
  an[i]=mean(runtime)
  print(paste("  runtime for numbers of length",number,"is",an[i]))
  i=i+1
}

## -----------------------------------------------------------------------------
tn = numbers*log(numbers)
mymodel = lm(an ~ tn)
plot(tn,an)
abline(mymodel,col = "orange")    # add regression line

## -----------------------------------------------------------------------------
antitheticV_mc = function(nsim = 1e5) {    # antithetic variate approach,nsim is the number of simulation.
  u = runif(nsim/2)
  v = 1 - u
  u = c(u, v)
  g = exp(u)
  mean(g)
}
simple_mc = function(nsim = 1e5) {    # simple monte carlo
  u = runif(nsim)
  g <- exp(u)
  mean(g)
}

n_sim = 1e4
MC1 <- MC2 <- numeric(n_sim)
for (i in 1:n_sim) {
  MC1[i] <- simple_mc(nsim=n_sim)
  MC2[i] <- antitheticV_mc(nsim=n_sim)
}

(var(MC1) - var(MC2))/var(MC1)

## -----------------------------------------------------------------------------
# clear all
rm(list = ls())

## -----------------------------------------------------------------------------
g = function (x) {
  x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
}
xseq = seq(1,10,0.01)
yg = g(xseq)
y1 = exp(-xseq)
y2 = dnorm(xseq)

# to draw the density functions:
ymax = max(c(yg,y1,y2))
plot(xseq,yg,type="l",xlim = c(1,10),ylim=c(0,ymax),main="functions comparison",xlab="x",ylab="y")
lines(xseq,y1,col="orange")
lines(xseq,y2,col="green")
legend("topright",legend = c("g","f1","f2"),col=c("black","orange","green"),lty=1)

# to draw the ratios:
ymax = max(c(yg/y1,yg/y2))
plot(xseq,yg/y1,type="l",col="orange",xlim = c(1,10),ylim=c(0,ymax),main="ratios comparison",xlab="x",ylab="g/f")
lines(xseq,yg/y2,col="green")
legend("topleft",legend = c("f1","f2"),col=c("orange","green"),lty=1)


## -----------------------------------------------------------------------------
m=10000
theta_hat = se = numeric(2)

x1=rexp(m,1)    # importance sampling for f1
fg1 = g(x1)/exp(-x1)
theta_hat[1] = mean(fg1);se[1] = sd(fg1)

x2=rnorm(m)   # importance sampling for f2
fg2 = g(x2)/dnorm(x2)
theta_hat[2] = mean(fg2);se[2] = sd(fg2)

t = rbind(theta_hat,se)
colnames(t) = c("f1","f2")
print(t)

## -----------------------------------------------------------------------------
M = 1e5;k = 5
r = M/k #replicates per stratum
N = 50
T2_k = numeric(k)
stra_est = numeric(N)
newfunc = function(x)  (1-exp(-1))/(1+x^2)*(x<1)*(x>0)
for (i in 1:N) {
  for(j in 1:k) {
    fk1 = -log( 1-(1-exp(-1))*(j-1)/k )
    fk2 = -log( 1-(1-exp(-1))*j/k )
    u = runif(M*(fk2-fk1))    # so there are M in total for one N.
    X = -log(exp(-fk1)-(1-exp(-1))*u/k )    # random numbers of 5*f3 in this subinterval
    T2_k[j] = mean( (1-exp(-1))/((1+X^2)*5) )
  }
  stra_est[i] = sum(T2_k)
}

mean(stra_est)
sd(stra_est)


## -----------------------------------------------------------------------------
# clear all
rm(list = ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

get_data = function(n,mu,sig){     # function for data generation
  x = rlnorm(n,mu,sig)
  return( x )
}

make_ucl = function(m,n,mu,sig,alpha){    # function for making confidence interval
  ucl = replicate(m,
                  expr = {
                    x = rlnorm(n,meanlog = mu,sdlog = sig)
                    y = log(x)
                    abs( sqrt(n)*(mean(y)-mu)/sd(y) )
                  }   )
  return( mean(ucl<qt(1-alpha/2,n-1)) )
}


## -----------------------------------------------------------------------------
n=50
mu=0
sig=2
m=1e4
alpha=0.05
set.seed(0)
res = make_ucl(m,n,mu,sig,alpha)    # we perform simulation here.
print(paste("The obtained empirical estimate of the confidence level is",res,"this is a good result."))

rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

count5test = function(x, y){    # function for count-five-test, from textbook
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

ftest = function(x,y,alpha){      # function for F-test
  ftp = var.test(x,y)$p.value
  return (as.integer(ftp<=alpha))
}

get_data = function(mu1,sig1,n1,mu2,sig2,n2){     # function for data generation
  x = rnorm(n1,mu1,sig1)
  y = rnorm(n2,mu2,sig2)
  return( list('x'=x,'y'=y) )
}

data_output = function(f_res,c5t_res,size_list){     # make a dataframe for result output
  mydf = as.data.frame( matrix( data = c(c5t_res,f_res),
                                nrow = length(size_list), ncol = 2,
                                byrow = FALSE
                                ) )
  colnames(mydf) = c("count5test","f_test")
  rownames(mydf) = paste0("n =",size_list)
  mydf
}



## -----------------------------------------------------------------------------
set.seed(0)
m=1e4
alpha0 = 0.055
n_list = c(10,50,100,500,1000) # from small to middle to large
ft_result = c5t_result = numeric(  length(n_list)  )

i=1
for(n in n_list){       # compute results for f_test and c5test for each sample size
  sims = replicate(m,
                   expr = {
                     mydata = get_data(mu1=0,sig1=1,n1=n,mu2=0,sig2=1.5,n2=n)
                     x = mydata$x; y = mydata$y
                     c('c5t'=count5test(x,y),'ft'=ftest(x,y,alpha = alpha0))
                   }  )
  results = rowMeans(sims)
  c5t_result[i] = results[[1]]; ft_result[i] = results[[2]]
  i = i+1
}

print(data_output(ft_result,c5t_result,n_list))
rm(list=ls())

## -----------------------------------------------------------------------------
# clear all
rm(list = ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

bootstrap_result_get = function(x,R){   # get bootstrap
  lambda_star = numeric(R)
  for(r in 1:R){
    xstar = sample(x,replace = T)
    lambda_star[r] = 1/mean(xstar)
  }
  return(lambda_star)
}

bootstrap_result_show = function(res,x){  # output
  est.samp = 1/mean(x)
  print(c(bias.boot = mean(res)-est.samp,
        se.boot = sd(res)))
}

## -----------------------------------------------------------------------------
# set.seed(2000)
library(boot)
mydata = aircondit$hours

R=1e4
bootstrapRes = bootstrap_result_get(mydata,R)
bootstrap_result_show(bootstrapRes,mydata)

## -----------------------------------------------------------------------------
x=mydata
f = function(data,inds){
  d = data[inds]
  return(1/mean(d))
}
# set.seed(2000)
bootmle=boot(x,f,R=1e4)
bootmle

## -----------------------------------------------------------------------------
f = function(data,inds){
  d = data[inds]
  return(mean(d))
}
# set.seed(2000)
boot.obj = boot(mydata,statistic = f,R=1e4)
print(boot.ci(boot.obj,type = c("norm","basic", "perc","bca")))

## -----------------------------------------------------------------------------
rm(list=ls())
library(boot)
set.seed(1)

n = 1e1
m = 1e4
mu=0;sig=1

boot.mean = function(x,i){  mean(x[i]) }
ci.norm=ci.basic=ci.perc=matrix(NA,m,2)
for(i in 1:m){
  mydata = rnorm(n,mu,sig)
  de <- boot(data=mydata,statistic=boot.mean, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}


## -----------------------------------------------------------------------------
myprint = function(ci.norm,ci.basic,ci.perc){
  print("for normal interval:")
  print(paste("  empirical coverage rates:",mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),";",
               "miss on left:",mean(ci.norm[,1]>mu),";",
               "miss on right:",mean(ci.norm[,2]<mu),
               sep=" "))
  
  print("for basic interval:")
  print(paste("  empirical coverage rates:",mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),";",
               "miss on left:",mean(ci.basic[,1]>mu),";",
               "miss on right:",mean(ci.basic[,2]<mu),
               sep=" "))
  
  print("for percent interval:")
  print(paste("  empirical coverage rates:",mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),";",
               "miss on left:",mean(ci.perc[,1]>mu),";",
               "miss on right:",mean(ci.perc[,2]<mu),
               sep=" "))
}
myprint(ci.norm,ci.basic,ci.perc)

mat = matrix(c(mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),mean(ci.norm[,1]>mu),mean(ci.norm[,2]<mu),
               mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),mean(ci.basic[,1]>mu),mean(ci.basic[,2]<mu),
               mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),mean(ci.perc[,1]>mu),mean(ci.perc[,2]<mu) ),
             3,3,
             byrow = TRUE)
rownames(mat) = c("for normal interval:","for basic interval:","for percent interval:")
colnames(mat) = c("empirical coverage rates","miss on left","miss on right")
library(kableExtra)
kable(as.data.frame(mat))


## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

get_theta = function(mydata){    # how we can get the parameter theta
  lam = eigen(cov(mydata))$values
  my_theta = lam[1]/sum(lam)
  
  return(my_theta)
}

get_jack = function(mytable){   # data analysis: obtain the jackknife estimates of bias and sd
  theta_Hat = get_theta(mytable)
  theta_Jack =theta_Jack2 = numeric(n)
  for(i in 1:n){
    theta_Jack[i] = get_theta(mytable[-i,])
  }

  bias_Jack = (n-1)*(mean(theta_Jack)-theta_Hat)
  se_Jack = sqrt((n-1)*mean((theta_Jack-mean(theta_Jack))^2)) 
  return(c(theta.hat=theta_Hat,bias.jack=bias_Jack,se.jack=se_Jack))
}


## -----------------------------------------------------------------------------
# to load the required data:
library(bootstrap)
n = nrow(scor)

# jackknife result output:
print(get_jack(scor))

# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())


# this function is to fit models on leave-two-out samples on n-fold cross validation and compare errors:
get_modErrs = function(magnetic,chemical){
  n <- length(magnetic)
  e1 <- e2 <- e3 <- e4 <- c() # errors for each model,initialized as empty
  for (k1 in 1:(n-1)) {
    for(k2 in (k1+1):n){
    y = magnetic[-c(k1,k2)]
    x = chemical[-c(k1,k2)]
    
    Model1 = lm(y ~ x)
    yhat11 = Model1$coef[1] + Model1$coef[2] * chemical[k1]
    yhat12 = Model1$coef[1] + Model1$coef[2] * chemical[k2]
    e1 = c(e1,c(magnetic[k1]-yhat11, magnetic[k2]-yhat12))
  
    Model2 = lm(y ~ x + I(x^2))
    yhat21 = Model2$coef[1] + Model2$coef[2] * chemical[k1] + Model2$coef[3] * chemical[k1]^2
    yhat22 = Model2$coef[1] + Model2$coef[2] * chemical[k2] + Model2$coef[3] * chemical[k2]^2
    e2 = c(e2,c(magnetic[k1] - yhat21, magnetic[k2]-yhat22))
    
    Model3 = lm(log(y) ~ x)
    logyhat31 = Model3$coef[1] + Model3$coef[2] * chemical[k1]
    yhat31 = exp(logyhat31)
    logyhat32 = Model3$coef[1] + Model3$coef[2] * chemical[k2]
    yhat32 = exp(logyhat32)
    e3 = c(e3,magnetic[k1] - yhat31,magnetic[k2] - yhat32)
    
    Model4 = lm(log(y) ~ log(x))
    logyhat41 = Model4$coef[1] + Model4$coef[2] * log(chemical[k1])
    yhat41 = exp(logyhat41)
    logyhat42 = Model4$coef[1] + Model4$coef[2] * log(chemical[k2])
    yhat42 = exp(logyhat42)
    e4 = c(e4, magnetic[k1] - yhat41,magnetic[k2] - yhat42)
    }
  }
  
  return( list(errs=c(Model1.error=mean(e1^2), Model2.error=mean(e2^2),Model3.error=mean(e3^2), Model4.error=mean(e4^2))) )
}


## -----------------------------------------------------------------------------
library(DAAG)
results = get_modErrs(ironslag$magnetic,ironslag$chemical)
print( results$errs )


## -----------------------------------------------------------------------------
lm(formula = magnetic ~ chemical + I(chemical^2) ,data = ironslag)

# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all at the start:
rm(list=ls())

get_data = function(){
  attach(chickwts)
  my_x = as.vector(weight[feed == "linseed"])
  my_y = as.vector(weight[feed == "sunflower"])
  detach(chickwts)
  return(list(x=my_x,y=my_y))
}

compare = function(x,y){  #data_analysis
  set.seed(12345)
  mydata = get_data()
  x = mydata$x; y = mydata$y
  
  n = length(x)  # x and y are equal in length
  R = 500
  z =c(x,y)
  reps = numeric(R)
  K = length(z)
    
  cor0 = cor(x,y,method = "spearman")
  for (i in 1:R) {
    k = sample(K,size=n,replace = F)
    x1=z[k];y1=z[-k]
    reps[i] = cor.test(x1,y1)$statistic
  }

  return( c(achieved_significance_level_of_permutation=mean(c(cor0,reps)>=cor0),
            cor_test_report_p_val=cor.test(x=x,y=y)$p.value) )
}


## -----------------------------------------------------------------------------
set.seed(12345)
mydata = get_data()
x = mydata$x; y = mydata$y
compare(x,y)

# remove all at the end:
rm(list=ls())

## ----fig.height=8, fig.width=8------------------------------------------------
# remove all before start:
rm(list=ls())

# pdf for standard laplace distribution
standard_d_laplace = function(x){
  return( 1/2*exp(-abs(x)) )
}

# function for a random Walk Metropolis
randomWalk.Metropolis = function(sigma, x0, N) {
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (standard_d_laplace(y) / standard_d_laplace(x[i-1])))
    x[i] = y else {
      x[i] = x[i-1]
      k = k + 1
    }
  }
  return(list(x=x, k=k))
}

# function for gelman-rubin method
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

mychains <- function(sigma, N, X1) {
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)

  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
      if (u[i] <= (standard_d_laplace(y) / standard_d_laplace(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
      }
  }
  return(x)
}

#gelman_rubin for monitoring
gelman_rubin_monitor = function(sigma,seed,drawplot=F){

  set.seed(seed)
  
  k <- 4 #number of chains to generate
  n <- 15000 #length of chains
  b <- 1000 #burn-in length
  
  #choose overdispersed initial values
  x0 <- rep(25,4)
  #generate the chains
  X <- matrix(0, nrow=k, ncol=n)
  for (i in 1:k){
    X[i, ] <- mychains(sigma, n, x0[i])
  }
  #compute diagnostic statistics
  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
  
  if(drawplot==T){
    #plot psi for the four chains
    par(mfrow=c(2,2))
    for (i in 1:k)
    plot(psi[i, (b+1):n], type="l",
    xlab=i, ylab=bquote(psi))
    par(mfrow=c(1,1)) #restore default
    #plot the sequence of R-hat statistics
    rhat <- rep(0, n)
    for (j in (b+1):n)
    rhat[j] <- Gelman.Rubin(psi[,1:j])
    plot(rhat[(b+1):n], type="l", xlab="", ylab="R",main=paste("sigma=",sigma))
    abline(h=1.2, lty=2,col="red")
  }
  
  return(psi)
}

## -----------------------------------------------------------------------------
set.seed(2022)
N = 15000
sigma = c(.05, .5, 2, 16)
x0 = 25
rw1 = randomWalk.Metropolis(sigma[1], x0, N)
rw2 = randomWalk.Metropolis(sigma[2], x0, N)
rw3 = randomWalk.Metropolis(sigma[3], x0, N)
rw4 = randomWalk.Metropolis(sigma[4], x0, N)

#number of candidate points rejected
Rejected = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
Accepted = matrix(1-Rejected/N,nrow=1,ncol=4)
rownames(Accepted) = "Accept rates"
colnames(Accepted) = paste("sigma",sigma)
print(Accepted)

## ----fig.height=8, fig.width=8------------------------------------------------
#plot
par(mfrow=c(2,2)) #display 4 graphs together
    rw = cbind(rw1$x, rw2$x, rw3$x, rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]) )
    }
par(mfrow = c(1,1))  # reset default


## -----------------------------------------------------------------------------
seed=0
monitor1=gelman_rubin_monitor(0.05,seed)
monitor2=gelman_rubin_monitor(0.5,seed)
monitor3=gelman_rubin_monitor(2,seed)
monitor4=gelman_rubin_monitor(16,seed)



## -----------------------------------------------------------------------------

# remove all at the end:
rm(list=ls())

## ----fig.height=8, fig.width=8------------------------------------------------
# remove all before start:
rm(list=ls())


#function for bivariate sample chain generation:
my_chain_generator = function(N=5000,burn=1000,rho=0.9,mu1=0,mu2=0,sigma1=1,sigma2=1,seed=0){
  set.seed(seed)
  
  s1 = sqrt(1-rho^2)*sigma1
  s2 = sqrt(1-rho^2)*sigma2
  
  X = matrix(0, N, 2) # to store the chain
  X[1, ] = c(mu1, mu2) # to initiate
  for (i in 2:N) {
    x2 = X[i-1, 2]
    m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] = rnorm(1, m1, s1)
    x1 = X[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] = rnorm(1, m2, s2)
  }
  b = burn + 1
  x = X[b:N, ]

  return(list(x=x, X=X))
}

# function for regression model fit and checking residuals for normality and constant variance
model_fit = function(x){
  myfit = lm(x[,2]~x[,1])
  print(myfit)
  
  par(mfrow=c(2,2),mar=rep(0,4))
  plot(myfit,cex=0.3) 
}


Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}


## -----------------------------------------------------------------------------
N=5000
mysample = my_chain_generator(N=N)
x = mysample$x
X = mysample$X

# to plot the generated sample after discarding a suitable burn-in sample:
plot(x, main="", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(x[,2]))


## -----------------------------------------------------------------------------
model_fit(x=x)

## ----fig.height=8, fig.width=8------------------------------------------------
  set.seed(0)
  
  k <- 2 #number of chains to generate
  n <- N #length of chains
  b <- 1000 #burn-in length
  
  X=t(X)
  #compute diagnostic statistics
  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
  #plot psi for the four chains
  par(mfrow=c(2,2),mar=rep(0,4))
  for (i in 1:k)
  plot(psi[i, (b+1):n], type="l",  xlab=i, ylab=bquote(psi))
  par(mfrow=c(1,1)) #restore default
  #plot the sequence of R-hat statistics
  rhat <- rep(0, n)
  for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
  abline(h=1.2, lty=2,col="red")

## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

# function for generation of X, whose distribution is based on our assignment, here is N(0,3)
getX = function(N,myseed=NULL){  # N is length of X
  return( rnorm(N,0,3) )
}

# function for generation of M
getM = function(x,am,alp,myseed=NULL){
  set.seed(myseed)
  N = length(x)
  return( am+alp*x+rnorm(N,0,1) )
}

# function for generation of Y
getY = function(m,x,ay,bet,gam,myseed=NULL){
  set.seed(myseed)
  N = length(m)
  return( ay+bet*m+gam*x+rnorm(N,0,1) )
}


## -----------------------------------------------------------------------------
library(mediation)

# function for test of mediation effects with package mediation:
mediate_test = function(myalpha=0,mybeta=0,mygamma=1,myseed=NULL){
  x = getX(N = 1e3,myseed = myseed)
  m = getM(x = x,am = 0,alp = myalpha,myseed = myseed)
  y = getY(m = m,x = x,ay = 0,bet = mybeta,gam = mygamma,myseed = myseed)
  dat = matrix(c(x,m,y),ncol=3)
  colnames(dat) = c("x","m","y")
  
  model_xm = lm(m~x)
  model_xy = lm(y~m+x)
  model_MedEffect = mediate(model_xm,model_xy,treat = 'x',mediator = 'm', sims = 1000,boot = T)
  
  print(summary(lm(y~m)))
  print(summary(model_MedEffect))
}


## -----------------------------------------------------------------------------
mediate_test(myalpha = 0,mybeta = 0)

## -----------------------------------------------------------------------------
mediate_test(myalpha = 0,mybeta = 1,myseed = 0)


## -----------------------------------------------------------------------------
mediate_test(myalpha = 1,mybeta = 0,myseed = 0)

# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())


# The following are the functions that are needed for this exercise:  

# expit function
expit = function(x){
  1/(1 + exp(-x))
}

# find alpha by solving the equation
find_alpha = function(N,b1,b2,b3,f0,seed=NULL){
  set.seed(seed)
  x1 = rpois(N,1)
  x2 = rexp(N,1)
  x3 = rbinom(N,1,0.5)
  
  my_pd = function(alpha,b1,b2,b3){
    temp = exp(-alpha-b1*x1-b2*x2-b3*x3)
    temp = 1/(1+temp)
    mean(temp) - f0
  }
  mysolution = uniroot(my_pd,c(-20,10),b1,b2,b3)
  myresult = (round(unlist(mysolution),5)[1:3])
  
  return(list(result=myresult, root=myresult[[1]]))
}


## -----------------------------------------------------------------------------
set.seed(0)
N = 1e6; b1 = 0; b2 = 1; b3 = -1
f0_list = c(1e-1,1e-2,1e-3,1e-4)
alpha_list = numeric(length = length(f0_list))


i=1
for (f0 in f0_list) {
  alpha_list[i] = find_alpha(N=N,b1=b1,b2=b2,b3=b3,f0=f0,seed=NULL)$root
  i=i+1
}

# the result of trying this function is as follows:
print("result output:")
names(alpha_list) = paste0(" f0=",f0_list)
print(alpha_list)


## -----------------------------------------------------------------------------
# because I think 4 points are too few for a scatter plot so I added some more points
f0_list = c(c(1e-4,1e-3,5*1e-3,1e-2,5*1e-2),seq(1e-1,0.99,0.05))
alpha_list = numeric(length = length(f0_list))
i=1
for (f0 in f0_list) {
  alpha_list[i] = find_alpha(N=N,b1=b1,b2=b2,b3=b3,f0=f0,seed=NULL)$root
  i=i+1
}

# the scatter plot is:
plot(f0_list,alpha_list,xlab = "f0",ylab = "alpha",main = "f0 vs alpha",col="orange")

# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# remove all before start:
rm(list=ls())

# function for the deviation of log-likelihood
log_likeh = function(lam,observs){
  myres = 0
  n = nrow(observs)
  for (i in 1:n) {
    ui = observs[i,1]
    vi = observs[i,2]
    myres = myres + log(exp(-lam*ui)-exp(-lam*vi))
  }
  return(myres)
}


# function for the deviation of log-likelihood, its root is MLE:
log_likeh_dev = function(lam,observs){
  myres = 0
  n = nrow(observs)
  for(i in 1:n){
    ui = observs[i,1]
    vi = observs[i,2]
    fac1 = exp(-lam*ui)
    fac2 = exp(-lam*vi)
    myres = myres+(-ui*fac1+vi*fac2)/(fac1-fac2)
  }
  return(myres)
}

# plot the above function to see vaguely where the MLE is and return the root:
getMLE_method1= function(){
  lamlist = seq(0.01,0.25,0.001)
  reslist = numeric(length = length(lamlist))
  for (i in 1:length(lamlist)) {
    reslist[i] = log_likeh_dev(lamlist[i],myobservs)
  }
  # take a peek:
  plot(lamlist,reslist,col="orange",
       xlab = "lambda",ylab = "loglikelihood_dev",
       main = "where the MLE vaguely is"
       )
  abline(h=0,col="red")
  
  res <- uniroot(log_likeh_dev,observs=myobservs,lower=0.03,upper=0.2)
  return(res$root)
}

## -----------------------------------------------------------------------------
# the observations of data:
myobservs = matrix(
  c(11,12, 8,9, 27,28, 13,14, 16,17, 0,1, 23,24, 10,11, 24,25, 2,3),
  ncol=2,byrow = T)

# get mle via method 1
mle1 = getMLE_method1()
print(paste("For the first method(maximizing likehood), the answer is ",mle1,"."))

## -----------------------------------------------------------------------------
method2 = function(lam0,lam1=5,obs1,obs2,N,eps){ # lam0: initial value; N: max_iterations
  lam2 = lam0

  i=1
  while( abs(lam1-lam2)>=eps ){
    lam1 = lam2
    lam2 = length(obs1)/sum(((lam1*obs1+1)*exp(-lam1*obs1)-(lam1*obs2+1)*exp(-lam1*obs2))/(lam1*(exp(-lam1*obs1)-exp(-lam1*obs2))))
    if(i==N) break
    i<-i+1
  }
  return(lam2)
}


em_sol = method2(lam0=0.5,
                 lam1=5,
                 obs1=myobservs[,1],
                 obs2=myobservs[,2],
                 N=1e5,
                 eps=1e-3)


print(paste("For the first method(maximizing likehood), the answer is ",em_sol,"."))



## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
this_vec=as.vector(c(1:5))
# this shows that x is a vector:
print(is.vector(this_vec))
# this shows return value of dim() when applied to a vector:
print(dim(this_vec))


## -----------------------------------------------------------------------------
x = matrix(1:12,3,4)

print(is.matrix(x))
print(is.array(x))


## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
mydf = as.data.frame(matrix(0,2,3))
# this is a data frame:
print(is.data.frame(mydf))
# and has these following attributes: 
print(attributes(mydf))

## -----------------------------------------------------------------------------
mydf = as.data.frame(matrix(c(1:3,'I Am','Another','Type',4:6,T,F,F,5.5,5.6,5.7),3,5))
print(mydf)
# as shown in this case, according to the convertion order, they are all into characters.
print(as.matrix(mydf))


## -----------------------------------------------------------------------------
# This dataframe has 4 cols and 0 row:
df_0row = data.frame(row.names=1:4)
print(dim(df_0row))

# This dataframe has 0 col and 2 rows:
df_0row = data.frame(col1=integer(length = 0),col2=numeric(length = 0))
print(dim(df_0row))

# This dataframe has 0 col and 0 row:
df_empty = data.frame()
print(dim(df_empty))


## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
# we first create a sample matrix for demonstration:
mymat = matrix(rnorm(12),3,4)
colnames(mymat) = c("col1","col2","col3","col4")
mydf = as.data.frame(mymat)  # my dataframe

# to compute sd of each column in a numeric data frame:
vapply(mydf, sd, numeric(1))


## -----------------------------------------------------------------------------
# we first create another sample matrix for demonstration:
mymat = matrix(rnorm(12),3,4)
colnames(mymat) = c("col1","col2","col3","col4")
mydf = cbind(as.data.frame(mymat),col5=c("Not","Numeric","here!"))  # my dataframe
print(mydf)
print("we can see that this dataframe does not contain numeric columns only!")

# to compute sd of each column in a numeric data frame:
vapply(mydf[vapply(mydf, is.numeric, logical(1))],
       sd, numeric(1))


## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp('../src/gibbsC.cpp')

## -----------------------------------------------------------------------------
# FIRST, WE OBTAIN THE SAMPLE USING PURE R:

# R function for gibbs sampler:
my_gibbs_xy_r = function(mu, sig, rho, init0, N,seed=NULL) {
  # gibbs sampler for data generation in R
  # mu,sig and rho being parameters for distribution,init0 being initials and N being length of chain
  set.seed(seed)
  s = sqrt(1 - rho^2) * sig
  X = Y = numeric(N)
  X[1] = init0[1]; Y[1] = init0[2]
  for (i in 2:N) {
    y = Y[i-1]
    m1 = mu[1] + rho * (y - mu[2]) * sig[1] / sig[2]
    X[i] = rnorm(1, m1, s[1])
    x = X[i]
    m2 = mu[2] + rho * (x - mu[1]) * sig[2] / sig[1]
    Y[i] = rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}
# data gerneration in R:
N = 1e3;b = 1000;rho = 0.9;mu = c(0, 0);sigma = c(1, 1)
my_chainR = my_gibbs_xy_r(mu, sigma, rho, c(0,0), N)


# NOW WE OBTAIN SAMPLE WITH Rcpp:

mychainC=gibbsC(0,0,1,1,0.9,0,0,1000)
mychainC=as.data.frame(mychainC)
colnames(mychainC) = c('x','y')

## -----------------------------------------------------------------------------
# NOW WE MAKE QQPLOT:
# qqplot of x_t vs x_t:
xt_r = my_chainR$X
xt_c = mychainC$x
qqplot(xt_r,xt_c)

# qqplot of y_t vs y_t:
yt_r = my_chainR$Y
yt_c = mychainC$y
qqplot(yt_r,yt_c)

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark('with Pure R'=my_gibbs_xy_r(mu, sigma, rho, c(0,0), N),
                     'with Rcpp'=gibbsC(0,0,1,1,0.9,0,0,100))
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
# remove all at the end:
rm(list=ls())

