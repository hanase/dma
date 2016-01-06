\name{dma}
\alias{dma}
\alias{dma.default}
\alias{print.dma}
\alias{plot.dma}
\alias{coef.dma}
\alias{makf4}
\alias{model.update3}
\alias{rm.Kalman}
\title{Dynamic model averaging for continuous outcomes}
\description{
Implemtent dynamic model averaging for continuous outcomes as described in 
Raftery, A.E., Karny, M., and Ettler, P. (2010). Online Prediction Under Model 
Uncertainty Via Dynamic Model Averaging: Application to a Cold Rolling Mill. Technometrics 52:52-66. Along with the values
described below, plot() creates a plot of the posterior model probabilities over time and 
model-averaged fitted values and print() returns model matrix and posterior
model probabilities.  
There are TT time points, K models, and d total covariates.}
\usage{
dma(x, y, models.which, lambda=0.99, gamma=0.99, 
 eps=.001/nrow(models.which), delay=1, initialperiod=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
TTxd matrix of system inputs
}
  \item{y}{
TT-vector of system outputs
}
  \item{models.which}{
Kxd matrix, with 1 row per model and 1 col per variable
		indicating whether that variable is in the model
		(the state theta is of dim (model.dim+1); the extra 1 for the intercept)
}
  \item{lambda}{
	parameter forgetting factor
}
  \item{gamma}{
	flatterning parameter for model updating
}
  \item{eps}{
regularization parameter for regularizing posterior model
		 model probabilities away from zero
}
  \item{delay}{
When y_t is controlled, only y_{t-delay-1} and before 
		 are available. This is determined by the machine.
		Note that delay as defined here corresponds to (k-1)
		 in the Ettler et al (2007, MixSim) paper.
		Thus we use the default delay=24, corresponding to k=25.
}
  \item{initialperiod}{
length of initial period. Performance is summarized with and
		 without the first initialperiod samples.
}
}

\value{
 \item{yhat.bymodel}{TTxK matrix whose tk element gives yhat for yt for model k}
 \item{yhat.ma}{TT vector whose t element gives the model-averaged yhat for yt}
 \item{pmp}{TTxK matrix whose tk element is the post prob of model k at t}
 \item{thetahat.ma}{TTx(nvar+1) matrix whose tk element is the model-averaged 
		 estimate of theta_{j-1} at t}
 \item{Vtheta.ma}{TTx(nvar+1) matrix whose tk element is the model-averaged
		 variance of thetahat_{j-1} at t}
 \item{mse.bymodel}{MSE for each model}
 \item{mse.ma}{MSE of model-averaged prediction}
 \item{mseinitialperiod.bymodel}{MSE for each model excluding the first initialperiod samples}
 \item{mseinitialperiod.ma}{MSE of model averaging excluding the first initialperiod samples}
 \item{model.forget}{forgetting factor for the model switching matrix}
}
\references{
Raftery, A.E., Karny, M., and Ettler, P. (2010). Online Prediction Under Model 
Uncertainty Via Dynamic Model Averaging: Application to a Cold Rolling Mill. Technometrics 52:52-66.}
\author{
Adrian Raftery, Tyler H. McCormick}
\examples{
#simulate some data to test
#first, static coefficients
coef<-c(1.8,3.4,-2,3,-2.8,3)
coefmat<-cbind(rep(coef[1],200),rep(coef[2],200),
            rep(coef[3],200),rep(coef[4],200),
            rep(coef[5],200),rep(coef[6],200))
#then, dynamic ones
coefmat<-cbind(coefmat,seq(1,2.45,length.out=nrow(coefmat)),
            seq(-.75,-2.75,length.out=nrow(coefmat)),
            c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)))
npar<-ncol(coefmat)-1
dat<-matrix(rnorm(200*(npar),0,1),200,(npar))
ydat<-rowSums((cbind(rep(1,nrow(dat)),dat))[1:100,]*coefmat[1:100,])
ydat<-c(ydat,rowSums((cbind(rep(1,nrow(dat)),dat)*coefmat)[-c(1:100),c(6:9)]))
mmat<-matrix(c(c(1,0,1,0,0,rep(1,(npar-7)),0,0),
            c(rep(0,(npar-4)),rep(1,4)),rep(1,npar)),3,npar,byrow=TRUE)
dma.test<-dma(dat,ydat,mmat,lambda=.99,gamma=.99,initialperiod=20)
plot(dma.test)
}
