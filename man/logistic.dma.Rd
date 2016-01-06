\name{logistic.dma}
\alias{logistic.dma}
\alias{logistic.dma.default}
\alias{print.logistic.dma}
\alias{plot.logistic.dma}
\alias{dlogr.init}
\alias{dlogr.predict}
\alias{dlogr.step}
\alias{laplace.fn}
\alias{tunemat.fn}
\title{
Dynamic model averaging for binary outcomes
}
\description{
Implemtent dynamic model averaging for continuous outcomes as described in 
McCormick, T.M., Raftery, A.E., Madigan, D. and Burd, R.S. (2011) "Dynamic Logistic Regression and 
Dynamic Model Averaging for Binary Classification." Biometrics, 66:1162-1173.  Along with the values
described below, plot() creates a plot of the posterior model probabilities over time and 
model-averaged fitted values (with smooth curve overlay) and print() returns model matrix and posterior
model probabilities.  There are K candidate
models, T time points, and d total covariates (including the intercept).
}
\usage{
logistic.dma(x, y, models.which, lambda=0.99, alpha=0.99,autotune=TRUE, 
 initmodelprobs=NULL,initialsamp=NULL)
}
\arguments{
  \item{x}{T by (d-1) matrix of observed covariates.  Note that a column of 1's is added
    automatically for the intercept.}
  \item{y}{T vector of binary responses}
\item{models.which}{K by (d-1) matrix defining models.  A 1 indicates a covariate is included
    in a particular model, a 0 if it is excluded.  Model averaging is done over all
    modeld specified in models.which.}
  \item{lambda}{scalar forgetting factor with each model}
  \item{alpha}{scalar forgetting factor for model averaging}
  \item{autotune}{ T/F indicates whether or not the automatic tuning procedure desribed in 
    McCormick et al. should be applied.  Default is true.}
    \item{initmodelprobs}{K vector of starting probabilities for model averaging.  If null (default),
    then use 1/K for each model.}
      \item{initialsamp}{scalar indicating how many observations to use for generating initial 
    values.  If null (default), then use the first 10 percent of observations.}
}

\value{
  \item{x}{T by (d-1) matrix of covariates}
  \item{y}{T by 1 vector of binary responses}
  \item{models.which}{K by (d-1) matrix of candidate models}
  \item{lambda}{scalar, tuning factor within models}
  \item{alpha}{scalar, tuning factor for model averaging}
  \item{autotune}{T/F, indicator of whether or not to use autotuning algorithm}
  \item{alpha.used}{T vector of alpha values used}
  \item{theta}{T vector of alpha values used}
  \item{vartheta}{K by T by d array of dynamic logistic regression variances for each model}
  \item{pmp}{K by T array of posterior model probabilities}
  \item{yhatdma}{T vector of model-averaged predictions}
  \item{yhatmodel}{K by T vector of fitted values for each model}
}
\references{
McCormick, T.M., Raftery, A.E., Madigan, D. and Burd, R.S. (2011) "Dynamic Logistic Regression and 
Dynamic Model Averaging for Binary Classification." Biometrics, 66:1162-1173.}
\author{
Tyler H. McCormick, David Madigan, Adrian Raftery
}
\examples{
#simulate some data to test
#first, static coefficients
coef<-c(.08,-.4,-.1)
coefmat<-cbind(rep(coef[1],200),rep(coef[2],200),rep(coef[3],200))
#then, dynamic ones
coefmat<-cbind(coefmat,seq(1,.45,length.out=nrow(coefmat)),
          seq(-.75,-.15,length.out=nrow(coefmat)),
          c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)))
npar<-ncol(coefmat)-1

#simulate data
dat<-matrix(rnorm(200*(npar),0,1),200,(npar))
ydat<-exp(rowSums((cbind(rep(1,nrow(dat)),dat))[1:100,]*coefmat[1:100,]))/
          (1+exp(rowSums(cbind(rep(1,nrow(dat)),dat)[1:100,]*coefmat[1:100,])))
y<-c(ydat,exp(rowSums(cbind(rep(1,nrow(dat)),dat)[-c(1:100),c(1,5,6)]*
               coefmat[-c(1:100),c(1,5,6)]))/
          (1+exp(rowSums(cbind(rep(1,nrow(dat)),dat)[-c(1:100),c(1,5,6)]*
               coefmat[-c(1:100),c(1,5,6)]))))
u <- runif (length(y))
y <- as.numeric (u < y)

#Consider three candidate models
mmat<-matrix(c(1,1,1,1,1,0,0,0,1,1,1,0,1,0,1),3,5,byrow=TRUE)

#Fit model and plot
#autotuning is turned off for this demonstration example
ldma.test<-logistic.dma(dat,y,mmat,lambda=.99,alpha=.99,autotune=FALSE)
plot(ldma.test)

}
