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
\alias{logdma.init}
\alias{logdma.update}
\alias{logdma.predict}
\alias{logdma.average}

\title{
Dynamic model averaging for binary outcomes
}
\description{
Implements dynamic model averaging for binary outcomes as described in 
McCormick et al. (2011, Biometrics).  It can be either performed for all data at once (using \code{logistic.dma}), or dynamically for one observation at a time (combining the remaining functions, see Example).
Along with the values described below, plot() creates a plot of the posterior model probabilities over time and 
model-averaged fitted values (with smooth curve overlay) and print() returns model matrix and posterior
model probabilities.  There are K candidate
models, T time points, and d total covariates (including the intercept).
}
\usage{
logistic.dma(x, y, models.which, lambda = 0.99, alpha = 0.99,autotune = TRUE, 
    initmodelprobs = NULL, initialsamp = NULL)
 
logdma.init(x, y, models.which)

logdma.predict(fit, newx)

logdma.update(fit, newx, newy, lambda = 0.99, autotune = TRUE)

logdma.average(fit, alpha = 0.99, initmodelprobs = NULL)

}
\arguments{
  \item{x}{T by (d-1) matrix of observed covariates.  Note that a column of 1's is added
    automatically for the intercept. In \code{logdma.init}, this matrix contains only the training set.}
  \item{y}{T vector of binary responses. In \code{logdma.init}, these correspond to the training set only.}
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
    \item{newx, newy}{Subset of \code{x} and \code{y} corresponding to new observations.}
    \item{fit}{List with estimation results that are outputs of functions \code{logdma.init}, \code{logdma.update} and \code{logdma.average}.}
}

\details{
The function \code{logistic.dma} is composed of three parts, which can be also used separately: First, the model is trained with a subset of the data  (function \code{logdma.init}), where the size of the training set is determined by \code{initialsamp}. Note that arguments \code{x} and \code{y} in \code{logdma.init} should contain the training subset only. Then, the estimation is updated with new observations (function \code{logdma.update}). Lastly, a dynamic model averaging is performed on the final estimates (function \code{logdma.average}). The updating, averaging and in addition predicting (\code{logdma.predict}) can be performed dynamically for one observation at a time, see Example below.
}


\value{
Functions \code{logistic.dma} and \code{logdma.average} return an object of class \code{logistic.dma}. Functions \code{logdma.init} and \code{logdma.update} return a list with estimation results which is a subset of the \code{logistic.dma} object. It has the following components:
  \item{x}{T by (d-1) matrix of covariates}
  \item{y}{T by 1 vector of binary responses}
  \item{models.which}{K by (d-1) matrix of candidate models}
  \item{lambda}{scalar, tuning factor within models}
  \item{alpha}{scalar, tuning factor for model averaging}
  \item{autotune}{T/F, indicator of whether or not to use autotuning algorithm}
  \item{alpha.used}{T vector of alpha values used}
  \item{theta}{K by T by d array of dynamic logistic regression estimates for each model}
  \item{vartheta}{K by T by d array of dynamic logistic regression variances for each model}
  \item{pmp}{K by T array of posterior model probabilities}
  \item{yhatdma}{T vector of model-averaged predictions}
  \item{yhatmodel}{K by T vector of fitted values for each model}
  Function \code{logdma.predict} returns a matrix with predictions corresponding to the \code{newx} covariates.
}
\references{
McCormick, T.M., Raftery, A.E., Madigan, D. and Burd, R.S. (2011) "Dynamic Logistic Regression and 
Dynamic Model Averaging for Binary Classification." Biometrics, 66:1162-1173.}
\author{
Tyler H. McCormick, David Madigan, Adrian Raftery

Sevvandi Kandanaarachchi and Hana Sevcikova implemented the "streaming" functionality, i.e. the original function was decomposed into standalone  parts that can be used separately for one observation at a time.
}
% Example takes too much time. It's tested in tests.
\examples{
# simulate some data to test
# first, static coefficients
coef <- c(.08,-.4,-.1)
coefmat <- cbind(rep(coef[1],200),rep(coef[2],200),rep(coef[3],200))
# then, dynamic ones
coefmat <- cbind(coefmat,seq(1,.45,length.out=nrow(coefmat)),
            seq(-.75,-.15,length.out=nrow(coefmat)),
            c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)))
npar <- ncol(coefmat)-1

# simulate data
set.seed(1234)
dat <- matrix(rnorm(200*(npar),0,1),200,(npar))
ydat <- exp(rowSums((cbind(rep(1,nrow(dat)),dat))[1:100,]*coefmat[1:100,]))/
          (1+exp(rowSums(cbind(rep(1,nrow(dat)),dat)[1:100,]*coefmat[1:100,])))
y <- c(ydat,exp(rowSums(cbind(rep(1,nrow(dat)),dat)[-c(1:100),c(1,5,6)]*
               coefmat[-c(1:100),c(1,5,6)]))/
          (1+exp(rowSums(cbind(rep(1,nrow(dat)),dat)[-c(1:100),c(1,5,6)]*
               coefmat[-c(1:100),c(1,5,6)]))))
u <- runif (length(y))
y <- as.numeric (u < y)

# Consider three candidate models
mmat <- matrix(c(1,1,1,1,1,0,0,0,1,1,1,0,1,0,1),3,5, byrow = TRUE)

# Fit model and plot
# autotuning is turned off for this demonstration example
ldma.test <- logistic.dma(dat, y, mmat, lambda = .99, alpha = .99, 
    autotune = FALSE, initialsamp = 20)
plot(ldma.test)

# Using DMA in a "streaming" mode
modl <- logdma.init(dat[1:20,], y[1:20], mmat)
yhat <- matrix(0, ncol=3, nrow=200)
for(i in 21:200){
  # if prediction is desired, use logdma.predict
  yhat[i,] <- logdma.predict(modl, dat[i,])
  # update
  modl <- logdma.update(modl, dat[i,], y[i], 
                lambda = .99, autotune = FALSE)
}
# the averaging step could be also done within the loop above
ldma.stream <- logdma.average(modl, alpha = .99)
plot(ldma.stream)
}
