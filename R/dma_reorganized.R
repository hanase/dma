InitializeDMALogistic <-
  function(x, y, models.which, lambda=0.99, alpha=0.99, 
           initmodelprobs=NULL) {
    
    #load packages
    #require(mnormt)
    #require(MASS)
    
    #K is number of candidate models
    #T is time
    #d is total number of covariates (including intercept)
    #inputs
    #x (d-1) by T matrix of observed covariates.  Note that a column of 1's is added
    #automatically for the intercept
    #y T vector of binary responses
    #models.which K by (d-1) matrix defining models.  A 1 indicates a covariate is included
    #in a particular model, a 0 if it is excluded.  Model averaging is done over all
    #modeld specified in models.which
    #lambda scalar forgetting factor with each model
    #alpha scalar forgetting factor for model averaging
    #autotune T/F indicates whether or not the automatic tuning procedure desribed in 
    #McCormick et al. should be applied.  Default is true.
    #initmodelprobs K vector of starting probabilities for model averaging.  If null (default),
    #then use 1/K for each model.
    #initialsamp scalar indicating how many observations to use for generating initial 
    #values.  If null (default), then use 10 percent.  
    
    K<-nrow(models.which)
    
    #set up arrays with output for each candidate model
    baytah<-array(dim=c(nrow(models.which),nrow(x),(ncol(x)+1)))
    varbaytah<-array(dim=c(nrow(models.which),nrow(x),(ncol(x)+1)))
    laplacemodel<-array(dim=c(nrow(models.which),nrow(x)))
    yhatmodel<-array(dim=c(nrow(models.which),nrow(x)))
    varcovar<-array(dim=c(nrow(models.which),(ncol(x)+1),(ncol(x)+1)))
    
    #how much of the sample should be used to generate initial values?
    # if(is.null(initialsamp)){initialsamp<-round(nrow(x)/10,0)}
    initialsamp <- nrow(x) ## SK - Use all data for training
    
    #dynamic logistic regression for each candidate model
    #this section could also be done in parallel
    for(mm in 1:nrow(models.which)){
      #data matrix for model mm
      xdat<-x[,c(which(models.which[mm,]==1)), drop=FALSE]
      xdat <- cbind(rep(1,dim(xdat)[1]),xdat)
      d <- dim(xdat)[2] 
      
      #generate inital values using glm and use them for the first initialsamp observations
      init.temp<-dlogr.init(xdat[1:initialsamp,],y[1:initialsamp])
      betahat.tm1 <- init.temp$BetaHat
      baytah[mm,1:initialsamp,c(1,1+which(models.which[mm,]==1))] <- matrix(rep(init.temp$BetaHat,initialsamp),initialsamp,length(init.temp$BetaHat),byrow=TRUE)
      varbaytah[mm,1:initialsamp,c(1,1+which(models.which[mm,]==1))] <- matrix(rep(diag(init.temp$VarBetaHat),initialsamp),initialsamp,length(init.temp$BetaHat),byrow=TRUE)
      varbetahat.tm1 <- init.temp$VarBetaHat ## SK - This stores the initial variance covariance matrix
      laplacemodel[mm,1:initialsamp]<-rep(0,initialsamp)
      yhatmodel[mm,1:initialsamp]<-exp(xdat[1:initialsamp,]%*%init.temp$BetaHat)/(1+exp(xdat[1:initialsamp,]%*%init.temp$BetaHat))
      
      sel.rows <- c(1,1+which(models.which[mm,]==1))
      varcovar[mm,sel.rows,sel.rows] <- varbetahat.tm1
    }
    est<-(list(x=x,y=y,models=models.which,lambda=lambda,alpha=alpha,
               varcov = varcovar,  
               theta=baytah,
               vartheta=varbaytah,
               yhatmodel=yhatmodel))
    
    
    #define as class to allow other commands later
    #class(est)<-"logistic.dma"
    
    return(est)
  }



PredictDMALogistic <- function(fit, models.which, newx){
  if(is.null(dim(newx))){
    newx <- as.matrix(newx,  nrow=1)
    dim(newx) <- c(1,length(newx))
  }
  
  newx <- cbind(rep(1,dim(newx)[1]), newx)
  num.mods <- dim(fit$theta)[1]
  last <- dim(fit$theta)[2]
  
  yhat <- matrix(0, nrow=nrow(newx), ncol=nrow(models.which))
  for(mm in 1:num.mods){
    sel.rows <- c(1,1+which(models.which[mm,]==1))
    xx <- newx[,sel.rows]
    BetaHat <- fit$theta[mm,last,sel.rows]
    dim(BetaHat) <- c(length(sel.rows),1)
    yhat[,mm] <-   dlogr.predict(xx,BetaHat)
  }
  return(yhat)
}



UpdateDMALogistic <- function(fit, models.which, newx, newy, autotune=TRUE){
  require("mnormt")
  require("MASS")
  #newx <- as.data.frame(newx)
  # x <- matrix(0, nrow=(dim(fit$x)[1]+dim(newx)[1]), ncol=dim(newx)[2])
  # x[1:nrow(fit$x),] <- fit$x
  # x[(nrow(fit$x)+1):dim(x)[1],] <- newx
  x <- rbind(fit$x,newx)
  rownames(x) <- NULL
  dimnames(x) <- NULL
  y <- c(fit$y, newy)
  if(is.null(dim(newx))){
    newx <- c(1, newx)
    klm <- 1
    newx <- as.matrix(newx)
    dim(newx) <- c(1,length(newx))
  }else{
    klm <- dim(newx)[1]
    newx <- cbind(rep(1,klm), newx)
  }
  last <- dim(fit$theta)[2]
  
  lambda <- fit$lambda
  
  
  
  big.len <- dim(fit$x)[1]+ klm
  
  K<-nrow(models.which)
  
  # set up arrays with output for each candidate model
  baytah<-array(dim=c(nrow(models.which),big.len,(ncol(fit$x)+1)))
  varbaytah<-array(dim=c(nrow(models.which),big.len,(ncol(fit$x)+1)))
  laplacemodel<-array(dim=c(nrow(models.which),big.len))
  yhatmodel<-array(dim=c(nrow(models.which),big.len))
  varcovar<-array(dim=c(nrow(models.which),(ncol(fit$x)+1),(ncol(fit$x)+1)))
  
  # Store fit values in the arrays
  len <- dim(fit$x)[1]
  baytah[,1:len,] <- fit$theta
  varbaytah[,1:len,] <- fit$vartheta
  if(is.null(fit$laplacemodel)){
    laplacemodel[,1:len]<-rep(0,len)
  }else{
    laplacemodel[,1:len]<-fit$laplacemodel
    
  }
  yhatmodel[,1:len] <- fit$yhatmodel
  
  
  # for(ll in 1:K){
  #   x1 <- cbind(rep(1,dim(fit$x)[1]), fit$x)
  #   sel.rows <- c(1,1+which(models.which[ll,]==1))
  #   
  #   yhatmodel[ll,1:len]<-apply(exp(x1[,sel.rows]*fit$theta[ll,,sel.rows]),1, sum)/apply((1+exp(x1[,sel.rows]*(fit$theta[ll,,sel.rows]))),1,sum)
  # }  

  
  for(mm in 1:K){
    sel.rows <- c(1,1+which(models.which[mm,]==1))
    xx <- newx[,sel.rows]
    d <- length(sel.rows)
    
    #matrix of possible combinations of lambda
    tune.mat <- tunemat.fn(lambda,1,d)  ## SK Commented it
    if(autotune==FALSE){tune.mat<-matrix(rep(lambda,d),nrow=1,ncol=d)}
    
    BetaHat <- fit$theta[mm,last,sel.rows]
    dim(BetaHat) <- c(length(sel.rows),1)
    varbetahat.tm1 <- fit$varcov[mm,sel.rows,sel.rows]
    step.tmp <- dlogr.step(xx,newy,BetaHat,varbetahat.tm1,tune.mat)
    
    
    
    #gather output
    baytah[mm,(len+1):big.len,c(1,1+which(models.which[mm,]==1))] <- step.tmp$betahat.t
    varbaytah[mm,(len+1):big.len,c(1,1+which(models.which[mm,]==1))] <- diag(step.tmp$varbetahat.t)
    laplacemodel[mm,(len+1):big.len]<-step.tmp$laplace.t 
    #compute fitted value
    yhatmodel[mm,(len+1):big.len]<-exp(xx%*%step.tmp$betahat.t)/(1+exp(xx%*%step.tmp$betahat.t))
    sel.rows <- c(1,1+which(models.which[mm,]==1))
    varcovar[mm,sel.rows,sel.rows] <- step.tmp$varbetahat.t
    
    #prepare for next obs
    #betahat.tm1 <- step.tmp$betahat.t
    #varbetahat.tm1 <- step.tmp$varbetahat.t  
  }
  
  
  est<-(list(x=x,y=y,models=models.which,lambda=fit$lambda,alpha=fit$alpha,
             theta=baytah,
             vartheta=varbaytah,
             yhatmodel=yhatmodel,
             varcov = varcovar,
             laplacemodel=laplacemodel))
  
  #define as class to allow other commands later
  #class(est)<-"logistic.dma"
  
  return(est)
}


DynamicModelAvg <- function(fit,models.which, initmodelprobs=NULL){
  laplacemodel <- fit$laplacemodel
  K<-nrow(models.which)
  x <- fit$x
  alpha <- fit$alpha
  yhatmodel <- fit$yhatmodel
  baytah <- fit$theta
  varbaytah <- fit$vartheta
  varcovar <- fit$varcov
  lambda <- fit$lambda
  
  
  if(sum(is.na(laplacemodel))>0|sum(laplacemodel==Inf)>0|sum(laplacemodel==-Inf)>0)
  {print("Warning: At least one laplace approximation is not well behaved.  This will likely lead to issues with posterior model probabilities.
         This is likely a computation issue")}
  
  #Dynamic model averaging
  #set up arrays for posterior model probabilities 
  dimnames <- list(paste("model", 1:K, sep = "_"), paste("time", 1:nrow(x), sep = "_"))
  pi.t.t <- array(dim = c(K, nrow(x)), dimnames = dimnames)
  pi.t.tm1 <- array(dim = c(K, nrow(x)), dimnames = dimnames)
  omega.tk <- array(dim = c(K, nrow(x)), dimnames = dimnames)
  
  #initial probability estimates
  #default is uniform or user-specified
  if(is.null(initmodelprobs)){
    pi.t.t[,1]<-1/K
    pi.t.tm1[,1]<-1/K
    omega.tk[,1]<-1/K}
  else{pi.t.t[,1]<-initmodelprobs
  pi.t.tm1[,1]<-initmodelprobs
  omega.tk[,1]<-initmodelprobs}
  
  #initialize model forgetting factor, alpha
  alpha.vec<-rep(NA,nrow(x))
  alpha.vec[1]<-alpha
  #can also include the option of always forgetting just a little by 
  #setting alpha.noforget<1
  alpha.noforget=1           
  #iterate for each subsequent
  for(t in 2:length(alpha.vec)){
    if(t>2){rm(alpha.tmp)}
    pred.forget<-sum(exp(laplacemodel[,t])*((pi.t.t[,t-1]^alpha)/sum(pi.t.t[,t-1]^alpha)))
    pred.noforget<-sum(exp(laplacemodel[,t])*((pi.t.t[,t-1]^alpha.noforget)/sum(pi.t.t[,t-1]^alpha.noforget)))
    alpha.vec[t]<-ifelse(pred.forget>=pred.noforget,alpha,alpha.noforget)      
    alpha.tmp<-alpha.vec[t]    
    for(k in 1:K){
      pi.t.tm1[k,t]<-(pi.t.t[k,t-1]^alpha.tmp)/sum(pi.t.t[,t-1]^alpha.tmp)
      omega.tk[k,t]<-(pi.t.tm1[k,t]*exp(laplacemodel[k,t]))
      #omega.tk[k,t]<-exp(laplace.mat[k,t])
    }
    for(k in 1:K){
      pi.t.t[k,t]<-omega.tk[k,t]/sum(omega.tk[,t])
      #if a model probability is within rounding error of zero,
      #the alg tends to get "stuck" at that value, so move 
      #slightly away to avoid this issue
      pi.t.t[k,t]<-ifelse(pi.t.t[k,t]>0,pi.t.t[k,t],.001)
      pi.t.t[k,t]<-ifelse(pi.t.t[k,t]<1,pi.t.t[k,t],.999)
    }
  }
  yhatdma=colSums(yhatmodel*pi.t.t)
  
  #outputs#
  #x time by (d-1) matrix of covariates
  #y time by 1 vector of binary responses
  #models.which K by (d-1) matrix of candidate models
  #lambda scalar, tuning factor within models
  #alpha scalar, tuning factor for model averaging 
  #autotune T/F, indicator of whether or not to use autotuning algorithm
  #alpha.used time-length vector of alpha values used
  #theta K by time by d array of dynamic logistic regression estimates for each model
  #vartheta K by time by d array of dynamic logistic regression variances for each model
  #pmp K by time array of posterior model probabilities
  #yhatbma time length vector of model-averaged predictions
  #yhatmodel K by time vector of fitted values for each model
  est<-(list(x=x,y=y,models=models.which,lambda=lambda,alpha=alpha,
             pmp=pi.t.t,
             alpha.used=alpha.vec,
             theta=baytah,
             vartheta=varbaytah,
             yhatdma=yhatdma,yhatmodel=yhatmodel,
             varcov=varcovar,
             laplacemodel=laplacemodel))
  
  #define outputs
  est$fitted.values<-yhatdma
  est$residuals<-y-yhatdma
  
  #define as class to allow other commands later
  class(est)<-"logistic.dma"
  
  return(est)
  
  }