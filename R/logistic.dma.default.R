logistic.dma.default <-
    function(x, y, models.which, lambda=0.99, alpha=0.99,autotune=TRUE, 
             initmodelprobs=NULL,initialsamp=NULL) {
        
        #K is number of candidate models
        #T is time
        #d is total number of covariates (including intercept)
        #inputs
        #########
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

      #how much of the sample should be used to generate initial values?
      if(is.null(initialsamp)){initialsamp<-round(nrow(x)/10,0)}
      if(initialsamp >= nrow(x)) stop("Initial sample must be smaller than the number of observations (", nrow(x), ").")
        
      # take the initial sample as training data
      trainX <- x[1:initialsamp,, drop = FALSE]
      trainY <- y[1:initialsamp]
 
      # train the dma logistic model
      mod <- logdma.init(trainX, trainY, models.which)
      
      # take the rest as test data  
      st <- initialsamp+1
      en <- dim(x)[1]
      testX <- x[st:en,, drop = FALSE]
      testY <- y[st:en]
      
      # update the model
      mod <- logdma.update(mod, testX, testY, lambda = lambda, autotune = autotune)
      
      # dynamic model averaging
      mod <- logdma.average(mod, alpha = alpha, initmodelprobs = initmodelprobs)
      
      # removing unwanted attributes(these were needed for updating the model)
      mod$varcov <- NULL
      mod$laplacemodel <- NULL

      class(mod)<-"logistic.dma"
      
      mod
}

logdma.init <- function(x, y, models.which) {
    # Prepare for estimation using training data
    #K is number of candidate models
    #T is time
    #d is total number of covariates (including intercept)
    #inputs
    #x (d-1) by T matrix of observed covariates (training data only).  Note that a column of 1's is added
    #automatically for the intercept
    #y T vector of binary responses (training data only)
    #models.which K by (d-1) matrix defining models.  A 1 indicates a covariate is included
    #in a particular model, a 0 if it is excluded.  Model averaging is done over all
    #modeld specified in models.which
    
    K<-nrow(models.which)
    nx <- nrow(x)
    
    #set up arrays with output for each candidate model
    baytah<-array(dim=c(K, nx, (ncol(x)+1)))
    varbaytah<-array(dim=c(K, nx, (ncol(x)+1)))
    yhatmodel<-array(dim=c(K, nx))
    varcovar<-array(dim=c(K, (ncol(x)+1),(ncol(x)+1)))
    
    #dynamic logistic regression for each candidate model
    #this section could also be done in parallel
    for(mm in 1:K){
        #data matrix for model mm
        xdat<-x[,c(which(models.which[mm,]==1)), drop=FALSE]
        xdat <- cbind(rep(1,dim(xdat)[1]),xdat)
        sel.rows <- c(1,1+which(models.which[mm,]==1))
        #generate inital values using glm
        init.temp <- dlogr.init(xdat, y)
        d <- length(init.temp$BetaHat)
        baytah[mm, , sel.rows] <- matrix(rep(init.temp$BetaHat, nx), nx, d, byrow=TRUE)
        varbaytah[mm, , sel.rows] <- matrix(rep(diag(init.temp$VarBetaHat), nx), nx, d, byrow=TRUE)
        varcovar[mm, sel.rows, sel.rows] <- init.temp$VarBetaHat ## SK - This stores the initial variance covariance matrix
        yhatmodel[mm, ] <- 1/(1 + exp( - xdat %*% init.temp$BetaHat))
    }
    return(list(x=x,y=y,models=models.which,
               varcov = varcovar,  
               theta=baytah,
               vartheta=varbaytah,
               yhatmodel=yhatmodel))
}



logdma.predict <- function(fit, newx){
    models.which <- fit$models
    if(is.null(dim(newx)))
        newx <- rbind(newx) # creates a matrix with one row
    
    newx <- cbind(rep(1,dim(newx)[1]), newx)
    num.mods <- dim(fit$theta)[1]
    last <- dim(fit$theta)[2]
    
    yhat <- matrix(0, nrow=nrow(newx), ncol=nrow(models.which))
    for(mm in 1:num.mods){
        sel.rows <- c(1,1+which(models.which[mm,]==1))
        xx <- newx[,sel.rows]
        BetaHat <- cbind(fit$theta[mm,last,sel.rows]) # creates a column from a vector
        yhat[,mm] <- dlogr.predict(xx,BetaHat)
    }
    return(yhat)
}



logdma.update <- function(fit, newx, newy, lambda = 0.99, autotune=TRUE){
    # Update fit for new observation(s)
    #autotune T/F indicates whether or not the automatic tuning procedure desribed in 
    #McCormick et al. should be applied.  Default is true.
    #lambda scalar tuning factor within models
    models.which <- fit$models
    if(is.null(dim(newx)))
        newx <- rbind(newx) # creates a matrix with one row
    x <- rbind(fit$x,newx)
    dimnames(x) <- NULL
    y <- c(fit$y, newy)
    
    newx <- cbind(1, newx)
    last <- dim(fit$theta)[2]
    
    big.len <- nrow(fit$x) + nrow(newx)
    K<-nrow(models.which)
    L <- ncol(fit$x)+1

    # Store fit values in the arrays
    len <- dim(fit$x)[1]
    update.index <- (len+1):big.len
    
    # set up arrays with output for each candidate model
    theta<-array(dim=c(K, big.len, L))
    vartheta<-array(dim=c(K, big.len, L))
    laplacemodel<-array(dim=c(K, big.len))
    yhatmodel<-array(dim=c(K, big.len))
    varcov<-array(dim=c(K, L, L))
    
    theta[,1:len,] <- fit$theta
    vartheta[,1:len,] <- fit$vartheta
    yhatmodel[,1:len] <- fit$yhatmodel
    
    if(is.null(fit$laplacemodel)){
        laplacemodel[,1:len]<-rep(0,len)
    }else{
        laplacemodel[,1:len]<-fit$laplacemodel
    }
    
    for(mm in 1:K){
        sel.rows <- c(1,1+which(models.which[mm,]==1))
        xx <- newx[,sel.rows, drop = FALSE]
        d <- length(sel.rows)
        #matrix of possible combinations of lambda
        tune.mat<- if(autotune==FALSE) matrix(rep(lambda,d),nrow=1,ncol=d) else tunemat.fn(lambda,1,d)

        BetaHat <- cbind(fit$theta[mm,last,sel.rows]) # creates a column from a vector
        varbetahat.tm1 <- fit$varcov[mm,sel.rows,sel.rows]

        # update
        for(i in 1:nrow(newx)) {
            x.t <- xx[i,]
            y.t <- newy[i]
            step.tmp <- dlogr.step(x.t,y.t,BetaHat,varbetahat.tm1,tune.mat)
            #compute fitted value
            yhatmodel[mm, update.index[i]] <- 1/(1 + exp(- x.t%*%step.tmp$betahat.t))
            # gather output
            theta[mm, update.index[i], sel.rows] <- BetaHat <- step.tmp$betahat.t
            vartheta[mm, update.index[i], sel.rows] <- diag(step.tmp$varbetahat.t)
            laplacemodel[mm, update.index[i]] <- step.tmp$laplace.t
            varbetahat.tm1 <- step.tmp$varbetahat.t
        }
        varcov[mm, sel.rows, sel.rows] <- step.tmp$varbetahat.t
    }
    est <- fit
    # update only items that were changed here
    for(item in c("x", "y", "lambda", "theta", "vartheta", "yhatmodel", "varcov", "laplacemodel"))
        est[[item]] <- get(item)
    
    return(est)
}


logdma.average <- function(fit, alpha = 0.99, initmodelprobs=NULL){
    # inputs
    #alpha scalar forgetting factor for model averaging
    #initmodelprobs K vector of starting probabilities for model averaging.  If null (default),
    #then use 1/K for each model.
    
    laplacemodel <- fit$laplacemodel
    K <- nrow(fit$models)
    x <- fit$x
    
    if(any(is.na(laplacemodel)) || any(laplacemodel==Inf) || any(laplacemodel==-Inf))
        warning("At least one laplace approximation is not well behaved. This will likely lead to issues with posterior model probabilities. This is likely a computation issue")

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
        omega.tk[,1]<-1/K
    }else{pi.t.t[,1]<-initmodelprobs
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
            omega.tk[k,t]<-pi.t.tm1[k,t]*max(exp(laplacemodel[k,t]), .Machine$double.xmin)
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
    yhatdma=colSums(fit$yhatmodel*pi.t.t)
    
    #outputs#
    #x time by (d-1) matrix of covariates
    #y time by 1 vector of binary responses
    #models.which K by (d-1) matrix of candidate models
    #lambda scalar, tuning factor within models
    #alpha scalar, forgetting factor for model averaging 
    #autotune T/F, indicator of whether or not to use autotuning algorithm
    #alpha.used time-length vector of alpha values used
    #theta K by time by d array of dynamic logistic regression estimates for each model
    #vartheta K by time by d array of dynamic logistic regression variances for each model
    #pmp K by time array of posterior model probabilities
    #yhatbma time length vector of model-averaged predictions
    #yhatmodel K by time vector of fitted values for each model
    
    #define outputs
    est <- fit
    est$yhatdma <- yhatdma
    est$pmp <- pi.t.t
    est$alpha <- alpha
    est$alpha.used <- alpha.vec
    est$fitted.values <- yhatdma
    est$residuals <- fit$y - yhatdma
    
    #define as class to allow other commands later
    class(est)<-"logistic.dma"
    
    return(est)
}

logistic.dma.default.static <-
function(x, y, models.which, lambda=0.99, alpha=0.99,autotune=TRUE, 
 initmodelprobs=NULL,initialsamp=NULL) {
  
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

#how much of the sample should be used to generate initial values?
if(is.null(initialsamp)){initialsamp<-round(nrow(x)/10,0)}

#dynamic logistic regression for each candidate model
#this section could also be done in parallel
for(mm in 1:nrow(models.which)){
	#data matrix for model mm
	xdat<-x[,c(which(models.which[mm,]==1)), drop=FALSE]
	xdat <- cbind(rep(1,dim(xdat)[1]),xdat)
	d <- dim(xdat)[2] 

#matrix of possible combinations of lambda
tune.mat<- if(autotune==FALSE) matrix(rep(lambda,d),nrow=1,ncol=d) else tunemat.fn(lambda,1,d)

#generate inital values using glm and use them for the first initialsamp observations
init.temp<-dlogr.init(xdat[1:initialsamp,],y[1:initialsamp])
betahat.tm1 <- init.temp$BetaHat
baytah[mm,1:initialsamp,c(1,1+which(models.which[mm,]==1))] <- matrix(rep(init.temp$BetaHat,initialsamp),initialsamp,length(init.temp$BetaHat),byrow=TRUE)
varbaytah[mm,1:initialsamp,c(1,1+which(models.which[mm,]==1))] <- matrix(rep(diag(init.temp$VarBetaHat),initialsamp),initialsamp,length(init.temp$BetaHat),byrow=TRUE)
varbetahat.tm1 <- init.temp$VarBetaHat
laplacemodel[mm,1:initialsamp]<-rep(0,initialsamp)
yhatmodel[mm,1:initialsamp]<-exp(xdat[1:initialsamp,]%*%init.temp$BetaHat)/(1+exp(xdat[1:initialsamp,]%*%init.temp$BetaHat))
bite.size <- 1
for (i in seq((initialsamp+1),nrow(xdat),by=bite.size)) {
  x.t <- (xdat[i:(i+bite.size-1),])
  y.t <- y[i:(i+bite.size-1)]
  #update
  step.tmp <- dlogr.step(x.t,y.t,betahat.tm1,varbetahat.tm1,tune.mat)
  #gather output
  baytah[mm,i:(i+bite.size-1),c(1,1+which(models.which[mm,]==1))] <- step.tmp$betahat.t
  varbaytah[mm,i:(i+bite.size-1),c(1,1+which(models.which[mm,]==1))] <- diag(step.tmp$varbetahat.t)
  laplacemodel[mm,i:(i+bite.size-1)]<-step.tmp$laplace.t 
  #compute fitted value
  yhatmodel[mm,i]<-exp(x.t%*%step.tmp$betahat.t)/(1+exp(x.t%*%step.tmp$betahat.t))
  #prepare for next obs
  betahat.tm1 <- step.tmp$betahat.t
  varbetahat.tm1 <- step.tmp$varbetahat.t  
} 
print(paste("Finished processing model",mm))  
}
if(sum(is.na(laplacemodel))>0|sum(laplacemodel==Inf)>0|sum(laplacemodel==-Inf)>0)
  {print("Warning: At least one laplace approximation is not well behaved.  This will likely lead to issues with posterior model probabilities.
         This is likely a computation issue")}
  
#Dynamic model averaging
#set up arrays for posterior model probabilities 
dimnames <- list(paste("model", 1:K, sep="_"), paste("time", 1:nrow(x), sep="_"))
pi.t.t<-array(dim=c(K,nrow(x)), dimnames=dimnames)
pi.t.tm1<-array(dim=c(K,nrow(x)),dimnames=dimnames)
omega.tk<-array(dim=c(K,nrow(x)),dimnames=dimnames)

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
                        yhatdma=yhatdma,yhatmodel=yhatmodel))
        
#define outputs
    est$fitted.values<-yhatdma
    est$residuals<-y-yhatdma
    
    #define as class to allow other commands later
    class(est)<-"logistic.dma"
    
    est
}

