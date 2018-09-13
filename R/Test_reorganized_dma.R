library("dma")
#simulate some data to test
#first, static coefficients
coef<-c(.08,-.4,-.1)
coefmat<-cbind(rep(coef[1],200),rep(coef[2],200),rep(coef[3],200))
#then, dynamic ones
coefmat<-cbind(coefmat,seq(1,.45,length.out=nrow(coefmat)),
               seq(-.75,-.15,length.out=nrow(coefmat)),
               c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)))
npar<-ncol(coefmat)-1

set.seed(1234)
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


modl <- InitializeDMALogistic(dat[1:20,],y[1:20],mmat,lambda=.99,alpha=.99,autotune=FALSE)
yhat <- matrix(0, ncol=3, nrow=200)
for(i in 21:200){
  yhat[i,] <- PredictDMALogistic(modl,mmat,dat[i,])
  modl <- UpdateDMALogistic(modl, mmat, dat[i,], y[i])
}
#yhat <- PredictDMALogistic(modl,mmat,dat[21:25,])
sum(abs(modl$yhatmodel[,21:200]-ldma.test$yhatmodel[,21:200]))
modl$yhatmodel[2,150:200]
ldma.test$yhatmodel[2,150:200]

mod.yhat <- t(modl$yhatmodel[,21:200])
mod.yhat[1:10,] - yhat[21:30,]

mod2<- DynamicModelAvg(modl,mmat)
sum(abs(mod2$yhatdma - ldma.test$yhatdma))
