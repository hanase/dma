context('Logistic DMA')

generate.data <- function() {
    #simulate some data to test
    coef <- c(.08,-.4,-.1)
    coefmat <- cbind(rep(coef[1],200),rep(coef[2],200),rep(coef[3],200))
    #then, dynamic ones
    coefmat <- cbind(coefmat,seq(1,.45,length.out=nrow(coefmat)),
                     seq(-.75,-.15,length.out=nrow(coefmat)),
                     c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)))
    npar <- ncol(coefmat)-1
    
    #simulate data
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
    
    #Consider three candidate models
    mmat <- matrix(c(1,1,1,1,1,0,0,0,1,1,1,0,1,0,1),3,5, byrow = TRUE)
    return(list(dat = dat, y = y, mmat = mmat))
}

run.logistic <- function(lambda = .99, alpha = .99, autotune = FALSE, initialsamp = 20) {
    d <- generate.data()
    # static mode
    ldma.stat <- logistic.dma(d$dat, d$y, d$mmat, lambda = lambda, alpha = alpha, 
                              autotune = autotune, initialsamp = initialsamp)
    
    # Using DMA in a "streaming" mode
    modl <- logdma.init(d$dat[1:20,], d$y[1:20], d$mmat)
    yhat <- matrix(0, ncol=3, nrow=200)
    for(i in 21:200){
        yhat[i,] <- logdma.predict(modl, d$dat[i,])
        modl <- logdma.update(modl, d$dat[i,], d$y[i], 
                              lambda = lambda, autotune = autotune)
    }
    ldma.stream <- logdma.average(modl, alpha = alpha)
    ldma.stream$laplacemodel <- NULL
    ldma.stream$varcov <- NULL
    return(list(stat = ldma.stat, stream = ldma.stream))
}

test_that('Example in logistic.dma without tuning works', {
    ldma <- run.logistic()
    expect_true(length(names(ldma$stat)) == 13)
    expect_equal(dim(ldma$stat$pmp), c(3,200))
    expect_false(any(is.na(ldma$stat$pmp)))
    expect_identical(ldma$stat[order(names(ldma$stat))], 
                     ldma$stream[order(names(ldma$stream))])
})

test_that('Example in logistic.dma with tuning works', {
    skip_on_cran()
    ldma <- run.logistic(autotune = TRUE)
    expect_true(length(names(ldma$stat)) == 13)
    expect_equal(dim(ldma$stat$pmp), c(3,200))
    expect_false(any(is.na(ldma$stat$pmp)))
    expect_identical(ldma$stat[order(names(ldma$stat))], 
                     ldma$stream[order(names(ldma$stream))])
})
