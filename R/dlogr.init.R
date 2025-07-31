dlogr.init <-
function (x,y) {
# function to get the initial estimate of beta and var(beta)
	MyData <- data.frame(x,y)
	MyModel <- glm(y~x-1,data=MyData,family=binomial(link=logit))
	VarBetaHat <- vcov(MyModel)
	BetaHat <- coefficients(MyModel)
	# if(regul) { 
	#     # regularization (Bayesian estimate)
	#     gamma <- diag(c(4, apply(x[,-1], 2, var))) # precision
	#     invVarBetaHat <- solve(VarBetaHat)
	#     VarBetaHat <- solve(gamma + invVarBetaHat)
	#     BetaHat <- VarBetaHat %*% invVarBetaHat %*% BetaHat
	# }
	rm(MyData)	
	return(list(BetaHat = BetaHat,VarBetaHat = VarBetaHat))
}

