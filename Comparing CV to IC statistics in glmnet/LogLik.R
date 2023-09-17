LogLik <- function(gfit, X, Y, n){
  # Seperate function to compute the Log Likelihood from the data
  # gfit = fitted estimators from YXBeb
  
  fun <- function(lambda){
    hatBeta = coef(gfit, s = lambda)
    hatMu = X%*%hatBeta[-1] + hatBeta[1]
    hatSigma = sqrt((1/(n-1)) * sum((Y - hatMu)^2))
    
    # the log likelihood
    return(sum(dnorm(Y, mean = hatMu, sd = hatSigma, log = TRUE)))
  }
  
  loglik <- apply(as.array(gfit$lambda), MARGIN = 1, FUN = fun)
  return(loglik)
}
