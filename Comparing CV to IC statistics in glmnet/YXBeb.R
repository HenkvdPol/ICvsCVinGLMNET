YXBeb <- function(n, p, q, rho, sd = 1, a = 0.1, b = 1, intercept = 0, alternate_beta = FALSE){
  # n = sample size, 
  # p = number of parameters, 
  # q = sparsity index, 
  # rho = correlation factor
  # sd = sd of eps ~ N(0,sd^2), 
  # a, b =  lower and upper bounds for the uniform distribution,
  #         to create the block matrix and multivariate normal matrix.
  # intercept = beta_0,
  # alternate_beta = binary value indicating if we use alternate non-zero
  #                  to zero elements as our true beta.
  
  # Install specific libraries here within the function such that
  # the forach() function works properly.
  library(Matrix)
  library(mvtnorm)
  
  
  # Make p even so that we can create a proper block matrix. Not really
  # necessary.
  p <- p + p%%2

  if(alternate_beta){
    # Create Beta that has alternating zero and non-zero components.
    Beta <-replace(c(runif(ceiling(p-(p /4)) , 
                           min = a,max = b) * sample(c(-1,1),
                                                     size = ceiling(p-(p/4)),replace = TRUE)),
                   list = 2*((ceiling(p-(p/4))/2) : 1),
                   values = rep(0, (ceiling(p-(p/4))/2)))
    
    Beta <- c(Beta, c(runif(floor(q/4),min = a,max = b) * sample(c(-1, 1), size = floor(q/4), replace = TRUE), 
                      rep(0, ((p/4) - floor(q/4)))))
  }
  else{
    # or Beta is some random parameter with q values of non-zero elements.
    Beta <- c(runif(q, 
                    min = a, 
                    max = b) * sample(c(-1, 1), 
                                      size = q, 
                                      replace = TRUE), 
              rep(0, p-q))
  }
  
  # Noise, generated from N(0, sd^2)
  Epsilon <- rnorm(n, mean = 0, sd = sd)
  
  # Block matrix that is the variance matrix of X.
  block_sigma <- matrix(c(1, rho, 
                          rho, 1), ncol = 2)
  Sigma <- as.matrix(bdiag(replicate(p/2,
                                     block_sigma, 
                                     simplify = FALSE)))
  
  # Input, generated from N_n((0,...,0), Sigma)
  X <- rmvnorm(n, mean = rep(0,p), sigma = Sigma)
  
  # Create the observation
  Y <- X %*% Beta + Epsilon + intercept
  
  # Return a list as result.
  out = list(Y, X, Epsilon, Beta)
  names(out) = c('Y', 'X', 'eps', 'Beta')
  return(out)
}
