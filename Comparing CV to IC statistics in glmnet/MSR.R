# Load all the libraries, just in case.
library(Matrix)
library(glmnet)
library(iterators)
library(parallel)
library(mvtnorm)
library(foreach)
library(doParallel)

# Import functions
# setwd("directorypath") Misschien nodig voor ALICE
source("YXBeb.R")
source("LogLik.R")
source("accuracy.R")

# Set seed value.
set.seed(123)

# Monte Carlo function.
compute <- function(n, p, q, rho = 0, sd, a = 0, b= 1, 
                    intercept = 0, nlambda, 
                    Lik = TRUE, incl_intercept = TRUE, 
                    alternate_beta = FALSE){
  # n = sample size,
  # p = dimension of the parameters,
  # q = sparsity index,
  # rho = correlation factor,
  # sd = sqrt(variance) of the noise,
  # a = lower bound of uniform distribution,
  # b = upper bound of uniform distribution,
  # intercept = Beta_0,
  # nlambda = number of lambda values to consider in glmnet,
  # Lik = Binary value to either use the LogLik() or the deviance as LR,
  # incl_intercept = include the intercept to be fitted in glmnet,
  # alternate_beta = choose Beta such that has alternating non-zero elements.
  
  
  # Generate random data
  df <- YXBeb(n, p, q ,rho, sd, a, b, intercept, alternate_beta)
  X_train <- df$X
  Y_train <- df$Y
  Beta <- df$Beta
  
  # Initial array for the results.
  TP_non_zero = c()
  TN_non_zero = c()
  T_sum_non_zero = c()
  P_shrunk = c()
  rmse = c()
  mae = c()
  med = c() # The median.
  
  ## add call to requested libraries:
  library(glmnet)
  
  # Fit the data
  fit <- cv.glmnet(X_train, Y_train, nlambda = nlambda, 
                 intercept = incl_intercept, 
                 lambda.min.ratio = 0.001)
  
  # Either use the Log Likelihood or the deviance for the test statistics.
  if(Lik) L <- LogLik(fit$glmnet.fit, X_train, Y_train, n) else L <- deviance(fit$glmnet.fit)
  
  # Get number of free parameters
  k <- fit$glmnet.fit$df + 1 
  
  # Compute the criterion statistics.
  aic <- 2*k - 2*L 
  aicc <- 2*k - 2*L + ((2*(k^2) + 2*k) / (n-k-1))
  bic <- k*log(n) - 2*L 
  # Not fully sure about the GIC but this is what I make up from the article and what I can do
  # from the article, the deviance should be scaled, but it is not really 
  # clear to me what it means with respect to the deviance.
  gic <- (1/n) * (L - (log(log(n))*log(p))*k) 
  
  # get position of best AIC, AICc, BIC and GIC:
  pos_aic <- which(aic == min(aic))
  pos_aicc <- which(aicc == min(aicc))
  pos_bic <- which(bic == min(bic)) 
  pos_gic <- which(gic == min(gic))
  
  # Extra test data set for testing to avoid bias.
  df_test <- YXBeb(n, p, q ,rho, sd, a, b, intercept, alternate_beta)
  X_test <- df_test$X
  Y_test <- df_test$Y
  
  # Compute the RMSE, MAE and median for the criterion statistics.
  for(i in c(pos_aic, pos_aicc, pos_bic, pos_gic)){
    # Estimated parameter and response.
    beta_i <- coef(fit, s = fit$lambda[i])
    hat_Y <- beta_i[1] + X_test %*% beta_i[-1]
    
    mae  = c(mae, mean(abs(hat_Y - Y_test)))
    rmse = c(rmse, sqrt(mean((hat_Y - Y_test)^2)))
    med = c(med, median(abs(hat_Y - Y_test)))
    
    # This is for the accuracy
    acc <- accuracy(Beta, beta_i[-1])
    TP_non_zero <- c(TP_non_zero, acc$TP)
    TN_non_zero <- c(TN_non_zero, acc$TN)
    T_sum_non_zero <- c(T_sum_non_zero, acc$T_sum)
    P_shrunk <- c(P_shrunk, acc$shrunk)
  }
  
  # Compute the estimators using CV using the default measure
  cvBeta <- coef(fit, s = fit$lambda.min)
  cv_Y <- X_test %*% cvBeta[-1] + cvBeta[1]

  # Accuracy of CV coefficients
  acc <- accuracy(Beta, cvBeta[-1])
  TP_non_zero <- c(TP_non_zero, acc$TP)
  TN_non_zero <- c(TN_non_zero, acc$TN)
  T_sum_non_zero <- c(T_sum_non_zero, acc$T_sum)
  P_shrunk <- c(P_shrunk, acc$shrunk)
  
  # CV with mae type.measure method:
  # Here I have to do cv.glmnet again because we use a different
  # measure.
  cvfit_mae <- cv.glmnet(X_train, Y_train, nlambda = nlambda, 
                         type.measure = 'mae',
                         lambda.min.ratio = 0.001)
  beta_cv_mae <- coef(cvfit_mae, 
                      s = cvfit_mae$lambda.min)
  cv_mae_Y <- X_test %*% beta_cv_mae[-1] + beta_cv_mae[1]
  
  # Accuracy of CV_mae
  acc <- accuracy(Beta, beta_cv_mae[-1])
  TP_non_zero <- c(TP_non_zero, acc$TP)
  TN_non_zero <- c(TN_non_zero, acc$TN)
  T_sum_non_zero <- c(T_sum_non_zero, acc$T_sum)
  P_shrunk <- c(P_shrunk, acc$shrunk)
  
  # Now computing prediction error, accuracy and median
  # for the 1se of both fits.
  cv_1se_Beta <- coef(fit, s = fit$lambda.1se)
  cv_mae_1se_Beta <- coef(cvfit_mae, s = cvfit_mae$lambda.1se)
  
  Y_1se <- X_test %*% cv_1se_Beta[-1] + cv_1se_Beta[1]
  Y_1se_mae <- X_test %*% cv_mae_1se_Beta[-1] + cv_mae_1se_Beta[1]
  
  acc <- accuracy(Beta, cv_1se_Beta[-1])
  TP_non_zero <- c(TP_non_zero, acc$TP)
  TN_non_zero <- c(TN_non_zero, acc$TN)
  T_sum_non_zero <- c(T_sum_non_zero, acc$T_sum)
  P_shrunk <- c(P_shrunk, acc$shrunk)
  
  acc <- accuracy(Beta, cv_mae_1se_Beta[-1])
  TP_non_zero <- c(TP_non_zero, acc$TP)
  TN_non_zero <- c(TN_non_zero, acc$TN)
  T_sum_non_zero <- c(T_sum_non_zero, acc$T_sum)
  P_shrunk <- c(P_shrunk, acc$shrunk)
  
  rmse <- c(rmse,
            sqrt(mean((cv_Y - Y_test)^2)),
            sqrt(mean((cv_mae_Y - Y_test)^2)),
            sqrt(mean((Y_1se - Y_test)^2)),
            sqrt(mean((Y_1se_mae - Y_test)^2)))
  
  mae <- c(mae,
           mean(abs(cv_Y - Y_test)),
           mean(abs(cv_mae_Y - Y_test)),
           mean(abs(Y_1se - Y_test)),
           mean(abs(Y_1se_mae - Y_test)))
  
  med <- c(med, 
           median(abs(cv_Y - Y_test)),
           median(abs(cv_mae_Y - Y_test)),
           median(abs(Y_1se - Y_test)),
           median(abs(Y_1se_mae - Y_test)))
  
  
  pos_cv <- which(fit$lambda == fit$lambda.min)
  pos_1se <- which(fit$lambda == fit$lambda.1se)
  pos_cv_mae <- which(cvfit_mae$lambda == cvfit_mae$lambda.min)
  pos_1se_mae <- which(cvfit_mae$lambda == cvfit_mae$lambda.min)
  
  
  pos <- c(pos_aic, pos_aicc, pos_bic, pos_gic, 
           pos_cv, pos_cv_mae, pos_1se, pos_1se_mae)
  
  # Return the result in a data frame. positional value is only interesting if we do it one time.
  result <- data.frame('RMSE' = rmse, 
                       'MAE' = mae, 
                       'median' = med,
                       'TP' = TP_non_zero, 
                       'TN' = TN_non_zero,
                       'T_sum' = T_sum_non_zero,
                       'shrunked' = P_shrunk,
                       'pos' = pos, 
                       row.names =  c('aic', 'aicc', 'bic', 'gic', 
                                      'cv_dev', 'cv_mae', '1se_dev', '1se_mae')
                       )
  return(result)
} 


# ---- Running the function -----

# Get the sum of the elements. Special function.
simplify <- function(x){Reduce('+', x) / M}

# Want to know the time
start <- proc.time()
# Make clusters to parallel
cl <- parallel::makeForkCluster(8)
doParallel::registerDoParallel(cl)

# For high to low dimensional setting I define a specific p for when n approx. p
# Where I take N=100.
P_100 = c(10*(1:7), (5*(1:11)+70), (50*(1:5)+100))
P_200 = c(20*(1:8), (5*(1:11)+170), (50*(1:5)+200))
P_1000 = 100*(1:10)

M <- 500 # Number of repeats.

# The Monte Carlo simulation.
g.x1 <-foreach(r = P_100) %:% foreach(1:M, 
                                      .packages = c('Matrix', 
                                                    'glmnet', 
                                                    'mvtnorm')) %dopar% compute(n=100, 
                                                                                p = r, 
                                                                                q = floor(0.7*r), 
                                                                                sd = .8*sqrt(r),
                                                                                a = 0.1,
                                                                                rho = .3, 
                                                                                intercept = 2, 
                                                                                nlambda = 1000)#,alternate_beta = TRUE)

# Stop clusturs.
parallel::stopCluster(cl)
print(proc.time() - start)

# Final product
f.x <- lapply(g.x1 , simplify)
# Print the file as a csv file so that it is easy to copy from ALICE
write.csv(f.x)
# Just to be sure, I also print the result as tables.
print(f.x)







