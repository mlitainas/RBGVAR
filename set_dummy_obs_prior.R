dummy_obs_prior <- function(vardata, lamda = 0.5, tau = 0, epsilon = 0.1, epsilonx = 0.5, epsilonxl = 0.5){
  

  lamdaP    = lamda;
  tauP      = tau;
  epsilonP  = epsilon;
  epsilonXP = epsilonx
  epsilonXLP= epsilonxl
  
vardata$number_of_non_endogenous_regressors
  X            = vardata$x_rhs 
  Y            = vardata$y_lhs
  L            = vardata$number_of_lags
  XL           = vardata$number_of_xlags
  t_lags       = vardata$number_of_observations
  n            = vardata$number_of_endogenous_variables
  det          = vardata$number_of_deterministic
  
  
  
  det <- vardata$number_of_deterministic
  
  nexo_cont<-vardata$number_of_contemporaneous_star_variables + vardata$number_of_exogenous_variables
  nexo_lag <- vardata$number_of_lagged_star_variables + vardata$number_of_lagged_exogenous_variables
  nexo <- vardata$number_of_non_endogenous_regressors
  
  if (nexo != det + nexo_cont + nexo_lag) {
    stop(
      paste0(
        "Internal error: the non-endogenous regressor counts ",
        "are inconsistent."
      ),
      call. = FALSE
    )
  }
  
  star_cont    = X[, (det+1):(det+nexo_cont), drop = F] 
  muP = apply(Y, 2, mean)
  
  s_prior       <-  rep(0,vardata$number_of_endogenous_variables)
  delta_prior   <-  rep(0,vardata$number_of_endogenous_variables)
  
  for (i in 1:vardata$number_of_endogenous_variables) {
    #Run a AR(1) and save the std of the residual for every endogenous variable
    xtemp          <-  dplyr::lag(Y[, i, drop = F])
    xtemp          <-  as.matrix(cbind("cons" = rep(1, nrow(xtemp)), xtemp)[-1,])
    ytemp          <-  as.matrix(Y[-c(1), i, drop = F])
    bols           <-  solve(crossprod(xtemp)) %*% crossprod(xtemp, ytemp)
    restemp        <-  ytemp - xtemp %*% bols
    if (bols[2] > 1){
      bols[2] <- 1
    }
    delta_prior[i] <-  bols[2]
    s_prior[i]     <-  sqrt(drop(crossprod(restemp)/(nrow(ytemp))))
  }
  
  deltaP = delta_prior
  sigmaP = s_prior
  
  s_prior       <- NULL#rep(0,vardata$number_of_endogenous*nexo_cont)
  delta_prior   <- NULL#  rep(0,vardata$number_of_endogenous*nexo_cont)
 
  y_dummy1 = x_dummy1 = y_dummy2 = x_dummy2 = NULL
 
   if (lamdaP > 0) {
    
    #Roughly speaking, the first block of dummies imposes prior beliefs on the
    # autoregressive coefficients.
    ar_prior  <- rbind( as.matrix(diag(sigmaP*deltaP)/lamdaP),
                        matrix(0,length(sigmaP)*(L-1), length(sigmaP)) )
    
    # the second block implements the prior for the covariance matrix and the
    cov_prior <- diag(sigmaP)
    
    # third block reflects the uninformative prior for the intercept
    # (is a very small number).
    det_prior <- NULL
    if (epsilonP >0) { # constant
      det_prior <- matrix(0,det, length(sigmaP))
    }
    
    exo_pr_cont_y <- NULL
    if (epsilonXP > 0){
      exo_pr_cont_y <- matrix(0,nexo_cont, length(sigmaP))
    }
    
    exo_pr_lag_y <-  NULL
    if (epsilonXLP > 0) {
      exo_pr_lag_y <- matrix(0,nexo_lag, length(sigmaP))
    }
    
    # the uninformative prior for the intercept (is a very small number)
    y_dummy1 <- rbind(det_prior,
                      exo_pr_cont_y,
                      exo_pr_lag_y,
                      ar_prior,
                      cov_prior)
    
    jp <- diag(1:L)
    
    # eq.5 2008(ecb) left matrix
    # autoregressive coefficients.
    
    ar_prior_x    <- cbind( matrix(0,(n*L),nexo) ,  jp %x% as.matrix(diag(sigmaP)/lamdaP) )
    #dim(ar_prior_x) # 33 23
    
    # covariance matrix
    cov_prior_x   <-  matrix(0, n, ((n*L) + nexo))
    #dim(cov_prior_x)
    
    # third block reflects the uninformative prior for the intercept
    det_prior_x <- NULL
    if (epsilonP > 0) {
      det_prior_x <- cbind(diag(det)*epsilonP, matrix(0, det, (n*L)+(nexo-det)))
    }
    
    exo_pr_cont_x <- NULL
    if (epsilonXP > 0) {
      exo_pr_cont_x <- cbind(matrix(0, nexo_cont, det),diag(nexo_cont) * epsilonXP, matrix(0, nexo_cont, (n * L) + (nexo - det - nexo_cont)))
    }
    
    exo_pr_lag_x <- NULL
    if (epsilonXLP > 0) {
      exo_pr_lag_x <- cbind(matrix(0, nexo_lag, det + nexo_cont),diag(nexo_lag) * epsilonXLP,matrix(0, nexo_lag, n * L) )
    }
    
    
    x_dummy1 = rbind(det_prior_x,
                     exo_pr_cont_x,
                     exo_pr_lag_x,
                     ar_prior_x,
                     cov_prior_x)
    
    
    
    #dim(x_dummy1);dim(y_dummy1)
  }
  
  
  
  # Get additional dummy matrices - see equation (9) of Banbura et al. 2007:
  
  if (tauP > 0) {
    y_dummy2 <- diag(deltaP * muP) / tauP
    x_dummy2 <- cbind( matrix(0, n, nexo), matrix(1, 1, L) %x% y_dummy2 )
  }
  
  y_dummy = rbind(y_dummy1, y_dummy2)
  x_dummy = rbind(x_dummy1, x_dummy2)
  
  # yd and xd are the dummy data. Append this to actual data
  Y0=rbind(Y,y_dummy)
  X0=rbind(X,x_dummy)
  
  dummy_obs = list("Y_dummy" = Y0, "X_dummy" = X0)
  return(dummy_obs)
  
}
