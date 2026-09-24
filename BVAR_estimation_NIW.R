BVAR_estimation_NIW <- function(data, 
                               lags = 4, 
                               const = TRUE, trend = FALSE, trend_qua = FALSE,
                               star = NULL, xlags = 0, det_controls = NULL,
                               ex = NULL, ex_lag = NULL, feedback = F,
                               reps = 1000, burn = 500,
                               lamda = 1, tau = 0, epsilon = 0.1 , epsilonxl=0.5, epsilonx = 0.5,
                               maxtries = 20000) {
  
  
  vardata <- vardata_form(
    data          = data,
    lags          = lags,
    const         = const,
    trend         = trend,
    trend_qua     = trend_qua,
    star          = star,
    xlags         = xlags,
    det_controls  = det_controls,
    ex            = ex,
    ex_lag        = ex_lag,
    feedback_only = feedback
  )
  

  # Extract model dimensions from the vardata structure
  n      <- vardata$number_of_endogenous_variables
  nexo   <- vardata$number_of_non_endogenous_regressors
  L      <- vardata$number_of_lags
  n_obs  <- vardata$number_of_observations
  
  # Number of retained posterior draws
  if (burn < 0L || reps <= burn) {
    stop("`reps` must be greater than `burn`, and `burn` cannot be negative.")
  }
  
  draws <- reps - burn
  t_lag = vardata$number_of_effective_observations
    
  # Construct the dummy-observation dataset
  dum <- dummy_obs_prior(
    vardata   = vardata,
    lamda     = lamda,
    tau       = tau,
    epsilon   = epsilon,
    epsilonx  = epsilonx,
    epsilonxl = epsilonxl
  )
  
  
  # conditional mean of the VAR coefficients
  X0 = dum$X_dummy; Y0 = dum$Y_dummy
  
  mstar= matrix(ols(x = X0,y = Y0))  #ols on the appended data
  ixx = chol2inv(chol(crossprod(X0))) # inv(X0'X0) to be used later in the Gibbs sampling algorithm
  
  sigma= diag(n); #starting value for sigma
  
  CM = array(0, list(n * L , n * L, draws))
  Sigma = array(0, list(n, n, draws))
  res = array(0, list(t_lag, n, draws))
  fit = array(0, list(t_lag, n, draws))
  beta_s = array(0, list((n*L + nexo),n, draws), dimnames=list( colnames(vardata$x_rhs) , vardata$names_of_endog_variables , NULL ))
  
  stab = 0
  i=1
  
  while (i <= reps) {
    
    if ((i %% 50) == 0) {
      print(paste0("Iteration ", i, " of ", reps, "."))
    }
    
    mxtries = 0
    stable_draw = FALSE
    
    while (!stable_draw && mxtries < maxtries) {
      
      # Draw beta conditional on the current sigma
      vstar = sigma %x% ixx
      beta = mstar + t(matrix(rnorm(n * (n * L + nexo)),nrow = 1) %*% chol(vstar))
      beta_pos = matrix(beta,nrow = n * L + nexo,ncol = n)
      
      # Construct the companion matrix
      CM_pot = rbind(t(beta_pos[(1 + nexo):nrow(beta_pos), , drop = FALSE]), cbind(diag(n * L - n),matrix(0, n * L - n, n)))
      eig = eigen(CM_pot, only.values = TRUE)$values
      
      if (max(Mod(eig)) < 1) {
        stable_draw = TRUE
      } else {
        stab = stab + 1
        mxtries = mxtries + 1
      }
    }
    
    if (!stable_draw) {
      stop("Could not obtain a stable draw within maxtries attempts.")
    }
    
    # Draw sigma only after accepting a stable beta
    e = Y0 - X0 %*% beta_pos
    scale = crossprod(e)
    
    sigma = MCMCpack::riwish(nrow(Y0),scale)
    
    # Retain draws after burn-in
    if (i > burn) {
      
      draw_number = i - burn
      beta_s[, , draw_number] = beta_pos
      CM[, , draw_number] = CM_pot
      Sigma[, , draw_number] = sigma
      fit[, , draw_number] = vardata$x_rhs %*% beta_pos
      res[, , draw_number] = vardata$y_lhs - fit[, , draw_number]
    }
    i = i + 1
  }
  
  arguments = list("tau" = tau, "lamda" = lamda, "epsilon" =  epsilon)
  bvar = list(beta_s, Sigma, CM, res, fit, stab, vardata, draws, arguments)
  names(bvar) = c("posterior","Sigma", "CM", "res", "yfit","unstable_draws", "vardata", "draws","priors")
  return(bvar)
  
  
}

