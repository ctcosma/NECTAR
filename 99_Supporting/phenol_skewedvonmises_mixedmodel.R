## For extending the skewed Von Mises approach to a mixed model framework (which could be especially useful for identifying active periods of multi-voltine species)
## Relies on characterization of skewed Von Mises distribution by Ahsanullah1 and Anis (2019) in https://www.aligarhjournalstatistics.com/issues/ajs-v-39-2019/2-mza-ma.pdf
## Note: Not used in Baiotto and Cosma (2026), and not thoroughly tested

# PDF from equation (1.3) - single component
pdf_skewed_vm <- function(theta, lambda, kappa, mu = 0) {
  if (lambda < -1 || lambda > 1) stop("lambda must be between -1 and 1")
  if (kappa < 0) stop("kappa must be >= 0")
  
  # Center theta around mu and wrap to [-π, π]
  theta_centered <- theta - mu
  theta_centered <- atan2(sin(theta_centered), cos(theta_centered))
  
  # Check bounds on centered theta
  if (any(theta_centered < -pi - 1e-10 | theta_centered > pi + 1e-10)) {
    stop("theta must be between -pi and pi")
  }
  
  I0_k <- besselI(kappa, nu = 0)
  pdf_val <- (exp(kappa * cos(theta_centered)) / (2 * pi * I0_k)) * (1 + lambda * sin(theta_centered))
  return(pdf_val)
}

# PDF for mixture of skewed von Mises distributions
pdf_mixture_skewed_vm <- function(theta, params_list, pi_weights) {
  # params_list: list of lists, each containing mu, kappa, lambda
  # pi_weights: vector of mixture weights (must sum to 1)
  
  if (abs(sum(pi_weights) - 1) > 1e-6) stop("pi_weights must sum to 1")
  if (any(pi_weights < 0)) stop("pi_weights must be non-negative")
  
  n_components <- length(params_list)
  pdf_val <- 0
  
  for (i in 1:n_components) {
    # Add weighted component PDF
    pdf_val <- pdf_val + pi_weights[i] * pdf_skewed_vm(theta, 
                                                       params_list[[i]]$lambda, 
                                                       params_list[[i]]$kappa,
                                                       params_list[[i]]$mu)
  }
  
  return(pdf_val)
}

# Fit mixture of skewed von Mises distributions using EM algorithm
fit_mixture_skewed_vm <- function(theta_data, n_components = 1, max_iter = 100, 
                                  tol = 1e-6, n_init = 5, control = list()) {
  theta_data <- as.numeric(theta_data)
  n <- length(theta_data)
  
  # Number of parameters: each component has 3 params (mu, kappa, lambda)
  # plus (n_components - 1) mixing weights
  n_params <- 3 * n_components + (n_components - 1)
  
  # Special case: single component
  if (n_components == 1) {
    single_fit <- fit_skewed_von_mises_single(theta_data, control)
    return(list(
      n_components = 1,
      params = list(list(mu = single_fit$mu, kappa = single_fit$kappa, lambda = single_fit$lambda)),
      pi_weights = 1,
      logLik = single_fit$logLik,
      n_params = 3,
      AIC = -2 * single_fit$logLik + 2 * 3,
      BIC = -2 * single_fit$logLik + log(n) * 3,
      convergence = single_fit$convergence
    ))
  }
  
  # Try multiple random initializations
  best_fit <- NULL
  best_loglik <- -Inf
  
  for (init in 1:n_init) {
    # Initialize parameters
    init_params <- initialize_mixture_params(theta_data, n_components)
    
    # Run EM algorithm
    fit <- tryCatch({
      em_mixture_skewed_vm(theta_data, init_params, max_iter, tol, control)
    }, error = function(e) {
      message(paste("Initialization", init, "failed:", e$message))
      NULL
    })
    
    if (!is.null(fit) && fit$logLik > best_loglik) {
      best_loglik <- fit$logLik
      best_fit <- fit
      
      # Early stopping if we have a good fit
      if (init >= 3 && best_fit$convergence == 0) {
        message(paste("Early stopping at initialization", init, "- good fit found"))
        break
      }
    }
  }
  
  if (is.null(best_fit)) {
    stop("EM algorithm failed for all initializations")
  }
  
  # Calculate information criteria
  best_fit$n_params <- n_params
  best_fit$AIC <- -2 * best_fit$logLik + 2 * n_params
  best_fit$BIC <- -2 * best_fit$logLik + log(n) * n_params
  
  return(best_fit)
}

# Initialize mixture parameters
initialize_mixture_params <- function(theta_data, n_components) {
  n <- length(theta_data)
  
  # Use k-means-like clustering on circular data
  # Randomly assign initial cluster memberships
  cluster_assign <- sample(1:n_components, n, replace = TRUE)
  
  params_list <- list()
  pi_weights <- rep(1/n_components, n_components)
  
  for (k in 1:n_components) {
    cluster_data <- theta_data[cluster_assign == k]
    
    if (length(cluster_data) < 3) {
      # If cluster too small, use random initialization
      params_list[[k]] <- list(
        mu = runif(1, -pi, pi),
        kappa = rexp(1, rate = 0.5),
        lambda = runif(1, -0.5, 0.5)
      )
    } else {
      # Estimate from cluster data
      circ_data <- circular(cluster_data)
      mu_k <- as.numeric(mean(circ_data))
      kappa_k <- est.kappa(circ_data)
      
      params_list[[k]] <- list(
        mu = mu_k,
        kappa = max(kappa_k, 0.1),
        lambda = runif(1, -0.3, 0.3)
      )
    }
    
    pi_weights[k] <- sum(cluster_assign == k) / n
  }
  
  list(params = params_list, pi_weights = pi_weights)
}

# EM algorithm for mixture model
em_mixture_skewed_vm <- function(theta_data, init_params, max_iter = 50, 
                                 tol = 1e-4, control = list()) {
  n <- length(theta_data)
  n_components <- length(init_params$params)
  
  params_list <- init_params$params
  pi_weights <- init_params$pi_weights
  
  loglik_old <- -Inf
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    # E-step: Calculate responsibilities
    responsibilities <- matrix(0, nrow = n, ncol = n_components)
    
    for (k in 1:n_components) {
      responsibilities[, k] <- pi_weights[k] * pdf_skewed_vm(
        theta_data, 
        params_list[[k]]$lambda,
        params_list[[k]]$kappa,
        params_list[[k]]$mu
      )
    }
    
    # Normalize responsibilities (with safeguard)
    resp_sum <- rowSums(responsibilities)
    resp_sum[resp_sum < 1e-300] <- 1e-300
    responsibilities <- responsibilities / resp_sum
    
    # Calculate log-likelihood
    loglik <- sum(log(resp_sum))
    
    # Check for improvement
    if (iter > 1 && (loglik < loglik_old - 1e-6)) {
      # Likelihood decreased - stop and return previous iteration
      message(paste("Likelihood decreased at iteration", iter, "- stopping"))
      converged <- TRUE
      break
    }
    
    # Check convergence
    if (abs(loglik - loglik_old) < tol) {
      converged <- TRUE
      break
    }
    
    loglik_old <- loglik
    
    # M-step: Update parameters
    for (k in 1:n_components) {
      # Update mixing weight
      n_k <- sum(responsibilities[, k])
      pi_weights[k] <- n_k / n
      
      # Skip update if component has very few points
      if (n_k < 3) {
        next
      }
      
      # Update component parameters using weighted MLE
      params_list[[k]] <- fit_weighted_skewed_vm(theta_data, responsibilities[, k], control)
    }
  }
  
  # Return results
  return(list(
    n_components = n_components,
    params = params_list,
    pi_weights = pi_weights,
    logLik = loglik_old,
    convergence = ifelse(converged, 0, 1),
    iterations = iter
  ))
}

# Fit weighted skewed von Mises (for M-step)
fit_weighted_skewed_vm <- function(theta_data, weights, control = list()) {
  theta_data <- as.numeric(theta_data)
  
  # Weighted mean direction
  weighted_sum_sin <- sum(weights * sin(theta_data))
  weighted_sum_cos <- sum(weights * cos(theta_data))
  mu0 <- atan2(weighted_sum_sin, weighted_sum_cos)
  
  # Initial kappa estimate
  R <- sqrt(weighted_sum_sin^2 + weighted_sum_cos^2) / sum(weights)
  kappa0 <- A1inv(R)
  
  # Negative weighted log-likelihood
  neg_loglik <- function(p) {
    mu     <- p[1]
    kappa  <- exp(p[2])
    lambda <- tanh(p[3])
    
    # Center data around mu
    u <- theta_data - mu
    u <- atan2(sin(u), cos(u))
    
    # Calculate weighted likelihood
    pdf_vals <- pdf_skewed_vm(u, lambda, kappa, mu = 0)
    pdf_vals[pdf_vals <= 1e-300] <- 1e-300
    
    return(-sum(weights * log(pdf_vals)))
  }
  
  # Try fewer lambda starting values for speed
  lambda_starts <- c(-0.3, 0, 0.3)
  
  best_opt <- NULL
  best_val <- Inf
  
  for (lambda_init in lambda_starts) {
    par0 <- c(mu0, log(max(kappa0, 0.1)), atanh(lambda_init))
    
    opt <- tryCatch({
      optim(
        par = par0,
        fn = neg_loglik,
        method = "BFGS",
        control = modifyList(list(reltol = 1e-6, maxit = 500), control)
      )
    }, error = function(e) NULL)
    
    if (!is.null(opt) && opt$value < best_val) {
      best_val <- opt$value
      best_opt <- opt
    }
  }
  
  if (is.null(best_opt)) {
    # Fallback to initial estimates
    return(list(mu = mu0, kappa = max(kappa0, 0.1), lambda = 0))
  }
  
  # Back-transform
  mu_hat <- atan2(sin(best_opt$par[1]), cos(best_opt$par[1]))
  kappa_hat <- exp(best_opt$par[2])
  lambda_hat <- tanh(best_opt$par[3])
  
  list(mu = mu_hat, kappa = kappa_hat, lambda = lambda_hat)
}

# Helper function for kappa estimation (inverse of A_1(kappa))
A1inv <- function(R) {
  if (R < 0.53) {
    return(2 * R + R^3 + (5 * R^5)/6)
  } else if (R < 0.85) {
    return(-0.4 + 1.39 * R + 0.43/(1 - R))
  } else {
    return(1/(R^3 - 4 * R^2 + 3 * R))
  }
}

# Fit single skewed von Mises distribution (original function)
fit_skewed_von_mises_single <- function(theta_data, control = list()) {
  theta_data <- as.numeric(theta_data)
  n <- length(theta_data)
  
  # Get initial parameter estimates from standard von Mises
  circ_data <- circular(theta_data)
  mu0 <- as.numeric(mean(circ_data))
  kappa0 <- est.kappa(circ_data)
  
  # Negative log-likelihood function
  neg_loglik <- function(p) {
    mu     <- p[1]
    kappa  <- exp(p[2])
    lambda <- tanh(p[3])
    
    # Center data around mu and wrap
    u <- theta_data - mu
    u <- atan2(sin(u), cos(u))
    
    # Calculate likelihood
    pdf_vals <- pdf_skewed_vm(u, lambda, kappa, mu = 0)
    pdf_vals[pdf_vals <= 1e-300] <- 1e-300
    
    return(-sum(log(pdf_vals)))
  }
  
  # Gradient function
  grad_fun <- function(p) {
    mu     <- p[1]
    kappa  <- exp(p[2])
    lambda <- tanh(p[3])
    
    u <- theta_data - mu
    u <- atan2(sin(u), cos(u))
    s <- sin(u)
    c <- cos(u)
    
    I0_k <- besselI(kappa, 0)
    
    # Gradient with respect to original parameters
    d_mu     <- kappa * sum(s) - lambda * sum(c / (1 + lambda * s))
    d_kappa  <- sum(c) - n * besselI(kappa, 1) / I0_k
    d_lambda <- sum(s / (1 + lambda * s))
    
    # Chain rule for transformed parameters
    grad <- c(
      -d_mu,
      -kappa * d_kappa,
      -(1 - lambda^2) * d_lambda
    )
    return(grad)
  }
  
  # Try multiple starting values for lambda to avoid local maxima
  lambda_starts <- seq(-0.8, 0.8, by = 0.4)
  
  best_opt <- NULL
  best_loglik <- -Inf
  
  for (lambda_init in lambda_starts) {
    par0 <- c(mu0, log(kappa0), atanh(lambda_init))
    
    opt <- tryCatch({
      optim(
        par = par0,
        fn = neg_loglik,
        gr = grad_fun,
        method = "BFGS",
        control = modifyList(list(reltol = 1e-10, maxit = 2000), control)
      )
    }, error = function(e) NULL)
    
    # Keep the best result
    if (!is.null(opt) && -opt$value > best_loglik) {
      best_loglik <- -opt$value
      best_opt <- opt
    }
  }
  
  if (is.null(best_opt)) {
    stop("Optimization failed for all starting values")
  }
  
  # Back-transform final estimates
  mu_hat     <- atan2(sin(best_opt$par[1]), cos(best_opt$par[1]))
  kappa_hat  <- exp(best_opt$par[2])
  lambda_hat <- tanh(best_opt$par[3])
  
  list(
    mu = mu_hat,
    kappa = kappa_hat,
    lambda = lambda_hat,
    logLik = best_loglik,
    convergence = best_opt$convergence
  )
}

# Select best model across different numbers of components
select_best_mixture_model <- function(theta_data, max_n_components = 2, 
                                      criterion = "BIC", control = list(),
                                      verbose = FALSE) {
  # criterion can be "AIC" or "BIC"
  
  n <- length(theta_data)
  models <- list()
  criteria_values <- numeric(max_n_components)
  
  for (k in 1:max_n_components) {
    if (verbose) message(paste("Fitting model with", k, "component(s)..."))
    
    model <- tryCatch({
      fit_mixture_skewed_vm(theta_data, n_components = k, control = control)
    }, error = function(e) {
      if (verbose) message(paste("Failed to fit", k, "component model:", e$message))
      NULL
    })
    
    if (!is.null(model)) {
      models[[k]] <- model
      criteria_values[k] <- model[[criterion]]
      if (verbose) {
        message(paste("  ", criterion, "=", round(criteria_values[k], 2), 
                      "| LogLik =", round(model$logLik, 2),
                      "| Iterations =", model$iterations))
      }
    } else {
      models[[k]] <- NULL
      criteria_values[k] <- Inf
    }
  }
  
  # Select best model (minimum AIC/BIC)
  best_k <- which.min(criteria_values)
  
  if (is.infinite(criteria_values[best_k])) {
    stop("All models failed to fit")
  }
  
  if (verbose) {
    message(paste("Best model has", best_k, "component(s) with", criterion, "=", 
                  round(criteria_values[best_k], 2)))
  }
  
  list(
    best_model = models[[best_k]],
    all_models = models,
    criteria_values = criteria_values,
    criterion_used = criterion,
    best_n_components = best_k
  )
}

# Calculate activity periods based on probability density threshold
# Returns a binary vector of length year_length indicating active days
calculate_activity_periods <- function(model, year_length = 365, 
                                       density_threshold = 0.90) {
  # Calculate probability density for each day of the year
  doy_seq <- 1:year_length
  theta_seq <- doy_to_radians(doy_seq)
  
  # Ensure theta is in valid range [-π, π]
  theta_seq <- atan2(sin(theta_seq), cos(theta_seq))
  
  if (model$n_components == 1) {
    pdf_vals <- pdf_skewed_vm(theta_seq, 
                              model$params[[1]]$lambda,
                              model$params[[1]]$kappa,
                              model$params[[1]]$mu)
  } else {
    pdf_vals <- pdf_mixture_skewed_vm(theta_seq, model$params, model$pi_weights)
  }
  
  # Normalize to sum to 1 (discrete approximation)
  pdf_vals <- pdf_vals / sum(pdf_vals)
  
  # Sort days by probability density (descending)
  sorted_indices <- order(pdf_vals, decreasing = TRUE)
  sorted_pdf <- pdf_vals[sorted_indices]
  
  # Find cumulative sum threshold
  cumsum_pdf <- cumsum(sorted_pdf)
  n_active_days <- which(cumsum_pdf >= density_threshold)[1]
  
  # Create binary activity vector
  active_days <- rep(0, year_length)
  active_days[sorted_indices[1:n_active_days]] <- 1
  
  return(active_days)
}

# Calculate summary statistics from activity periods
summarize_activity_periods <- function(activity_vector) {
  # Find runs of 1s (active periods)
  rle_result <- rle(activity_vector)
  
  active_runs <- which(rle_result$values == 1)
  n_periods <- length(active_runs)
  
  if (n_periods == 0) {
    return(list(
      n_periods = 0,
      total_active_days = 0,
      periods = data.frame()
    ))
  }
  
  # Calculate start and end of each period
  run_ends <- cumsum(rle_result$lengths)
  run_starts <- c(1, run_ends[-length(run_ends)] + 1)
  
  periods_df <- data.frame(
    period_id = 1:n_periods,
    start_doy = run_starts[active_runs],
    end_doy = run_ends[active_runs],
    length = rle_result$lengths[active_runs]
  )
  
  list(
    n_periods = n_periods,
    total_active_days = sum(activity_vector),
    periods = periods_df
  )
}

# Wrapper function with mixture model selection
phenol_wrapper_mixture_skewed <- function(sp_data, 
                                          sp_name, 
                                          doy_col = "DOY", 
                                          max_n_components = 2,
                                          model_selection_criterion = "BIC",
                                          density_threshold = 0.90,
                                          year_length = 365,
                                          min_obs = 10,
                                          verbose = FALSE) {
  
  # Convert DOY to theta
  sp_theta <- doy_to_radians(sp_data[[doy_col]])
  
  # Return NA structure if not enough observations
  if (length(sp_theta) <= min_obs) {
    result <- list(
      species = sp_name,
      n_components = NA,
      activity_vector = rep(NA, year_length),
      n_periods = NA,
      total_active_days = NA
    )
    return(result)
  }
  
  # Select best mixture model
  model_selection <- tryCatch({
    select_best_mixture_model(sp_theta, 
                              max_n_components = max_n_components,
                              criterion = model_selection_criterion,
                              verbose = verbose)
  }, error = function(e) {
    if (verbose) message(paste("Model selection failed for", sp_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(model_selection)) {
    result <- list(
      species = sp_name,
      n_components = NA,
      activity_vector = rep(NA, year_length),
      n_periods = NA,
      total_active_days = NA
    )
    return(result)
  }
  
  best_model <- model_selection$best_model
  
  # Calculate activity periods
  activity_vector <- calculate_activity_periods(best_model, 
                                                year_length = year_length,
                                                density_threshold = density_threshold)
  
  # Summarize periods
  period_summary <- summarize_activity_periods(activity_vector)
  
  # Compile results
  result <- list(
    species = sp_name,
    n_components = best_model$n_components,
    n_observations = length(sp_theta),
    model_logLik = best_model$logLik,
    model_AIC = best_model$AIC,
    model_BIC = best_model$BIC,
    activity_vector = activity_vector,
    n_periods = period_summary$n_periods,
    total_active_days = period_summary$total_active_days,
    period_details = period_summary$periods,
    model_params = best_model$params,
    pi_weights = best_model$pi_weights
  )
  
  return(result)
}

# Helper function: convert DOY to radians (centered at 0, range [-π, π])
doy_to_radians <- function(doy, year_length = 365) {
  # Convert to [0, 2π] then shift to [-π, π]
  theta <- (doy - 1) / year_length * 2 * pi
  # Shift to center at 0: [0, 2π] -> [-π, π]
  theta <- theta - pi
  return(theta)
}

# Helper function: convert radians to DOY
radians_to_doy <- function(theta, year_length = 365) {
  # Shift from [-π, π] back to [0, 2π]
  theta_shifted <- theta + pi
  # Convert to DOY
  doy <- (theta_shifted / (2 * pi)) * year_length + 1
  # Wrap to [1, year_length]
  doy <- ((doy - 1) %% year_length) + 1
  return(doy)
}
