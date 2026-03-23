#### Circular pheno functions ####
# Function to convert day of year to radians
doy_to_radians <- function(doy, year_length = 365) {
  theta <- 2 * pi * (doy - 1) / year_length
  return(theta)
}

# Function to convert radians back to day of year
radians_to_doy <- function(theta, year_length = 365) {
  theta <- theta %% (2 * pi)
  doy <- (theta * year_length) / (2 * pi) + 1
  return(doy)
}

# PDF from equation (1.3)
pdf_skewed_vm <- function(theta, lambda, kappa) {
  if (lambda < -0.99 || lambda > 0.99) stop("lambda must be between -1 and 1")
  if (kappa < 0) stop("kappa must be >= 0")
  if (any(theta < -pi | theta > pi)) stop("theta must be between -pi and pi")
  
  # Handle extreme kappa values that cause numerical issues
  if (kappa > 500) {
    warning("kappa too large, capping at 500")
    kappa <- 500
  }
  
  I0_k <- besselI(kappa, nu = 0)
  pdf_val <- (exp(kappa * cos(theta)) / (2 * pi * I0_k)) * (1 + lambda * sin(theta))
  
  # Safety check
  if (any(pdf_val <= 0)) {
    warning("Non-positive PDF values detected")
    pdf_val[pdf_val <= 0] <- 1e-300
  }
  
  return(pdf_val)
}

# CDF from equation (2.3)
## NOTE - kappa from term 2 omitted in paper's equation
# Pre-compute Bessel values once
cdf_skewed_vm <- function(theta, kappa, lambda, bessel_cache = NULL) {
  if (is.null(bessel_cache)) {
    n_terms <- 1000
    j_seq <- 1:n_terms
    I0_k <- besselI(kappa, 0)
    I_j <- besselI(kappa, nu = j_seq)
    bessel_cache <- list(I0_k = I0_k, I_j = I_j, j_seq = j_seq)
  }
  
  # Vectorized for multiple theta
  theta <- as.vector(theta)
  I0_k <- bessel_cache$I0_k
  I_j <- bessel_cache$I_j
  j_seq <- bessel_cache$j_seq
  
  # Vectorized bessel sum
  bessel_sum <- sapply(theta, function(t) {
    sum((I_j / I0_k) * (sin(j_seq * t) / j_seq))
  })
  
  term1 <- (1 / (2 * pi)) * ((pi + theta) + 2 * bessel_sum)
  term2 <- (lambda / (2 * pi * kappa * I0_k)) * (exp(-kappa) - exp(kappa * cos(theta)))
  
  cdf_val <- term1 + term2
  pmax(0, pmin(1, cdf_val))
}

fit_skewed_von_mises <- function(theta_data, control = list()) {
  tryCatch({
    theta_data <- as.numeric(theta_data)
    n <- length(theta_data)
    
    # Get initial parameter estimates from standard von Mises
    circ_data <- circular::circular(theta_data)
    mu0 <- as.numeric(circular::mean(circ_data))
    kappa0 <- circular::est.kappa(circ_data)
    
    # Cap kappa at 500
    kappa0 <- min(kappa0, 500)
    
    # Validate initial parameters
    if (!is.finite(kappa0) || kappa0 <= 0) {
      return(NULL)
    }
    
    if (!is.finite(mu0)) {
      return(NULL)
    }
    
    # Negative log-likelihood function
    neg_loglik <- function(p) {
      mu     <- p[1]
      kappa  <- exp(p[2])
      kappa  <- min(kappa, 500)  # Cap kappa to match pdf_skewed_vm
      lambda <- 0.99 * tanh(p[3])  # Constrain to (-0.99, 0.99) to avoid negative probability density
      
      # Center data around mu and wrap
      u <- theta_data - mu
      u <- atan2(sin(u), cos(u))
      
      # Calculate likelihood
      pdf_vals <- pdf_skewed_vm(u, lambda, kappa)
      
      # CHANGE: Check for invalid PDF values before taking log
      if (any(!is.finite(pdf_vals)) || any(pdf_vals <= 0)) {
        return(1e10)
      }
      
      pdf_vals[pdf_vals <= 1e-300] <- 1e-300
      
      return(-sum(log(pdf_vals)))
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
          gr = NULL,
          method = "Nelder-Mead",
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
      return(NULL)
    }
    
    # Back-transform final estimates
    mu_hat     <- atan2(sin(best_opt$par[1]), cos(best_opt$par[1]))
    kappa_hat  <- exp(best_opt$par[2])
    kappa_hat  <- min(kappa_hat, 500)  # Cap final estimate
    lambda_hat <- 0.99 * tanh(best_opt$par[3])  # Match constraint above
    
    list(
      mu = mu_hat,
      kappa = kappa_hat,
      lambda = lambda_hat,
      logLik = best_loglik,
      convergence = best_opt$convergence
    )
  }, error = function(e) {
    return(NULL)
  })
}

# Quantile function for sine skewed von Mises distribution
qskewed_vonmises <- function(p, mu, kappa, lambda, n_terms = 1000) {
  # Pre-compute Bessel values once
  j_seq <- 1:n_terms
  I0_k <- besselI(kappa, 0)
  I_j <- besselI(kappa, nu = j_seq)
  bessel_cache <- list(I0_k = I0_k, I_j = I_j, j_seq = j_seq)
  
  find_quantile <- function(prob) {
    f <- function(theta) {
      cdf_skewed_vm(theta, kappa, lambda, bessel_cache) - prob
    }
    result <- uniroot(f, lower = -pi, upper = pi, tol = 1e-6)
    result$root
  }
  
  quantiles <- sapply(p, find_quantile)
  (quantiles + mu) %% (2 * pi)
}

# Updated bootstrap function for sine skewed von Mises
bootstrap_circular_quantiles_skewed <- function(theta_data, quantiles = c(0.05, 0.5, 0.95), 
                                                n_boot = 100, year_length = 365, 
                                                species_name = NULL, include_ci = FALSE,
                                                ci_level = 0.95, n_terms = 1000) {
  
  # Function for bootstrap sampling
  boot_quantiles <- function(data, indices) {
    # Resample the original data
    boot_sample_data <- data[indices]
    
    # Fit sine skewed von Mises distribution to bootstrap sample
    svm_params <- fit_skewed_von_mises(boot_sample_data)
    
    # Calculate quantiles from the fitted distribution
    quants <- qskewed_vonmises(quantiles, 
                               mu = svm_params$mu, 
                               kappa = svm_params$kappa, 
                               lambda = svm_params$lambda,
                               n_terms = n_terms)
    
    return(as.numeric(quants))
  }
  
  # Perform bootstrap
  boot_results <- boot::boot(data = theta_data, statistic = boot_quantiles, R = n_boot)
  
  # Create named list of results
  result_list <- list()
  
  # Add species name if provided
  if (!is.null(species_name)) {
    result_list$species <- species_name
  }
  
  # Add quantile estimates and optionally CIs
  for (i in 1:length(quantiles)) {
    q_name <- paste0("q", sprintf("%02d", quantiles[i] * 100))
    
    # Point estimate
    result_list[[q_name]] <- radians_to_doy(boot_results$t0[i], year_length)
    
    # Add confidence intervals if requested
    if (include_ci) {
      # Calculate CI using percentile method
      ci <- boot::boot.ci(boot_results, index = i, type = "perc", conf = ci_level)
      
      # Extract CI bounds and convert to DOY
      ci_lower <- radians_to_doy(ci$percent[4], year_length)
      ci_upper <- radians_to_doy(ci$percent[5], year_length)
      
      # Add CI to results
      result_list[[paste0(q_name, "_ci_lower")]] <- ci_lower
      result_list[[paste0(q_name, "_ci_upper")]] <- ci_upper
    }
  }
  
  return(result_list)
}
