## Estimate phenometrics using a standard Von Mises distribution.
## Note: Not used in Baiotto and Cosma (2026)

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

# Function to fit von Mises distribution and estimate parameters
fit_von_mises <- function(theta_data) {
  # Convert to circular object
  circ_data <- circular(theta_data)
  
  # Estimate von Mises parameters using maximum likelihood
  # mu: mean direction, kappa: concentration parameter
  mu <- mean(circ_data)
  kappa <- est.kappa(circ_data)
  
  return(list(mu = as.numeric(mu), kappa = kappa))
}

# Bootstrap function for quantile estimation
bootstrap_circular_quantiles <- function(theta_data, quantiles = c(0.05, 0.5, 0.95), 
                                         n_boot = 1000, year_length = 365, 
                                         species_name = NULL, include_ci = FALSE,
                                         ci_level = 0.95) {
  
  # Function for bootstrap sampling
  boot_quantiles <- function(data, indices) {
    # Resample the original data
    boot_sample_data <- data[indices]
    
    # Fit von Mises distribution to bootstrap sample
    vm_params <- fit_von_mises(boot_sample_data)
    
    # Calculate quantiles directly from the fitted von Mises distribution
    quants <- qvonmises(quantiles, mu = vm_params$mu, kappa = vm_params$kappa)
    
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
      ci_lower <- radians_to_doy(ci$percent[4], year_length)  # Lower bound
      ci_upper <- radians_to_doy(ci$percent[5], year_length)  # Upper bound
      
      # Add CI to results
      result_list[[paste0(q_name, "_ci_lower")]] <- ci_lower
      result_list[[paste0(q_name, "_ci_upper")]] <- ci_upper
    }
  }
  
  return(result_list)
}

phenol_wrapper_circular <- function(sp_data, sp_name, doy_col = "DOY", quantiles = c(0.05, 0.95), iterations = 1000, include_ci = FALSE) {
  
  # Convert DOY to theta
  sp_theta <- doy_to_radians(sp_data[[doy_col]])
  
  # Return named NA structure if not more than 15 phenol observations for the species
  if (length(sp_theta) <= 15) {
    # Create named list with NA values for each quantile
    result <- list(species = sp_name)
    for (q in quantiles) {
      q_name <- paste0("q", sprintf("%02d", q * 100))
      result[[q_name]] <- NA
    }
    return(result)
  }
  
  # Nonparametric bootstrap quantile estimates
  boot_quantiles <- bootstrap_circular_quantiles(sp_theta, 
                                                 species_name = sp_name, 
                                                 quantiles = quantiles, 
                                                 n_boot = iterations, 
                                                 year_length = 365,
                                                 include_ci = include_ci)
  
  return(boot_quantiles)
}
