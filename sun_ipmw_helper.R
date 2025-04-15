

#' Scale a variable ignoring missing values
#' @param data Dataframe
#' @param x Column name to scale
#' @return Scaled vector
scale_rmna <- function(data, x) {
  mean <- mean(data[[x]], na.rm = TRUE)
  sd <- sd(data[[x]], na.rm = TRUE)
  (data[[x]] - mean) / sd
}

#' Scale all variables in a dataset
#' @param data Dataframe with variables to scale
#' @param vars Vector of variable names to scale
#' @return Dataframe with scaled variables
scale_dataset <- function(data, vars) {
  result <- data
  for (v in vars) {
    result[[v]] <- scale_rmna(data, v)
  }
  return(result)
}

#' Create pattern indicators
#' @param data Dataframe
#' @param vars Variables to check for missingness
#' @return Dataframe with pattern indicators
create_pattern_indicators <- function(data, vars) {
  
  # Step 1: Create a subset of the data containing only the specified variables
  subset_data <- data[, vars, drop = FALSE]
  
  # Step 2: Create a binary matrix where TRUE means the value is missing (NA)
  is_missing <- is.na(subset_data)
  
  # Step 3: Convert each row to a unique pattern identifier
  # - First convert each row to a string representation (e.g., "TFFT")
  # - This creates a factor with one level per unique missingness pattern
  pattern_strings <- apply(is_missing, 1, function(row) paste(as.integer(row), collapse = ""))
  
  # Step 4: Convert the factor to numeric values starting from 1
  # - Pattern 1 means no missing values (all FALSE)
  # - Other patterns get assigned 2, 3, 4, etc.
  
  # First create a factor with the pattern strings
  pattern_factor <- factor(pattern_strings)
  
  # Get the "no missingness" pattern (all zeros)
  no_missing_pattern <- paste(rep(0, length(vars)), collapse = "")
  
  # Reorder the levels so the "no missingness" pattern is first
  if(no_missing_pattern %in% levels(pattern_factor)) {
    new_levels <- c(no_missing_pattern, 
                    levels(pattern_factor)[levels(pattern_factor) != no_missing_pattern])
    pattern_factor <- factor(pattern_factor, levels = new_levels)
  }
  
  # Convert to numeric
  pattern_numeric <- as.numeric(pattern_factor)
  
  # Step 5: Add the pattern indicator to the original dataset
  data$M <- pattern_numeric
  
  # Return the modified dataset
  return(data)
}

#' Prepare data for JAGS
#' @param data Dataframe with M pattern indicator
#' @param vars Variables in model
#' @return List formatted for JAGS
prepare_jags_data <- function(data, vars) {
  # Sort so complete cases (M=1) are first
  sorted <- data %>% arrange(M)
  
  # Create list for JAGS
  dat <- list(R = sorted$M)
  
  # Create matrix of variables, replacing NA with 0 (will be ignored in model)
  dat$L <- as.matrix(sorted[, vars])
  #dat$L[is.na(dat$L)] <- 0  # Replace NA with 0 instead of -9999
  #MM: changed back
  dat$L[is.na(dat$L)] <- -9999
  
  dat$N <- nrow(sorted)
  dat$f <- rep(1, dat$N)
  dat$Nc <- sum(sorted$M == 1)
  dat$onesc <- rep(1, dat$Nc)
  dat$c <- 10^-8
  
  return(dat)
}

#' Generate random initial values
#' @param pattern_counts Vector with counts of variables observed in each pattern
#' @return Random initial values for JAGS
generate_initial_values <- function(pattern_counts) {
  ints <- runif(length(pattern_counts), min = -4, max = -1)
  lim <- 0.06
  
  result <- c()
  for (i in 1:length(pattern_counts)) {
    result <- c(result, ints[i], runif(pattern_counts[i], min = -lim, max = lim))
  }
  
  return(result)
}



#' Create JAGS model for missingness patterns
#' @param num_patterns Number of patterns
#' @param vars_per_pattern List with variables available for each pattern
#' @return JAGS model as text
create_jags_model <- function(num_patterns, vars_per_pattern) {
  # Build JAGS model text with better handling of missing values
  model_text <- "
  model {
    for(i in 1:N) {
      f[i] ~ dbern(pmiss[i,R[i]]) # f = 1 for all obs 
      "
  
  # Add pattern probability equations
  param_idx <- 1
  for (p in 2:num_patterns) {
    model_text <- paste0(model_text, "\n      logit(pmiss[i, ", p, "]) <- g[", param_idx, "]")
    param_idx <- param_idx + 1
    
    # Only include variables if they're observed in this pattern
    if (length(vars_per_pattern[[p]]) > 0) {
      # Add coefficients for observed variables
      for (v_idx in vars_per_pattern[[p]]) {
        # Add a check for missing values (though we've replaced them with 0)
        model_text <- paste0(model_text, " + g[", param_idx, "]*L[i,", v_idx, "]")
        param_idx <- param_idx + 1
      }
    }
    # If no variables are observed, we've already added the intercept only
  }
  
  # Add calculation for pattern 1 probability
  model_text <- paste0( model_text, "
      
      # Calculate pattern 1 probability
      pmiss[i,1] <- 1" )
  
  # Subtract all other pattern probabilities
  for (p in 2:num_patterns) {
    model_text <- paste0(model_text, " - pmiss[i,", p, "]")
  }
  
  model_text <- paste0(model_text, "
    }
    
    # constraint
    for (j in 1:Nc) {
      onesc[j] ~ dbern(C[j]) 
      C[j] <- step(pmiss[j,1]-c)
    }
    
    # priors
    for(k in 1:", param_idx - 1, ") {
      g[k] ~ dnorm(0,1/100) # diffuse prior
    }
  }")
  
  return(model_text)
}

#' Run Bayesian estimation for missingness model 
#' @param data Dataset with pattern indicators
#' @param vars Variables for model
#' @return JAGS fit object
run_missingness_model <- function(data, vars) {
  
  vars = sort(vars)
  
  if (!"M" %in% names(data)) {
    working_data <- create_pattern_indicators(data, vars)
  } else {
    working_data = data
  }
  
  # # Make a copy to avoid modifying the original - save for future use in pkg
  # working_data <- data
  
  # Ensure all variables exist, adding dummy columns if needed
  for (v in vars) {
    if (!v %in% names(working_data)) {
      working_data[[v]] <- NA
    }
  }
  
  # Scale variables
  scaled_data <- scale_dataset(working_data, vars)
  
  # Get number of patterns
  num_patterns <- length(unique(scaled_data$M))
  message(paste("Number of patterns:", num_patterns))
  
  # Determine which variables are available for each pattern
  vars_per_pattern <- list()
  for (p in 1:num_patterns) {
    pattern_data <- scaled_data[scaled_data$M == p, ]
    
    # Check which variables have non-missing values for this pattern
    if (nrow(pattern_data) > 0) {
      observed_vars <- sapply(vars, function(v) {
        !all(is.na(pattern_data[[v]]))
      })
      #bm: added sort here
      vars_per_pattern[[p]] <- sort(which(observed_vars))
    } else {
      vars_per_pattern[[p]] <- integer(0)
    }
    
    message(paste("Pattern", p, "has", length(vars_per_pattern[[p]]), 
                  "observed variables:", 
                  ifelse(length(vars_per_pattern[[p]]) > 0, 
                         paste(vars_per_pattern[[p]], collapse=", "), 
                         "none (intercept-only model)")))
  }
  
  
  jags_data <- prepare_jags_data(scaled_data, vars)
  
  # Count parameters per pattern (intercept + coefficients for observed variables)
  var_counts <- sapply(2:num_patterns, function(p) {
    # Each pattern has at least an intercept
    1 + length(vars_per_pattern[[p]])
  })
  
  # Calculate total parameters
  total_params <- sum(var_counts)
  message(paste("Total parameters:", total_params))
  
  # Create initial values for 3 chains
  # generalized initialization fn
  # generates one SET of coefficients for each chain
  initialvals2 <- function(numcoeff) {
    # numcoeff: a numeric vector where each element represents the number of slope coefficients for that group.
    n_groups <- length(numcoeff)              # Determine the number of groups/patterns.
    
    #intercepts <- runif(n_groups, min = -4, max = -1)  # Generate one intercept per group.
    # *** If jags throws "invalid parent values", try making the intercepts more negative
    #  (maybe need to make them more negative as the number of patterns increases, to keep the complete-case probability from getting too low?)
    intercepts <- runif(n_groups, min = -5, max = -3)  # Initialize an empty vector for storing the initial values.
    lim <- 0.06                               # Limit for the slope coefficients.
    
    init_values <- c()  # Initialize an empty vector for storing the initial values.
    
    # Loop over each group and append its intercept and its slope coefficients.
    for (i in seq_along(numcoeff)) {
      group_vals <- c(intercepts[i], runif(numcoeff[i], min = -lim, max = lim))
      init_values <- c(init_values, group_vals)
    }
    
    return(init_values)
  }
  # calculate the number of slope coeffs for each pattern M>1
  nCoeff = unlist( lapply(vars_per_pattern[2: length(vars_per_pattern)], length) )
  # init has one set of start vals per chain
  nchains = 3  #TEMP: LATER USE AN ARG PASSED TO WRAPPER FN
  init <- lapply(1:nchains, function(i) list(g = initialvals2(nCoeff)))
  
  
  # Create JAGS model
  jags_model_text <- create_jags_model(num_patterns, vars_per_pattern)
  
  # Print model for debugging
  message("JAGS model:")
  cat(jags_model_text)
  
  # Create a temporary file for the JAGS model
  model_file <- tempfile()
  writeLines(jags_model_text, model_file)
  
  # Run JAGS with try-catch to provide better error messages
  jags_fit <- tryCatch({
    R2jags::jags(
      data = jags_data,
      inits = init,
      n.chains = 3,
      parameters.to.save = 'g',
      n.iter = 1000,
      n.burnin = 500,
      n.thin = 1,
      model.file = model_file
    )
  }, error = function(e) {
    message("JAGS error: ", e$message)
    message("Trying a simpler approach...")
    
    # If error occurs, try a simpler model with stronger priors
    simpler_model <- "
    model {
      for(i in 1:N) {
        f[i] ~ dbern(pmiss[i,R[i]])
        
        # Fixed probability for each pattern
        for (p in 2:16) {
          pmiss[i,p] <- theta[p-1]
        }
        
        # Complete cases probability
        pmiss[i,1] <- max(0.01, 1 - sum(theta[]))
      }
      
      # Constraint
      for (j in 1:Nc) {
        onesc[j] ~ dbern(C[j])
        C[j] <- step(pmiss[j,1]-c)
      }
      
      # Dirichlet prior on pattern probabilities
      theta[1:15] ~ ddirch(alpha[])
    }
    "
    writeLines(simpler_model, model_file)
    
    # Create alpha parameter for Dirichlet (all 1's for uniform)
    jags_data$alpha <- rep(1, num_patterns-1)
    
    # New init values
    init_simple <- list(
      list(theta = rep(1/(num_patterns*2), num_patterns-1)),
      list(theta = rep(1/(num_patterns*2), num_patterns-1)),
      list(theta = rep(1/(num_patterns*2), num_patterns-1))
    )
    
    # Try the simpler model
    tryCatch({
      R2jags::jags(
        data = jags_data,
        inits = init_simple,
        n.chains = 3,
        parameters.to.save = 'theta',
        n.iter = 1000,
        n.burnin = 500, 
        n.thin = 1,
        model.file = model_file
      )
    }, error = function(e2) {
      message("Second attempt also failed: ", e2$message)
      return(NULL)
    })
  })
  
  # Clean up the temporary file
  file.remove(model_file)
  
  if (is.null(jags_fit)) {
    stop("JAGS modeling failed after multiple attempts")
  }
  
  jags_results = list(
    fit = jags_fit,
    scaled_data = scaled_data,
    vars_per_pattern = vars_per_pattern,
    model_text = jags_model_text
  )
  
  return(list(
    fit = jags_fit,
    scaled_data = scaled_data,
    vars_per_pattern = vars_per_pattern,
    model_text = jags_model_text
  ))
}


#' Calculate IPMW weights using JAGS results
#' @param jags_results Results from run_missingness_model
#' @param data Original dataset
#' @param use_posterior_draws If TRUE, draws parameters from the posterior samples themselves. If FALSE, just uses posterior medians.
#' @return Dataset with IPMW weights
calculate_ipmw_weights <- function(jags_results, data, use_posterior_draws = TRUE) {
  
  # Extract parameter estimates (medians)
  g_estimates <- jags_results$fit$BUGSoutput$median$g
  
  # g_draws has one row per draw and one column per parameter to be estimated
  # NOTE for generalization: this breaks if there is only 1 parameter to be estimate, bc sims.matrix has a "deviance" col and a "g" col
  if (use_posterior_draws) {
    g_draws <- jags_results$fit$BUGSoutput$sims.matrix
    keepers = grep("^g", colnames(g_draws))
    g_draws <- g_draws[, keepers]
    
    # catch the case where there's only 1 parameter; in that case g_draws is a vector
    if ( length(keepers) == 1 ) g_draws = matrix(g_draws, ncol=1)
  }
  

  # Get complete cases from scaled data
  cc_scaled <- jags_results$scaled_data %>%
    filter(M == 1)
  
  # Get variables used in model
  # MM: this only works if data, as passed to this fn, *only* contains M variables and the vars in PS model
  vars <- names(cc_scaled)[names(cc_scaled) %in% names(data) & 
                             !names(cc_scaled) %in% c("id", "M", paste0("M", 1:16))]
  
  
  # Get number of patterns
  num_patterns <- length(unique(jags_results$scaled_data$M))
  vars_per_pattern <- jags_results$vars_per_pattern
  
  # Print variable information
  message("Variables in model:")
  print(vars)
  message("Parameter estimates:")
  print(g_estimates)
  
  # Initialize pattern probabilities
  pattern_probs <- matrix(0, nrow = nrow(cc_scaled), ncol = num_patterns - 1)
  
  # For each non-complete pattern, calculate probability
  # Parameter index counter
  param_idx <- 1
  
  #browser()
  for (p in 2:num_patterns) {
    message(paste("Pattern", p, "variable indices:", paste(vars_per_pattern[[p]], collapse = ", ")))
    
    # Compute linear predictor
    if (use_posterior_draws) {
      num_draws <- nrow(g_draws)
      num_obs <- nrow(cc_scaled)
      
      # Get intercept
      intercept <- as.numeric(g_draws[, param_idx])
      param_idx <- param_idx + 1
      
      # Initialize linear predictor matrix
      # Repeat intercept into a matrix with dimensions: [num_draws × num_obs]
      # note: this matrix will look blank if you print it, but that's only how R displays it. View() works.
      # lp_mat has one row per draw and one col per *observation* in the CC dataset
      lp_mat <- matrix(intercept, nrow = num_draws, ncol = num_obs)  # broadcasts intercept over columns
      
      # Add covariate terms
      if (length(vars_per_pattern[[p]]) > 0) {
        for (v_idx in vars_per_pattern[[p]]) {
          var_name <- vars[v_idx]
          # recall that g_draws has one row per draw and one column per parameter to be estimated
          coef <- g_draws[, param_idx]
          param_idx <- param_idx + 1
          
          X <- matrix(rep(cc_scaled[[var_name]], each = num_draws), nrow = num_draws)
          lp_mat <- lp_mat + coef * X
        }
      }
      
      # Convert to probabilities and take mean over posterior draws
      prob_mat <- plogis(lp_mat)
      pattern_probs[, p - 1] <- colMeans(prob_mat)
      
      # Also store mean plogis in the data frame for reference
      cc_scaled[[paste0("p", p)]] <- colMeans(prob_mat)
      
      message(paste("Pattern", p, "probability summary (draw-averaged):"))
      print(summary(cc_scaled[[paste0("p", p)]]))
      
    }

    
    #@@ haven't tested this yet
    if (!use_posterior_draws) {
      # Use median point estimates
      intercept <- g_estimates[param_idx]
      param_idx <- param_idx + 1
      
      lp <- rep(intercept, nrow(cc_scaled))
      
      if (length(vars_per_pattern[[p]]) > 0) {
        for (v_idx in vars_per_pattern[[p]]) {
          var_name <- vars[v_idx]
          coef <- g_estimates[param_idx]
          param_idx <- param_idx + 1
          
          lp <- lp + coef * cc_scaled[[var_name]]
        }
      }
      
      prob_vec <- plogis(lp)
      pattern_probs[, p - 1] <- prob_vec
      cc_scaled[[paste0("p", p)]] <- prob_vec
      
      message(paste("Pattern", p, "probability summary (from medians):"))
      print(summary(prob_vec))
    }
  }
      
  
  # Calculate p1 (complete case probability)
  # If all pattern probabilities are close to 1, something's wrong
  # Check if this is the case
  if (all(colMeans(pattern_probs) > 0.9)) {
    message("WARNING: All pattern probabilities are very high - likely a calculation error.")
    message("Falling back to empirical frequencies.")
    
    # Use empirical frequencies instead
    emp_freqs <- prop.table(table(jags_results$scaled_data$M))
    for (p in 2:num_patterns) {
      cc_scaled[[paste0("p", p)]] <- emp_freqs[as.character(p)]
    }
    
    # Calculate p1
    cc_scaled$p1 <- emp_freqs["1"]
  } else {
    # Calculate p1 as 1 minus sum of other probabilities
    cc_scaled$p1 <- 1 - rowSums(pattern_probs)
  }
  
  # ### TEMP DEBUGGING - look at rows with negative p1
  # #browser()
  # summary(cc_scaled$p1)
  # x = cc_scaled %>% filter(p1 < 0)
  # ###
  
  
  # Calculate IPMW
  mnum <- mean(data$M == 1, na.rm = TRUE)
  cc_scaled$ipmw <- mnum / cc_scaled$p1
  
  # Trim extreme weights
  if (any(cc_scaled$ipmw > 10)) {
    message("Some weights are large. Trimming at 99th percentile.")
    cc_scaled$ipmw <- pmin(cc_scaled$ipmw, quantile(cc_scaled$ipmw, 0.99))
  }
  
  # Print weight summary
  message("IPMW weight summary:")
  print(summary(cc_scaled$ipmw))
  
  # Merge with original data
  if ("id" %in% names(data)) {
    cc_result <- cc_scaled %>%
      select(id, ipmw, p1) %>%
      right_join(data %>% filter(M == 1), by = "id")
  } else {
    cc_data <- data %>% filter(M == 1)
    cc_data$ipmw <- cc_scaled$ipmw
    cc_data$p1 <- cc_scaled$p1
    cc_result <- cc_data
  }
  
  return(cc_result)
}



#' Fit outcome model using IPMW weights
#' @param weighted_data Data with IPMW weights
#' @param outcome Outcome variable name
#' @param exposure Exposure variable name
#' @param adjustment Variables to adjust for
#' @return Model output and results
fit_outcome_model <- function(weighted_data, outcome, exposure, adjustment = NULL) {
  # Create outcome model formula with adjustment variables and interactions
  if (!is.null(adjustment) && length(adjustment) > 0) {
    # Create all possible interactions
    all_vars <- c(exposure, adjustment)
    interaction_terms <- c()
    
    for (i in 2:length(all_vars)) {
      combos <- combn(all_vars, i, simplify = FALSE)
      for (combo in combos) {
        interaction_terms <- c(interaction_terms, paste(combo, collapse = ":"))
      }
    }
    
    # Combine main effects and interactions
    formula_terms <- c(all_vars, interaction_terms)
    outcome_formula <- as.formula(paste(outcome, "~", paste(formula_terms, collapse = " + ")))
  } else {
    # No adjustment - just exposure
    outcome_formula <- as.formula(paste(outcome, "~", exposure))
  }
  
  # Fit linear model with ipmw weights
  model <- lm(outcome_formula, data = weighted_data, weights = weighted_data$ipmw)
  
  # Get robust standard errors (HC0)
  robust_se <- sqrt(diag(sandwich::vcovHC(model, type = "HC0")))
  
  # Extract coefficient for exposure
  effect_estimate <- coef(model)[[exposure]]
  effect_se <- robust_se[names(coef(model)) == exposure]
  
  # Create results
  results <- tibble(
    estimate = effect_estimate,
    se = effect_se,
    lcl = effect_estimate - 1.96 * effect_se,
    ucl = effect_estimate + 1.96 * effect_se
  )
  
  return(list(
    model = model,
    results = results,
    robust_se = robust_se
  ))
}



# May be out of date: currently running the guts of this
#' Main function to implement IPMW analysis
#' @param data Dataset
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param adjustment Adjustment variables
#' @return Full analysis results
ipmw_analysis <- function(data, outcome, exposure, adjustment = NULL) {
  # Ensure ID exists
  if (!"id" %in% names(data)) {
    data$id <- 1:nrow(data)
  }
  
  # Identify variables to be used in the analysis
  analysis_vars <- unique(c(exposure, outcome, adjustment))
  
  # Add pattern indicators
  data_with_patterns <- create_pattern_indicators(data, analysis_vars)
  
  # Print pattern distribution
  message("Distribution of missingness patterns:")
  print(table(data_with_patterns$M))
  
  # Run Bayesian model for missingness patterns
  message("Running Bayesian model for missingness patterns...")
  jags_results <- run_missingness_model(data_with_patterns, analysis_vars)
  
  # Calculate IPMW weights
  message("Calculating IPMW weights...")
  weighted_data <- calculate_ipmw_weights(jags_results, data_with_patterns)
  
  message("IPMW weight summary:")
  print(summary(weighted_data$ipmw))
  
  # Fit outcome model
  message("Fitting outcome model...")
  model_results <- fit_outcome_model(weighted_data, outcome, exposure, adjustment)
  
  message("IPMW analysis results:")
  print(model_results$results)
  
  # Return complete results
  return(list(
    ipmw_results = model_results$results,
    model = model_results$model,
    weighted_data = weighted_data,
    jags_fit = jags_results$fit,
    missingness_patterns = table(data_with_patterns$M)
  ))
}

# Example usage:
# results <- ipmw_analysis(
#   data = my_data,
#   outcome = "D",
#   exposure = "A",
#   adjustment = c("B", "C")
# )
