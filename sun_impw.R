

rm(list=ls()) 

# PRELIMINARIES  -------------------------------------------------

library(tidyverse)
library(mice)
library(R2jags)
library(sandwich)
library(lmtest)
library(tibble)

options(scipen=999)

# also need to source both helper_IAI.R and helper_for_draft_ipmw.R


# GENERATE DIFFERENT DATASETS  -------------------------------------------------

# ~ My own DAG 12B and variations -------------------------------------------------



sim_obj = sim_data(.p = data.frame(N=10000,
                                   coef_of_interest = "A",
                                   dag_name = "12B" ) )

# # WORKS (only 1 incomplete var)
# # if all vars incomplete, results in error: JAGS error: Error in node f[1]
# # Invalid parent values
# outcome = "B"
# exposure = "A1"
# adjustment = c("C", "D")  
# #

# # 2025-04-14: Works now!!
# # still 2 missingness patterns, but expect NOT to have all weights = 1
# # in particular, D predicts RC
# # here, the issue is that complete cases get extreme wts (362!)
# outcome = "B1"
# exposure = "A1"
# adjustment = c("C", "D")


# # WORKS
# # only 1 variable, and it's incomplete
# data = sim_obj$du %>% select(A)
# outcome = "A"
# exposure = NULL
# adjustment = NULL

# # WORKS
# # 2 vars: one complete and one incomplete
# data = sim_obj$du %>% select(A, B1)
# outcome = "A"
# exposure = NULL
# adjustment = NULL

# # WORKS
# # 2 vars: BOTH incomplete
# data = sim_obj$du %>% select(A, B)
# outcome = "B"
# exposure = "A"
# adjustment = NULL

# # WORKS
# # 3 vars: all incomplete
# data = sim_obj$du %>% select(A, B, C)
# outcome = "B"
# exposure = "A"
# adjustment = "C"

# WORKS
# 3 vars: all incomplete, but D instead of C
# data = sim_obj$du %>% select(A, B, D)
# outcome = "B"
# exposure = "A"
# adjustment = "D"

# ****BREAKS!!
# 4 incomplete vars
data = sim_obj$du %>% select(A, B, C, D)
outcome = "B"
exposure = "A"
adjustment = c("C", "D")





# # ~ Ross' "simple example" DAG -------------------------------------------------
# 
# # YES!! WITH THE RIGHT INITIAL VALS, AGREES WITH HERS.
# 
# set.seed(13)
# 
# # Full data
# gendata <- function(n){
#   tibble(
#     id = c(1:n),
#     z = rbinom(n, size=1, prob=0.5), #binary confounder
#     x = rbinom(n, size=1, prob=0.2 + 0.2*z), #binary exposure
#     py = 1/(1+exp(-1*(log(0.1/0.9) + log(2)*z))), #probability of y under no exposure
#     y0 = rbinom(n, size=1, prob=py), #binary potential outcome under no exposure
#     y1 = rbinom(n, size=1, prob=py + .1), #binary potential outcomes under exposure
#     y = x*y1 + (1-x)*y0) #observed binary outcome
# }
# 
# full <- gendata(20000)
# apply(full, 2, mean)
# 
# 
# add_missing <- function(data,parm){
# 
#   pR2 <-  with(data,plogis(parm[1] + parm[2]*z  + parm[3]*x  + parm[4]*y)) #probability of being in pattern 2 (observed: Z, X)
#   pR3 <-  with(data,plogis(parm[5] + parm[6]*z  + parm[7]*x  + parm[8]*y)) #probability of being in pattern 3 (observed: X)
#   pR4 <-  with(data,plogis(parm[9] + parm[10]*z + parm[11]*x + parm[12]*y)) #probability of being in pattern 4 (observed: X, Y)
#   pR1 <- 1 - pR2 - pR3 - pR4 #probability of being in pattern 1
#   probs <- cbind(pR1,pR2,pR3,pR4)
#   colnames(probs) <- NULL
# 
#   # Generate R
#   data %>%
#     mutate(R = Hmisc::rMultinom(probs,1), #multinomial draw using the probabilities for each pattern
#            z = ifelse(R %in% c(3,4),NA,z),
#            y = ifelse(R %in% c(2,3),NA,y),
#            R1 = ifelse(R==1,1,0),
#            R2 = ifelse(R==2,1,0),
#            R3 = ifelse(R==3,1,0),
#            R4 = ifelse(R==4,1,0))
# }
# 
# # Set parameters of missingness models
# g20 <- -0.945 # selected to produce 25% pattern 2
# g21 <- log(.5)
# g22 <- log(1.7)
# g23 <- 0  # MM: by MAR, since this is the coef for Y, which is missing for this pattern
# g30 <- -1.815 # selected to produce 15% pattern 3
# g31 <- 0
# g32 <- log(1.3)
# g33 <- 0
# g40 <- -2.41 # selected to produce 10% pattern 4
# g41 <- 0
# g42 <- log(2)
# g43 <- log(0.8)
# gams <- c(g20, g21, g22, g23,
#           g30, g31, g32, g33,
#           g40, g41, g42, g43)
# 
# withmiss <- add_missing(full, gams)
# prop.table(table(withmiss$R))
# 
# # relabel vars to match what we'll use below
# data = withmiss %>% select(id, z, x, y) %>%
#   rename(Z1 = z, A1 = x, B1 = y)
# 
# data = as.data.frame(data)
# 
# 
# 
# 
# exposure = "A1"
# outcome = "B1"
# adjustment = "Z1" # different from Ross: she estimated ATE





# CLAUDE FROM ROSS BAYESIAN CODE  -------------------------------------------------

# run guts of wrapper
# Ensure ID exists
if (!"id" %in% names(data)) {
  data$id <- 1:nrow(data)
}

# Identify variables to be used in the analysis
analysis_vars <- unique(c(adjustment, exposure, outcome))

# Add pattern indicators
data_with_patterns <- create_pattern_indicators(data, analysis_vars)

# Print pattern distribution
message("Distribution of missingness patterns:")
print(table(data_with_patterns$M))
# MM: fit the true model to check for extreme coefs (for 12B)
# glm( I(M == 2) ~ A1 + B1 + C1 + D1, data = data_with_patterns)

# # Run Bayesian model for missingness patterns
# message("Running Bayesian model for missingness patterns...")
# # debug(run_missingness_model)
# jags_results <- run_missingness_model(data_with_patterns, analysis_vars)
# jags_results$fit
# #bm

# ~ ALTERNATE: Run guts of run_missingness_model step by step ----
# Create pattern indicators if not already present
# important: alphabetize them to keep a consistent order with make_pattern_indicators below
vars = sort(analysis_vars)

if (!"M" %in% names(data_with_patterns)) {
  working_data <- create_pattern_indicators(data_with_patterns, vars)
} else {
  working_data = data_with_patterns
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

# so far agrees with Ross code when using their dataset, though patterns are in different order


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
# MM: much less dispersed than in "simple example"
# init <- list(
#   list(g = runif(total_params, -0.1, 0.1)),
#   list(g = runif(total_params, -0.1, 0.1)),
#   list(g = runif(total_params, -0.1, 0.1))
# )
# TEMP: USE INITIAL VALS FROM SIMPLE EXAMPLE (MORE DISPERSED)
# THIS MADE A HUGE DIFFERENCE!
# initialvals <- function(numcoeff){
#   ints <- runif(3, min=-4,max=-1)  # intercepts
#   lim <- .06
#   c(ints[1],runif(numcoeff[1], min=-lim,max=lim),
#     ints[2],runif(numcoeff[2], min=-lim,max=lim),
#     ints[3],runif(numcoeff[3], min=-lim,max=lim))
# }
# init <- list(list(g=initialvals(c(2,1,2))),
#              list(g=initialvals(c(2,1,2))),
#              list(g=initialvals(c(2,1,2))))
# END OF ROSS' CODE
# generalized initialization fn
# generates one SET of coefficients for each chain
initialvals2 <- function(numcoeff) {
  # numcoeff: a numeric vector where each element represents the number of slope coefficients for that group.
  n_groups <- length(numcoeff)              # Determine the number of groups/patterns.
  #bm: TEMP DEBUGGING
  #intercepts <- runif(n_groups, min = -4, max = -1)  # Generate one intercept per group.
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
    #n.chains = 5,  # MM increased (original:3 ) to try to improve convergence
    parameters.to.save = 'g',
    n.iter = 1000,
    #n.iter = 10000,  # MM edited to try to improve convergence
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
## END GUTS OF run_missingness_model

# sanity check
jags_results$fit


# MM: reduce dataset to only vars in IPMW model and M indicator
all_vars = names(data_with_patterns)
m_vars <- all_vars[startsWith(all_vars, "M")]
data_with_patterns2 = data_with_patterns %>% select( all_of(c( "id", analysis_vars, m_vars ) ) )
names(data_with_patterns2)
# "id" "C"  "D"  "A1" "B1" "M" 

# Calculate IPMW weights
message("Calculating IPMW weights...")
weighted_data <- calculate_ipmw_weights(jags_results, data_with_patterns2, use_posterior_draws = TRUE)
#weighted_data <- calculate_ipmw_weights(jags_results, data_with_patterns2, use_posterior_draws = FALSE)

message("IPMW weight summary:")
print(summary(weighted_data$ipmw))

# Fit outcome model
message("Fitting outcome model...")
model_results <- fit_outcome_model(weighted_data, outcome, exposure, adjustment)

message("IPMW analysis results:")
print(model_results$results)


#bm: indeed, seems to work if at least 1 var is complete, but not if all are incomplete, yayyy
#next: try to dx the issue with having all vars incomplete. maybe start with a dataset that has only 1 variables and it's incomplete?





# # PACKAGE  -------------------------------------------------
# 
# library(NMMIPW)
# 
# nmm_fit(data = du, O = cbind(du$RA, du$RB, du$RC), formula = B ~ A * C, func = lm)
# 
# 
# debug(nmm_fit)
# undebug(nmm_fit)
# 
# 
# # their example
# n = 100
# X = rnorm(n, 0, 1)
# Y = rnorm(n, 1 * X, 1)
# O1 = rbinom(n, 1, 1/(1 + exp(-1 - 0.5 * X)))
# O2 = rbinom(n, 1, 1/(1 + exp(+0.5 + 1 * Y)))
# O = cbind(O1, O2)
# df <- data.frame(Y = Y, X = X)
# fit <- nmm_fit(data = df, O = O, formula = Y ~ X, funct = lm)
# summary(fit)
# 
# 
# 
# # CLAUDE  -------------------------------------------------
# 
# # problem was that optimizers aren't finding the min -- not sure why
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# #
# ##########################################
# 
# library(tidyverse)
# library(geepack)
# library(mice)
# library(Hmisc)
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# # With direct adjustment and interaction terms
# #
# ##########################################
# 
# library(tidyverse)
# library(geepack)
# library(mice)
# library(Hmisc)
# library(Formula)
# 
# 
# 
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# #
# ##########################################
# 
# 
# ### TEMP ONLY
# data = du
# outcome = "B"
# exposure = "A"
# adjustment = c("C", "D")
# ###
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# # With direct adjustment and interaction terms
# #
# ##########################################
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# # With direct adjustment and interaction terms
# #
# ##########################################
# 
# ###########################################
# #
# # Flexible IPMW Implementation for Custom Datasets
# # Based on Sun & Tchetgen Tchetgen (2018)
# # With direct adjustment and conditional effects
# #
# ##########################################
# 
# library(tidyverse)
# library(mice)
# library(Hmisc)
# 
# 
# # Function to implement IPMW for any dataset with variables A1, B1, C1, D1 (or subset)
# # Requires specifying:
# # - data: your dataset
# # - outcome: the outcome variable name (e.g., "D1")
# # - exposure: the exposure variable name (e.g., "A1")
# # - adjustment: vector of variable names to adjust for in outcome model (e.g., c("B1", "C1"))
# ipmw_analysis <- function(data, outcome, exposure, adjustment = NULL, init_value = 0) {
#   
#   # Create a copy of the original data to work with
#   original_data <- data
#   
#   # Identify variables to be used in the analysis
#   analysis_vars <- unique(c(exposure, outcome, adjustment))
#   
#   # Extract only needed variables
#   working_data <- original_data %>%
#     select(all_of(analysis_vars))
#   
#   working_data$id = 1:nrow(working_data)
#   
#   # Determine missingness patterns
#   message("Analyzing missingness patterns...")
#   md_pattern <- mice::md.pattern(working_data[, analysis_vars])
#   print(md_pattern)
#   
#   # Identify complete cases and create pattern indicator
#   working_data <- working_data %>%
#     mutate(M = 0)
#   
#   # Create indicators for each pattern
#   num_patterns <- nrow(md_pattern) - 1  # Subtract 1 for the summary row
#   
#   # Loop through patterns and assign pattern numbers
#   # Pattern 1 is complete cases
#   for(p in 1:num_patterns) {
#     pattern_mask <- rep(TRUE, nrow(working_data))
#     
#     # For each variable, check if it should be missing in this pattern
#     for(v in 1:length(analysis_vars)) {
#       var_name <- analysis_vars[v]
#       # If md_pattern shows 0 for this variable in this pattern, it should be missing
#       if(md_pattern[p, v] == 0) {
#         pattern_mask <- pattern_mask & is.na(working_data[[var_name]])
#       } else {
#         pattern_mask <- pattern_mask & !is.na(working_data[[var_name]])
#       }
#     }
#     
#     # Assign pattern number
#     working_data$M[pattern_mask] <- p
#   }
#   
#   # Create indicators for each pattern (M1, M2, etc.)
#   for(p in 1:num_patterns) {
#     working_data[[paste0("M", p)]] <- ifelse(working_data$M == p, 1, 0)
#   }
#   
#   # Display pattern distribution
#   pattern_freq <- table(working_data$M)
#   message("Distribution of missingness patterns:")
#   print(pattern_freq)
#   
# 
#   # Create outcome model formula with adjustment variables and interaction terms
#   if(!is.null(adjustment) && length(adjustment) > 0) {
#     # Create main effects terms
#     main_effects <- paste(c(exposure, adjustment), collapse = " + ")
#     
#     # Create interaction terms for adjustment variables
#     if(length(adjustment) > 1) {
#       # Get all pairwise and higher-order interactions among adjustment variables
#       adjustment_interactions <- c()
#       for(i in 2:length(adjustment)) {
#         interact_terms <- combn(adjustment, i, simplify = FALSE)
#         for(terms in interact_terms) {
#           adjustment_interactions <- c(adjustment_interactions, paste(terms, collapse = ":"))
#         }
#       }
#       
#       # Add interaction with exposure
#       exposure_interactions <- paste(exposure, adjustment, sep = ":")
#       
#       # Combine all terms
#       all_terms <- c(main_effects, adjustment_interactions, exposure_interactions)
#       outcome_formula <- as.formula(paste(outcome, "~", paste(all_terms, collapse = " + ")))
#     } else {
#       # Just one adjustment variable, so only interaction with exposure
#       outcome_formula <- as.formula(paste(outcome, "~", exposure, "+", adjustment, "+", 
#                                           exposure, ":", adjustment))
#     }
#   } else {
#     # No adjustment variables
#     outcome_formula <- as.formula(paste(outcome, "~", exposure))
#   }
#   
#   
#   # ~ UMLE for IPMW if there are at least 2 patterns
#   if(num_patterns > 1) {
#     message("Implementing UMLE for IPMW...")
#     
#     # Build the negative log-likelihood function based on observed patterns
#     umleLogL <- function(g, data) {
#       # Initialize a matrix to hold pattern probabilities (one column per non-complete pattern)
#       pattern_probs <- matrix(NA, nrow = nrow(data), ncol = num_patterns - 1)
#       
#       # Parameter index counter
#       param_idx <- 1
#       
#       # For each non-complete pattern, build a probability model
#       for(p in 2:num_patterns) {
#         # Identify which variables are observed in this pattern
#         observed_in_pattern <- colnames(md_pattern)[md_pattern[p, ] == 1]
#         observed_in_pattern <- observed_in_pattern[observed_in_pattern %in% analysis_vars]
#         
#         # Create formula terms for this pattern's model
#         formula_terms <- c()
#         for(var in observed_in_pattern) {
#           # Add term if it's observed in this pattern
#           if(var %in% names(data)) {
#             formula_terms <- c(formula_terms, paste0("g[", param_idx, "]*", var))
#             param_idx <- param_idx + 1
#           }
#         }
#         
#         # Create the linear predictor with intercept
#         lp <- paste0("g[", param_idx, "]", 
#                      ifelse(length(formula_terms) > 0, paste0(" + ", paste(formula_terms, collapse = " + ")), ""))
#         param_idx <- param_idx + 1
#         
#         # Calculate probability using logistic function
#         pattern_probs[, p-1] <- with(data, plogis(eval(parse(text = lp))))
#       }
#       
#       # Calculate sum of probabilities
#       sump <- rowSums(pattern_probs, na.rm = TRUE)
#       # Ensure sump is less than 1 to avoid negative probabilities for complete cases
#       sump <- pmin(sump, 0.999)
#       
#       
#       # Log-likelihood for each observation
#       ilogL <- rep(NA, nrow(data))
#       
#       # Complete cases
#       ilogL[data$M == 1] <- log(1 - sump[data$M == 1])
#       
#       # Non-complete cases
#       for(p in 2:num_patterns) {
#         ilogL[data$M == p] <- log(pattern_probs[data$M == p, p-1])
#       }
#       
#       # Return negative log-likelihood
#       return(-1 * sum(ilogL, na.rm = TRUE))
#     }
#     
#     # Estimate number of parameters needed
#     total_params <- 0
#     for(p in 2:num_patterns) {
#       observed_in_pattern <- colnames(md_pattern)[md_pattern[p, ] == 1]
#       observed_in_pattern <- observed_in_pattern[observed_in_pattern %in% analysis_vars]
#       total_params <- total_params + length(observed_in_pattern) + 1  # +1 for intercept
#     }
#     
#     # Initial values: all zeros
#     #init <- rep(0, total_params)
#     init <- rep(init_value, total_params)
#     
#     # *** Optimize to find UMLE estimates
#     umlefit <- tryCatch({
#       nlm(umleLogL, init, data = working_data)
#     }, error = function(e) {
#       message("Error in UMLE optimization: ", e$message)
#       return(NULL)
#     })
#     
#     # c.f. neglogL for initial values:
#     #init = c(-2, 0, -1, 0, -1)
#     umleLogL(init, working_data)
#     
#     if(!is.null(umlefit)) {
#       gumle <- umlefit$estimate
#       
#       # Calculate probability of each pattern for complete cases
#       cc_umle <- working_data %>%
#         filter(M == 1)
#       
#       # Calculate pattern probabilities
#       param_idx <- 1
#       pattern_probs <- matrix(NA, nrow = nrow(cc_umle), ncol = num_patterns - 1)
#       
#       for(p in 2:num_patterns) {
#         observed_in_pattern <- colnames(md_pattern)[md_pattern[p, ] == 1]
#         observed_in_pattern <- observed_in_pattern[observed_in_pattern %in% analysis_vars]
#         
#         formula_terms <- c()
#         for(var in observed_in_pattern) {
#           if(var %in% names(cc_umle)) {
#             formula_terms <- c(formula_terms, paste0("gumle[", param_idx, "]*", var))
#             param_idx <- param_idx + 1
#           }
#         }
#         
#         lp <- paste0("gumle[", param_idx, "]", 
#                      ifelse(length(formula_terms) > 0, paste0(" + ", paste(formula_terms, collapse = " + ")), ""))
#         param_idx <- param_idx + 1
#         
#         pattern_probs[, p-1] <- with(cc_umle, plogis(eval(parse(text = lp))))
#       }
#       
#       # Calculate probability of being a complete case
#       cc_umle$p1 <- 1 - rowSums(pattern_probs)
#       # Ensure probabilities are positive
#       cc_umle$p1 <- pmax(cc_umle$p1, 0.001)
#       
#       # Obtain marginal probability of complete cases
#       mnum <- mean(working_data$M1)
#       
#       # Create IPMW weights
#       cc_umle <- cc_umle %>%
#         mutate(ipmw = mnum / p1)
#       
#       message("IPMW summary:")
#       print(summary(cc_umle$ipmw))
#       
#       # ***Fit weighted outcome model using geeglm
#       uoutmod <- tryCatch({
#         geeglm(outcome_formula, 
#                data = cc_umle, 
#                weights = cc_umle$ipmw, 
#                id = cc_umle$id, 
#                corstr = "independence")
#       }, error = function(e) {
#         message("Error in UMLE outcome model: ", e$message)
#         return(NULL)
#       })
#       
#       # ***added
#       # Calculate simple treatment effect and return results
#       # Calculate treatment effect and return results
#       if(!is.null(uoutmod)) {
#         # Extract the coefficient for the exposure variable
#         treatment_effect <- coef(uoutmod)[[exposure]]
#         se <- coef(summary(uoutmod))[exposure, 2]
#         
#         results_umle <- tibble(
#           estimate = treatment_effect,
#           se = se,
#           lcl = treatment_effect - 1.96 * se,
#           ucl = treatment_effect + 1.96 * se
#         )
#         
#         message("UMLE-IPMW analysis results:")
#         print(results_umle)
#         
#         # Return results
#         return(list(
#           cc_results = if(exists("results_cc")) results_cc else NULL,
#           ipmw_results = results_umle,
#           ipmw_model = uoutmod,
#           cc_model = cc_out,
#           weighted_data = cc_umle,
#           missingness_patterns = table(working_data$M),
#           umle_parameters = gumle
#         ))
#       }
#       
#   }  # end if(!is.null(umlefit))
#         
# } # end if(num_patterns > 1)
#   
# 
#   
# }
#   
# res = ipmw_analysis(
#    data = du,
#    outcome = "B",
#    exposure = "A",
#    adjustment = c("C", "D") )
# 
# res$ipmw_results
# 
# 
# ##### sanity checks:
# 
# # compare to existing results for manually implemented IPMW for monotone MAR (should agree)
# sim_obj = sim_data(.p = data.frame(N=10000,
#                                    coef_of_interest = "A",
#                                    dag_name = "6A" ) )
# 
# 
# du = sim_obj$du
# 
# # with init values of 0
# # res = ipmw_analysis(
# #   data = du,
# #   outcome = "B",
# #   exposure = "A",
# #   adjustment = c("C") )
# # 
# # res$ipmw_results
# # res$umle_parameters
# 
# res = ipmw_analysis(
#   data = du,
#   outcome = "B",
#   exposure = "A",
#   adjustment = c("C"),
#   init_value = -2)
# 
# res$ipmw_results
# res$umle_parameters
# 
# # with these init values, it comes up with more reasonable umle_parameters, but still overestimates the CATE. hmmmm.
# 
