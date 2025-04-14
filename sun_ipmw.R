

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

# WORKS (only 1 incomplete var)
# if all vars incomplete, results in error: JAGS error: Error in node f[1]
# Invalid parent values
outcome = "B1"
exposure = "A"
adjustment = c("C1", "D1")

# WORKS
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

# # **NOW WORKS WITH LOWER START VALUES FOR INTERCEPTS
# # 4 incomplete vars
# outcome = "B"
# exposure = "A"
# adjustment = c("C", "D")

data = sim_obj$du %>% select( all_of( c(exposure, outcome, adjustment) ) )




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





# RUN EVERYTHING  -------------------------------------------------

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
## END GUTS OF run_missingness_model

# sanity check
jags_results$fit


# Reduce dataset to only vars in IPMW model and M indicator
all_vars = names(data_with_patterns)
m_vars <- all_vars[startsWith(all_vars, "M")]
data_with_patterns2 = data_with_patterns %>% select( all_of(c( "id", analysis_vars, m_vars ) ) )
names(data_with_patterns2)

# Calculate IPMW weights
message("Calculating IPMW weights...")
# using posterior medians themselves
weighted_data <- calculate_ipmw_weights(jags_results, data_with_patterns2, use_posterior_draws = TRUE)
# using just the posterior medians (as in both Sun's and Ross' implementations)
#weighted_data <- calculate_ipmw_weights(jags_results, data_with_patterns2, use_posterior_draws = FALSE)

message("IPMW weight summary:")
print(summary(weighted_data$ipmw))

# Fit outcome model
message("Fitting outcome model...")
model_results <- fit_outcome_model(weighted_data, outcome, exposure, adjustment)

message("IPMW analysis results:")
print(model_results$results)
