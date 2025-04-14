

rm(list=ls()) 


# PRELIMINARIES ---------------------------------------------------------------

to.load = c("here",
            "tidyverse",
            "mice",
            "R2jags",
            "sandwich",
            "lmtest",
            "tibble")

# load within installation if needed
for (pkg in to.load) {
  
  cat( paste("\nAbout to try loading package", pkg) )
  
  tryCatch({
    # eval below needed because library() will otherwise be confused
    # https://www.mitchelloharawild.com/blog/loading-r-packages-in-a-loop/
    eval( bquote( library( .(pkg) ) ) )
  }, error = function(err) {
    install.packages(pkg)
  })
}


# set working directories
setwd( here() )
source("sun_ipmw_helper.R")

# if simulating using IAI code, also need this
setwd("/Users/mmathur/Dropbox/Personal computer/Independent studies/*Inchoate/2024/2024-10-20 - IAI (Incomplete auxiliaries in imputation)/Simulation study/Code")
source("helper_IAI.R")



# GENERATE DIFFERENT DATASETS  -------------------------------------------------

# ~ My own DAG 12B and variations -------------------------------------------------

sim_obj = sim_data(.p = data.frame(N=10000,
                                   coef_of_interest = "A",
                                   dag_name = "12B" ) )

# # WORKS (only 1 incomplete var)
# # if all vars incomplete, results in error: JAGS error: Error in node f[1]
# # Invalid parent values
# outcome = "B1"
# exposure = "A"
# adjustment = c("C1", "D1")

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
head(data)


# ~ Ross' "simple example" DAG -------------------------------------------------
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
# sanity check: fit the true model to check for extreme coefs (for 12B)
# glm( I(M == 2) ~ A1 + B1 + C1 + D1, data = data_with_patterns)

# # Run Bayesian model for missingness patterns
# message("Running Bayesian model for missingness patterns...")
# # debug(run_missingness_model)
jags_results <- run_missingness_model(data_with_patterns, analysis_vars)
# sanity check
#@@ for more general implementation, should add checks for convergence (e.g., Rhat, n.eff, etc.)
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
