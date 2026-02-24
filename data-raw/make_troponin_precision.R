# Generate troponin_precision dataset for valytics package
# High-sensitivity cardiac troponin I precision study
#
# Design: 6 concentration levels × 5 days × 2 runs × 2 replicates = 120 observations
# CV pattern follows hyperbolic model: CV = sqrt(a^2 + (b/x)^2)
# with a = 3% (asymptotic CV) and b = 25 (controls low-end CV)

set.seed(20240201)

# Concentration levels (ng/L) - spanning clinical range
levels <- c(5, 10, 25, 50, 100, 500)
level_names <- paste0("L", 1:6)

# True CV at each level using hyperbolic model
# CV = sqrt(3^2 + (25/x)^2)
# This gives: L1(5)=15.8%, L2(10)=8.5%, L3(25)=4.2%, L4(50)=3.4%, L5(100)=3.1%, L6(500)=3.0%
a <- 3
b <- 25
true_cv <- sqrt(a^2 + (b / levels)^2)

# Design parameters
n_days <- 5
n_runs <- 2
n_reps <- 2

# Generate data
data_list <- list()

for (i in seq_along(levels)) {
  conc <- levels[i]
  cv <- true_cv[i] / 100  # Convert to proportion
  sd_total <- conc * cv
  
  # Partition variance: 40% between-day, 20% between-run, 40% within-run
  var_total <- sd_total^2
  var_day <- 0.40 * var_total
  var_run <- 0.20 * var_total
  var_error <- 0.40 * var_total
  
  for (d in 1:n_days) {
    day_effect <- rnorm(1, 0, sqrt(var_day))
    
    for (r in 1:n_runs) {
      run_effect <- rnorm(1, 0, sqrt(var_run))
      
      for (rep in 1:n_reps) {
        error <- rnorm(1, 0, sqrt(var_error))
        
        # Measured value = true + day effect + run effect + error
        measured <- conc + day_effect + run_effect + error
        
        # Ensure positive (truncate at 0.1 for very low values)
        measured <- max(0.1, measured)
        
        # Round to realistic precision
        if (conc < 20) {
          measured <- round(measured, 2)
        } else if (conc < 100) {
          measured <- round(measured, 1)
        } else {
          measured <- round(measured, 0)
        }
        
        data_list[[length(data_list) + 1]] <- data.frame(
          level = level_names[i],
          day = paste0("D", d),
          run = paste0("R", r),
          replicate = rep,
          value = measured,
          target = conc
        )
      }
    }
  }
}

troponin_precision <- do.call(rbind, data_list)
troponin_precision$level <- factor(troponin_precision$level, levels = level_names)
troponin_precision$day <- factor(troponin_precision$day)
troponin_precision$run <- factor(troponin_precision$run)
rownames(troponin_precision) <- NULL

# Verify the data
cat("Dataset structure:\n")
str(troponin_precision)
cat("\nSummary by level:\n")
aggregate(value ~ level + target, data = troponin_precision, 
          FUN = function(x) c(mean = mean(x), sd = sd(x), cv = 100 * sd(x) / mean(x)))

# Save
usethis::use_data(troponin_precision, overwrite = TRUE)