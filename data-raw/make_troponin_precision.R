# Data Generation Script for Troponin Precision Dataset
# ======================================================
# This script generates a synthetic multi-level precision dataset suitable
# for demonstrating precision_study() with multiple concentration levels
# and precision_profile() for functional sensitivity estimation.
#
# The data mimics a high-sensitivity cardiac troponin I (hs-cTnI) assay
# precision evaluation across the clinically relevant measurement range.
#
# Run this script to regenerate the .rda file in data/

set.seed(20240206)  # Reproducible generation

# Troponin Precision Dataset ----
# Context: Precision evaluation of a high-sensitivity cardiac troponin I assay
# across 6 concentration levels spanning the analytical measurement interval.
#
# Realistic characteristics:
# - Concentration range: 5-500 ng/L (near LoD to elevated AMI levels)
# - CV decreases with concentration (hyperbolic relationship)
# - Asymptotic CV ~3% at high concentrations
# - Higher CV (~8-10%) near the limit of detection
# - Design: 5 days x 2 runs x 2 replicates per level (EP05-style)
# -----------------------------------------------------------------------------

# Define concentration levels (ng/L)
# These represent typical QC material levels for hs-cTnI assays
conc_levels <- c(5, 10, 25, 50, 100, 500)
n_levels <- length(conc_levels)

# True precision parameters (hyperbolic model: CV = sqrt(a^2 + (b/x)^2))
# a = asymptotic CV at high concentrations (~3%)
# b = concentration-dependent component
true_a <- 3.0
true_b <- 25

# Experimental design parameters
n_days <- 5
n_runs <- 2
n_reps <- 2

# Generate data
troponin_precision <- data.frame()

for (i in seq_along(conc_levels)) {
  conc <- conc_levels[i]
  
  # Calculate true CV for this concentration level
  true_cv_pct <- sqrt(true_a^2 + (true_b / conc)^2)
  true_sd <- conc * true_cv_pct / 100
  
  # Variance components (as fractions of total variance)
  # Between-day variance ~40% of total
  # Between-run variance ~20% of total
  # Within-run (error) variance ~40% of total
  total_var <- true_sd^2
  var_day <- 0.40 * total_var
  var_run <- 0.20 * total_var
  var_error <- 0.40 * total_var
  
  # Generate day effects
  day_effects <- rnorm(n_days, mean = 0, sd = sqrt(var_day))
  
  # Generate data for each day/run/replicate combination
  for (d in 1:n_days) {
    # Generate run effects (nested within day)
    run_effects <- rnorm(n_runs, mean = 0, sd = sqrt(var_run))
    
    for (r in 1:n_runs) {
      # Generate replicates
      for (rep in 1:n_reps) {
        # Measurement = true concentration + day effect + run effect + error
        error <- rnorm(1, mean = 0, sd = sqrt(var_error))
        measurement <- conc + day_effects[d] + run_effects[r] + error
        
        # Ensure positive values (realistic for concentration data)
        measurement <- max(measurement, 0.1)
        
        # Add to data frame
        row <- data.frame(
          level = factor(paste0("L", i), levels = paste0("L", 1:n_levels)),
          concentration = conc,
          day = factor(d),
          run = factor(r),
          replicate = factor(rep),
          value = round(measurement, 2)
        )
        
        troponin_precision <- rbind(troponin_precision, row)
      }
    }
  }
}

# Add level labels with concentration info
level_labels <- paste0("L", 1:n_levels, " (", conc_levels, " ng/L)")
troponin_precision$level_label <- factor(
  troponin_precision$level,
  levels = paste0("L", 1:n_levels),
  labels = level_labels
)

# Reorder columns
troponin_precision <- troponin_precision[, c("level", "level_label", "concentration",
                                             "day", "run", "replicate", "value")]

# Reset row names
rownames(troponin_precision) <- NULL

# Save dataset
usethis::use_data(troponin_precision, overwrite = TRUE)

# Print summary
cat("Dataset created: troponin_precision\n")
cat("Dimensions:", nrow(troponin_precision), "observations x",
    ncol(troponin_precision), "variables\n")
cat("Concentration levels:", paste(conc_levels, collapse = ", "), "ng/L\n")
cat("Design:", n_days, "days x", n_runs, "runs x", n_reps, "replicates per level\n")
cat("\nExpected CV by level:\n")
for (i in seq_along(conc_levels)) {
  expected_cv <- sqrt(true_a^2 + (true_b / conc_levels[i])^2)
  cat(sprintf("  Level %d (%d ng/L): %.1f%%\n", i, conc_levels[i], expected_cv))
}