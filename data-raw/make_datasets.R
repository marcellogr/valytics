# Data Generation Script for valytics Package ----
# This script generates synthetic but realistic datasets for method comparison
# examples. The data mimics patterns seen in real clinical laboratory validation
# studies without using actual patient data.
#
# Run this script to regenerate the .rda files in data/

set.seed(20240115)  # Reproducible generation

# Dataset 1: glucose_methods ----
# Context: Comparison of two glucose analyzers in a clinical laboratory.
# Reference method: Established hexokinase-based analyzer (method_a)
# New method: Point-of-care glucose meter (method_b)
#
# Realistic characteristics:
# - Range: 50-350 mg/dL (covers hypoglycemia to hyperglycemia)
# - Small positive bias in POC meter (~3-5 mg/dL)
# - Slight proportional error (POC reads ~2% higher at high concentrations)
# - Random error: CV ~3% for reference, ~5% for POC
# -----------------------------------------------------------------------------

n_glucose <- 60

# True glucose concentrations (underlying latent values)
# Mix of normal and pathological samples
true_glucose <- c(
  runif(15, 50, 80),    # Hypoglycemic range
  runif(30, 80, 140),   # Normal/prediabetic range
  runif(15, 140, 350)   # Diabetic range
)

# Reference method: low CV, no bias
method_a <- true_glucose * rnorm(n_glucose, mean = 1.0, sd = 0.025)

# POC method: higher CV, small constant + proportional bias
method_b <- true_glucose * rnorm(n_glucose, mean = 1.02, sd = 0.045) +
            rnorm(n_glucose, mean = 3, sd = 1.5)

# Round to realistic precision
method_a <- round(method_a, 0)
method_b <- round(method_b, 0)

# Create data frame
glucose_methods <- data.frame(
  sample_id = sprintf("GLU%03d", seq_len(n_glucose)),
  reference = method_a,
  poc_meter = method_b
)

# Shuffle to remove concentration ordering
glucose_methods <- glucose_methods[sample(n_glucose), ]
row.names(glucose_methods) <- NULL


# Dataset 2: creatinine_serum ----
# Context: Comparison of enzymatic vs Jaffe creatinine methods
# This is a classic example where Jaffe method has known positive interference
# from proteins and other substances.
#
# Realistic characteristics:
# - Range: 0.4-8.0 mg/dL (normal to severe kidney disease)
# - Jaffe shows positive bias, especially at low concentrations
# - Bias decreases proportionally at high concentrations
# - Some outliers due to interfering substances (bilirubin, ketones)
# -----------------------------------------------------------------------------

n_creat <- 80

# True creatinine concentrations
# Weighted toward normal range but including CKD patients
true_creat <- c(
  runif(40, 0.5, 1.2),   # Normal range
  runif(25, 1.2, 3.0),   # Mild-moderate CKD
  runif(15, 3.0, 8.0)    # Severe CKD
)

# Enzymatic method (reference): more specific, lower CV
enzymatic <- true_creat * rnorm(n_creat, mean = 1.0, sd = 0.03)

# Jaffe method: positive interference, higher at low concentrations
# Bias model: ~0.2 mg/dL constant + decreasing proportional effect
jaffe_bias <- 0.15 + 0.08 * exp(-true_creat / 2)  # Decreasing bias at high values
jaffe <- true_creat * rnorm(n_creat, mean = 1.0, sd = 0.04) +
         jaffe_bias + rnorm(n_creat, mean = 0, sd = 0.05)

# Add a few outliers (high bilirubin or hemolysis interference)
outlier_idx <- sample(n_creat, 4)
jaffe[outlier_idx] <- jaffe[outlier_idx] + runif(4, 0.3, 0.6)

# Round to realistic precision
enzymatic <- round(enzymatic, 2)
jaffe <- round(jaffe, 2)

# Ensure no negative values
enzymatic <- pmax(enzymatic, 0.3)
jaffe <- pmax(jaffe, 0.3)

# Create data frame
creatinine_serum <- data.frame(
  sample_id = sprintf("CREAT%03d", seq_len(n_creat)),
  enzymatic = enzymatic,
  jaffe = jaffe
)

# Shuffle
creatinine_serum <- creatinine_serum[sample(n_creat), ]
row.names(creatinine_serum) <- NULL


# Dataset 3: troponin_cardiac ----
# Context: High-sensitivity cardiac troponin I assay comparison
# Two different hs-cTnI platforms from different manufacturers
#
# Realistic characteristics:
# - Range: 2-5000 ng/L (near LOD to acute MI levels)
# - Log-normal distribution (most values low, few very high)
# - Proportional difference between methods (no 1:1 equivalence)
# - Higher CV at low concentrations (near detection limit)
# -----------------------------------------------------------------------------

n_trop <- 50

# Log-normal distribution for true troponin values
true_trop <- exp(rnorm(n_trop, mean = 3.5, sd = 1.8))
true_trop <- pmin(true_trop, 5000)  # Cap at 5000 ng/L
true_trop <- pmax(true_trop, 2)     # Floor at 2 ng/L (near LOD)

# Platform A (Abbott)
# CV increases at low concentrations
cv_a <- 0.05 + 0.10 * exp(-true_trop / 20)
platform_a <- true_trop * rnorm(n_trop, mean = 1.0, sd = cv_a)

# Platform B (Roche) - systematically reads ~15% lower
cv_b <- 0.06 + 0.12 * exp(-true_trop / 20)
platform_b <- true_trop * 0.85 * rnorm(n_trop, mean = 1.0, sd = cv_b)

# Round based on magnitude (clinical reporting conventions)
round_trop <- function(x) {
  ifelse(x < 50, round(x, 1),
         ifelse(x < 500, round(x, 0),
                round(x / 10) * 10))
}

platform_a <- round_trop(platform_a)
platform_b <- round_trop(platform_b)

# Ensure minimum values
platform_a <- pmax(platform_a, 2)
platform_b <- pmax(platform_b, 2)

# Create data frame
troponin_cardiac <- data.frame(
  sample_id = sprintf("TROP%03d", seq_len(n_trop)),
  platform_a = platform_a,
  platform_b = platform_b
)

# Shuffle
troponin_cardiac <- troponin_cardiac[sample(n_trop), ]
row.names(troponin_cardiac) <- NULL


# Save datasets ----

usethis::use_data(glucose_methods, overwrite = TRUE)
usethis::use_data(creatinine_serum, overwrite = TRUE)
usethis::use_data(troponin_cardiac, overwrite = TRUE)

cat("Datasets created:\n")
cat("- glucose_methods:", nrow(glucose_methods), "observations\n")
cat("- creatinine_serum:", nrow(creatinine_serum), "observations\n")
cat("- troponin_cardiac:", nrow(troponin_cardiac), "observations\n")
