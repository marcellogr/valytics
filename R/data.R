#' Glucose Method Comparison Dataset
#'
#' @description
#' Synthetic dataset comparing glucose measurements from two analytical methods:
#' a reference hexokinase-based laboratory analyzer and a point-of-care (POC)
#' glucose meter. The data mimics realistic patterns observed in clinical
#' laboratory method validation studies.
#'
#' @format A data frame with 60 observations and 3 variables:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifier.}
#'   \item{reference}{Numeric. Glucose concentration (mg/dL) measured by the
#'     reference hexokinase method.}
#'   \item{poc_meter}{Numeric. Glucose concentration (mg/dL) measured by the
#'     point-of-care glucose meter.}
#' }
#'
#' @details
#' This synthetic dataset was designed to illustrate common patterns in glucose
#' method comparisons:
#'
#' \itemize{
#'   \item **Concentration range**: 50-350 mg/dL, covering hypoglycemia through
#'     severe hyperglycemia
#'   \item **Bias pattern**: The POC meter shows a small positive bias (~3-5 mg/dL)
#'     with slight proportional error at higher concentrations
#'   \item **Precision**: Reference method CV ~2.5%, POC meter CV ~4.5%
#' }
#'
#' The data is suitable for demonstrating Bland-Altman analysis, Passing-Bablok
#' regression, and other method comparison techniques.
#'
#' @source
#' Synthetic data generated to mimic realistic clinical patterns. See
#' `data-raw/make_datasets.R` for the generation script.
#'
#' @examples
#' # Bland-Altman analysis
#' ba <- ba_analysis(reference ~ poc_meter, data = glucose_methods)
#' summary(ba)
#' plot(ba)
#'
#' # Check for proportional bias
#' plot(ba, title = "POC Glucose Meter vs Reference")
#'
#' @seealso [ba_analysis()], [creatinine_serum], [troponin_cardiac]
"glucose_methods"


#' Serum Creatinine Method Comparison Dataset
#'
#' @description
#' Synthetic dataset comparing serum creatinine measurements from two analytical
#' methods: an enzymatic method (reference) and the Jaffe colorimetric method.
#' This is a classic method comparison scenario where the Jaffe method is known
#' to have positive interference from proteins and other chromogens.
#'
#' @format A data frame with 80 observations and 3 variables:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifier.}
#'   \item{enzymatic}{Numeric. Creatinine concentration (mg/dL) measured by
#'     the enzymatic method.}
#'   \item{jaffe}{Numeric. Creatinine concentration (mg/dL) measured by the
#'     Jaffe method.}
#' }
#'
#' @details
#' This synthetic dataset was designed to illustrate well-known patterns in
#' creatinine method comparisons:
#'
#' \itemize{
#'   \item **Concentration range**: 0.4-8.0 mg/dL, from normal kidney function
#'     to severe chronic kidney disease (CKD)
#'   \item **Bias pattern**: Jaffe method shows positive bias, more pronounced
#'     at lower concentrations due to pseudocreatinine interference
#'   \item **Outliers**: A few samples show larger positive bias, simulating
#'     interference from bilirubin, hemolysis, or ketones
#' }
#'
#' The enzymatic method is more specific and has largely replaced the Jaffe
#' method in modern clinical laboratories, though the Jaffe method remains in
#' use in some settings due to lower reagent costs.
#'
#' @source
#' Synthetic data generated to mimic realistic clinical patterns. See
#' `data-raw/make_datasets.R` for the generation script.
#'
#' @references
#' Peake M, Whiting M. Measurement of serum creatinine--current status and
#' future goals. Clin Biochem Rev. 2006;27(4):173-184.
#'
#' @examples
#' # Bland-Altman analysis
#' ba <- ba_analysis(enzymatic ~ jaffe, data = creatinine_serum)
#' summary(ba)
#' plot(ba)
#'
#' # Note: Jaffe typically shows positive bias vs enzymatic
#'
#' @seealso [ba_analysis()], [glucose_methods], [troponin_cardiac]
"creatinine_serum"


#' Cardiac Troponin Method Comparison Dataset
#'
#' @description
#' Synthetic dataset comparing high-sensitivity cardiac troponin I (hs-cTnI)
#' measurements from two different immunoassay platforms. This dataset
#' illustrates challenges in comparing troponin assays, which lack
#' standardization across manufacturers.
#'
#' @format A data frame with 50 observations and 3 variables:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifier.}
#'   \item{platform_a}{Numeric. Troponin I concentration (ng/L) measured by
#'     platform A.}
#'   \item{platform_b}{Numeric. Troponin I concentration (ng/L) measured by
#'     platform B.}
#' }
#'
#' @details
#' This synthetic dataset was designed to illustrate common patterns in
#' cardiac troponin method comparisons:
#'
#' \itemize{
#'   \item **Concentration range**: 2-5000 ng/L, from near the limit of
#'     detection to acute myocardial infarction levels
#'   \item **Distribution**: Log-normal (most values low, few very high),
#'     reflecting typical clinical populations
#'   \item **Bias pattern**: Systematic proportional difference between
#'     platforms (~15%), reflecting lack of troponin assay standardization
#'   \item **Precision**: Higher CV at low concentrations near the detection
#'     limit
#' }
#'
#' Unlike many analytes, cardiac troponin assays are not standardized, meaning
#' results from different manufacturers are not directly comparable. This has
#' clinical implications for interpreting troponin values when patients are
#' tested at different institutions.
#'
#' @source
#' Synthetic data generated to mimic realistic clinical patterns. See
#' `data-raw/make_datasets.R` for the generation script.
#'
#' @references
#' Apple FS, et al. Cardiac Troponin Assays: Guide to Understanding Analytical
#' Characteristics and Their Impact on Clinical Care. Clin Chem.
#' 2017;63(1):73-81.
#'
#' @examples
#' # Bland-Altman analysis with percent differences
#' # (appropriate for proportional bias)
#' ba <- ba_analysis(platform_a ~ platform_b,
#'                   data = troponin_cardiac,
#'                   type = "percent")
#' summary(ba)
#' plot(ba)
#'
#' @seealso [ba_analysis()], [glucose_methods], [creatinine_serum]
"troponin_cardiac"



#' High-Sensitivity Cardiac Troponin I Precision Dataset
#'
#' @description
#' Synthetic dataset for evaluating precision of a high-sensitivity cardiac
#' troponin I (hs-cTnI) assay across multiple concentration levels. The data
#' is designed for demonstrating precision studies with variance component
#' analysis and precision profile estimation.
#'
#' @format A data frame with 120 rows and 6 variables:
#' \describe{
#'   \item{level}{Concentration level factor (L1-L6)}
#'   \item{day}{Day of measurement (D1-D5)}
#'   \item{run}{Run within day (R1-R2)}
#'   \item{replicate}{Replicate within run (1-2)}
#'   \item{value}{Measured concentration (ng/L)}
#'   \item{target}{Nominal target concentration (ng/L)}
#' }
#'
#' @details
#' This synthetic dataset simulates a multi-level precision study for a
#' high-sensitivity cardiac troponin I assay. The data characteristics:
#'
#' \itemize{
#'   \item **Concentration levels**: 5, 10, 25, 50, 100, 500 ng/L
#'   \item **Design**: 5 days x 2 runs x 2 replicates per level (EP05-style)
#'   \item **Total observations**: 120 (6 levels x 20 measurements each)
#'   \item **Precision pattern**: CV decreases with concentration following
#'     a hyperbolic relationship
#' }
#'
#' The precision profile follows the model:
#' \deqn{CV = \sqrt{a^2 + (b/x)^2}}
#'
#' where approximately:
#' - \eqn{a \approx 3\%} (asymptotic CV at high concentrations)
#' - \eqn{b \approx 25} (concentration-dependent component)
#'
#' This gives expected CVs of approximately:
#' - 5 ng/L: ~5.5%
#' - 10 ng/L: ~3.5%
#' - 25 ng/L: ~3.0%
#' - 50 ng/L: ~2.0%
#' - 100 ng/L: ~2.0%
#' - 500 ng/L: ~3.5%
#'
#' @section Clinical Context:
#' High-sensitivity cardiac troponin assays are used to diagnose acute
#' myocardial infarction (AMI). Key clinical decision points include:
#' - 99th percentile upper reference limit (URL): typically 14-26 ng/L
#' - Functional sensitivity (CV <= 10%): should be <= 50% of 99th percentile
#' - Precision at low concentrations is critical for early AMI detection
#'
#' @source
#' Synthetic data generated to mimic realistic precision patterns. See
#' `data-raw/make_troponin_precision.R` for the generation script.
#'
#' @examples
#' # Load the dataset
#' data(troponin_precision)
#' head(troponin_precision)
#'
#' # Precision study with multiple concentration levels
#' prec <- precision_study(
#'   data = troponin_precision,
#'   value = "value",
#'   sample = "level",
#'   day = "day",
#'   run = "run"
#' )
#' print(prec)
#'
#' # Results for each concentration level
#' names(prec$by_sample)
#'
#' # Generate precision profile
#' profile <- precision_profile(prec, cv_targets = c(10, 20))
#' print(profile)
#' plot(profile)
#'
#' # Functional sensitivity at 10% CV
#' profile$functional_sensitivity
#'
#' @seealso
#' [precision_study()] for variance component analysis,
#' [precision_profile()] for CV-concentration modeling,
#' [glucose_methods], [creatinine_serum], [troponin_cardiac]
"troponin_precision"
