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
