# valytics 0.1.0

Initial CRAN release.
 
## Features

* `ba_analysis()`: Bland-Altman analysis for method comparison
  - Vector and formula interfaces
  - Absolute and percent difference types
  - Confidence intervals for bias and limits of agreement (Bland & Altman 1999)
  - `print()`, `summary()`, and `plot()` methods
  
* `autoplot.ba_analysis()`: ggplot2-style plotting

* Example datasets for method comparison:
  - `glucose_methods`: POC glucose meter vs laboratory analyzer (n=60)
  - `creatinine_serum`: Enzymatic vs Jaffe creatinine methods (n=80)
  - `troponin_cardiac`: Two high-sensitivity troponin platforms (n=50)