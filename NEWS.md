# sdmTMB

# sdmTMB 0.0.6.9005

* Fixed bug in predictions with standard errors where one(?)
  parameter (a breakpoint parameter) would be passed in at its initial
  instead of MLE value.

# sdmTMB 0.0.6.9004

* Fixed bug with predictions on new data in models with break points

* Overhauled cross validation function. The function now:
    * uses Eric's weights hack so it can also be used for forecasting
    * initializes subsequent folds at the MLE of the first fold for considerable speed increases
    * works in parallel if a future plan initialized; see examples
    
* Added threshold parameters to the print method

* Added forecasting example with the weights hack

* Fixed bug in linear break point models

#  sdmTMB 0.0.6.9002

* Fixed GAM predictions with all 0s in new data.

* Add linear and logistic threshold models. #17

# sdmTMB 0.0.5.9000

* Added parsing of mgcv formulas for splines. #16

* Added ability to predict with standard errors at the population level. This
  helps with making marginal-effect plots. #15

* Added optimization options to aid convergence. Also added
  `run_extra_optimization()` to run these on already fit models. Default is
  for no extra optimization.

* Added binomial likelihood to cross validation. Git hash `ee3f3ba`.

* Started keeping track of news in `NEWS.md`.
