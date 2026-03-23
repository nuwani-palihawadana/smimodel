## Resubmission

This is a resubmission following feedback from CRAN.

* No published references are currently available for the implemented method; therefore none were added to DESCRIPTION.
* Added @return tags for exported print S3 methods.
* The vignette is pre-computed due to dependency on the Gurobi commercial MIP solver; therefore, the code in vignette is non-executable. 
* Examples requiring Gurobi are wrapped in \dontrun{} as they depend on the external Gurobi Optimization solver and cannot be run without a valid license.
* Made informational console messages suppressible.

Additionally, while revising the package, I identified and fixed two minor bugs:
* Fixed compatibility of output from model_smimodel() as input to cb_cvforecast() and bb_cvforecast().
* Fixed an error when using h = 1 in cb_cvforecast().