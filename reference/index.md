# Package index

## Estimating SMI model

Functions to estimate a SMI model with or without penalty parameter
tuning.

- [`model_smimodel()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md)
  : Sparse Multiple Index (SMI) Models
- [`greedy_smimodel()`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy_smimodel.md)
  : SMI model estimation through a greedy search for penalty parameters

## Estimating other benchmark models

Functions to fit some benchmark comparison methods.

- [`model_backward()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md)
  : Nonparametric Additive Model with Backward Elimination
- [`model_gaim()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md)
  : Groupwise Additive Index Models (GAIM)
- [`model_ppr()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_ppr.md)
  : Projection Pursuit Regression (PPR) models
- [`model_gam()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gam.md)
  : Generalised Additive Models
- [`model_lm()`](https://nuwani-palihawadana.github.io/smimodel/reference/model_lm.md)
  : Linear Regression models

## Augment methods

Obtain residuals and fitted values of the fitted models.

- [`augment(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.smimodel.md)
  :

  Augment function for class `smimodel`

- [`augment(`*`<backward>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.backward.md)
  :

  Augment function for class `backward`

- [`augment(`*`<gaimFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.gaimFit.md)
  :

  Augment function for class `gaimFit`

- [`augment(`*`<pprFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.pprFit.md)
  :

  Augment function for class `pprFit`

- [`augment(`*`<gamFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.gamFit.md)
  :

  Augment function for class `gamFit`

- [`augment(`*`<lmFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.lmFit.md)
  :

  Augment function for class `lmFit`

## Predict methods

Obtain residuals and fitted values of the fitted models.

- [`predict(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.smimodel.md)
  :

  Obtaining forecasts on a test set from a fitted `smimodel`

- [`predict(`*`<backward>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.backward.md)
  :

  Obtaining forecasts on a test set from a fitted `backward`

- [`predict(`*`<gaimFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.gaimFit.md)
  :

  Obtaining forecasts on a test set from a fitted `gaimFit`

- [`predict(`*`<pprFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.pprFit.md)
  :

  Obtaining forecasts on a test set from a fitted `pprFit`

- [`predict(`*`<gamFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.gamFit.md)
  :

  Obtaining forecasts on a test set from a fitted `gamFit`

- [`predict(`*`<lmFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.lmFit.md)
  :

  Obtaining forecasts on a test set from a fitted `lmFit`

- [`predict(`*`<cgaim>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.cgaim.md)
  : Predictions from a fitted CGAIM object - copied from
  cgaim:::predict.cgaim() and modified

## Point estimate accuracy measures

Calculate point estimate accuracy measures.

- [`MSE()`](https://nuwani-palihawadana.github.io/smimodel/reference/point_measures.md)
  [`MAE()`](https://nuwani-palihawadana.github.io/smimodel/reference/point_measures.md)
  [`point_measures`](https://nuwani-palihawadana.github.io/smimodel/reference/point_measures.md)
  : Point estimate accuracy measures

## Prediction Intervals

Functions to construct and evaluate prediction intervals in time series
forecasting problems.

- [`bb_cvforecast()`](https://nuwani-palihawadana.github.io/smimodel/reference/bb_cvforecast.md)
  : Single season block bootstrap prediction intervals through time
  series cross-validation forecasting
- [`cb_cvforecast()`](https://nuwani-palihawadana.github.io/smimodel/reference/cb_cvforecast.md)
  : Conformal bootstrap prediction intervals through time series
  cross-validation forecasting
- [`avgCoverage()`](https://nuwani-palihawadana.github.io/smimodel/reference/avgCoverage.md)
  : Calculate interval forecast coverage
- [`avgWidth()`](https://nuwani-palihawadana.github.io/smimodel/reference/avgWidth.md)
  : Calculate interval forecast width

## Other user-facing functions

Other exported functions that can be called by users.

- [`autoplot(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/autoplot.smimodel.md)
  :

  Plot estimated smooths from a fitted `smimodel`

- [`residuals(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/residuals.smimodel.md)
  :

  Extract residuals from a fitted `smimodel`

- [`lag_matrix()`](https://nuwani-palihawadana.github.io/smimodel/reference/lag_matrix.md)
  : Function for adding lags of time series variables

- [`print(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/print.smimodel.md)
  :

  Printing a `smimodel` object

- [`print(`*`<smimodelFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/print.smimodelFit.md)
  :

  Printing a `smimodelFit` object

- [`print(`*`<backward>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/print.backward.md)
  :

  Printing a `backward` object

- [`print(`*`<gaimFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/print.gaimFit.md)
  :

  Printing a `gaimFit` object

- [`print(`*`<pprFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/print.pprFit.md)
  :

  Printing a `pprFit` object

- [`forecast(`*`<smimodel>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/forecast.smimodel.md)
  : Forecasting using SMI models

- [`forecast(`*`<backward>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/forecast.backward.md)
  : Forecasting using nonparametric additive models with backward
  elimination

- [`forecast(`*`<gaimFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/forecast.gaimFit.md)
  : Forecasting using GAIMs

- [`forecast(`*`<pprFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/forecast.pprFit.md)
  : Forecasting using PPR models

- [`forecast(`*`<gamFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/forecast.gamFit.md)
  : Forecasting using GAMs

- [`leadlagMat()`](https://nuwani-palihawadana.github.io/smimodel/reference/leadlagMat.md)
  : Create lags or leads of a matrix

## Other non-user-facing functions

Other internal functions that are not exported.

- [`smimodel.fit()`](https://nuwani-palihawadana.github.io/smimodel/reference/smimodel.fit.md)
  : SMI model estimation

- [`new_smimodelFit()`](https://nuwani-palihawadana.github.io/smimodel/reference/new_smimodelFit.md)
  :

  Constructor function for the class `smimodelFit`

- [`update_smimodelFit()`](https://nuwani-palihawadana.github.io/smimodel/reference/update_smimodelFit.md)
  :

  Updating a `smimodelFit`

- [`make_smimodelFit()`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)
  :

  Converting a fitted `gam` object to a `smimodelFit` object

- [`inner_update()`](https://nuwani-palihawadana.github.io/smimodel/reference/inner_update.md)
  : Updating index coefficients and non-linear functions iteratively

- [`init_alpha()`](https://nuwani-palihawadana.github.io/smimodel/reference/init_alpha.md)
  : Initialising index coefficients

- [`update_alpha()`](https://nuwani-palihawadana.github.io/smimodel/reference/update_alpha.md)
  : Updating index coefficients using MIP

- [`greedy.fit()`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy.fit.md)
  : Greedy search for tuning penalty parameters

- [`tune_smimodel()`](https://nuwani-palihawadana.github.io/smimodel/reference/tune_smimodel.md)
  : SMI model with a given penalty parameter combination

- [`normalise_alpha()`](https://nuwani-palihawadana.github.io/smimodel/reference/normalise_alpha.md)
  : Scaling index coefficient vectors to have unit norm

- [`loss()`](https://nuwani-palihawadana.github.io/smimodel/reference/loss.md)
  : Calculating the loss of the MIP used to estimate a SMI model

- [`allpred_index()`](https://nuwani-palihawadana.github.io/smimodel/reference/allpred_index.md)
  : Constructing index coefficient vectors with all predictors in each
  index

- [`split_index()`](https://nuwani-palihawadana.github.io/smimodel/reference/split_index.md)
  : Splitting predictors into multiple indices

- [`scaling()`](https://nuwani-palihawadana.github.io/smimodel/reference/scaling.md)
  : Scale data

- [`unscaling()`](https://nuwani-palihawadana.github.io/smimodel/reference/unscaling.md)
  :

  Unscale a fitted `smimodel`

- [`augment(`*`<smimodelFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/augment.smimodelFit.md)
  :

  Augment function for class `smimodelFit`

- [`predict(`*`<smimodelFit>`*`)`](https://nuwani-palihawadana.github.io/smimodel/reference/predict.smimodelFit.md)
  :

  Obtaining forecasts on a test set from a `smimodelFit`

- [`predict_gam()`](https://nuwani-palihawadana.github.io/smimodel/reference/predict_gam.md)
  :

  Obtaining recursive forecasts on a test set from a fitted
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

- [`eliminate()`](https://nuwani-palihawadana.github.io/smimodel/reference/eliminate.md)
  : Eliminate a variable and fit a nonparametric additive model

- [`blockBootstrap()`](https://nuwani-palihawadana.github.io/smimodel/reference/blockBootstrap.md)
  : Futures through single season block bootstrapping

- [`residBootstrap()`](https://nuwani-palihawadana.github.io/smimodel/reference/residBootstrap.md)
  : Generate multiple single season block bootstrap series

- [`seasonBootstrap()`](https://nuwani-palihawadana.github.io/smimodel/reference/seasonBootstrap.md)
  : Single season block bootstrap

- [`randomBlock()`](https://nuwani-palihawadana.github.io/smimodel/reference/randomBlock.md)
  : Randomly sampling a block

- [`possibleFutures_smimodel()`](https://nuwani-palihawadana.github.io/smimodel/reference/possibleFutures_smimodel.md)
  :

  Possible future sample paths (multi-step) from `smimodel` residuals

- [`possibleFutures_benchmark()`](https://nuwani-palihawadana.github.io/smimodel/reference/possibleFutures_benchmark.md)
  : Possible future sample paths (multi-step) from residuals of a fitted
  benchmark model

- [`prep_newdata()`](https://nuwani-palihawadana.github.io/smimodel/reference/prep_newdata.md)
  : Prepare a data set for recursive forecasting

- [`remove_lags()`](https://nuwani-palihawadana.github.io/smimodel/reference/remove_lags.md)
  : Remove actual values from a data set for recursive forecasting

- [`truncate_vars()`](https://nuwani-palihawadana.github.io/smimodel/reference/truncate_vars.md)
  : Truncating predictors to be in the in-sample range

## Package metadata

Package information.

- [`smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/smimodel-package.md)
  [`smimodel-package`](https://nuwani-palihawadana.github.io/smimodel/reference/smimodel-package.md)
  : smimodel: Sparse Multiple Index Models for Nonparametric Forecasting
