# smimodel: Sparse Multiple Index Models for Nonparametric Forecasting

Implements a general algorithm for estimating Sparse Multiple Index
(SMI) models for nonparametric forecasting and prediction. Estimation of
SMI models requires the Gurobi mixed integer programming (MIP) solver
via the gurobi R package. To use this functionality, the Gurobi
Optimizer must be installed, and a valid license obtained and activated
from <https://www.gurobi.com>. The gurobi R package must then be
installed and configured following the instructions at
<https://support.gurobi.com/hc/en-us/articles/14462206790033-How-do-I-install-Gurobi-for-R>.
The package also includes functions for fitting nonparametric additive
models with backward elimination, group-wise additive index models, and
projection pursuit regression models as benchmark comparison methods. In
addition, it provides tools for generating prediction intervals to
quantify uncertainty in point forecasts produced by the SMI model and
benchmark models, using the classical block bootstrap and a new method
called conformal bootstrap, which integrates block bootstrap with split
conformal prediction.

## See also

Useful links:

- <https://github.com/nuwani-palihawadana/smimodel>

- <https://nuwani-palihawadana.github.io/smimodel/>

- Report bugs at
  <https://github.com/nuwani-palihawadana/smimodel/issues>

## Author

**Maintainer**: Nuwani Palihawadana <nuwanipalihawadana@gmail.com>
([ORCID](https://orcid.org/0009-0008-6395-7797)) \[copyright holder\]

Other contributors:

- Xiaoqian Wang <Xiaoqian.Wang@amss.ac.cn>
  ([ORCID](https://orcid.org/0000-0003-4827-496X)) \[contributor\]
