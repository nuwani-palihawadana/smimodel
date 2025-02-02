#' @keywords package
#' @aliases smimodel-package
"_PACKAGE"

#' @importFrom cgaim cgaim g s
#' @importFrom conformalForecast coverage lagmatrix width
#' @importFrom dplyr arrange bind_cols bind_rows distinct filter group_split lag mutate mutate_at pull rename row_number select
#' @importFrom future multisession plan
#' @importFrom furrr future_map
#' @importFrom generics augment
#' @importFrom ggplot2 autoplot
#' @importFrom gratia derivatives draw
#' @importFrom gtools permutations
#' @importFrom Matrix bdiag
#' @importFrom methods as
#' @importFrom purrr map set_names
#' @importFrom ROI Q_objective L_constraint OP V_bound ROI_solve 
#' @importFrom stats as.formula as.ts end frequency gaussian lm model.frame na.omit predict ppr quantile runif sd start time ts tsp var window 
#' @importFrom tibble as_tibble is_tibble tibble
#' @importFrom tidyr drop_na
#' @importFrom tidyselect all_of
#' @importFrom tsibble as_tsibble index is_tsibble key
NULL 

# Generics to re-export

#' @export
generics::augment
#' @export
generics::forecast
#' @export
ggplot2::autoplot
