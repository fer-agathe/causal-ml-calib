# ______________________________________________________________________________
# Nuisance Outcome Estimators (g-models)
# ______________________________________________________________________________

#' Estimate Outcome Models Across Folds
#'
#' Fits chosen outcome model learners $g(d, X) = \mathbb{E}[Y \mid D=d, X]$ on each
#' training fold using cross-fitting.
#'
#' @param data A tibble containing the outcome `y`, binary treatment `d`, and
#'   covariate columns \eqn{X_1, \dots, X_p}.
#' @param learner_g_name String specifying the outcome model learner. Options:
#'   `"ols"`, `"lasso", `"rf"`, `"lgbm"`.
#' @param folds Integer vector of fold indices assigned to each row of `data`.
#' @param equilibrium Logical. If `TRUE`, also fits an equilibrium outcome
#'   regression model \eqn{g^{\text{eq}}(X) = \mathbb{E}[Y \mid X]} ignoring
#'   treatment assignment. Default is `FALSE`.
#'
#' @returns A list of length K (number of unique folds). Each element contains the
#'   fitted models returned by the specific estimation helper function.
#'
#' @md
estim_g <- function(data,
                    learner_g_name,
                    folds,
                    equilibrium = FALSE) {

  estim_g_helper <- match.fun(paste0("estim_g_", learner_g_name))

  K <- sort(unique(folds))
  fitted_g <- vector(mode = "list", length = length(K))

  for (k in K) {
    ind_train <- which(folds != k)
    fitted_g[[k]] <- estim_g_helper(
      x = data |> dplyr::select(-"y", -!!"d") |> dplyr::slice(ind_train),
      y = data$y[ind_train],
      d = data$d[ind_train],
      equilibrium = equilibrium
    )
  }

  fitted_g
}

# OLS----

#' OLS Learner for Outcome Models
#'
#' Fits separate Ordinary Least Squares (OLS) regressions for control (\eqn{D=0})
#' and treated (\eqn{D=1}) units, and optionally for the equilibrium model.
#'
#' @param y Numeric vector of outcomes.
#' @param x Matrix or data frame of covariates.
#' @param d Vector of binary treatment indicators (`0` or `1`).
#' @param equilibrium Logical. If `TRUE`, fits a pooled OLS model on all units.
#'
#' @returns A list containing:
#' * fitted_g_0: Fitted \code{\link[stats]{lm}} model for untreated units (\eqn{D=0}).
#' * fitted_g_1: Fitted \code{\link[stats]{lm}} model for treated units (\eqn{D=1}).
#' * fitted_g_eqb: Fitted \code{\link[stats]{lm}} model for all units, or `NULL`
#'   if `equilibrium = FALSE`.
#'
#' @md
estim_g_ols <- function(y,
                        x,
                        d,
                        equilibrium = FALSE) {
  ind_d_0 <- which(d == 0)

  tb_0 <- cbind(x[ind_d_0, , drop = FALSE], y = y[ind_d_0])
  tb_1 <- cbind(x[-ind_d_0, , drop = FALSE], y = y[-ind_d_0])

  fitted_g_0 <- stats::lm(y ~ ., data = tb_0)
  fitted_g_1 <- stats::lm(y ~ ., data = tb_1)

  if (isTRUE(equilibrium)) {
    tb <- cbind(x, y = y)
    fitted_g_eqb <- stats::lm(y ~ ., data = tb)
  } else {
    fitted_g_eqb <- NULL
  }

  list(
    fitted_g_0 = fitted_g_0,
    fitted_g_1 = fitted_g_1,
    fitted_g_eqb = fitted_g_eqb
  )
}


#' Predict Function for OLS Outcome Models
#'
#' @param object A fitted \code{\link[stats]{lm}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted outcomes.
#' @md
predict_g_ols <- function(object, new_data) {
  predict(object, newdata = new_data) |> as.vector()
}

# Lasso----

#' Lasso Learner for Outcome Models
#'
#' Fits Lasso linear regressions via \code{\link[hdm]{rlasso}} separately for
#' control (\eqn{D=0}) and treated (\eqn{D=1}) units, and optionally for the
#' equilibrium model.
#'
#' @param y Numeric vector of outcomes.
#' @param x Matrix or data frame of covariates.
#' @param d Vector of binary treatment indicators (`0` or `1`).
#' @param equilibrium Logical. If `TRUE`, fits a pooled Lasso model on all units.
#'
#' @returns A list containing:
#' * fitted_g_0: Fitted \code{\link[hdm]{rlasso}} model for untreated units (\eqn{D=0}).
#' * fitted_g_1: Fitted \code{\link[hdm]{rlasso}} model for treated units (\eqn{D=1}).
#' * fitted_g_eqb: Fitted \code{\link[hdm]{rlasso}} model for all units, or
#'   `NULL` if `equilibrium = FALSE`.
#'
#' @md
estim_g_lasso <- function(y,
                          x,
                          d,
                          equilibrium = FALSE) {
  ind_d_0 <- which(d == 0)

  fitted_g_0 <- hdm::rlasso(x = as.matrix(x[ind_d_0, , drop = FALSE]), y = y[ind_d_0])
  fitted_g_1 <- hdm::rlasso(x = as.matrix(x[-ind_d_0, , drop = FALSE]), y = y[-ind_d_0])

  if (isTRUE(equilibrium)) {
    fitted_g_eqb <- hdm::rlasso(x = as.matrix(x), y = y)
  } else {
    fitted_g_eqb <- NULL
  }

  list(
    fitted_g_0 = fitted_g_0,
    fitted_g_1 = fitted_g_1,
    fitted_g_eqb = fitted_g_eqb
  )
}


#' Predict Function for Lasso Outcome Models
#'
#' @param object A fitted \code{\link[hdm]{rlasso}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted outcomes.
#' @md
predict_g_lasso <- function(object, new_data) {
  predict(object, newdata = as.matrix(new_data)) |> as.vector()
}

# Random Forest----

#' Random Forest Learner for Outcome Models
#'
#' Fits Random Forest regression models via \code{\link[randomForest]{randomForest}}
#' separately for control (\eqn{D=0}) and treated (\eqn{D=1}) units, and
#' optionally for the equilibrium model.
#'
#' @param y Numeric vector of outcomes.
#' @param x Matrix or data frame of covariates.
#' @param d Vector of binary treatment indicators (`0` or `1`).
#' @param equilibrium Logical. If `TRUE`, fits a pooled Random Forest model on
#'   all units.
#'
#' @returns A list containing:
#' * fitted_g_0: Fitted \code{\link[randomForest]{randomForest}} model for
#'   untreated units (\eqn{D=0}).
#' * fitted_g_1: Fitted \code{\link[randomForest]{randomForest}} model for
#'   treated units (\eqn{D=1}).
#' * fitted_g_eqb: Fitted \code{\link[randomForest]{randomForest}} model for
#'   all units, or `NULL` if `equilibrium = FALSE`.
#'
#' @export
estim_g_rf <- function(y,
                       x,
                       d,
                       equilibrium = FALSE) {
  ind_d_0 <- which(d == 0)

  fitted_g_0 <- randomForest::randomForest(x = x[ind_d_0, , drop = FALSE], y = y[ind_d_0])
  fitted_g_1 <- randomForest::randomForest(x = x[-ind_d_0, , drop = FALSE], y = y[-ind_d_0])

  if (isTRUE(equilibrium)) {
    fitted_g_eqb <- randomForest::randomForest(x = x, y = y)
  } else {
    fitted_g_eqb <- NULL
  }

  list(
    fitted_g_0 = fitted_g_0,
    fitted_g_1 = fitted_g_1,
    fitted_g_eqb = fitted_g_eqb
  )
}


#' Predict Function for Random Forest Outcome Models
#'
#' @param object A fitted \code{\link[randomForest]{randomForest}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted outcomes.
#' @md
predict_g_rf <- function(object, new_data) {
  predict(object, newdata = new_data) |> as.vector()
}

# Boosting----

#' LightGBM Learner for Outcome Models
#'
#' Fits LightGBM gradient boosting regression models via \code{\link[lightgbm]{lgb.train}}
#' separately for control (\eqn{D=0}) and treated (\eqn{D=1}) units, and
#' optionally for the equilibrium model.
#'
#' @param y Numeric vector of outcomes.
#' @param x Matrix or data frame of covariates.
#' @param d Vector of binary treatment indicators (`0` or `1`).
#' @param equilibrium Logical. If `TRUE`, fits a pooled LightGBM model on all units.
#' @param boosting_type Boosting type. Default is `"gbdt"`.
#' @param num_leaves Maximum tree leaves for base learners. Default is `31`.
#' @param max_depth Maximum tree depth. Default is `-1` (no limit).
#' @param learning_rate Boosting learning rate. Default is `0.1`.
#' @param n_estimators Number of boosting iterations. Default is `100`.
#' @param subsample_for_bin Number of samples for constructing bins. Default is `200000`.
#' @param min_split_gain Minimum loss reduction required for split. Default is `0.0`.
#' @param min_child_weight Minimum sum of instance weight in a leaf. Default is `0.001`.
#' @param min_child_samples Minimum number of data points in a leaf. Default is `20`.
#' @param subsample Subsample ratio of training instances. Default is `1.0`.
#' @param subsample_freq Frequency of subsample. Default is `0`.
#' @param colsample_bytree Subsample ratio of columns per tree. Default is `1.0`.
#' @param reg_alpha L1 regularization term on weights. Default is `0.0`.
#' @param reg_lambda L2 regularization term on weights. Default is `0.0`.
#' @param random_state Optional seed for reproducibility.
#' @param n_jobs Number of parallel threads.
#'
#' @returns A list containing:
#' * fitted_g_0: Fitted \code{\link[lightgbm]{lgb.Booster}} model for untreated
#'   units (\eqn{D=0}).
#' * fitted_g_1: Fitted \code{\link[lightgbm]{lgb.Booster}} model for treated
#'   units (\eqn{D=1}).
#' * fitted_g_eqb: Fitted \code{\link[lightgbm]{lgb.Booster}} model for all
#'   units, or `NULL` if `equilibrium = FALSE`.
#'
#' @md
estim_g_lgbm <- function(y,
                         x,
                         d,
                         equilibrium = FALSE,
                         boosting_type = "gbdt",
                         num_leaves = 31,
                         max_depth = -1,
                         learning_rate = 0.1,
                         n_estimators = 100,
                         subsample_for_bin = 200000,
                         min_split_gain = 0.0,
                         min_child_weight = 0.001,
                         min_child_samples = 20,
                         subsample = 1.0,
                         subsample_freq = 0,
                         colsample_bytree = 1.0,
                         reg_alpha = 0.0,
                         reg_lambda = 0.0,
                         random_state = NULL,
                         n_jobs = NULL) {

  params <- list(
    boosting = boosting_type,
    objective = "regression",
    num_leaves = num_leaves,
    max_depth = max_depth,
    learning_rate = learning_rate,
    subsample_for_bin = subsample_for_bin,
    min_split_gain = min_split_gain,
    min_child_weight = min_child_weight,
    min_data_in_leaf = min_child_samples,
    bagging_fraction = subsample,
    bagging_freq = subsample_freq,
    feature_fraction = colsample_bytree,
    lambda_l1 = reg_alpha,
    lambda_l2 = reg_lambda
  )

  if (!is.null(random_state)) params$seed <- random_state
  if (!is.null(n_jobs)) params$num_threads <- n_jobs

  ind_d_0 <- which(d == 0)

  X0 <- as.matrix(x[ind_d_0, , drop = FALSE])
  y0 <- y[ind_d_0]
  X1 <- as.matrix(x[-ind_d_0, , drop = FALSE])
  y1 <- y[-ind_d_0]

  dtrain0 <- lightgbm::lgb.Dataset(data = X0, label = y0)
  dtrain1 <- lightgbm::lgb.Dataset(data = X1, label = y1)

  fitted_g_0 <- lightgbm::lgb.train(
    params = params,
    data = dtrain0,
    nrounds = n_estimators,
    verbose = -1
  )

  fitted_g_1 <- lightgbm::lgb.train(
    params = params,
    data = dtrain1,
    nrounds = n_estimators,
    verbose = -1
  )

  if (isTRUE(equilibrium)) {
    X <- as.matrix(x)
    dtrain <- lightgbm::lgb.Dataset(data = X, label = y)

    fitted_g_eqb <- lightgbm::lgb.train(
      params = params,
      data = dtrain,
      nrounds = n_estimators,
      verbose = -1
    )
  } else {
    fitted_g_eqb <- NULL
  }

  list(
    fitted_g_0 = fitted_g_0,
    fitted_g_1 = fitted_g_1,
    fitted_g_eqb = fitted_g_eqb
  )
}


#' Predict Function for LightGBM Outcome Models
#'
#' @param object A fitted \code{\link[lightgbm]{lgb.Booster}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted outcomes.
#' @md
predict_g_lgbm <- function(object, new_data) {
  preds <- predict(object, newdata = as.matrix(new_data))
  as.vector(preds)
}
