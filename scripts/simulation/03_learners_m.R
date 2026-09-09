# ______________________________________________________________________________
# Nuisance Propensity Score Estimators (m-models)
# ______________________________________________________________________________

#' Estimate Propensity Score Models Across Folds
#'
#' Fits chosen propensity score model learners \eqn{m(X) = \Pr(D = 1 \mid X)}
#' on each training fold using cross-fitting.
#'
#' @param data A tibble containing the binary treatment
#'   column `d`, outcome column `y`, and covariate columns \eqn{X_1, \dots, X_p}.
#' @param learner_m_name String specifying the propensity score model learner.
#'   Options: `"logit"`, `"rf_classif"`, `"rf_reg"`, `"lgbm"`.
#' @param folds Integer vector of fold indices assigned to each row of `data}.
#'
#' @returns A list of length K (number of unique folds). Each element contains the
#'   fitted propensity score model returned by the specific estimation helper function.
#'
#' @md
estim_m <- function(data,
                    learner_m_name,
                    folds) {

  estim_m_helper <- match.fun(paste0("estim_m_", learner_m_name))

  K <- sort(unique(folds))
  fitted_m <- vector(mode = "list", length = length(K))

  for (k in K) {
    ind_train <- which(folds != k)
    fitted_m[[k]] <- estim_m_helper(
      d = data$d[ind_train],
      x = data |> dplyr::select(-"y", -!!"d") |> dplyr::slice(ind_train)
    )
  }

  fitted_m
}

# Logit----


#' Logistic Regression Learner for Propensity Score
#'
#' Fits a logistic regression model \eqn{\Pr(D = 1 \mid X) = \text{logit}^{-1}(X \beta)}
#' via \code{\link[stats]{glm}}.
#'
#' @param d Integer or numeric vector of binary treatment indicators (`0` or `1`).
#' @param x Matrix or data frame of covariates.
#'
#' @returns A fitted \code{\link[stats]{glm}} object with binomial family link.
#'
#' @md
estim_m_logit <- function(d, x) {
  tb <- cbind(x, d = d)
  stats::glm(d ~ ., data = tb, family = "binomial")
}


#' Predict Function for Logistic Regression Propensity Score
#'
#' @param object A fitted \code{\link[stats]{glm}} logistic model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted propensity scores bounded in \eqn{(0, 1)}.
#'
#' @md
predict_m_logit <- function(object, new_data) {
  predict(object, newdata = new_data, type = "response") |> as.vector()
}

# Random Forest Classifier----

#' Random Forest Classifier Learner for Propensity Score
#'
#' Fits a classification Random Forest model via \code{\link[randomForest]{randomForest}}
#' treating treatment status $D$ as a categorical factor.
#'
#' @param d Integer, numeric vector, or factor of binary treatment indicators (`0` or `1`).
#' @param x Matrix or data frame of covariates.
#'
#' @returns A fitted classification \code{\link[randomForest]{randomForest}} model.
#'
#' @md
estim_m_rf_classif <- function(d, x) {
  if (!is.factor(d)) d <- as.factor(as.character(d))
  randomForest::randomForest(x = x, y = d)
}


#' Predict Function for Random Forest Classifier Propensity Score
#'
#' @param object A fitted classification \code{\link[randomForest]{randomForest}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted class-1 probabilities.
#'
#' @md
predict_m_rf_classif <- function(object, new_data) {
  preds <- randomForest:::predict.randomForest(object, newdata = new_data, type = "prob")[, 2]
  as.vector(preds)
}

# Random Forest Regression----

#' Random Forest Regression Learner for Propensity Score
#'
#' Fits a regression Random Forest model via \code{\link[randomForest]{randomForest}}
#' treating treatment status $D$ as a continuous response in \eqn{\{0, 1\}}.
#'
#' @param d Numeric vector of binary treatment indicators (`0` or `1`).
#' @param x Matrix or data frame of covariates.
#'
#' @returns A fitted regression \code{\link[randomForest]{randomForest}} model.
#'
#' @md
estim_m_rf_reg <- function(d, x) {
  randomForest::randomForest(x = x, y = as.numeric(d))
}


#' Predict Function for Random Forest Regression Propensity Score
#'
#' @param object A fitted regression \code{\link[randomForest]{randomForest}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted propensity scores.
#'
#' @export
predict_m_rf_reg <- function(object, new_data) {
  predict(object, newdata = new_data, type = "response") |> as.vector()
}


# Boosting----

#' LightGBM Learner for Propensity Score
#'
#' Fits a gradient boosted binary classification model via \code{\link[lightgbm]{lgb.train}}.
#'
#' @param d Vector of binary treatment indicators (`0` or `1`).
#' @param x Matrix or data frame of covariates.
#' @param boosting_type Boosting type. Default is `"gbdt"`.
#' @param num_leaves Maximum tree leaves for base learners. Default is `31`.
#' @param max_depth Maximum tree depth. Default is `-1` (no limit).
#' @param learning_rate Boosting learning rate. Default is `0.1}.
#' @param n_estimators Number of boosting iterations. Default is `100}.
#' @param subsample_for_bin Number of samples for constructing bins. Default is `200000`.
#' @param min_split_gain Minimum loss reduction required for split. Default is `0.0`.
#' @param min_child_weight Minimum sum of instance weight in a leaf. Default is `0.001`.
#' @param min_child_samples Minimum number of data points in a leaf. Default is `20`.
#' @param subsample Subsample ratio of training instances. Default is `1.0`.
#' @param subsample_freq Frequency of subsample. Default is `0}.
#' @param colsample_bytree Subsample ratio of columns per tree. Default is `1.0`.
#' @param reg_alpha L1 regularization term on weights. Default is `0.0`.
#' @param reg_lambda L2 regularization term on weights. Default is `0.0`.
#' @param random_state Optional seed for reproducibility.
#' @param n_jobs Number of parallel threads.
#'
#' @returns A fitted \code{\link[lightgbm]{lgb.Booster}} binary classification model.
#'
#' @export
estim_m_lgbm <- function(d, x,
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

  X <- as.matrix(x)
  y <- as.numeric(d)

  dtrain <- lightgbm::lgb.Dataset(data = X, label = y)

  params <- list(
    boosting = boosting_type,
    objective = "binary",
    metric = "binary_logloss",
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

  lightgbm::lgb.train(
    params = params,
    data = dtrain,
    nrounds = n_estimators,
    verbose = -1
  )
}


#' Predict Function for LightGBM Propensity Score
#'
#' @param object A fitted \code{\link[lightgbm]{lgb.Booster}} model.
#' @param new_data Data frame or matrix of new covariate values.
#'
#' @returns Numeric vector of predicted propensity score probabilities.
#'
#' @export
predict_m_lgbm <- function(object, new_data) {
  preds <- predict(object, newdata = as.matrix(new_data))
  as.vector(preds)
}
