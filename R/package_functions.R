
#' Quantile-normalization of expression data
#'
#'
#' @param seqFISH_expression A cell by gene expression matrix of testing data
#' @param scRNAseq_expression A cell by gene expression matrix of training data data
#' @param normalize_scRNAseq Whether to first quantile-normalize the scRNA-seq data, default TRUE
#'
#'
#' @return Returns a list of the normalized matrices
#'
#'
#' @export
#'
#'


quantileNormalize <- function(
  seqFISH_expression,
  scRNAseq_expression,
  normalize_scRNAseq = TRUE
){

  seqFISH_genes <- colnames(seqFISH_expression)
  seqFISH_cells <- rownames(seqFISH_expression)
  scRNAseq_genes <- colnames(scRNAseq_expression)
  scRNAseq_cells <- rownames(scRNAseq_expression)

  if(normalize_scRNAseq) scRNAseq_expression <- preprocessCore::normalize.quantiles(scRNAseq_expression)

  scRNAseq_target <- preprocessCore::normalize.quantiles.determine.target(t(scRNAseq_expression))
  seqFISH_expression <- preprocessCore::normalize.quantiles.use.target(seqFISH_expression, scRNAseq_target)

  seqFISH_expression <- round(seqFISH_expression, 4)
  scRNAseq_expression <- round(scRNAseq_expression, 4)

  colnames(seqFISH_expression) <- seqFISH_genes
  rownames(seqFISH_expression) <- seqFISH_cells
  colnames(scRNAseq_expression) <- scRNAseq_genes
  rownames(scRNAseq_expression) <- scRNAseq_cells

  result <- list(seqFISH_expression = seqFISH_expression, scRNAseq_expression = scRNAseq_expression)

  return(result)

}


#' Train a multi-response elastic-net regression for prediction using
#' cross-validation to select lambda
#'
#'
#' @param training_x A cell by gene expression matrix of training predictor data
#' @param training_y A cell by gene expression matrix of training outcome data
#' @param prediction_x A cell by gene expression matrix of testing predictor data
#' @param lambda_mse Whether to select the lambda that minimizes the cross-validation, default TRUE
#' @param family family parameter passed on to the glmnet function, default 'mgaussian'
#' @param ... additional parameters to pass on to the glmnet function
#'
#'
#' @return Returns a list of the predicted outcome and model fit glmnet object
#'
#'
#' @export
#'
#'

predictExpression <- function(
  training_x,
  training_y,
  prediction_x,
  lambda_mse = TRUE,
  family = 'mgaussian',
  ...
){

  cat('computing lambda \n')

  cv <- glmnet::cv.glmnet(
    x = training_x,
    y = training_y,
    family = 'mgaussian',
    parallel = TRUE,
    ...
  )

  if(lambda_mse)
    best_lambda <- cv$lambda.min
  else
    best_lambda <- cv$lambda.se1

  cat('training model \n')

  en_fit <- glmnet::glmnet(
    x = training_x,
    y = training_y,
    lambda = best_lambda,
    family = 'mgaussian',
    ...
  )

  cat('predicting expression \n')

  outcome <- glmnet::predict.glmnet(
    object = en_fit,
    s = best_lambda,
    newx = prediction_x
  )

  if(family == 'mgaussian') outcome <- outcome[,,1]

  result <- list(outcome = outcome, fit = en_fit)

  return(result)

}


#' Determine estimable genes using a local polynomial regression
#' on model components
#'
#'
#' @param outcome A cell by gene estimated expression matrix
#' @param fit Model fit glmnet object
#' @param quantile_threshold Percentile at which to threshold estimate standard deviation, default 75\%
#'
#'
#' @return Returns a character vector of the estimable genes
#'
#'
#' @export
#'


determineEstimableGenes <- function(
  outcome,
  fit,
  quantile_threshold = .75
){

  predicted_genes <- colnames(outcome)
  l1norm <- unlist(lapply(fit$beta, FUN = function(x){sum(abs(x))}))
  stdev <- apply(outcome, 2, stats::sd)
  local_mean <- stats::predict(stats::loess(stdev ~ l1norm))
  threshold <- stats::quantile(stdev, probs = quantile_threshold)
  estimable_genes <- predicted_genes[stdev > local_mean & stdev >= threshold]

  return(estimable_genes)

}

