#' Calculate the Continuous Boyce Index
#'
#' Measures how far the ratio of predicted to expected presences rises with
#' suitability. A well calibrated model gives a value near 1, a random one near 0.
#'
#' @param pos_scores Numeric suitability scores at presence points. NA is dropped.
#' @param all_scores Numeric suitability scores at all points, presence and
#'   background together. NA is dropped.
#' @param window_width Width of the moving window, as a proportion of the score
#'   range.
#' @param n_bins Number of window positions to evaluate.
#' @return A single number in [-1, 1], the Spearman correlation between window
#'   position and the predicted-to-expected ratio. Returns 0 when there are no
#'   presence scores, when all scores are equal, or when the correlation is
#'   undefined.
#' @export
calculate_cbi <- function(pos_scores, all_scores, window_width = 0.1, n_bins = 100) {
  pos_scores <- pos_scores[!is.na(pos_scores)]
  all_scores <- all_scores[!is.na(all_scores)]
  if (length(pos_scores) == 0) return(0.0)
  
  min_score <- min(all_scores, na.rm = TRUE)
  max_score <- max(all_scores, na.rm = TRUE)
  score_range <- max_score - min_score
  
  if (score_range < 1e-9) return(0.0)
  
  eval_points <- seq(min_score, max_score, length.out = n_bins)
  half_width <- score_range * window_width / 2
  
  f_obs <- numeric(n_bins)
  f_exp <- numeric(n_bins)
  
  for (i in seq_along(eval_points)) {
    p <- eval_points[i]
    low <- p - half_width
    high <- p + half_width
    
    n_pos <- sum(pos_scores >= low & pos_scores <= high)
    n_all <- sum(all_scores >= low & all_scores <= high)
    
    f_obs[i] <- n_pos / length(pos_scores)
    f_exp[i] <- n_all / length(all_scores)
  }
  
  valid <- f_exp > 0
  if (!any(valid)) return(0.0)
  
  p_e <- f_obs[valid] / f_exp[valid]
  
  correlation <- stats::cor(eval_points[valid], p_e, method = "spearman")
  
  if (is.na(correlation)) return(0.0)
  return(as.numeric(correlation))
}

#' Calculate discrimination and calibration metrics
#'
#' @param scores_pos Numeric suitability scores at presence points. NA is dropped.
#' @param scores_neg Numeric suitability scores at background or absence points.
#'   NA is dropped.
#' @return A named list: `cbi`, `auc_roc`, `auc_prg`, `tss`, `ba` and `cor`, each a
#'   single number. Scores need only be on a common scale within one call, since
#'   every metric except `cor` depends on the ranking alone. When either class is
#'   empty the list is filled with the no-skill values.
#' @export
calculate_classifier_metrics <- function(scores_pos, scores_neg) {
  scores_pos <- scores_pos[!is.na(scores_pos)]
  scores_neg <- scores_neg[!is.na(scores_neg)]
  scores_all <- c(scores_pos, scores_neg)
  cbi <- calculate_cbi(scores_pos, scores_all)
  
  n_pos <- length(scores_pos)
  n_neg <- length(scores_neg)
  
  if (n_pos == 0 || n_neg == 0) {
    return(list(
      cbi = cbi,
      auc_roc = 0.5,
      auc_prg = 0.0,
      tss = 0.0,
      ba = 0.5,
      cor = 0.0
    ))
  }
  
  # AUC-ROC (Mann-Whitney U)
  all_labels <- c(rep(1, n_pos), rep(0, n_neg))
  ranks <- rank(scores_all, ties.method = "average")
  pos_ranks <- ranks[1:n_pos]
  auc_roc <- (sum(pos_ranks) - (n_pos * (n_pos + 1) / 2)) / (n_pos * n_neg)

  # AUC-PRG (Precision-Recall Gain)
  ord <- order(scores_all, decreasing = TRUE)
  sorted_labels <- all_labels[ord]
  sorted_scores <- scores_all[ord]
  
  tp <- cumsum(sorted_labels)
  fp <- cumsum(1 - sorted_labels)
  
  # PrecGain = 1 - (fp/tp) / (n_neg/n_pos), RecGain = 1 - (fn/tp) / (n_neg/n_pos)
  fn <- n_pos - tp
  prec_gain <- 1 - (fp / tp) / (n_neg / n_pos)
  rec_gain <- 1 - (fn / tp) / (n_neg / n_pos)
  
  # Keep only non-negative gains: precision-recall-gain space is defined so a
  # random classifier sits at 0, and negative values fall outside the unit square.
  valid_idx <- which(prec_gain >= 0 & rec_gain >= 0)
  if (length(valid_idx) > 0) {
    pg <- c(0, prec_gain[valid_idx])
    rg <- c(0, rec_gain[valid_idx])
    ord_prg <- order(rg)
    pg <- pg[ord_prg]
    rg <- rg[ord_prg]
    # Trapezoidal rule over the gain curve.
    auc_prg <- sum(diff(rg) * (pg[-1] + pg[-length(pg)]) / 2)
  } else {
    auc_prg <- 0.0
  }
  
  # Evaluate sensitivity and specificity only at score boundaries. A threshold has
  # to put every point sharing a score on the same side of it, so a sweep that
  # splits a run of tied scores is not a threshold any model could apply. Without
  # this, tied or constant scores reach TSS = 1, because presences come first in
  # scores_all and are counted before the absences they are tied with.
  keep <- c(sorted_scores[-length(sorted_scores)] != sorted_scores[-1], TRUE)
  sens <- tp[keep] / n_pos
  spec <- (n_neg - fp[keep]) / n_neg
  tss <- max(sens + spec - 1)
  ba <- max((sens + spec) / 2)
  
  # A degenerate model can return one constant score for every point. Pearson
  # correlation is undefined at zero variance and warns, so report no association.
  cor_val <- if (stats::sd(scores_all) < 1e-12) 0 else
    stats::cor(all_labels, scores_all, method = "pearson")
  
  return(list(
    cbi = as.numeric(cbi),
    auc_roc = as.numeric(auc_roc),
    auc_prg = as.numeric(auc_prg),
    tss = as.numeric(tss),
    ba = as.numeric(ba),
    cor = as.numeric(cor_val)
  ))
}
