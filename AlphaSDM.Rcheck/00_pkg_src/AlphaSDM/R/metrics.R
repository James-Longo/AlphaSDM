#' Calculate Continuous Boyce Index (CBI)
#' 
#' @param pos_scores Scores for presence points
#' @param all_scores Scores for all points (presence + background)
#' @param window_width Width of the moving window (proportion of score range)
#' @param n_bins Number of bins for evaluation
#' @export
calculate_cbi <- function(pos_scores, all_scores, window_width = 0.1, n_bins = 100) {
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
  
  # Spearman correlation
  correlation <- stats::cor(eval_points[valid], p_e, method = "spearman")
  
  if (is.na(correlation)) return(0.0)
  return(as.numeric(correlation))
}

#' Calculate Classifier Metrics
#' 
#' @param scores_pos Scores for presence points
#' @param scores_neg Scores for background/absence points
#' @export
calculate_classifier_metrics <- function(scores_pos, scores_neg) {
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
  
  # 1. AUC-ROC (Mann-Whitney U)
  all_labels <- c(rep(1, n_pos), rep(0, n_neg))
  ranks <- rank(scores_all, ties.method = "average")
  pos_ranks <- ranks[1:n_pos]
  auc_roc <- (sum(pos_ranks) - (n_pos * (n_pos + 1) / 2)) / (n_pos * n_neg)
  
  # 2. AUC-PRG (Precision-Recall Gain)
  ord <- order(scores_all, decreasing = TRUE)
  sorted_labels <- all_labels[ord]
  
  tp <- cumsum(sorted_labels)
  fp <- cumsum(1 - sorted_labels)
  
  # Formula: PrecGain = 1 - (fp/tp) / (n_neg/n_pos), RecGain = 1 - (fn/tp) / (n_neg/n_pos)
  fn <- n_pos - tp
  prec_gain <- 1 - (fp / tp) / (n_neg / n_pos)
  rec_gain <- 1 - (fn / tp) / (n_neg / n_pos)
  
  # Filter for non-negative gains (random classifier is 0)
  valid_idx <- which(prec_gain >= 0 & rec_gain >= 0)
  if (length(valid_idx) > 0) {
    # Add origin (0,0)
    pg <- c(0, prec_gain[valid_idx])
    rg <- c(0, rec_gain[valid_idx])
    # Sort for integration
    ord_prg <- order(rg)
    pg <- pg[ord_prg]
    rg <- rg[ord_prg]
    # Trapezoidal rule
    auc_prg <- sum(diff(rg) * (pg[-1] + pg[-length(pg)]) / 2)
  } else {
    auc_prg <- 0.0
  }
  
  # 3. TSS and Balanced Accuracy
  sens <- tp / n_pos
  spec <- (n_neg - fp) / n_neg
  tss <- max(sens + spec - 1)
  ba <- max((sens + spec) / 2)
  
  # 4. Pearson Correlation (COR)
  cor_val <- stats::cor(all_labels, scores_all, method = "pearson")
  
  return(list(
    cbi = as.numeric(cbi),
    auc_roc = as.numeric(auc_roc),
    auc_prg = as.numeric(auc_prg),
    tss = as.numeric(tss),
    ba = as.numeric(ba),
    cor = as.numeric(cor_val)
  ))
}
