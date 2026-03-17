#' Analyze Embeddings
#'
#' This function calculates the similarity between each location's embedding and the species centroid (mean embedding).
#'
#' @param df Data frame containing columns A00 to A63.
#' @param method Mapping method: 'centroid' or 'ridge'.
#' @return A list containing the results.
#' @keywords internal
analyze_embeddings <- function(df, method = "centroid", gee_project = NULL, cv = FALSE, scale = 10) {
  # 1. GEE Auth
  ensure_gee_authenticated(project = gee_project)

  # 2. Preparation (Sample embeddings from GEE)
  message("Extracting embeddings for analysis...")
  prep <- prep_training_data_gee(df, class_property = "present", scale = scale)
  sampled_fc <- prep$fc

  # 3. Training
  message(sprintf("Training %s model...", method))
  model_res <- train_gee_model(sampled_fc, method, class_property = "present")

  # 4. Results
  # For centroid, we return the mean embedding
  mean_emb <- model_res$weights
  
  # Calculate scores locally for the input data if possible
  emb_cols <- sprintf("A%02d", 0:63)
  emb_mat <- as.matrix(df[, emb_cols])
  
  dot_products <- NULL
  if (!model_res$is_classifier) {
    dot_products <- as.numeric((emb_mat %*% model_res$weights) + model_res$intercept)
  }

  return(list(
    mean_embedding = mean_emb,
    dot_products = dot_products,
    data = df,
    method = method,
    is_classifier = model_res$is_classifier
  ))
}
