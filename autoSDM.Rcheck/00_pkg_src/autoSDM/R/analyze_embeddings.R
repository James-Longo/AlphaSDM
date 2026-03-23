#' Analyze Embeddings
#'
#' This function calculates the similarity between each location's embedding and the species centroid (mean embedding).
#'
#' @param df Data frame containing columns A00 to A63.
#' @param method Mapping method: 'centroid' or 'ridge'.
#' @return A list containing the results.
#' @keywords internal
analyze_embeddings <- function(df, method, gee_project = NULL, cv = FALSE, scale) {
  if (is.null(method)) stop("Parameter 'method' (e.g., 'centroid' or 'ridge') must be provided.")
  if (is.null(scale)) stop("Parameter 'scale' (resolution in meters) must be provided.")
  # 1. GEE Auth
  ensure_gee_authenticated(project = gee_project)

  # 2. Preparation (Upload and sample on GEE)
  prep <- prep_training_data_gee(df, class_property = "present", scale = scale)
  sampled_fc <- prep$fc

  # 3. Training
  message(sprintf("Training %s model...", method))
  model_res <- train_gee_model(sampled_fc, method, class_property = "present")

  # 4. Score on GEE (EXCLUDING raw embedding download)
  mean_emb <- model_res$weights
  dot_products <- NULL
  if (!model_res$is_classifier) {
    message("Calculating scores on GEE...")
    scored_fc <- predict_points_gee(sampled_fc, model_res)
    
    # Retrieve only the score property
    scores_list <- scored_fc$aggregate_array("score")$getInfo()
    dot_products <- as.numeric(scores_list)
  }

  return(list(
    mean_embedding = mean_emb,
    dot_products = dot_products,
    data = df,
    method = method,
    is_classifier = model_res$is_classifier
  ))
}
