params <- function(methods, ...) {
  defaults <- list(
    n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
    shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
    svm_type = "EPSILON_SVR", svm_kernel = "RBF", svm_cost = 10, svm_gamma = 0.05,
    maxent_beta = 1, maxent_features = "auto",
    knn_k = NULL, knn_search_method = NULL, knn_metric = NULL
  )
  do.call(build_method_params, c(list(methods = methods), modifyList(defaults, list(...))))
}

test_that("only the requested methods get an entry", {
  expect_named(params(c("rf", "maxent")), c("rf", "maxent"))
  expect_named(params("svm"), "svm")
})

test_that("random forest gets depth the shared defaults withhold", {
  rf <- params("rf")$rf
  expect_equal(rf$numberOfTrees, 500L)
  expect_equal(rf$variablesPerSplit, 8L)
  expect_equal(rf$minLeafPopulation, 1L)
  expect_null(rf$maxNodes)
})

test_that("boosting stays shallow with small steps", {
  gbt <- params("gbt")$gbt
  expect_equal(gbt$numberOfTrees, 150L)
  expect_equal(gbt$shrinkage, 0.05)
  expect_equal(gbt$maxNodes, 6L)
})

test_that("an explicit argument always beats the per-method override", {
  # Each override is conditional on the argument still holding its default, so a
  # caller who names a value must get that value back.
  expect_equal(params("rf",  n_trees = 42L)$rf$numberOfTrees, 42L)
  expect_equal(params("gbt", n_trees = 42L)$gbt$numberOfTrees, 42L)
  expect_equal(params("gbt", shrinkage = 0.2)$gbt$shrinkage, 0.2)
  expect_equal(params("rf",  variables_per_split = 3L)$rf$variablesPerSplit, 3L)
  expect_equal(params("rf",  max_nodes = 12L)$rf$maxNodes, 12L)
  expect_equal(params("rf",  min_leaf_population = 9L)$rf$minLeafPopulation, 9L)
})

test_that("libsvm arguments are passed through", {
  svm <- params("svm", svm_cost = 3, svm_gamma = 0.1)$svm
  expect_equal(svm$svmType, "EPSILON_SVR")
  expect_equal(svm$kernelType, "RBF")
  expect_equal(svm$cost, 3)
  expect_equal(svm$gamma, 0.1)
})

test_that("kNN tuning arguments apply only when supplied", {
  expect_null(params("knn")$knn$k)
  expect_null(params("knn")$knn$searchMethod)
  tuned <- params("knn", knn_k = 25, knn_search_method = "COVER_TREE", knn_metric = "MAHALANOBIS")$knn
  expect_equal(tuned$k, 25L)
  expect_equal(tuned$searchMethod, "COVER_TREE")
  expect_equal(tuned$metric, "MAHALANOBIS")
})

test_that("maxent feature classes switch off autoFeature", {
  expect_equal(params("maxent")$maxent$betaMultiplier, 1)
  expect_null(params("maxent")$maxent$autoFeature)
  lqh <- params("maxent", maxent_features = "LQH")$maxent
  expect_false(lqh$autoFeature)
  expect_true(lqh$quadratic)
  expect_true(lqh$hinge)
  expect_false(lqh$product)
})

test_that("evaluate_models and generate_map cannot drift apart", {
  # Both entry points build their parameters here, so the model that gets
  # cross-validated is the model that gets mapped. Guard that invariant.
  formals_of <- function(f) names(formals(f))
  shared <- c("n_trees", "min_leaf_population", "bag_fraction", "shrinkage",
              "max_nodes", "variables_per_split", "svm_type", "svm_kernel",
              "svm_cost", "svm_gamma", "maxent_beta", "maxent_features",
              "knn_k", "knn_search_method", "knn_metric")
  expect_true(all(shared %in% formals_of(evaluate_models)))
  expect_true(all(shared %in% formals_of(generate_map)))
})
