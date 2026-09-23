test_that("each model starts from Earth Engine's defaults plus AlphaSDM's", {
  expect_identical(build_gee_clf_params("rf"), list(numberOfTrees = 500L))
  expect_identical(build_gee_clf_params("gbt"), list(numberOfTrees = 150L))
  expect_identical(build_gee_clf_params("knn"), list(k = 15L))
  # No AlphaSDM defaults: Earth Engine's apply untouched.
  expect_identical(build_gee_clf_params("svm"), list())
  expect_identical(build_gee_clf_params("maxent"), list())
  expect_identical(build_gee_clf_params("cart"), list())
})

test_that("user settings override defaults and reach only their model", {
  expect_equal(build_gee_clf_params("gbt", list(shrinkage = 0.01))$shrinkage, 0.01)
  expect_equal(build_gee_clf_params("gbt", list(shrinkage = 0.01))$numberOfTrees, 150L)
  # Values equal to an old default are honoured, not replaced.
  expect_identical(build_gee_clf_params("rf", list(numberOfTrees = 100))$numberOfTrees, 100L)
  s <- method_settings(c("svm", "rf", "gbt"), list(gbt = list(shrinkage = 0.01)))
  expect_null(s$rf)
  expect_equal(s$gbt$shrinkage, 0.01)
})

test_that("integer settings are coerced and a linear SVM drops gamma", {
  expect_type(build_gee_clf_params("rf", list(numberOfTrees = 300))$numberOfTrees, "integer")
  expect_type(build_gee_clf_params("knn", list(k = 25))$k, "integer")
  expect_null(build_gee_clf_params("svm", list(kernelType = "LINEAR", gamma = 0.1))$gamma)
  expect_true(svm_is_regression(list(svmType = "EPSILON_SVR")))
  expect_false(svm_is_regression(list()))   # Earth Engine's default is C_SVC
})

test_that("settings for a model that is not being fitted are an error", {
  expect_error(method_settings(c("svm", "rf"), list(gbm = list(shrinkage = 0.1))), "gbm")
  expect_error(method_settings("rf", list(list(numberOfTrees = 10))), "named list")
})

test_that("evaluate_models and generate_map take the same settings", {
  expect_true("params" %in% names(formals(evaluate_models)))
  expect_true("params" %in% names(formals(generate_map)))
})
