# Tests for classification metrics

test_that("calc_accuracy computes correctly", {
  actual    <- factor(c("A", "B", "A", "B", "A"))
  predicted <- factor(c("A", "B", "A", "A", "A"))

  # 4 out of 5 correct
  expect_equal(calc_accuracy(actual, predicted), 0.8)
})


test_that("calc_accuracy handles perfect prediction", {
  actual <- factor(c("A", "B", "A", "B"))
  expect_equal(calc_accuracy(actual, actual), 1.0)
})


test_that("calc_sensitivity computes correctly", {
  # TP=2, FN=1 → sensitivity = 2/3
  actual    <- factor(c("pos", "pos", "pos", "neg", "neg"))
  predicted <- factor(c("pos", "pos", "neg", "neg", "neg"))

  expect_equal(calc_sensitivity(actual, predicted, "pos"), 2/3)
})


test_that("calc_specificity computes correctly", {
  # TN=2, FP=0 → specificity = 1.0
  actual    <- factor(c("pos", "pos", "pos", "neg", "neg"))
  predicted <- factor(c("pos", "pos", "neg", "neg", "neg"))

  expect_equal(calc_specificity(actual, predicted, "pos"), 1.0)
})


test_that("calc_precision computes correctly", {
  # TP=2, FP=0 → precision = 1.0
  actual    <- factor(c("pos", "pos", "pos", "neg", "neg"))
  predicted <- factor(c("pos", "pos", "neg", "neg", "neg"))

  expect_equal(calc_precision(actual, predicted, "pos"), 1.0)
})


test_that("calc_f_measure computes correctly", {
  actual    <- factor(c("pos", "pos", "pos", "neg", "neg"))
  predicted <- factor(c("pos", "pos", "neg", "neg", "neg"))

  prec <- 1.0
  sens <- 2/3
  expected_f1 <- 2 * prec * sens / (prec + sens)
  expect_equal(calc_f_measure(actual, predicted, "pos"), expected_f1)
})


test_that("calc_cohens_kappa computes correctly", {
  # Perfect agreement → kappa = 1
  actual <- factor(c("A", "B", "A", "B"))
  expect_equal(calc_cohens_kappa(actual, actual), 1.0)

  # Random agreement → kappa ≈ 0
  set.seed(42)
  n <- 1000
  actual_big <- factor(sample(c("A", "B"), n, replace = TRUE))
  predicted_big <- factor(sample(c("A", "B"), n, replace = TRUE))
  kappa <- calc_cohens_kappa(actual_big, predicted_big)
  expect_true(abs(kappa) < 0.15)
})


test_that("calc_all_metrics returns all expected metrics", {
  actual    <- factor(c("pos", "pos", "neg", "neg"))
  predicted <- factor(c("pos", "pos", "neg", "neg"))

  metrics <- calc_all_metrics(actual, predicted, "pos")

  expect_type(metrics, "list")
  expected_names <- c("Accuracy", "Sensitivity", "Specificity",
                      "Precision", "F_measure", "Cohens_kappa")
  expect_true(all(expected_names %in% names(metrics)))

  # Perfect prediction → all should be 1.0
  expect_equal(metrics$Accuracy, 1.0)
  expect_equal(metrics$Sensitivity, 1.0)
  expect_equal(metrics$Specificity, 1.0)
  expect_equal(metrics$Precision, 1.0)
  expect_equal(metrics$F_measure, 1.0)
  expect_equal(metrics$Cohens_kappa, 1.0)
})


test_that("metrics handle edge case: no positive predictions", {
  actual    <- factor(c("pos", "pos", "neg", "neg"))
  predicted <- factor(c("neg", "neg", "neg", "neg"))

  # Precision undefined (no positive predictions) → NA
  expect_true(is.na(calc_precision(actual, predicted, "pos")))
  # Sensitivity = 0/2 = 0
  expect_equal(calc_sensitivity(actual, predicted, "pos"), 0.0)
  # Specificity = 2/2 = 1
  expect_equal(calc_specificity(actual, predicted, "pos"), 1.0)
})
