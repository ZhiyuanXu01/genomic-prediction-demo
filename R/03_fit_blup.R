library(here)
library(ggplot2)

# -------------------------------------------------------------------------
# Check required objects from the relationship-matrix stage
# -------------------------------------------------------------------------

required_objects <- c(
  "Y",
  "A_pedigree",
  "G_genomic"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Missing required objects: ",
      paste(missing_objects, collapse = ", "),
      ". Run R/02_relationship_matrices.R first."
    )
  )
}

# Environment 1 standardized grain yield.
y <- as.numeric(Y[, 1])
names(y) <- rownames(Y)

stopifnot(
  !anyNA(y),
  !anyDuplicated(names(y)),
  identical(names(y), rownames(A_pedigree)),
  identical(names(y), rownames(G_genomic))
)

# -------------------------------------------------------------------------
# Generic kernel-BLUP fitting function
# -------------------------------------------------------------------------

fit_kernel_blup <- function(
    y,
    K,
    training_ids,
    testing_ids,
    model_name
) {
  if (
    anyDuplicated(training_ids) ||
    anyDuplicated(testing_ids)
  ) {
    stop("Training and testing IDs must not contain duplicates.")
  }
  
  if (
    length(model_name) != 1 ||
    is.na(model_name) ||
    !nzchar(model_name)
  ) {
    stop("model_name must be one non-empty value.")
  }
  
  if (is.null(names(y))) {
    stop("y must be a named numeric vector.")
  }
  
  individual_ids <- names(y)
  
  if (
    nrow(K) != length(y) ||
    ncol(K) != length(y)
  ) {
    stop("K dimensions must match the length of y.")
  }
  
  if (
    !identical(rownames(K), individual_ids) ||
    !identical(colnames(K), individual_ids)
  ) {
    stop("The row and column names of K must match names(y).")
  }
  
  if (length(intersect(training_ids, testing_ids)) > 0) {
    stop("Training and testing IDs must not overlap.")
  }
  
  if (
    !setequal(
      c(training_ids, testing_ids),
      individual_ids
    )
  ) {
    stop(
      "Training and testing IDs must partition all individuals."
    )
  }
  
  if (max(abs(K - t(K))) > 1e-10) {
    stop("K must be symmetric.")
  }
  
  # Mask the testing phenotypes.
  y_masked <- y
  y_masked[testing_ids] <- NA_real_
  
  stopifnot(
    sum(!is.na(y_masked)) == length(training_ids),
    all(!is.na(y_masked[training_ids])),
    all(is.na(y_masked[testing_ids]))
  )
  
  # Explicit incidence matrix for individual genetic effects.
  Z <- diag(length(individual_ids))
  
  rownames(Z) <- individual_ids
  colnames(Z) <- individual_ids
  
  # Intercept-only fixed-effect design.
  X <- matrix(
    1,
    nrow = length(individual_ids),
    ncol = 1,
    dimnames = list(
      individual_ids,
      "Intercept"
    )
  )
  
  fitted_model <- rrBLUP::mixed.solve(
    y = y_masked,
    Z = Z,
    K = K,
    X = X,
    method = "REML"
  )
  
  intercept <- unname(fitted_model$beta[1])
  
  estimated_breeding_values <- setNames(
    as.numeric(fitted_model$u),
    individual_ids
  )
  
  predicted_phenotypes <- intercept +
    estimated_breeding_values
  
  testing_predictions <- data.frame(
    individual_id = testing_ids,
    observed = unname(y[testing_ids]),
    predicted = unname(
      predicted_phenotypes[testing_ids]
    ),
    estimated_breeding_value = unname(
      estimated_breeding_values[testing_ids]
    ),
    model = model_name,
    stringsAsFactors = FALSE
  )
  
  prediction_errors <- (
    testing_predictions$observed -
      testing_predictions$predicted
  )
  
  metrics <- data.frame(
    model = model_name,
    n_training = length(training_ids),
    n_testing = length(testing_ids),
    predictive_correlation = cor(
      testing_predictions$observed,
      testing_predictions$predicted
    ),
    rmse = sqrt(mean(prediction_errors^2)),
    mae = mean(abs(prediction_errors)),
    genetic_variance = fitted_model$Vu,
    residual_variance = fitted_model$Ve,
    variance_ratio = (
      fitted_model$Vu /
        (fitted_model$Vu + fitted_model$Ve)
    ),
    lambda = fitted_model$Ve / fitted_model$Vu,
    stringsAsFactors = FALSE
  )
  
  list(
    model = fitted_model,
    predictions = testing_predictions,
    metrics = metrics,
    intercept = intercept,
    breeding_values = estimated_breeding_values
  )
}

# -------------------------------------------------------------------------
# Fixed 80/20 train-test split
# -------------------------------------------------------------------------

fixed_split_seed <- 20260810

set.seed(fixed_split_seed)

number_testing <- round(
  0.20 * length(y)
)

testing_ids <- sample(
  names(y),
  size = number_testing,
  replace = FALSE
)

training_ids <- setdiff(
  names(y),
  testing_ids
)

stopifnot(
  length(training_ids) == length(y)-round(0.20 * length(y)),
  length(testing_ids) == round(0.20 * length(y)),
  length(intersect(training_ids, testing_ids)) == 0
)

split_assignments <- data.frame(
  individual_id = names(y),
  dataset = ifelse(
    names(y) %in% testing_ids,
    "testing",
    "training"
  ),
  seed = fixed_split_seed,
  stringsAsFactors = FALSE
)

write.csv(
  split_assignments,
  here(
    "results",
    "fixed_split_assignments.csv"
  ),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Fit PBLUP and GBLUP
# -------------------------------------------------------------------------

pblup_fixed <- fit_kernel_blup(
  y = y,
  K = A_pedigree,
  training_ids = training_ids,
  testing_ids = testing_ids,
  model_name = "PBLUP"
)

gblup_fixed <- fit_kernel_blup(
  y = y,
  K = G_genomic,
  training_ids = training_ids,
  testing_ids = testing_ids,
  model_name = "GBLUP"
)

baseline_metrics <- rbind(
  pblup_fixed$metrics,
  gblup_fixed$metrics
)

baseline_predictions <- rbind(
  pblup_fixed$predictions,
  gblup_fixed$predictions
)

stopifnot(
  nrow(baseline_metrics) == 2,
  nrow(baseline_predictions) ==
    2 * length(testing_ids),
  all(is.finite(
    baseline_metrics$predictive_correlation
  )),
  all(is.finite(baseline_metrics$rmse)),
  all(is.finite(baseline_metrics$mae))
)

write.csv(
  baseline_metrics,
  here(
    "results",
    "baseline_model_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  baseline_predictions,
  here(
    "results",
    "baseline_model_predictions.csv"
  ),
  row.names = FALSE
)

cat("\n--- Fixed-split model comparison ---\n")
print(baseline_metrics)

cat(
  "\nSaved fixed split and baseline model results.\n"
)

# -------------------------------------------------------------------------
# Visualise fixed-split predictions
# -------------------------------------------------------------------------

baseline_predictions$model <- factor(
  baseline_predictions$model,
  levels = c(
    "PBLUP",
    "GBLUP"
  )
)

baseline_metrics$model <- factor(
  baseline_metrics$model,
  levels = c(
    "PBLUP",
    "GBLUP"
  )
)

metric_labels <- data.frame(
  model = baseline_metrics$model,
  label = sprintf(
    "r = %.3f\nRMSE = %.3f\nMAE = %.3f",
    baseline_metrics$predictive_correlation,
    baseline_metrics$rmse,
    baseline_metrics$mae
  )
)

common_axis_limits <- range(
  c(
    baseline_predictions$observed,
    baseline_predictions$predicted
  ),
  finite = TRUE
)

baseline_prediction_plot <- ggplot(
  baseline_predictions,
  aes(
    x = observed,
    y = predicted
  )
) +
  geom_point(
    color = "#2b5c8f",
    alpha = 0.70,
    size = 1.8
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.7
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "#d95f02",
    linewidth = 0.9
  ) +
  geom_text(
    data = metric_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.10,
    vjust = 1.15,
    size = 3.7
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  coord_equal(
    xlim = common_axis_limits,
    ylim = common_axis_limits
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Fixed-Split Prediction Performance",
    subtitle = paste0(
      "Environment 1 grain yield; ",
      length(training_ids),
      " training and ",
      length(testing_ids),
      " testing lines"
    ),
    x = "Observed standardized grain yield",
    y = "Predicted standardized grain yield",
    caption = paste0(
      "Single random 80/20 split (seed ",
      fixed_split_seed,
      "); shown for pipeline validation, ",
      "not final model ranking."
    )
  )

ggsave(
  filename = here(
    "figures",
    "baseline_observed_vs_predicted.png"
  ),
  plot = baseline_prediction_plot,
  width = 9,
  height = 5,
  dpi = 300
)

cat(
  "Saved plot to: ",
  "figures/baseline_observed_vs_predicted.png\n"
)
