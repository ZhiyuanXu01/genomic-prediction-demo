library(here)
library(ggplot2)

# -------------------------------------------------------------------------
# Check objects created by previous stages
# -------------------------------------------------------------------------

required_objects <- c(
  "y",
  "A_pedigree",
  "G_genomic",
  "fit_kernel_blup"
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
      ". Run R/02_relationship_matrices.R and ",
      "R/03_fit_blup.R first."
    )
  )
}

# -------------------------------------------------------------------------
# Repeated-fold generator
# -------------------------------------------------------------------------

make_repeated_folds <- function(
  individual_ids,
  n_folds = 5L,
  n_repeats = 5L,
  seed = 20260811L
) {
  if (
    is.null(individual_ids) ||
      length(individual_ids) == 0
  ) {
    stop("individual_ids must not be empty.")
  }

  if (anyDuplicated(individual_ids)) {
    stop("individual_ids must be unique.")
  }

  if (
    n_folds < 2 ||
      n_repeats < 1 ||
      length(individual_ids) < n_folds
  ) {
    stop("Invalid numbers of folds or repeats.")
  }

  set.seed(seed)

  assignment_list <- lapply(
    seq_len(n_repeats),
    function(repeat_number) {
      shuffled_ids <- sample(
        individual_ids,
        replace = FALSE
      )

      fold_labels <- rep(
        seq_len(n_folds),
        length.out = length(individual_ids)
      )

      data.frame(
        individual_id = shuffled_ids,
        `repeat` = repeat_number,
        fold = fold_labels,
        seed = seed,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  assignments <- do.call(
    rbind,
    assignment_list
  )

  rownames(assignments) <- NULL

  assignments
}

# -------------------------------------------------------------------------
# Generate 5-fold cross-validation repeated five times
# -------------------------------------------------------------------------

cv_seed <- 20260811L
n_folds <- 5L
n_repeats <- 5L

cv_fold_assignments <- make_repeated_folds(
  individual_ids = names(y),
  n_folds = n_folds,
  n_repeats = n_repeats,
  seed = cv_seed
)

cv_fold_assignments <- cv_fold_assignments[
  order(
    cv_fold_assignments$`repeat`,
    cv_fold_assignments$fold,
    cv_fold_assignments$individual_id
  ),
]

rownames(cv_fold_assignments) <- NULL

# -------------------------------------------------------------------------
# Validate fold assignments
# -------------------------------------------------------------------------

stopifnot(
  nrow(cv_fold_assignments) ==
    length(y) * n_repeats,
  setequal(
    cv_fold_assignments$individual_id,
    names(y)
  ),
  setequal(
    cv_fold_assignments$`repeat`,
    seq_len(n_repeats)
  ),
  setequal(
    cv_fold_assignments$fold,
    seq_len(n_folds)
  )
)

individual_repeat_counts <- table(
  cv_fold_assignments$individual_id,
  cv_fold_assignments$`repeat`
)

stopifnot(
  all(individual_repeat_counts == 1)
)

fold_size_table <- as.data.frame(
  table(
    `repeat` = cv_fold_assignments$`repeat`,
    fold = cv_fold_assignments$fold
  )
)

names(fold_size_table)[3] <- "n_testing"

expected_minimum_size <- floor(
  length(y) / n_folds
)

expected_maximum_size <- ceiling(
  length(y) / n_folds
)

stopifnot(
  all(
    fold_size_table$n_testing %in%
      expected_minimum_size:expected_maximum_size
  )
)

# -------------------------------------------------------------------------
# Check phenotype balance across folds
# -------------------------------------------------------------------------

fold_groups <- split(
  cv_fold_assignments,
  interaction(
    cv_fold_assignments$`repeat`,
    cv_fold_assignments$fold,
    drop = TRUE
  )
)

cv_fold_diagnostics <- do.call(
  rbind,
  lapply(
    fold_groups,
    function(fold_data) {
      fold_phenotypes <- y[
        fold_data$individual_id
      ]

      data.frame(
        `repeat` = fold_data$`repeat`[1],
        fold = fold_data$fold[1],
        n_testing = length(fold_phenotypes),
        phenotype_mean = mean(fold_phenotypes),
        phenotype_sd = sd(fold_phenotypes),
        phenotype_minimum = min(fold_phenotypes),
        phenotype_maximum = max(fold_phenotypes),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

cv_fold_diagnostics <- cv_fold_diagnostics[
  order(
    cv_fold_diagnostics$`repeat`,
    cv_fold_diagnostics$fold
  ),
]

rownames(cv_fold_diagnostics) <- NULL

stopifnot(
  nrow(cv_fold_diagnostics) ==
    n_folds * n_repeats,
  all(cv_fold_diagnostics$n_testing %in% c(floor(length(y) / n_folds), ceiling(length(y) / n_folds))),
  all(is.finite(cv_fold_diagnostics$phenotype_mean)),
  all(is.finite(cv_fold_diagnostics$phenotype_sd))
)

# -------------------------------------------------------------------------
# Save fold definitions and diagnostics
# -------------------------------------------------------------------------

write.csv(
  cv_fold_assignments,
  here(
    "results",
    "repeated_cv_fold_assignments.csv"
  ),
  row.names = FALSE
)

write.csv(
  cv_fold_diagnostics,
  here(
    "results",
    "repeated_cv_fold_diagnostics.csv"
  ),
  row.names = FALSE
)

cat("\n--- Repeated cross-validation fold sizes ---\n")
print(fold_size_table)

cat("\n--- Fold phenotype diagnostics ---\n")
print(cv_fold_diagnostics)

cat(
  "\nSaved repeated CV fold assignments and diagnostics.\n"
)

# -------------------------------------------------------------------------
# Fit PBLUP and GBLUP across all repeated CV folds
# -------------------------------------------------------------------------

cv_model_kernels <- list(
  PBLUP = A_pedigree,
  GBLUP = G_genomic
)

number_cv_fits <- (
  n_repeats *
    n_folds *
    length(cv_model_kernels)
)

cv_metric_list <- vector(
  "list",
  number_cv_fits
)

cv_prediction_list <- vector(
  "list",
  number_cv_fits
)

fit_index <- 0L

for (repeat_number in seq_len(n_repeats)) {
  for (fold_number in seq_len(n_folds)) {
    testing_ids_cv <- cv_fold_assignments$individual_id[
      cv_fold_assignments[["repeat"]] ==
        repeat_number &
        cv_fold_assignments$fold ==
        fold_number
    ]
    
    training_ids_cv <- setdiff(
      names(y),
      testing_ids_cv
    )
    
    stopifnot(
      length(testing_ids_cv) %in% c(floor(length(y) / n_folds), ceiling(length(y) / n_folds)),
      length(training_ids_cv) +
        length(testing_ids_cv) ==
        length(y),
      length(
        intersect(
          training_ids_cv,
          testing_ids_cv
        )
      ) == 0
    )
    
    for (model_name in names(cv_model_kernels)) {
      fit_index <- fit_index + 1L
      
      cat(
        sprintf(
          "\nFitting %s: repeat %d, fold %d (%d/%d)\n",
          model_name,
          repeat_number,
          fold_number,
          fit_index,
          number_cv_fits
        )
      )
      
      fitted_cv_model <- fit_kernel_blup(
        y = y,
        K = cv_model_kernels[[model_name]],
        training_ids = training_ids_cv,
        testing_ids = testing_ids_cv,
        model_name = model_name
      )
      
      fold_metrics <- fitted_cv_model$metrics
      fold_metrics[["repeat"]] <- repeat_number
      fold_metrics$fold <- fold_number
      
      fold_metrics <- fold_metrics[
        ,
        c(
          "repeat",
          "fold",
          setdiff(
            names(fold_metrics),
            c("repeat", "fold")
          )
        )
      ]
      
      fold_predictions <- fitted_cv_model$predictions
      fold_predictions[["repeat"]] <- repeat_number
      fold_predictions$fold <- fold_number
      
      fold_predictions <- fold_predictions[
        ,
        c(
          "repeat",
          "fold",
          setdiff(
            names(fold_predictions),
            c("repeat", "fold")
          )
        )
      ]
      
      cv_metric_list[[fit_index]] <- fold_metrics
      cv_prediction_list[[fit_index]] <- fold_predictions
    }
  }
}

stopifnot(
  fit_index == number_cv_fits
)

repeated_cv_metrics <- do.call(
  rbind,
  cv_metric_list
)

repeated_cv_predictions <- do.call(
  rbind,
  cv_prediction_list
)

rownames(repeated_cv_metrics) <- NULL
rownames(repeated_cv_predictions) <- NULL

repeated_cv_metrics <- repeated_cv_metrics[
  order(
    repeated_cv_metrics[["repeat"]],
    repeated_cv_metrics$fold,
    match(
      repeated_cv_metrics$model,
      names(cv_model_kernels)
    )
  ),
]

repeated_cv_predictions <- repeated_cv_predictions[
  order(
    repeated_cv_predictions[["repeat"]],
    repeated_cv_predictions$fold,
    match(
      repeated_cv_predictions$model,
      names(cv_model_kernels)
    ),
    repeated_cv_predictions$individual_id
  ),
]

rownames(repeated_cv_metrics) <- NULL
rownames(repeated_cv_predictions) <- NULL

# -------------------------------------------------------------------------
# Validate repeated CV results
# -------------------------------------------------------------------------

metric_keys <- paste(
  repeated_cv_metrics[["repeat"]],
  repeated_cv_metrics$fold,
  repeated_cv_metrics$model,
  sep = "::"
)

prediction_keys <- paste(
  repeated_cv_predictions$individual_id,
  repeated_cv_predictions[["repeat"]],
  repeated_cv_predictions$model,
  sep = "::"
)

assignment_keys <- paste(
  cv_fold_assignments$individual_id,
  cv_fold_assignments[["repeat"]],
  cv_fold_assignments$fold,
  sep = "::"
)

prediction_assignment_keys <- paste(
  repeated_cv_predictions$individual_id,
  repeated_cv_predictions[["repeat"]],
  repeated_cv_predictions$fold,
  sep = "::"
)

observed_from_y <- unname(
  y[repeated_cv_predictions$individual_id]
)

stopifnot(
  nrow(repeated_cv_metrics) ==
    n_repeats *
    n_folds *
    length(cv_model_kernels),
  nrow(repeated_cv_predictions) ==
    length(y) *
    n_repeats *
    length(cv_model_kernels),
  !anyDuplicated(metric_keys),
  !anyDuplicated(prediction_keys),
  all(
    prediction_assignment_keys %in%
      assignment_keys
  ),
  setequal(
    repeated_cv_metrics$model,
    names(cv_model_kernels)
  ),
  setequal(
    repeated_cv_predictions$model,
    names(cv_model_kernels)
  ),
  all(
    repeated_cv_metrics$n_training +
      repeated_cv_metrics$n_testing ==
      length(y)
  ),
  all(
    repeated_cv_metrics$n_testing %in%
      c(floor(length(y) / n_folds), ceiling(length(y) / n_folds))
  ),
  all(
    is.finite(
      repeated_cv_metrics$predictive_correlation
    )
  ),
  all(is.finite(repeated_cv_metrics$rmse)),
  all(is.finite(repeated_cv_metrics$mae)),
  all(
    is.finite(
      repeated_cv_metrics$genetic_variance
    )
  ),
  all(
    is.finite(
      repeated_cv_metrics$residual_variance
    )
  ),
  all(
    is.finite(
      repeated_cv_predictions$predicted
    )
  ),
  isTRUE(
    all.equal(
      repeated_cv_predictions$observed,
      observed_from_y,
      check.attributes = FALSE
    )
  )
)

# -------------------------------------------------------------------------
# Save repeated CV model results
# -------------------------------------------------------------------------

write.csv(
  repeated_cv_metrics,
  here(
    "results",
    "repeated_cv_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  repeated_cv_predictions,
  here(
    "results",
    "repeated_cv_predictions.csv"
  ),
  row.names = FALSE
)

cat("\n--- Repeated CV metrics ---\n")
print(repeated_cv_metrics)

cat(
  "\nSaved 50 model fits and their predictions.\n"
)

# -------------------------------------------------------------------------
# Summarise repeated CV performance
# -------------------------------------------------------------------------

primary_cv_metrics <- c(
  "predictive_correlation",
  "rmse",
  "mae"
)

cv_summary_list <- list()
summary_index <- 0L

for (model_name in names(cv_model_kernels)) {
  model_results <- repeated_cv_metrics[
    repeated_cv_metrics$model == model_name,
  ]
  
  for (metric_name in primary_cv_metrics) {
    summary_index <- summary_index + 1L
    
    metric_values <- model_results[[metric_name]]
    
    cv_summary_list[[summary_index]] <- data.frame(
      model = model_name,
      metric = metric_name,
      n_splits = length(metric_values),
      mean = mean(metric_values),
      standard_deviation = sd(metric_values),
      median = median(metric_values),
      minimum = min(metric_values),
      maximum = max(metric_values),
      stringsAsFactors = FALSE
    )
  }
}

repeated_cv_summary <- do.call(
  rbind,
  cv_summary_list
)

rownames(repeated_cv_summary) <- NULL

# -------------------------------------------------------------------------
# Calculate paired fold differences
# -------------------------------------------------------------------------

pblup_fold_metrics <- repeated_cv_metrics[
  repeated_cv_metrics$model == "PBLUP",
  c(
    "repeat",
    "fold",
    primary_cv_metrics
  )
]

gblup_fold_metrics <- repeated_cv_metrics[
  repeated_cv_metrics$model == "GBLUP",
  c(
    "repeat",
    "fold",
    primary_cv_metrics
  )
]

paired_fold_metrics <- merge(
  pblup_fold_metrics,
  gblup_fold_metrics,
  by = c("repeat", "fold"),
  suffixes = c("_PBLUP", "_GBLUP"),
  sort = TRUE
)

paired_cv_differences <- data.frame(
  `repeat` = paired_fold_metrics[["repeat"]],
  fold = paired_fold_metrics$fold,
  correlation_difference = (
    paired_fold_metrics$predictive_correlation_GBLUP -
      paired_fold_metrics$predictive_correlation_PBLUP
  ),
  rmse_difference = (
    paired_fold_metrics$rmse_GBLUP -
      paired_fold_metrics$rmse_PBLUP
  ),
  mae_difference = (
    paired_fold_metrics$mae_GBLUP -
      paired_fold_metrics$mae_PBLUP
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

paired_cv_summary <- data.frame(
  metric = primary_cv_metrics,
  difference_definition = "GBLUP minus PBLUP",
  mean_difference = c(
    mean(
      paired_cv_differences$
        correlation_difference
    ),
    mean(
      paired_cv_differences$
        rmse_difference
    ),
    mean(
      paired_cv_differences$
        mae_difference
    )
  ),
  standard_deviation = c(
    sd(
      paired_cv_differences$
        correlation_difference
    ),
    sd(
      paired_cv_differences$
        rmse_difference
    ),
    sd(
      paired_cv_differences$
        mae_difference
    )
  ),
  gblup_better_splits = c(
    sum(
      paired_cv_differences$
        correlation_difference > 0
    ),
    sum(
      paired_cv_differences$
        rmse_difference < 0
    ),
    sum(
      paired_cv_differences$
        mae_difference < 0
    )
  ),
  total_splits = nrow(paired_cv_differences),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(repeated_cv_summary) ==
    length(primary_cv_metrics) *
    length(cv_model_kernels),
  all(repeated_cv_summary$n_splits == 25),
  nrow(paired_cv_differences) == 25,
  nrow(paired_cv_summary) == 3,
  all(
    is.finite(
      paired_cv_summary$mean_difference
    )
  )
)

write.csv(
  repeated_cv_summary,
  here(
    "results",
    "repeated_cv_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_cv_differences,
  here(
    "results",
    "repeated_cv_paired_differences.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_cv_summary,
  here(
    "results",
    "repeated_cv_paired_summary.csv"
  ),
  row.names = FALSE
)

cat("\n--- Repeated CV summary ---\n")
print(repeated_cv_summary)

cat("\n--- Paired model comparison ---\n")
print(paired_cv_summary)

# -------------------------------------------------------------------------
# Visualise paired repeated CV performance
# -------------------------------------------------------------------------

metric_display_labels <- c(
  predictive_correlation =
    "Predictive correlation\n(higher is better)",
  rmse =
    "RMSE\n(lower is better)",
  mae =
    "MAE\n(lower is better)"
)

cv_plot_data <- do.call(
  rbind,
  lapply(
    primary_cv_metrics,
    function(metric_name) {
      data.frame(
        `repeat` =
          repeated_cv_metrics[["repeat"]],
        fold =
          repeated_cv_metrics$fold,
        model =
          repeated_cv_metrics$model,
        metric =
          metric_display_labels[[metric_name]],
        value =
          repeated_cv_metrics[[metric_name]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

cv_plot_data$split_id <- paste(
  cv_plot_data[["repeat"]],
  cv_plot_data$fold,
  sep = "-"
)

cv_plot_data$model <- factor(
  cv_plot_data$model,
  levels = c(
    "PBLUP",
    "GBLUP"
  )
)

cv_plot_data$metric <- factor(
  cv_plot_data$metric,
  levels = unname(metric_display_labels)
)

repeated_cv_plot <- ggplot(
  cv_plot_data,
  aes(
    x = model,
    y = value
  )
) +
  geom_line(
    aes(group = split_id),
    color = "grey70",
    linewidth = 0.4,
    alpha = 0.65
  ) +
  geom_boxplot(
    aes(fill = model),
    width = 0.45,
    alpha = 0.20,
    outlier.shape = NA,
    color = "grey35",
    linewidth = 0.6
  ) +
  geom_point(
    aes(color = model),
    size = 1.8,
    alpha = 0.80
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_color_manual(
    values = c(
      PBLUP = "#2b5c8f",
      GBLUP = "#d95f02"
    )
  ) +
  scale_fill_manual(
    values = c(
      PBLUP = "#2b5c8f",
      GBLUP = "#d95f02"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Repeated Cross-Validation Performance",
    subtitle = paste0(
      "Five repeats of five-fold cross-validation; ",
      "lines connect models evaluated on the same split"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Points represent 25 repeat-fold combinations. ",
      "Higher correlation and lower error indicate ",
      "better prediction."
    )
  )

ggsave(
  filename = here(
    "figures",
    "repeated_cv_model_comparison.png"
  ),
  plot = repeated_cv_plot,
  width = 11,
  height = 5,
  dpi = 300
)

cat(
  "\nSaved repeated CV summaries and comparison plot.\n"
)
