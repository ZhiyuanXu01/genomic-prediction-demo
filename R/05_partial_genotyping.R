library(here)
library(ggplot2)

# -------------------------------------------------------------------------
# Check objects created by previous stages
# -------------------------------------------------------------------------

required_objects <- c(
  "y",
  "cv_fold_assignments"
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
      ". Run the preceding analysis scripts first."
    )
  )
}

stopifnot(
  setequal(
    cv_fold_assignments$individual_id,
    names(y)
  ),
  all(
    table(
      cv_fold_assignments$individual_id,
      cv_fold_assignments[["repeat"]]
    ) == 1
  )
)

# -------------------------------------------------------------------------
# Define a fixed 50% genotyping scenario
# -------------------------------------------------------------------------

partial_genotyping_seed <- 20260812L
genotyped_proportion <- 0.50

number_genotyped <- round(
  length(y) * genotyped_proportion
)

set.seed(partial_genotyping_seed)

genotyped_ids <- sample(
  names(y),
  size = number_genotyped,
  replace = FALSE
)

partial_genotyping_assignments <- data.frame(
  individual_id = names(y),
  genotyped = names(y) %in% genotyped_ids,
  genotyping_status = ifelse(
    names(y) %in% genotyped_ids,
    "genotyped",
    "not_genotyped"
  ),
  selection_proportion = genotyped_proportion,
  selection_seed = partial_genotyping_seed,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(partial_genotyping_assignments) ==
    length(y),
  !anyDuplicated(
    partial_genotyping_assignments$individual_id
  ),
  sum(
    partial_genotyping_assignments$genotyped
  ) == number_genotyped,
  sum(
    !partial_genotyping_assignments$genotyped
  ) == length(y) - number_genotyped
)

# -------------------------------------------------------------------------
# Check representation within every CV testing fold
# -------------------------------------------------------------------------

partial_cv_membership <- merge(
  cv_fold_assignments[
    ,
    c(
      "individual_id",
      "repeat",
      "fold"
    )
  ],
  partial_genotyping_assignments[
    ,
    c(
      "individual_id",
      "genotyped",
      "genotyping_status"
    )
  ],
  by = "individual_id",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  nrow(partial_cv_membership) ==
    nrow(cv_fold_assignments),
  !anyNA(partial_cv_membership$genotyped),
  !anyNA(
    partial_cv_membership$genotyping_status
  )
)

partial_cv_membership$genotyping_status <- factor(
  partial_cv_membership$genotyping_status,
  levels = c(
    "genotyped",
    "not_genotyped"
  )
)

genotyping_count_table <- table(
  partial_cv_membership[["repeat"]],
  partial_cv_membership$fold,
  partial_cv_membership$genotyping_status
)

partial_genotyping_fold_diagnostics <- as.data.frame(
  genotyping_count_table,
  stringsAsFactors = FALSE
)

names(
  partial_genotyping_fold_diagnostics
) <- c(
  "repeat",
  "fold",
  "genotyping_status",
  "n_testing"
)

partial_genotyping_fold_diagnostics[["repeat"]] <-
  as.integer(
    partial_genotyping_fold_diagnostics[["repeat"]]
  )

partial_genotyping_fold_diagnostics$fold <-
  as.integer(
    partial_genotyping_fold_diagnostics$fold
  )

partial_genotyping_fold_diagnostics$n_testing <-
  as.integer(
    partial_genotyping_fold_diagnostics$n_testing
  )

partial_genotyping_fold_diagnostics <-
  partial_genotyping_fold_diagnostics[
    order(
      partial_genotyping_fold_diagnostics[["repeat"]],
      partial_genotyping_fold_diagnostics$fold,
      partial_genotyping_fold_diagnostics$
        genotyping_status
    ),
  ]

rownames(
  partial_genotyping_fold_diagnostics
) <- NULL

number_folds <- length(
  unique(cv_fold_assignments$fold)
)

expected_minimum_fold_size <- floor(
  length(y) / number_folds
)

expected_maximum_fold_size <- ceiling(
  length(y) / number_folds
)

testing_fold_totals <- tapply(
  partial_genotyping_fold_diagnostics$n_testing,
  list(
    partial_genotyping_fold_diagnostics[["repeat"]],
    partial_genotyping_fold_diagnostics$fold
  ),
  sum
)

stopifnot(
  nrow(partial_genotyping_fold_diagnostics) ==
    length(
      unique(
        cv_fold_assignments[["repeat"]]
      )
    ) *
    number_folds *
    2,
  all(
    partial_genotyping_fold_diagnostics$n_testing >
      0
  ),
  all(
    testing_fold_totals %in%
      expected_minimum_fold_size:
      expected_maximum_fold_size
  )
)

# -------------------------------------------------------------------------
# Save the fixed partial-genotyping scenario
# -------------------------------------------------------------------------

write.csv(
  partial_genotyping_assignments,
  here(
    "results",
    "partial_genotyping_assignments.csv"
  ),
  row.names = FALSE
)

write.csv(
  partial_genotyping_fold_diagnostics,
  here(
    "results",
    "partial_genotyping_fold_diagnostics.csv"
  ),
  row.names = FALSE
)

cat("\n--- Partial-genotyping scenario ---\n")
cat(
  "Genotyped individuals:",
  sum(partial_genotyping_assignments$genotyped),
  "\n"
)
cat(
  "Non-genotyped individuals:",
  sum(!partial_genotyping_assignments$genotyped),
  "\n"
)

cat("\n--- Testing-fold representation ---\n")
print(partial_genotyping_fold_diagnostics)

cat(
  "\nSaved partial-genotyping assignments ",
  "and fold diagnostics.\n"
)

# -------------------------------------------------------------------------
# Construct G22 using only genotyped individuals
# -------------------------------------------------------------------------

matrix_required_objects <- c(
  "M",
  "A_pedigree",
  "align_G_to_A"
)

missing_matrix_objects <- matrix_required_objects[
  !vapply(
    matrix_required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_matrix_objects) > 0) {
  stop(
    paste0(
      "Missing matrix objects: ",
      paste(missing_matrix_objects, collapse = ", "),
      ". Run R/02_relationship_matrices.R first."
    )
  )
}

# Preserve the common ordering used by y and A.
genotyped_ids <- names(y)[
  names(y) %in% genotyped_ids
]

not_genotyped_ids <- setdiff(
  names(y),
  genotyped_ids
)

stopifnot(
  length(genotyped_ids) == number_genotyped,
  length(not_genotyped_ids) ==
    length(y) - number_genotyped,
  identical(
    rownames(M),
    names(y)
  ),
  identical(
    rownames(A_pedigree),
    names(y)
  )
)

M_genotyped <- M[
  genotyped_ids,
  ,
  drop = FALSE
]

# Estimate marker frequencies using only the genotyped subset.
partial_marker_frequency <- colMeans(
  M_genotyped
)

partial_marker_maf <- pmin(
  partial_marker_frequency,
  1 - partial_marker_frequency
)

partial_maf_threshold <- 0.05

keep_partial_maf <- (
  partial_marker_maf >=
    partial_maf_threshold
)

M_genotyped_filtered <- M_genotyped[
  ,
  keep_partial_maf,
  drop = FALSE
]

# The binary marker states represent the two homozygous
# classes in these highly inbred wheat lines.
M_genotyped_rrblup <- (
  2 * M_genotyped_filtered - 1
)

stopifnot(
  nrow(M_genotyped_rrblup) ==
    number_genotyped,
  ncol(M_genotyped_rrblup) > 0,
  !anyNA(M_genotyped_rrblup),
  all(
    M_genotyped_rrblup %in%
      c(-1, 1)
  )
)

G22_raw <- rrBLUP::A.mat(
  X = M_genotyped_rrblup,
  min.MAF = NULL,
  max.missing = 0,
  impute.method = "mean",
  shrink = FALSE
)

A22 <- A_pedigree[
  genotyped_ids,
  genotyped_ids,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(G22_raw),
    genotyped_ids
  ),
  identical(
    colnames(G22_raw),
    genotyped_ids
  ),
  identical(
    rownames(A22),
    genotyped_ids
  ),
  identical(
    colnames(A22),
    genotyped_ids
  )
)

# Align G22 to the pedigree scale using the same
# genotyped subset in both matrices.
partial_alignment <- align_G_to_A(
  G = G22_raw,
  A = A22
)

G22_aligned <- partial_alignment$G
partial_alignment_alpha <- partial_alignment$alpha
partial_alignment_beta <- partial_alignment$beta

G22_raw_eigenvalues <- eigen(
  G22_raw,
  symmetric = TRUE,
  only.values = TRUE
)$values

G22_aligned_eigenvalues <- eigen(
  G22_aligned,
  symmetric = TRUE,
  only.values = TRUE
)$values

stopifnot(
  max(
    abs(
      G22_aligned -
        t(G22_aligned)
    )
  ) < 1e-10,
  min(G22_aligned_eigenvalues) >
    1e-8,
  abs(
    mean(diag(G22_aligned)) -
      mean(diag(A22))
  ) < 1e-10,
  abs(
    mean(
      G22_aligned[
        lower.tri(G22_aligned)
      ]
    ) -
      mean(
        A22[
          lower.tri(A22)
        ]
      )
  ) < 1e-10
)

# -------------------------------------------------------------------------
# Construct H inverse and H
# -------------------------------------------------------------------------

A_inverse <- solve(
  A_pedigree
)

A22_inverse <- solve(
  A22
)

G22_inverse <- solve(
  G22_aligned
)

H_inverse <- A_inverse

H_inverse[
  genotyped_ids,
  genotyped_ids
] <- (
  H_inverse[
    genotyped_ids,
    genotyped_ids
  ] +
    G22_inverse -
    A22_inverse
)

# Remove negligible numerical asymmetry.
H_inverse <- (
  H_inverse + t(H_inverse)
) / 2

H_relationship_raw <- solve(
  H_inverse
)

H_asymmetry_before_symmetrization <- max(
  abs(
    H_relationship_raw -
      t(H_relationship_raw)
  )
)

H_relationship <- (
  H_relationship_raw +
    t(H_relationship_raw)
) / 2

H_asymmetry_after_symmetrization <- max(
  abs(
    H_relationship -
      t(H_relationship)
  )
)

rownames(H_relationship) <- names(y)
colnames(H_relationship) <- names(y)

H_inverse_eigenvalues <- eigen(
  H_inverse,
  symmetric = TRUE,
  only.values = TRUE
)$values

H_eigenvalues <- eigen(
  H_relationship,
  symmetric = TRUE,
  only.values = TRUE
)$values

H_inverse_residual <- max(
  abs(
    H_inverse %*%
      H_relationship -
      diag(length(y))
  )
)

stopifnot(
  identical(
    rownames(H_relationship),
    names(y)
  ),
  identical(
    colnames(H_relationship),
    names(y)
  ),
  all(is.finite(H_inverse)),
  all(is.finite(H_relationship)),
  max(
    abs(
      H_inverse -
        t(H_inverse)
    )
  ) < 1e-10,
  max(
    abs(
      H_relationship -
        t(H_relationship)
    )
  ) < 1e-10,
  min(H_inverse_eigenvalues) >
    1e-8,
  min(H_eigenvalues) >
    1e-8,
  H_inverse_residual < 1e-8
)

# -------------------------------------------------------------------------
# Save matrix diagnostics
# -------------------------------------------------------------------------

partial_matrix_diagnostics <- data.frame(
  metric = c(
    "individuals_total",
    "individuals_genotyped",
    "individuals_not_genotyped",
    "markers_total",
    "markers_retained_partial_maf_0.05",
    "markers_removed_partial_maf_0.05",
    "G22_raw_mean_diagonal",
    "G22_raw_mean_off_diagonal",
    "G22_raw_minimum_eigenvalue",
    "G22_alignment_alpha",
    "G22_alignment_beta",
    "G22_aligned_mean_diagonal",
    "G22_aligned_mean_off_diagonal",
    "G22_aligned_minimum_eigenvalue",
    "H_inverse_minimum_eigenvalue",
    "H_minimum_eigenvalue",
    "H_maximum_eigenvalue",
    "H_asymmetry_before_symmetrization",
    "H_asymmetry_after_symmetrization",
    "H_inverse_reconstruction_residual"
  ),
  value = c(
    length(y),
    length(genotyped_ids),
    length(not_genotyped_ids),
    ncol(M),
    sum(keep_partial_maf),
    sum(!keep_partial_maf),
    mean(diag(G22_raw)),
    mean(
      G22_raw[
        lower.tri(G22_raw)
      ]
    ),
    min(G22_raw_eigenvalues),
    partial_alignment_alpha,
    partial_alignment_beta,
    mean(diag(G22_aligned)),
    mean(
      G22_aligned[
        lower.tri(G22_aligned)
      ]
    ),
    min(G22_aligned_eigenvalues),
    min(H_inverse_eigenvalues),
    min(H_eigenvalues),
    max(H_eigenvalues),
    H_asymmetry_before_symmetrization,
    H_asymmetry_after_symmetrization,
    H_inverse_residual
  ),
  stringsAsFactors = FALSE
)

write.csv(
  partial_matrix_diagnostics,
  here(
    "results",
    "partial_genotyping_matrix_diagnostics.csv"
  ),
  row.names = FALSE
)

cat("\n--- Partial-genotyping matrix diagnostics ---\n")
print(partial_matrix_diagnostics)

cat(
  "\nConstructed G22, H inverse, and H successfully.\n"
)

# -------------------------------------------------------------------------
# Fit ssGBLUP across the repeated CV folds
# -------------------------------------------------------------------------

ssgblup_required_objects <- c(
  "fit_kernel_blup",
  "repeated_cv_metrics",
  "repeated_cv_predictions"
)

missing_ssgblup_objects <- ssgblup_required_objects[
  !vapply(
    ssgblup_required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_ssgblup_objects) > 0) {
  stop(
    paste0(
      "Missing cross-validation objects: ",
      paste(missing_ssgblup_objects, collapse = ", "),
      ". Run R/03_fit_blup.R and ",
      "R/04_cross_validation.R first."
    )
  )
}

cv_repeat_values <- sort(
  unique(
    cv_fold_assignments[["repeat"]]
  )
)

cv_fold_values <- sort(
  unique(
    cv_fold_assignments$fold
  )
)

number_ssgblup_fits <- (
  length(cv_repeat_values) *
    length(cv_fold_values)
)

ssgblup_metric_list <- vector(
  "list",
  number_ssgblup_fits
)

ssgblup_prediction_list <- vector(
  "list",
  number_ssgblup_fits
)

ssgblup_fit_index <- 0L

for (repeat_number in cv_repeat_values) {
  for (fold_number in cv_fold_values) {
    ssgblup_fit_index <-
      ssgblup_fit_index + 1L
    
    testing_ids_ssgblup <-
      cv_fold_assignments$individual_id[
        cv_fold_assignments[["repeat"]] ==
          repeat_number &
          cv_fold_assignments$fold ==
          fold_number
      ]
    
    training_ids_ssgblup <- setdiff(
      names(y),
      testing_ids_ssgblup
    )
    
    stopifnot(
      length(
        intersect(
          training_ids_ssgblup,
          testing_ids_ssgblup
        )
      ) == 0,
      setequal(
        c(
          training_ids_ssgblup,
          testing_ids_ssgblup
        ),
        names(y)
      )
    )
    
    cat(
      sprintf(
        paste0(
          "\nFitting ssGBLUP: repeat %d, ",
          "fold %d (%d/%d)\n"
        ),
        repeat_number,
        fold_number,
        ssgblup_fit_index,
        number_ssgblup_fits
      )
    )
    
    fitted_ssgblup <- fit_kernel_blup(
      y = y,
      K = H_relationship,
      training_ids = training_ids_ssgblup,
      testing_ids = testing_ids_ssgblup,
      model_name = "ssGBLUP"
    )
    
    ssgblup_fold_metrics <-
      fitted_ssgblup$metrics
    
    ssgblup_fold_metrics[["repeat"]] <-
      repeat_number
    
    ssgblup_fold_metrics$fold <-
      fold_number
    
    ssgblup_fold_metrics <-
      ssgblup_fold_metrics[
        ,
        c(
          "repeat",
          "fold",
          setdiff(
            names(ssgblup_fold_metrics),
            c("repeat", "fold")
          )
        )
      ]
    
    ssgblup_fold_predictions <-
      fitted_ssgblup$predictions
    
    ssgblup_fold_predictions[["repeat"]] <-
      repeat_number
    
    ssgblup_fold_predictions$fold <-
      fold_number
    
    ssgblup_fold_predictions <-
      ssgblup_fold_predictions[
        ,
        c(
          "repeat",
          "fold",
          setdiff(
            names(ssgblup_fold_predictions),
            c("repeat", "fold")
          )
        )
      ]
    
    ssgblup_metric_list[[
      ssgblup_fit_index
    ]] <- ssgblup_fold_metrics
    
    ssgblup_prediction_list[[
      ssgblup_fit_index
    ]] <- ssgblup_fold_predictions
  }
}

stopifnot(
  ssgblup_fit_index ==
    number_ssgblup_fits
)

ssgblup_cv_metrics <- do.call(
  rbind,
  ssgblup_metric_list
)

ssgblup_cv_predictions <- do.call(
  rbind,
  ssgblup_prediction_list
)

rownames(ssgblup_cv_metrics) <- NULL
rownames(ssgblup_cv_predictions) <- NULL

# -------------------------------------------------------------------------
# Combine ssGBLUP with the existing paired PBLUP results
# -------------------------------------------------------------------------

pblup_partial_metrics <-
  repeated_cv_metrics[
    repeated_cv_metrics$model == "PBLUP",
  ]

pblup_partial_predictions <-
  repeated_cv_predictions[
    repeated_cv_predictions$model == "PBLUP",
  ]

partial_genotyping_cv_metrics <- rbind(
  pblup_partial_metrics,
  ssgblup_cv_metrics
)

partial_genotyping_cv_predictions <- rbind(
  pblup_partial_predictions,
  ssgblup_cv_predictions
)

partial_model_order <- c(
  "PBLUP",
  "ssGBLUP"
)

partial_genotyping_cv_metrics <-
  partial_genotyping_cv_metrics[
    order(
      partial_genotyping_cv_metrics[["repeat"]],
      partial_genotyping_cv_metrics$fold,
      match(
        partial_genotyping_cv_metrics$model,
        partial_model_order
      )
    ),
  ]

partial_genotyping_cv_predictions <-
  partial_genotyping_cv_predictions[
    order(
      partial_genotyping_cv_predictions[["repeat"]],
      partial_genotyping_cv_predictions$fold,
      match(
        partial_genotyping_cv_predictions$model,
        partial_model_order
      ),
      partial_genotyping_cv_predictions$
        individual_id
    ),
  ]

rownames(partial_genotyping_cv_metrics) <- NULL
rownames(partial_genotyping_cv_predictions) <- NULL

# Attach the fixed genotyping status to every prediction.
genotyping_status_lookup <- setNames(
  partial_genotyping_assignments$
    genotyping_status,
  partial_genotyping_assignments$
    individual_id
)

partial_genotyping_cv_predictions$
  genotyping_status <- unname(
    genotyping_status_lookup[
      partial_genotyping_cv_predictions$
        individual_id
    ]
  )

# -------------------------------------------------------------------------
# Calculate performance by genotyping status
# -------------------------------------------------------------------------

status_prediction_groups <- split(
  partial_genotyping_cv_predictions,
  list(
    partial_genotyping_cv_predictions[["repeat"]],
    partial_genotyping_cv_predictions$fold,
    partial_genotyping_cv_predictions$model,
    partial_genotyping_cv_predictions$
      genotyping_status
  ),
  drop = TRUE
)

partial_genotyping_status_metrics <- do.call(
  rbind,
  lapply(
    status_prediction_groups,
    function(group_predictions) {
      prediction_errors <- (
        group_predictions$observed -
          group_predictions$predicted
      )
      
      data.frame(
        `repeat` =
          group_predictions[["repeat"]][1],
        fold =
          group_predictions$fold[1],
        model =
          group_predictions$model[1],
        genotyping_status =
          group_predictions$
          genotyping_status[1],
        n_testing =
          nrow(group_predictions),
        predictive_correlation = cor(
          group_predictions$observed,
          group_predictions$predicted
        ),
        rmse = sqrt(
          mean(prediction_errors^2)
        ),
        mae = mean(
          abs(prediction_errors)
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

partial_genotyping_status_metrics <-
  partial_genotyping_status_metrics[
    order(
      partial_genotyping_status_metrics[["repeat"]],
      partial_genotyping_status_metrics$fold,
      match(
        partial_genotyping_status_metrics$model,
        partial_model_order
      ),
      match(
        partial_genotyping_status_metrics$
          genotyping_status,
        c(
          "genotyped",
          "not_genotyped"
        )
      )
    ),
  ]

rownames(
  partial_genotyping_status_metrics
) <- NULL

# -------------------------------------------------------------------------
# Validate the paired evaluation
# -------------------------------------------------------------------------

overall_metric_keys <- paste(
  partial_genotyping_cv_metrics[["repeat"]],
  partial_genotyping_cv_metrics$fold,
  partial_genotyping_cv_metrics$model,
  sep = "::"
)

prediction_keys <- paste(
  partial_genotyping_cv_predictions$
    individual_id,
  partial_genotyping_cv_predictions[["repeat"]],
  partial_genotyping_cv_predictions$model,
  sep = "::"
)

status_metric_keys <- paste(
  partial_genotyping_status_metrics[["repeat"]],
  partial_genotyping_status_metrics$fold,
  partial_genotyping_status_metrics$model,
  partial_genotyping_status_metrics$
    genotyping_status,
  sep = "::"
)

observed_from_y <- unname(
  y[
    partial_genotyping_cv_predictions$
      individual_id
  ]
)

stopifnot(
  nrow(ssgblup_cv_metrics) ==
    length(cv_repeat_values) *
    length(cv_fold_values),
  nrow(ssgblup_cv_predictions) ==
    length(y) *
    length(cv_repeat_values),
  nrow(partial_genotyping_cv_metrics) ==
    length(cv_repeat_values) *
    length(cv_fold_values) *
    length(partial_model_order),
  nrow(partial_genotyping_cv_predictions) ==
    length(y) *
    length(cv_repeat_values) *
    length(partial_model_order),
  nrow(partial_genotyping_status_metrics) ==
    length(cv_repeat_values) *
    length(cv_fold_values) *
    length(partial_model_order) *
    length(
      unique(
        partial_genotyping_assignments$
          genotyping_status
      )
    ),
  !anyDuplicated(overall_metric_keys),
  !anyDuplicated(prediction_keys),
  !anyDuplicated(status_metric_keys),
  !anyNA(
    partial_genotyping_cv_predictions$
      genotyping_status
  ),
  setequal(
    partial_genotyping_cv_metrics$model,
    partial_model_order
  ),
  setequal(
    partial_genotyping_cv_predictions$model,
    partial_model_order
  ),
  setequal(
    partial_genotyping_status_metrics$
      genotyping_status,
    c(
      "genotyped",
      "not_genotyped"
    )
  ),
  all(
    is.finite(
      partial_genotyping_cv_metrics$
        predictive_correlation
    )
  ),
  all(
    is.finite(
      partial_genotyping_cv_metrics$rmse
    )
  ),
  all(
    is.finite(
      partial_genotyping_cv_metrics$mae
    )
  ),
  all(
    is.finite(
      partial_genotyping_status_metrics$
        predictive_correlation
    )
  ),
  all(
    is.finite(
      partial_genotyping_status_metrics$rmse
    )
  ),
  all(
    is.finite(
      partial_genotyping_status_metrics$mae
    )
  ),
  isTRUE(
    all.equal(
      partial_genotyping_cv_predictions$
        observed,
      observed_from_y,
      check.attributes = FALSE
    )
  )
)

# -------------------------------------------------------------------------
# Save ssGBLUP cross-validation results
# -------------------------------------------------------------------------

write.csv(
  partial_genotyping_cv_metrics,
  here(
    "results",
    "partial_genotyping_cv_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  partial_genotyping_cv_predictions,
  here(
    "results",
    "partial_genotyping_cv_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  partial_genotyping_status_metrics,
  here(
    "results",
    "partial_genotyping_status_metrics.csv"
  ),
  row.names = FALSE
)

cat("\n--- Partial-genotyping CV metrics ---\n")
print(partial_genotyping_cv_metrics)

cat(
  "\nSaved PBLUP and ssGBLUP results ",
  "for the partial-genotyping scenario.\n"
)

# -------------------------------------------------------------------------
# Summarise performance under partial genotyping
# -------------------------------------------------------------------------

partial_primary_metrics <- c(
  "predictive_correlation",
  "rmse",
  "mae"
)

partial_evaluation_sets <- list(
  overall = partial_genotyping_cv_metrics,
  genotyped =
    partial_genotyping_status_metrics[
      partial_genotyping_status_metrics$
        genotyping_status == "genotyped",
    ],
  not_genotyped =
    partial_genotyping_status_metrics[
      partial_genotyping_status_metrics$
        genotyping_status == "not_genotyped",
    ]
)

partial_summary_list <- list()
partial_summary_index <- 0L

for (
  evaluation_name in
  names(partial_evaluation_sets)
) {
  evaluation_results <-
    partial_evaluation_sets[[evaluation_name]]
  
  for (model_name in partial_model_order) {
    model_results <- evaluation_results[
      evaluation_results$model == model_name,
    ]
    
    for (metric_name in partial_primary_metrics) {
      partial_summary_index <-
        partial_summary_index + 1L
      
      metric_values <-
        model_results[[metric_name]]
      
      partial_summary_list[[
        partial_summary_index
      ]] <- data.frame(
        evaluation_group = evaluation_name,
        model = model_name,
        metric = metric_name,
        n_splits = length(metric_values),
        total_testing_predictions = sum(
          model_results$n_testing
        ),
        mean = mean(metric_values),
        standard_deviation = sd(metric_values),
        median = median(metric_values),
        minimum = min(metric_values),
        maximum = max(metric_values),
        stringsAsFactors = FALSE
      )
    }
  }
}

partial_genotyping_cv_summary <- do.call(
  rbind,
  partial_summary_list
)

rownames(partial_genotyping_cv_summary) <- NULL

# -------------------------------------------------------------------------
# Calculate paired ssGBLUP-minus-PBLUP differences
# -------------------------------------------------------------------------

partial_paired_list <- list()
partial_paired_index <- 0L

for (
  evaluation_name in
  names(partial_evaluation_sets)
) {
  evaluation_results <-
    partial_evaluation_sets[[evaluation_name]]
  
  pblup_evaluation <- evaluation_results[
    evaluation_results$model == "PBLUP",
    c(
      "repeat",
      "fold",
      partial_primary_metrics
    )
  ]
  
  ssgblup_evaluation <- evaluation_results[
    evaluation_results$model == "ssGBLUP",
    c(
      "repeat",
      "fold",
      partial_primary_metrics
    )
  ]
  
  paired_evaluation <- merge(
    pblup_evaluation,
    ssgblup_evaluation,
    by = c("repeat", "fold"),
    suffixes = c(
      "_PBLUP",
      "_ssGBLUP"
    ),
    sort = TRUE
  )
  
  partial_paired_index <-
    partial_paired_index + 1L
  
  partial_paired_list[[
    partial_paired_index
  ]] <- data.frame(
    evaluation_group = evaluation_name,
    `repeat` =
      paired_evaluation[["repeat"]],
    fold =
      paired_evaluation$fold,
    correlation_difference = (
      paired_evaluation$
        predictive_correlation_ssGBLUP -
        paired_evaluation$
        predictive_correlation_PBLUP
    ),
    rmse_difference = (
      paired_evaluation$rmse_ssGBLUP -
        paired_evaluation$rmse_PBLUP
    ),
    mae_difference = (
      paired_evaluation$mae_ssGBLUP -
        paired_evaluation$mae_PBLUP
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

partial_genotyping_paired_differences <-
  do.call(
    rbind,
    partial_paired_list
  )

rownames(
  partial_genotyping_paired_differences
) <- NULL

difference_columns <- c(
  predictive_correlation =
    "correlation_difference",
  rmse =
    "rmse_difference",
  mae =
    "mae_difference"
)

partial_paired_summary_list <- list()
partial_paired_summary_index <- 0L

for (
  evaluation_name in
  names(partial_evaluation_sets)
) {
  evaluation_differences <-
    partial_genotyping_paired_differences[
      partial_genotyping_paired_differences$
        evaluation_group == evaluation_name,
    ]
  
  for (metric_name in names(difference_columns)) {
    partial_paired_summary_index <-
      partial_paired_summary_index + 1L
    
    difference_values <-
      evaluation_differences[[
          difference_columns[[metric_name]]
        ]]
    
    ssgblup_better <- if (
      metric_name == "predictive_correlation"
    ) {
      difference_values > 0
    } else {
      difference_values < 0
    }
    
    partial_paired_summary_list[[
      partial_paired_summary_index
    ]] <- data.frame(
      evaluation_group = evaluation_name,
      metric = metric_name,
      difference_definition =
        "ssGBLUP minus PBLUP",
      n_splits = length(difference_values),
      mean_difference =
        mean(difference_values),
      standard_deviation =
        sd(difference_values),
      ssgblup_better_splits =
        sum(ssgblup_better),
      stringsAsFactors = FALSE
    )
  }
}

partial_genotyping_paired_summary <- do.call(
  rbind,
  partial_paired_summary_list
)

rownames(
  partial_genotyping_paired_summary
) <- NULL

# -------------------------------------------------------------------------
# Validate and save summaries
# -------------------------------------------------------------------------

partial_difference_keys <- paste(
  partial_genotyping_paired_differences$
    evaluation_group,
  partial_genotyping_paired_differences[[
    "repeat"
  ]],
  partial_genotyping_paired_differences$fold,
  sep = "::"
)

stopifnot(
  nrow(partial_genotyping_cv_summary) ==
    length(partial_evaluation_sets) *
    length(partial_model_order) *
    length(partial_primary_metrics),
  nrow(
    partial_genotyping_paired_differences
  ) ==
    length(partial_evaluation_sets) *
    number_ssgblup_fits,
  nrow(partial_genotyping_paired_summary) ==
    length(partial_evaluation_sets) *
    length(partial_primary_metrics),
  all(
    partial_genotyping_cv_summary$n_splits ==
      number_ssgblup_fits
  ),
  all(
    partial_genotyping_paired_summary$n_splits ==
      number_ssgblup_fits
  ),
  !anyDuplicated(partial_difference_keys),
  all(
    is.finite(
      partial_genotyping_cv_summary$mean
    )
  ),
  all(
    is.finite(
      partial_genotyping_paired_summary$
        mean_difference
    )
  )
)

write.csv(
  partial_genotyping_cv_summary,
  here(
    "results",
    "partial_genotyping_cv_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  partial_genotyping_paired_differences,
  here(
    "results",
    "partial_genotyping_paired_differences.csv"
  ),
  row.names = FALSE
)

write.csv(
  partial_genotyping_paired_summary,
  here(
    "results",
    "partial_genotyping_paired_summary.csv"
  ),
  row.names = FALSE
)

cat("\n--- Partial-genotyping summary ---\n")
print(partial_genotyping_cv_summary)

cat("\n--- Paired ssGBLUP comparison ---\n")
print(partial_genotyping_paired_summary)

# -------------------------------------------------------------------------
# Plot performance by test-line genotyping status
# -------------------------------------------------------------------------

partial_metric_labels <- c(
  predictive_correlation =
    "Predictive correlation\n(higher is better)",
  rmse =
    "RMSE\n(lower is better)",
  mae =
    "MAE\n(lower is better)"
)

partial_status_labels <- c(
  genotyped =
    "Genotyped testing individuals",
  not_genotyped =
    "Non-genotyped testing individuals"
)

partial_status_plot_data <- do.call(
  rbind,
  lapply(
    partial_primary_metrics,
    function(metric_name) {
      data.frame(
        `repeat` =
          partial_genotyping_status_metrics[[
            "repeat"
          ]],
        fold =
          partial_genotyping_status_metrics$fold,
        model =
          partial_genotyping_status_metrics$model,
        genotyping_status =
          partial_status_labels[
            partial_genotyping_status_metrics$
              genotyping_status
          ],
        metric =
          partial_metric_labels[[metric_name]],
        value =
          partial_genotyping_status_metrics[[
            metric_name
          ]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

partial_status_plot_data$split_id <- paste(
  partial_status_plot_data[["repeat"]],
  partial_status_plot_data$fold,
  sep = "-"
)

partial_status_plot_data$model <- factor(
  partial_status_plot_data$model,
  levels = partial_model_order
)

partial_status_plot_data$
  genotyping_status <- factor(
    partial_status_plot_data$
      genotyping_status,
    levels = unname(partial_status_labels)
  )

partial_status_plot_data$metric <- factor(
  partial_status_plot_data$metric,
  levels = unname(partial_metric_labels)
)

partial_genotyping_plot <- ggplot(
  partial_status_plot_data,
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
    size = 1.7,
    alpha = 0.80
  ) +
  facet_wrap(
    vars(
      genotyping_status,
      metric
    ),
    scales = "free_y",
    nrow = 2
  ) +
  scale_color_manual(
    values = c(
      PBLUP = "#2b5c8f",
      ssGBLUP = "#7b3294"
    )
  ) +
  scale_fill_manual(
    values = c(
      PBLUP = "#2b5c8f",
      ssGBLUP = "#7b3294"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title =
      "PBLUP and ssGBLUP Under 50% Genotyping",
    subtitle = paste0(
      "Performance evaluated separately for ",
      "genotyped and non-genotyped testing lines"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Points represent 25 repeat-fold combinations; ",
      "lines connect models evaluated on the same split. ",
      "Repeated folds overlap, so summaries are descriptive."
    )
  )

ggsave(
  filename = here(
    "figures",
    "partial_genotyping_model_comparison.png"
  ),
  plot = partial_genotyping_plot,
  width = 12,
  height = 8,
  dpi = 300
)

cat(
  "\nSaved partial-genotyping summaries ",
  "and comparison plot.\n"
)
