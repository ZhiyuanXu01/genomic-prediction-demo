library(here)

# -------------------------------------------------------------------------
# Load completed analysis summaries
# -------------------------------------------------------------------------

required_result_files <- c(
  "repeated_cv_metrics.csv",
  "repeated_cv_summary.csv",
  "repeated_cv_paired_summary.csv",
  "partial_genotyping_cv_summary.csv",
  "partial_genotyping_paired_summary.csv"
)

required_result_paths <- file.path(
  here("results"),
  required_result_files
)

missing_result_files <- required_result_files[
  !file.exists(required_result_paths)
]

if (length(missing_result_files) > 0) {
  stop(
    paste0(
      "Missing result files: ",
      paste(
        missing_result_files,
        collapse = ", "
      ),
      ". Run the preceding analysis scripts first."
    )
  )
}

read_result <- function(file_name) {
  read.csv(
    here(
      "results",
      file_name
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

repeated_cv_metrics_final <- read_result(
  "repeated_cv_metrics.csv"
)

repeated_cv_summary_final <- read_result(
  "repeated_cv_summary.csv"
)

repeated_cv_paired_summary_final <- read_result(
  "repeated_cv_paired_summary.csv"
)

partial_cv_summary_final <- read_result(
  "partial_genotyping_cv_summary.csv"
)

partial_paired_summary_final <- read_result(
  "partial_genotyping_paired_summary.csv"
)

# -------------------------------------------------------------------------
# Combine model-performance summaries
# -------------------------------------------------------------------------

full_prediction_records <- aggregate(
  n_testing ~ model,
  data = repeated_cv_metrics_final,
  FUN = sum
)

full_summary_long <- repeated_cv_summary_final

full_summary_long$scenario <-
  "Full genotyping"

full_summary_long$evaluation_group <-
  "overall"

full_summary_long$prediction_records <-
  full_prediction_records$n_testing[
    match(
      full_summary_long$model,
      full_prediction_records$model
    )
  ]

full_summary_long <- full_summary_long[
  ,
  c(
    "scenario",
    "evaluation_group",
    "model",
    "metric",
    "n_splits",
    "prediction_records",
    "mean",
    "standard_deviation",
    "median",
    "minimum",
    "maximum"
  )
]

partial_summary_long <-
  partial_cv_summary_final

partial_summary_long$scenario <-
  "50% genotyping"

partial_summary_long$prediction_records <-
  partial_summary_long$
  total_testing_predictions

partial_summary_long <- partial_summary_long[
  ,
  c(
    "scenario",
    "evaluation_group",
    "model",
    "metric",
    "n_splits",
    "prediction_records",
    "mean",
    "standard_deviation",
    "median",
    "minimum",
    "maximum"
  )
]

combined_summary_long <- rbind(
  full_summary_long,
  partial_summary_long
)

primary_metrics <- c(
  "predictive_correlation",
  "rmse",
  "mae"
)

desired_performance_rows <- data.frame(
  scenario = c(
    "Full genotyping",
    "Full genotyping",
    rep(
      "50% genotyping",
      6
    )
  ),
  evaluation_group = c(
    "overall",
    "overall",
    "overall",
    "overall",
    "genotyped",
    "genotyped",
    "not_genotyped",
    "not_genotyped"
  ),
  model = c(
    "PBLUP",
    "GBLUP",
    "PBLUP",
    "ssGBLUP",
    "PBLUP",
    "ssGBLUP",
    "PBLUP",
    "ssGBLUP"
  ),
  stringsAsFactors = FALSE
)

make_performance_row <- function(row_number) {
  row_definition <-
    desired_performance_rows[
      row_number,
    ]
  
  selected_summary <- combined_summary_long[
    combined_summary_long$scenario ==
      row_definition$scenario &
      combined_summary_long$evaluation_group ==
      row_definition$evaluation_group &
      combined_summary_long$model ==
      row_definition$model,
  ]
  
  metric_order <- match(
    primary_metrics,
    selected_summary$metric
  )
  
  stopifnot(
    nrow(selected_summary) ==
      length(primary_metrics),
    !anyNA(metric_order),
    length(
      unique(
        selected_summary$n_splits
      )
    ) == 1,
    length(
      unique(
        selected_summary$
          prediction_records
      )
    ) == 1
  )
  
  selected_summary <-
    selected_summary[
      metric_order,
    ]
  
  data.frame(
    scenario =
      row_definition$scenario,
    evaluation_group =
      row_definition$evaluation_group,
    model =
      row_definition$model,
    n_splits =
      unique(
        selected_summary$n_splits
      ),
    prediction_records =
      unique(
        selected_summary$
          prediction_records
      ),
    mean_predictive_correlation =
      selected_summary$mean[1],
    sd_predictive_correlation =
      selected_summary$
      standard_deviation[1],
    mean_rmse =
      selected_summary$mean[2],
    sd_rmse =
      selected_summary$
      standard_deviation[2],
    mean_mae =
      selected_summary$mean[3],
    sd_mae =
      selected_summary$
      standard_deviation[3],
    stringsAsFactors = FALSE
  )
}

final_model_performance <- do.call(
  rbind,
  lapply(
    seq_len(
      nrow(
        desired_performance_rows
      )
    ),
    make_performance_row
  )
)

rownames(final_model_performance) <- NULL

# -------------------------------------------------------------------------
# Combine paired model-improvement summaries
# -------------------------------------------------------------------------

full_paired_long <- data.frame(
  scenario = rep(
    "Full genotyping",
    nrow(
      repeated_cv_paired_summary_final
    )
  ),
  evaluation_group = "overall",
  comparison =
    repeated_cv_paired_summary_final$
    difference_definition,
  metric =
    repeated_cv_paired_summary_final$
    metric,
  n_splits =
    repeated_cv_paired_summary_final$
    total_splits,
  mean_difference =
    repeated_cv_paired_summary_final$
    mean_difference,
  standard_deviation =
    repeated_cv_paired_summary_final$
    standard_deviation,
  better_splits =
    repeated_cv_paired_summary_final$
    gblup_better_splits,
  stringsAsFactors = FALSE
)

partial_paired_long <- data.frame(
  scenario = rep(
    "50% genotyping",
    nrow(
      partial_paired_summary_final
    )
  ),
  evaluation_group =
    partial_paired_summary_final$
    evaluation_group,
  comparison =
    partial_paired_summary_final$
    difference_definition,
  metric =
    partial_paired_summary_final$
    metric,
  n_splits =
    partial_paired_summary_final$
    n_splits,
  mean_difference =
    partial_paired_summary_final$
    mean_difference,
  standard_deviation =
    partial_paired_summary_final$
    standard_deviation,
  better_splits =
    partial_paired_summary_final$
    ssgblup_better_splits,
  stringsAsFactors = FALSE
)

combined_paired_long <- rbind(
  full_paired_long,
  partial_paired_long
)

desired_paired_rows <- data.frame(
  scenario = c(
    "Full genotyping",
    rep(
      "50% genotyping",
      3
    )
  ),
  evaluation_group = c(
    "overall",
    "overall",
    "genotyped",
    "not_genotyped"
  ),
  comparison = c(
    "GBLUP minus PBLUP",
    rep(
      "ssGBLUP minus PBLUP",
      3
    )
  ),
  stringsAsFactors = FALSE
)

make_paired_row <- function(row_number) {
  row_definition <-
    desired_paired_rows[
      row_number,
    ]
  
  selected_summary <- combined_paired_long[
    combined_paired_long$scenario ==
      row_definition$scenario &
      combined_paired_long$evaluation_group ==
      row_definition$evaluation_group &
      combined_paired_long$comparison ==
      row_definition$comparison,
  ]
  
  metric_order <- match(
    primary_metrics,
    selected_summary$metric
  )
  
  stopifnot(
    nrow(selected_summary) ==
      length(primary_metrics),
    !anyNA(metric_order),
    length(
      unique(
        selected_summary$n_splits
      )
    ) == 1
  )
  
  selected_summary <-
    selected_summary[
      metric_order,
    ]
  
  data.frame(
    scenario =
      row_definition$scenario,
    evaluation_group =
      row_definition$evaluation_group,
    comparison =
      row_definition$comparison,
    n_splits =
      unique(
        selected_summary$n_splits
      ),
    mean_correlation_difference =
      selected_summary$
      mean_difference[1],
    sd_correlation_difference =
      selected_summary$
      standard_deviation[1],
    correlation_better_splits =
      selected_summary$
      better_splits[1],
    mean_rmse_difference =
      selected_summary$
      mean_difference[2],
    sd_rmse_difference =
      selected_summary$
      standard_deviation[2],
    rmse_better_splits =
      selected_summary$
      better_splits[2],
    mean_mae_difference =
      selected_summary$
      mean_difference[3],
    sd_mae_difference =
      selected_summary$
      standard_deviation[3],
    mae_better_splits =
      selected_summary$
      better_splits[3],
    stringsAsFactors = FALSE
  )
}

final_paired_improvements <- do.call(
  rbind,
  lapply(
    seq_len(
      nrow(
        desired_paired_rows
      )
    ),
    make_paired_row
  )
)

rownames(final_paired_improvements) <- NULL

# -------------------------------------------------------------------------
# Validate and save final tables
# -------------------------------------------------------------------------

stopifnot(
  nrow(final_model_performance) == 8,
  nrow(final_paired_improvements) == 4,
  !anyDuplicated(
    final_model_performance[
      ,
      c(
        "scenario",
        "evaluation_group",
        "model"
      )
    ]
  ),
  !anyDuplicated(
    final_paired_improvements[
      ,
      c(
        "scenario",
        "evaluation_group",
        "comparison"
      )
    ]
  ),
  all(
    final_model_performance$n_splits ==
      25
  ),
  all(
    final_paired_improvements$n_splits ==
      25
  ),
  all(
    is.finite(
      as.matrix(
        final_model_performance[
          ,
          6:ncol(
            final_model_performance
          )
        ]
      )
    )
  ),
  all(
    is.finite(
      as.matrix(
        final_paired_improvements[
          ,
          5:ncol(
            final_paired_improvements
          )
        ]
      )
    )
  ),
  all(
    final_paired_improvements$
      correlation_better_splits <=
      final_paired_improvements$n_splits
  ),
  all(
    final_paired_improvements$
      rmse_better_splits <=
      final_paired_improvements$n_splits
  ),
  all(
    final_paired_improvements$
      mae_better_splits <=
      final_paired_improvements$n_splits
  )
)

write.csv(
  final_model_performance,
  here(
    "results",
    "final_model_performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_paired_improvements,
  here(
    "results",
    "final_paired_improvements.csv"
  ),
  row.names = FALSE
)

cat("\n--- Final model performance ---\n")
print(final_model_performance)

cat("\n--- Final paired improvements ---\n")
print(final_paired_improvements)

cat(
  "\nSaved final project-level summary tables.\n"
)

