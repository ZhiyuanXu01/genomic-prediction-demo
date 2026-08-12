library(BGLR)
library(here)
library(ggplot2)

# -------------------------------------------------------------------------
# Load the public CIMMYT wheat dataset
# -------------------------------------------------------------------------

data(wheat)

Y <- wheat.Y
M <- wheat.X
A_pedigree <- wheat.A

# Preserve the CIMMYT line identifiers.
stopifnot(
  nrow(Y) == nrow(M),
  nrow(Y) == nrow(A_pedigree),
  nrow(A_pedigree) ==
    ncol(A_pedigree),
  identical(
    rownames(Y),
    rownames(A_pedigree)
  ),
  identical(
    rownames(A_pedigree),
    colnames(A_pedigree)
  ),
  !anyDuplicated(
    rownames(Y)
  )
)

line_ids <- rownames(Y)

rownames(M) <- line_ids

# Environment 1 standardized grain yield.
y <- as.numeric(
  Y[, 1]
)

names(y) <- line_ids

stopifnot(
  identical(
    names(y),
    rownames(M)
  ),
  identical(
    rownames(M),
    rownames(A_pedigree)
  ),
  !anyNA(y),
  !anyNA(M),
  all(is.finite(y))
)

# -------------------------------------------------------------------------
# Inspect dimensions, marker coding, and missingness
# -------------------------------------------------------------------------

marker_values <- sort(
  unique(
    as.vector(M)
  )
)

stopifnot(
  identical(
    marker_values,
    c(0, 1)
  )
)

cat("\n--- Dataset dimensions ---\n")
cat(
  "Phenotype matrix:",
  paste(
    dim(Y),
    collapse = " x "
  ),
  "\n"
)
cat(
  "Marker matrix:",
  paste(
    dim(M),
    collapse = " x "
  ),
  "\n"
)
cat(
  "Pedigree relationship matrix:",
  paste(
    dim(A_pedigree),
    collapse = " x "
  ),
  "\n"
)

cat("\n--- Data quality diagnostics ---\n")
cat(
  "Marker values:",
  paste(
    marker_values,
    collapse = ", "
  ),
  "\n"
)
cat(
  "Missing Environment 1 phenotypes:",
  sum(is.na(y)),
  "\n"
)
cat(
  "Missing marker scores:",
  sum(is.na(M)),
  "\n"
)
cat(
  "Identifiers aligned:",
  identical(
    names(y),
    rownames(M)
  ) &&
    identical(
      rownames(M),
      rownames(A_pedigree)
    ),
  "\n"
)

cat("\n--- Environment 1 phenotype summary ---\n")
print(
  summary(y)
)

# -------------------------------------------------------------------------
# Save the dataset summary
# -------------------------------------------------------------------------

data_summary <- data.frame(
  Metric = c(
    "Total Lines",
    "Total Markers",
    "Phenotype Environments",
    "Marker Unique Values",
    "Missing Phenotypes (E1)",
    "Missing Markers",
    "E1 Phenotype Mean",
    "E1 Phenotype SD"
  ),
  Value = c(
    nrow(Y),
    ncol(M),
    ncol(Y),
    paste(
      marker_values,
      collapse = "/"
    ),
    sum(is.na(y)),
    sum(is.na(M)),
    round(
      mean(y),
      4
    ),
    round(
      sd(y),
      4
    )
  ),
  stringsAsFactors = FALSE
)

dir.create(
  here("results"),
  showWarnings = FALSE
)

write.csv(
  data_summary,
  here(
    "results",
    "data_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSaved summary to: results/data_summary.csv\n"
)

# -------------------------------------------------------------------------
# Visualise the phenotype distribution
# -------------------------------------------------------------------------

dir.create(
  here("figures"),
  showWarnings = FALSE
)

histogram_bins <- 30L

phenotype_bin_width <- (
  max(y) - min(y)
) / histogram_bins

phenotype_plot_data <- data.frame(
  phenotype = y
)

phenotype_distribution_plot <- ggplot(
  phenotype_plot_data,
  aes(
    x = phenotype
  )
) +
  geom_histogram(
    bins = histogram_bins,
    fill = "#2b5c8f",
    color = "white",
    alpha = 0.80
  ) +
  geom_density(
    aes(
      y = after_stat(count) *
        phenotype_bin_width
    ),
    color = "#d95f02",
    linewidth = 1
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title =
      "Distribution of Grain Yield (Environment 1)",
    subtitle = paste0(
      "CIMMYT wheat dataset (n = ",
      length(y),
      ")"
    ),
    x = "Standardized grain yield",
    y = "Frequency"
  )

ggsave(
  filename = here(
    "figures",
    "phenotype_distribution.png"
  ),
  plot = phenotype_distribution_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat(
  "Saved plot to: figures/phenotype_distribution.png\n"
)