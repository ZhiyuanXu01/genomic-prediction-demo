library(BGLR)
library(here)
library(ggplot2)

data(wheat)

Y <- wheat.Y
M <- wheat.X
A_pedigree <- wheat.A

# Preserve CIMMYT line identifiers.
stopifnot(
  nrow(Y) == nrow(M),
  nrow(Y) == nrow(A_pedigree),
  identical(rownames(Y), rownames(A_pedigree)),
  identical(rownames(A_pedigree), colnames(A_pedigree))
)

rownames(M) <- rownames(Y)

# Binary marker frequency and minor allele frequency.
marker_frequency <- colMeans(M)
marker_maf <- pmin(marker_frequency, 1 - marker_frequency)

cat("--- Marker frequency diagnostics ---\n")
cat("Number of markers:", ncol(M), "\n")
cat("Monomorphic markers:", sum(marker_maf == 0), "\n")
cat("Markers with MAF < 0.01:", sum(marker_maf < 0.01), "\n")
cat("Markers with MAF < 0.05:", sum(marker_maf < 0.05), "\n")

cat("\nMAF summary:\n")
print(summary(marker_maf))

cat("\nPedigree A diagonal summary:\n")
print(summary(diag(A_pedigree)))

cat("\nPedigree A off-diagonal summary:\n")
print(summary(A_pedigree[lower.tri(A_pedigree)]))

# -------------------------------------------------------------------------
# Construct the genomic relationship matrix
# -------------------------------------------------------------------------

maf_threshold <- 0.05
keep_maf <- marker_maf >= maf_threshold
M_filtered <- M[, keep_maf, drop = FALSE]

cat("\n--- Marker filtering ---\n")
cat("MAF threshold:", maf_threshold, "\n")
cat("Markers retained:", ncol(M_filtered), "\n")
cat("Markers removed:", ncol(M) - ncol(M_filtered), "\n")

stopifnot(
  ncol(M_filtered) == ncol(M) - sum(marker_maf < maf_threshold),
  !anyNA(M_filtered),
  all(M_filtered %in% c(0, 1))
  )

# The historical wheat lines are highly inbred. The binary DArT states are
# therefore represented as the two homozygous marker classes required by
# rrBLUP::A.mat(): 0 -> -1 and 1 -> +1.

M_rrblup <- 2 * M_filtered - 1

stopifnot(all(M_rrblup %in% c(-1,1)))

G_genomic <- rrBLUP::A.mat(
  X = M_rrblup,
  min.MAF = NULL,
  max.missing = 0,
  impute.method = "mean",
  shrink = FALSE
)

# -------------------------------------------------------------------------
# Validate dimensions, identifiers, symmetry, and eigenvalues
# -------------------------------------------------------------------------

stopifnot(
  nrow(G_genomic) == nrow(A_pedigree),
  ncol(G_genomic) == ncol(A_pedigree),
  identical(rownames(G_genomic), rownames(A_pedigree)),
  identical(colnames(G_genomic), colnames(A_pedigree))
)

maximum_asymmetry <- max(abs(G_genomic - t(G_genomic)))

eigenvalues_G <- eigen(
  G_genomic,
  symmetric = TRUE,
  only.values = TRUE
)$values

eigen_tolerance <- 1e-8

negative_eigenvalues <- sum(eigenvalues_G < -eigen_tolerance)
near_zero_eigenvalues <- sum(abs(eigenvalues_G) <= eigen_tolerance)

cat("\n--- Genomic G diagnostics ---\n")
cat("Dimensions:", paste(dim(G_genomic), collapse = " x "), "\n")
cat("Maximum asymmetry:", maximum_asymmetry, "\n")
cat("Minimum eigenvalue:", min(eigenvalues_G), "\n")
cat("Negative eigenvalues:", negative_eigenvalues, "\n")
cat("Near-zero eigenvalues:", near_zero_eigenvalues, "\n")

cat("\nGenomic G diagonal summary:\n")
print(summary(diag(G_genomic)))

cat("\nGenomic G off-diagonal summary:\n")
print(summary(G_genomic[lower.tri(G_genomic)]))

# -------------------------------------------------------------------------
# Compare pedigree and genomic relationships
# -------------------------------------------------------------------------

lower_A <- A_pedigree[lower.tri(A_pedigree)]
lower_G <- G_genomic[lower.tri(G_genomic)]

pearson_A_G <- cor(lower_A, lower_G, method = "pearson")
spearman_A_G <- cor(lower_A, lower_G, method = "spearman")

cat("\n--- Comparison of A and G ---\n")
cat("Pearson correlation:", pearson_A_G, "\n")
cat("Spearman correlation:", spearman_A_G, "\n")

relationship_summary <- data.frame(
  metric = c(
    "markers_total",
    "markers_retained_maf_0.05",
    "markers_removed_maf_0.05",
    "A_mean_diagonal",
    "A_mean_off_diagonal",
    "G_mean_diagonal",
    "G_mean_off_diagonal",
    "G_maximum_asymmetry",
    "G_minimum_eigenvalue",
    "G_negative_eigenvalues",
    "G_near_zero_eigenvalues",
    "A_G_pearson_off_diagonal",
    "A_G_spearman_off_diagonal"
  ),
  value = c(
    ncol(M),
    ncol(M_filtered),
    ncol(M) - ncol(M_filtered),
    mean(diag(A_pedigree)),
    mean(lower_A),
    mean(diag(G_genomic)),
    mean(lower_G),
    maximum_asymmetry,
    min(eigenvalues_G),
    negative_eigenvalues,
    near_zero_eigenvalues,
    pearson_A_G,
    spearman_A_G
  )
)

# -------------------------------------------------------------------------
# Align G to the pedigree relationship scale
# -------------------------------------------------------------------------

align_G_to_A <- function(G, A) {
  stopifnot(
    identical(dim(G), dim(A)),
    identical(rownames(G), rownames(A)),
    identical(colnames(G), colnames(A))
  )
  
  mean_diag_G <- mean(diag(G))
  mean_off_G <- mean(G[lower.tri(G)])
  
  mean_diag_A <- mean(diag(A))
  mean_off_A <- mean(A[lower.tri(A)])
  
  beta <- (
    mean_diag_A - mean_off_A
  ) / (
    mean_diag_G - mean_off_G
  )
  
  alpha <- mean_off_A - beta * mean_off_G
  
  G_aligned <- alpha + beta * G
  
  list(
    G = G_aligned,
    alpha = alpha,
    beta = beta
  )
}

alignment_full <- align_G_to_A(
  G = G_genomic,
  A = A_pedigree
)

G_aligned_full <- alignment_full$G
alignment_alpha <- alignment_full$alpha
alignment_beta <- alignment_full$beta

eigenvalues_G_aligned <- eigen(
  G_aligned_full,
  symmetric = TRUE,
  only.values = TRUE
)$values

aligned_maximum_asymmetry <- max(
  abs(G_aligned_full - t(G_aligned_full))
)

cat("\n--- Alignment of G to A ---\n")
cat("Alpha:", alignment_alpha, "\n")
cat("Beta:", alignment_beta, "\n")
cat(
  "Aligned G mean diagonal:",
  mean(diag(G_aligned_full)),
  "\n"
)
cat(
  "Aligned G mean off-diagonal:",
  mean(G_aligned_full[lower.tri(G_aligned_full)]),
  "\n"
)
cat(
  "Aligned G minimum eigenvalue:",
  min(eigenvalues_G_aligned),
  "\n"
)
cat(
  "Aligned G maximum asymmetry:",
  aligned_maximum_asymmetry,
  "\n"
)

stopifnot(
  aligned_maximum_asymmetry < 1e-10,
  min(eigenvalues_G_aligned) > -1e-8,
  abs(
    mean(diag(G_aligned_full)) -
      mean(diag(A_pedigree))
  ) < 1e-10,
  abs(
    mean(G_aligned_full[lower.tri(G_aligned_full)]) -
      mean(A_pedigree[lower.tri(A_pedigree)])
  ) < 1e-10
)

relationship_summary <- rbind(
  relationship_summary,
  data.frame(
    metric = c(
      "G_alignment_alpha",
      "G_alignment_beta",
      "G_aligned_mean_diagonal",
      "G_aligned_mean_off_diagonal",
      "G_aligned_minimum_eigenvalue",
      "G_aligned_maximum_asymmetry"
    ),
    value = c(
      alignment_alpha,
      alignment_beta,
      mean(diag(G_aligned_full)),
      mean(G_aligned_full[lower.tri(G_aligned_full)]),
      min(eigenvalues_G_aligned),
      aligned_maximum_asymmetry
    )
  )
)

write.csv(
  relationship_summary,
  here("results", "relationship_matrix_summary.csv"),
  row.names = FALSE
)

cat("\nSaved summary to: results/relationship_matrix_summary.csv\n")

# -------------------------------------------------------------------------
# Visualise pedigree and genomic pairwise relationships
# -------------------------------------------------------------------------

lower_G_aligned <- G_aligned_full[
  lower.tri(G_aligned_full)
]

set.seed(20260810)

number_of_pairs_to_plot <- min(
  50000L,
  length(lower_A)
)

plot_indices <- sample.int(
  length(lower_A),
  size = number_of_pairs_to_plot,
  replace = FALSE
)

relationship_plot_data <- rbind(
  data.frame(
    pedigree_relationship = lower_A[plot_indices],
    genomic_relationship = lower_G[plot_indices],
    genomic_matrix = "Raw G"
  ),
  data.frame(
    pedigree_relationship = lower_A[plot_indices],
    genomic_relationship = lower_G_aligned[plot_indices],
    genomic_matrix = "G aligned to A"
  )
)

relationship_plot <- ggplot(
  relationship_plot_data,
  aes(
    x = pedigree_relationship,
    y = genomic_relationship
  )
) +
  geom_point(
    color = "#2b5c8f",
    alpha = 0.10,
    size = 0.45
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.6
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "#d95f02",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ genomic_matrix,
    nrow = 1
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Pedigree and Genomic Pairwise Relationships",
    subtitle = paste0(
      "CIMMYT wheat lines; ",
      format(number_of_pairs_to_plot, big.mark = ","),
      " sampled pairs"
    ),
    x = "Pedigree relationship (A)",
    y = "Genomic relationship",
    caption = paste0(
      "Correlations were calculated using all ",
      format(length(lower_A), big.mark = ","),
      " pairwise relationships."
    )
  )

ggsave(
  filename = here(
    "figures",
    "A_vs_G_relationship.png"
  ),
  plot = relationship_plot,
  width = 10,
  height = 5,
  dpi = 300
)

cat(
  "Saved plot to: figures/A_vs_G_relationship.png\n"
)
