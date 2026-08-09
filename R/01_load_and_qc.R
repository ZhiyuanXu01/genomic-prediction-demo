library(BGLR)
library(here)
library(ggplot2)

# 1
data(wheat)

Y <- wheat.Y
M <- wheat.X
A <- wheat.A

# 2

cat("Phenotype dimensions:", dim(Y), "\n")
cat("Marker dimensions:", dim(M), "\n")
cat("Pedigree relationship dimensions:", dim(A), "\n")

cat("Marker values:\n")
print(sort(unique(as.vector(M))))

cat("Missing phenotypes:", sum(is.na(Y)), "\n")
cat("Missing marker scores:", sum(is.na(M)), "\n")

# 3
y <- Y[,1]

print(summary(y))

# 4
line_ids <- sprintf("Line_%03d", 1:599)

names(y) <- line_ids
rownames(M) <- line_ids
rownames(A) <- line_ids
colnames(A) <- line_ids

check_ids <- identical(names(y), rownames(M)) && identical(rownames(M), rownames(A))

cat("\n--- ID Order Alignment ---\n")
cat("Are ID orders in y, M, and A completely identical?:", check_ids, "\n")

stopifnot(check_ids)

# 5
summary_df <- data.frame(
  Metric=c(
    "Total Lines",
    "Total Markers",
    "Phenotype Environments",
    "Marker Unique Values",
    "Missing Phenotypes (E1)",
    "Missing Markers",
    "E1 Phenotype Mean",
    "E1 Phenotype SD"
  ),
  Value=c(
    nrow(Y),
    ncol(M),
    ncol(Y),
    paste(sort(unique(as.vector(M))), collapse = '/'),
    sum(is.na(y)),
    sum(is.na(M)),
    round(mean(y, na.rm = TRUE), 4),
    round(sd(y, na.rm = TRUE), 4)
  )
)

dir.create(here("results"), showWarnings = FALSE)
write.csv(summary_df, here("results", "data_summary.csv"), row.names = FALSE)
cat("\nSaved summary to: results/data_summary.csv\n")

# 6
dir.create(here("figures"), showWarnings = FALSE)

p <- ggplot(data.frame(y=y), aes(x=y))+
  geom_histogram(bins=30, fill = "#2b5c8f", color = "white", alpha = 0.8) +
  geom_density(aes(y=after_stat(count)*(max(y, na.rm = TRUE)-min(y, na.rm = TRUE))/30), color = "#d95f02", linewidth = 1) +
  theme_minimal(base_size=12) +
  labs(
    title = "Distribution of Grain Yield (Environment 1)",
    subtitle = paste0("CIMMYT Wheat Dataset (n = ", nrow(y), ")"),
    x = "Standardized Grain Yield",
    y = "Frequency"
  )

ggsave(
  filename=(here("figures", "phenotype_distribution.png")),
  plot=p,
  width=7,
  height=5,
  dpi=300
)

cat("Saved plot to: figures/phenotype_distribution.png\n")
