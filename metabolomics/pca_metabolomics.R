## pca_metabolomics.R ----------------------------------------------------------
## Principal component analysis of targeted metabolomics data.
##
## Input : a CSV with metabolites in rows and samples in columns, the first
##         column holding the compound names.
## Output: an .xlsx with PC1/PC2 scores per sample, the loadings, and the
##         variance explained by each component.
##
## Values are log2-transformed, then centred and scaled to unit variance before
## decomposition, so that abundant and scarce metabolites contribute equally.
##
## Usage: run the script, pick the matrix in the dialog, set the group labels.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("installing openxlsx ...")
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

library(openxlsx)

## The input matrix is chosen from a file dialog when the script runs.
message("select the metabolomics matrix (CSV)")
INFILE  <- file.choose()

OUTFILE <- "PCA_iBAT_meta.xlsx"

## Group label for each sample, matched by a pattern in the sample name.
GROUPS <- c("Vc-F1" = "_VcF1_", "dVc-F1" = "_dVcF1_")

N_PC <- 2                      # components to report


## ---- read -------------------------------------------------------------------

mat <- read.csv(INFILE, row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

stopifnot(all(is.finite(mat)), all(mat >= 0))

## samples in rows, metabolites in columns
X <- t(log2(mat + 1))


## ---- pca --------------------------------------------------------------------

pca <- prcomp(X, center = TRUE, scale. = TRUE)

var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

## assign groups by pattern; samples matching nothing are left blank
group <- rep("", nrow(X))
for (g in names(GROUPS)) group[grepl(GROUPS[[g]], rownames(X), fixed = TRUE)] <- g

scores <- data.frame(
  Sample = rownames(X),
  round(pca$x[, seq_len(N_PC), drop = FALSE], 6),
  Group  = group,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
rownames(scores) <- NULL

loadings <- data.frame(
  Compound = rownames(pca$rotation),
  round(pca$rotation[, seq_len(N_PC), drop = FALSE], 6),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
rownames(loadings) <- NULL

variance <- data.frame(
  PC                     = paste0("PC", seq_along(var_pct)),
  `Variance explained %` = round(var_pct, 2),
  `Cumulative %`         = round(cumsum(var_pct), 2),
  check.names            = FALSE
)


## ---- write ------------------------------------------------------------------

wb    <- createWorkbook()
arial <- createStyle(fontName = "Arial", fontSize = 10)

for (nm in c("scores", "loadings", "variance")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, get(nm), startRow = 1, startCol = 1)
  addStyle(wb, nm, arial,
           rows = 1:(nrow(get(nm)) + 1), cols = 1:6, gridExpand = TRUE)
}

saveWorkbook(wb, OUTFILE, overwrite = TRUE)

cat(sprintf("%d samples, %d metabolites\n", nrow(X), ncol(X)))
cat(sprintf("PC1 %.1f%%, PC2 %.1f%%\n", var_pct[1], var_pct[2]))
cat("written to", OUTFILE, "\n")
