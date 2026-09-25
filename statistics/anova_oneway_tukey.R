## anova_oneway_tukey.R --------------------------------------------------------
## One-way ANOVA followed by Tukey's HSD post hoc test (all pairwise comparisons).
##
## Use this when every group is to be compared with every other group. When all
## groups are compared against a single control instead, use
## anova_oneway_dunnett.R.
##
## The example data are Figure S8K as published with the paper.
##
## Usage: edit SHEET and `groups`, then run.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------
## Installed automatically on first run.

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("installing openxlsx ...")
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

library(openxlsx)

OUTFILE <- "FigureS8K.xlsx"
SHEET   <- "FigureS8K"

TEST_NAME <- "One-way ANOVA with Tukey's HSD post hoc test"


## ---- data -------------------------------------------------------------------
## Group names become the column headers and the post hoc labels, so write them
## exactly as in the figure legend. Groups may differ in size.
##
## Group names may contain hyphens: comparison labels are built from the names
## themselves rather than by splitting what TukeyHSD returns.

groups <- list(
  "CTL-F1"    = c(36.6, 36.7, 35.4, 35.7, 35.5, 35.5, 35.7),
  "DEX-F1"    = c(35.1, 33.7, 34.6, 34.0, 34.2),
  "DEX-F1+Vc" = c(35.2, 34.8, 34.5, 34.3, 34.2, 35.4),
  "DEX+Vc-F1" = c(38.4, 38.5, 38.8, 38.0, 38.6)
)


## ---- helpers ----------------------------------------------------------------

star <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## Tukey t statistics are derived from the ANOVA mean square error; TukeyHSD
## reports differences and adjusted p values but not t.
tukey_t <- function(groups, fit, pairs) {
  mse <- summary(fit)[[1]][["Mean Sq"]][2]
  vapply(seq_len(nrow(pairs)), function(i) {
    a <- groups[[ pairs[i, 1] ]]
    b <- groups[[ pairs[i, 2] ]]
    round((mean(a) - mean(b)) /
            sqrt(mse * (1 / length(a) + 1 / length(b))), 3)
  }, numeric(1))
}

## TukeyHSD returns one row per pair in the order (2-1), (3-1), ..., (k-1).
## Rebuilding that order from the group names avoids splitting labels on "-".
pair_order <- function(nms) {
  k <- length(nms)
  do.call(rbind, unlist(lapply(seq_len(k - 1), function(i)
    lapply((i + 1):k, function(j) c(nms[j], nms[i]))), recursive = FALSE))
}


## ---- run --------------------------------------------------------------------

long <- data.frame(
  value = unlist(groups, use.names = FALSE),
  group = factor(rep(names(groups), lengths(groups)), levels = names(groups)),
  stringsAsFactors = FALSE
)

fit <- aov(value ~ group, data = long)
a   <- summary(fit)[[1]]

## raw values, one column per group
n_max <- max(lengths(groups))
raw   <- as.data.frame(
  lapply(groups, function(v) c(v, rep(NA, n_max - length(v)))),
  check.names = FALSE
)

## omnibus block, label column plus value column
omni <- data.frame(
  c("statistics", "Test", "Groups", "n per group", "p", "Signif"),
  c("",
    TEST_NAME,
    paste(names(groups), collapse = ", "),
    paste(lengths(groups), collapse = ", "),
    signif(a[["Pr(>F)"]][1], 4),
    star(a[["Pr(>F)"]][1])),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

## post hoc: all pairwise comparisons
tk    <- TukeyHSD(fit)$group
pairs <- pair_order(names(groups))

post <- data.frame(
  `Post hoc`  = paste(pairs[, 1], "\u2212", pairs[, 2]),
  Difference  = round(unname(tk[, "diff"]), 4),
  t           = tukey_t(groups, fit, pairs),
  `p adj`     = signif(unname(tk[, "p adj"]), 4),
  Signif      = vapply(tk[, "p adj"], star, character(1)),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
rownames(post) <- NULL

r_stats <- n_max + 2                 # raw header + values, then statistics
r_post  <- r_stats + nrow(omni) + 1

wb    <- createWorkbook()
arial <- createStyle(fontName = "Arial", fontSize = 10)
addWorksheet(wb, SHEET)

writeData(wb, SHEET, raw,  startRow = 1,       startCol = 1)
writeData(wb, SHEET, omni, startRow = r_stats, startCol = 1, colNames = FALSE)
writeData(wb, SHEET, "Post hoc: Tukey's HSD test (all pairwise comparisons)",
          startRow = r_post, startCol = 1)
writeData(wb, SHEET, post, startRow = r_post + 1, startCol = 1)

addStyle(wb, SHEET, arial, rows = 1:40, cols = 1:10, gridExpand = TRUE)
saveWorkbook(wb, OUTFILE, overwrite = TRUE)

print(post)
cat(sprintf("ANOVA p = %s  %s\n", omni[[2]][5], omni[[2]][6]))
cat("written to", OUTFILE, "sheet", SHEET, "\n")
