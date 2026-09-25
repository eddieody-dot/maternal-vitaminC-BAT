## ttest.R ---------------------------------------------------------------------
## Two-tailed unpaired t-test, with Welch correction applied when an F test
## indicates unequal variances (F test p < 0.05).
##
## The example data below are the first two time points of Figure 1B, as
## published in the source data accompanying the paper.
##
## Usage: replace the contents of `panels` with your own values and run.
## -----------------------------------------------------------------------------

library(openxlsx)

OUTFILE <- "ttest_results.xlsx"

TEST_NAME <- "Two-tailed unpaired t-test (Welch correction where F test p < 0.05)"


## ---- data -------------------------------------------------------------------
## One entry per panel: x1 is the control group, x2 the comparison group.
## Figure1B_E0.5  -> equal variances, Student
## Figure1B_E11.5 -> unequal variances, Welch correction applied

panels <- list(

  "Figure1B_E0.5" = list(
    x1 = c(79.47, 91.00, 80.73),
    x2 = c(82.93, 89.97, 92.63)
  ),

  "Figure1B_E11.5" = list(
    x1 = c(77.54, 89.10, 71.50),
    x2 = c( 3.10,  1.01,  2.56)
  )

)


## ---- helpers ----------------------------------------------------------------

fmt_p <- function(p) if (p < 1e-4) "<0.0001" else formatC(p, digits = 4, format = "g")

fmt_n <- function(x, digits = 3) formatC(round(x, digits), format = "g")

star <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## Returns the statistics block: two columns, no header row.
ttest_block <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]

  f_p       <- var.test(x1, x2)$p.value
  use_equal <- f_p >= 0.05
  res       <- t.test(x1, x2, var.equal = use_equal)

  data.frame(
    c("statistics", "Test", "n", "F test p", "t", "df", "p", "Signif"),
    c("",
      TEST_NAME,
      paste0(length(x1), ", ", length(x2)),
      fmt_p(f_p),
      fmt_n(unname(res$statistic)),
      fmt_n(unname(res$parameter)),
      fmt_p(res$p.value),
      star(res$p.value)),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}


## ---- run --------------------------------------------------------------------

wb    <- createWorkbook()
arial <- createStyle(fontName = "Arial", fontSize = 10)

for (panel in names(panels)) {

  x1 <- panels[[panel]]$x1
  x2 <- panels[[panel]]$x2

  n_max <- max(length(x1), length(x2))
  raw   <- data.frame(
    x1 = c(x1, rep(NA, n_max - length(x1))),
    x2 = c(x2, rep(NA, n_max - length(x2))),
    check.names = FALSE
  )

  stats <- ttest_block(x1, x2)

  addWorksheet(wb, panel)
  writeData(wb, panel, raw,   startRow = 1, startCol = 1)
  writeData(wb, panel, stats, startRow = 1, startCol = 4, colNames = FALSE)
  addStyle(wb, panel, arial, rows = 1:20, cols = 1:6, gridExpand = TRUE)

  cat(sprintf("%-10s  %s  p = %s  %s\n",
              panel, stats[[2]][4], stats[[2]][7], stats[[2]][8]))
}

saveWorkbook(wb, OUTFILE, overwrite = TRUE)
cat("written to", OUTFILE, "\n")
