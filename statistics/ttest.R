## ttest.R ---------------------------------------------------------------------
## Two-tailed unpaired t-test, with Welch correction applied when an F test
## indicates unequal variances (F test p < 0.05).
##
## Produces one worksheet: raw values on top, statistics table underneath,
## matching the source-data layout of Figure 1B.
##
## The example data are Figure 1B as published with the paper.
##
## Usage: edit SHEET, GROUP1, GROUP2 and `levels` / `x1` / `x2`, then run.
## -----------------------------------------------------------------------------

library(openxlsx)

OUTFILE    <- "Figure1B.xlsx"
SHEET      <- "Figure1B"
LEVEL_NAME <- "day"                 # first column header: day / week / gene ...
GROUP1     <- "Vc-Dams"
GROUP2     <- "dVc-Dams"

TEST_NAME  <- "Two-tailed unpaired t-test (Welch correction where F test p < 0.05)"


## ---- data -------------------------------------------------------------------
## One entry per level, in the order they should appear.

levels <- c("E0.5", "E6.5", "E11.5", "E13.5", "E18.5", "PPD7")

x1 <- list(                          # GROUP1
  c(79.47,  91.00, 80.73),
  c(73.57, 108.33, 73.70),
  c(77.54,  89.10, 71.50),
  c(58.60,  68.87, 66.83),
  c(41.30,  41.10, 35.40),
  c(72.13, 107.53, 84.73)
)

x2 <- list(                          # GROUP2
  c(82.93, 89.97, 92.63),
  c(91.97, 74.57, 77.07),
  c( 3.10,  1.01,  2.56),
  c( 2.53,  1.20,  3.83),
  c( 1.73,  4.77,  2.73),
  c(73.57, 80.90, 77.83)
)


## ---- helpers ----------------------------------------------------------------

star <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## One row of the statistics table.
ttest_row <- function(level, a, b) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]

  f_p       <- var.test(a, b)$p.value
  use_equal <- f_p >= 0.05
  res       <- t.test(a, b, var.equal = use_equal)

  data.frame(
    level      = level,
    n          = paste0(length(a), ", ", length(b)),
    `F test p` = round(f_p, 4),
    t          = round(unname(res$statistic), 3),
    df         = round(unname(res$parameter), 3),
    p          = signif(res$p.value, 4),
    Signif     = star(res$p.value),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}


## ---- run --------------------------------------------------------------------

stats <- do.call(rbind, Map(ttest_row, levels, x1, x2))
names(stats)[1] <- LEVEL_NAME
rownames(stats) <- NULL

## raw block: one row per level, GROUP1 replicates then GROUP2 replicates
raw <- cbind(
  setNames(data.frame(levels, stringsAsFactors = FALSE), LEVEL_NAME),
  do.call(rbind, lapply(x1, function(v) as.data.frame(t(v)))),
  do.call(rbind, lapply(x2, function(v) as.data.frame(t(v))))
)

n1 <- length(x1[[1]])
n2 <- length(x2[[1]])
raw_header <- c(LEVEL_NAME, GROUP1, rep("", n1 - 1), GROUP2, rep("", n2 - 1))

row_stats <- nrow(raw) + 3          # blank row between raw block and statistics

wb    <- createWorkbook()
arial <- createStyle(fontName = "Arial", fontSize = 10)
addWorksheet(wb, SHEET)

writeData(wb, SHEET, t(raw_header), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, SHEET, raw,           startRow = 2, startCol = 1, colNames = FALSE)

writeData(wb, SHEET, "statistics",  startRow = row_stats,     startCol = 1)
writeData(wb, SHEET, t(c("Test", TEST_NAME)),
          startRow = row_stats + 1, startCol = 1, colNames = FALSE)
writeData(wb, SHEET, stats,         startRow = row_stats + 2, startCol = 1)

addStyle(wb, SHEET, arial, rows = 1:30, cols = 1:10, gridExpand = TRUE)

saveWorkbook(wb, OUTFILE, overwrite = TRUE)

print(stats)
cat("written to", OUTFILE, "sheet", SHEET, "\n")
