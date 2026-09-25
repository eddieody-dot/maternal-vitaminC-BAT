## anova_oneway.R --------------------------------------------------------------
## One-way ANOVA followed by Dunnett's post hoc test (each group vs one control).
##
## Produces one worksheet holding one block per dataset, laid out side by side,
## matching the source-data layout of Figure 5K.
##
## The example data are Figure 5K as published with the paper.
##
## Usage: edit SHEET and `blocks`, then run.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------
## Installed automatically on first run.

for (pkg in c("openxlsx", "DescTools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("installing ", pkg, " ...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(openxlsx)
library(DescTools)

OUTFILE <- "Figure5K.xlsx"
SHEET   <- "Figure5K"

TEST_NAME <- "One-way ANOVA with Dunnett's post hoc test (vs control)"


## ---- data -------------------------------------------------------------------
## One entry per block, left to right. The first group listed is the control.
## Group names become the column headers and the post hoc labels, so write them
## exactly as in the figure legend.
##
## Note: do not put an ASCII hyphen inside a group name - DescTools separates
## comparison labels with "-". Use the minus sign (U+2212) as below.

blocks <- list(

  list(
    "Vector (+Vc)"  = c(0.215, 0.213, 0.348, 0.291),
    "Vector (−Vc)"  = c(1.127, 1.191, 1.142, 1.294),
    "shEsrrg (+Vc)" = c(0.802, 0.647, 0.350, 0.338),
    "shKdm6b (+Vc)" = c(2.517, 1.615, 1.303, 1.568)
  ),

  list(
    "Vector (+Vc)"  = c(2.045, 1.843, 1.528, 1.044),
    "Vector (−Vc)"  = c(0.019, 0.094, 0.036, 0.061),
    "shEsrrg (+Vc)" = c(0.006, 0.018, 0.019, 0.022),
    "shKdm6b (+Vc)" = c(0.119, 0.522, 0.448, 0.474)
  )

)

GAP <- 1                 # blank columns between blocks


## ---- helpers ----------------------------------------------------------------

star <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## Dunnett t statistics are derived from the ANOVA mean square error, which is
## the quantity DescTools uses internally but does not return.
dunnett_t <- function(groups, fit) {
  mse   <- summary(fit)[[1]][["Mean Sq"]][2]
  ctrl  <- groups[[1]]
  vapply(groups[-1], function(v)
    round((mean(v) - mean(ctrl)) / sqrt(mse * (1 / length(v) + 1 / length(ctrl))), 3),
    numeric(1))
}


## ---- run --------------------------------------------------------------------

wb    <- createWorkbook()
arial <- createStyle(fontName = "Arial", fontSize = 10)
addWorksheet(wb, SHEET)

col <- 1

for (groups in blocks) {

  control <- names(groups)[1]

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
    c("statistics", "Test", "Groups", "Control", "n per group",
      "F", "df1", "df2", "p", "Signif"),
    c("",
      TEST_NAME,
      paste(names(groups), collapse = ", "),
      control,
      paste(lengths(groups), collapse = ", "),
      round(a[["F value"]][1], 3),
      a[["Df"]][1],
      a[["Df"]][2],
      signif(a[["Pr(>F)"]][1], 4),
      star(a[["Pr(>F)"]][1])),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  ## post hoc: t from the ANOVA MSE, adjusted p from DescTools
  d    <- DunnettTest(value ~ group, data = long, control = control)[[1]]
  padj <- d[, "pval"]

  post <- data.frame(
    `Post hoc` = paste(names(groups)[-1], "\u2212", control),
    t          = unname(dunnett_t(groups, fit)),
    `p adj`    = signif(unname(padj), 4),
    Signif     = vapply(padj, star, character(1)),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  rownames(post) <- NULL

  r_stats <- n_max + 2              # raw header + values, then statistics
  r_post  <- r_stats + nrow(omni) + 1

  writeData(wb, SHEET, raw,  startRow = 1,       startCol = col)
  writeData(wb, SHEET, omni, startRow = r_stats, startCol = col, colNames = FALSE)
  writeData(wb, SHEET,
            paste0("Post hoc: Dunnett's test (vs ", control, ")"),
            startRow = r_post, startCol = col)
  writeData(wb, SHEET, post, startRow = r_post + 1, startCol = col)

  col <- col + length(groups) + GAP

  cat(sprintf("F = %s  p = %s  %s\n",
              omni[[2]][6], omni[[2]][9], omni[[2]][10]))
}

addStyle(wb, SHEET, arial, rows = 1:30, cols = 1:20, gridExpand = TRUE)
saveWorkbook(wb, OUTFILE, overwrite = TRUE)
cat("written to", OUTFILE, "sheet", SHEET, "\n")
