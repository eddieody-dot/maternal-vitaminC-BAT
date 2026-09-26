## anova_twoway_tukey.R -------------------------------------------------------
## Two-way ANOVA followed by Tukey's HSD post hoc test.
##
## Factors are reported generically as "Factor 1" and "Factor 2" with their
## levels in brackets, so the output never has to commit to calling something a
## genotype or a treatment. Post hoc comparisons are written in a fixed order;
## differences and t values are signed to match that order.
##
## The example data are Figure S9E as published with the paper.
##
## Usage: edit SHEET, PREFIX, F1, F2 and `cells`, then run.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("installing openxlsx ...")
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

library(openxlsx)

OUTFILE <- "FigureS9E.xlsx"
SHEET   <- "FigureS9E"
PREFIX  <- "iBAT"        # prepended to the cell names in the output

TEST_NAME <- "Two-way ANOVA with Tukey's HSD post hoc test"


## ---- data -------------------------------------------------------------------
## One entry per panel. Cell names are "<f1 level>_<f2 level>"; the underscore
## separates the two factors and must not appear inside a level name.

F1 <- c("Vc", "dVc")            # factor 1 levels, in order
F2 <- c("RT", "Cold")           # factor 2 levels, in order

cells <- list(
  "Vc_RT"    = c(0.9315, 0.8395, 1.2291),
  "dVc_RT"   = c(1.1232, 1.3732, 1.0479),
  "Vc_Cold"  = c(1.3357, 0.9846, 1.7262),
  "dVc_Cold" = c(1.9966, 0.8631, 1.2991)
)


## ---- helpers ----------------------------------------------------------------

star <- function(p) {
  if (is.na(p)) "" else
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

## Fixed comparison order: within each level of factor 2, then within each level
## of factor 1, then the two diagonals.
pair_order <- function(f1, f2) list(
  c(paste(f1[1], f2[1], sep = "_"), paste(f1[2], f2[1], sep = "_")),
  c(paste(f1[1], f2[2], sep = "_"), paste(f1[2], f2[2], sep = "_")),
  c(paste(f1[1], f2[2], sep = "_"), paste(f1[1], f2[1], sep = "_")),
  c(paste(f1[2], f2[2], sep = "_"), paste(f1[2], f2[1], sep = "_")),
  c(paste(f1[1], f2[2], sep = "_"), paste(f1[2], f2[1], sep = "_")),
  c(paste(f1[2], f2[2], sep = "_"), paste(f1[1], f2[1], sep = "_"))
)


## ---- run --------------------------------------------------------------------

wb     <- createWorkbook()
arial  <- createStyle(fontName = "Arial", fontSize = 10)
sheet  <- SHEET
prefix <- PREFIX

parts <- do.call(rbind, strsplit(names(cells), "_", fixed = TRUE))
long  <- data.frame(
  value = unlist(cells, use.names = FALSE),
  f1    = factor(rep(parts[, 1], lengths(cells)), levels = F1),
  f2    = factor(rep(parts[, 2], lengths(cells)), levels = F2),
  cell  = factor(rep(names(cells), lengths(cells)), levels = names(cells)),
  stringsAsFactors = FALSE
)

fit <- aov(value ~ f1 * f2, data = long)
a   <- summary(fit)[[1]]

main <- data.frame(
  Source    = c(sprintf("Factor 1 (%s)", paste(F1, collapse = "/")),
                sprintf("Factor 2 (%s)", paste(F2, collapse = "/")),
                "Interaction", "Residuals"),
  Df        = a[["Df"]],
  `Sum Sq`  = round(a[["Sum Sq"]],  4),
  `Mean Sq` = round(a[["Mean Sq"]], 4),
  F         = c(round(a[["F value"]][1:3], 4), NA),
  p         = c(signif(a[["Pr(>F)"]][1:3], 4), NA),
  Signif    = c(vapply(a[["Pr(>F)"]][1:3], star, character(1)), ""),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

## post hoc on the four cells, then reordered and re-signed
tk  <- TukeyHSD(aov(value ~ cell, data = long))$cell
mse <- a[["Mean Sq"]][4]

get_pair <- function(a1, b1) {
  k1 <- paste(a1, b1, sep = "-"); k2 <- paste(b1, a1, sep = "-")
  if (k1 %in% rownames(tk)) return(list(d =  tk[k1, "diff"], p = tk[k1, "p adj"]))
  if (k2 %in% rownames(tk)) return(list(d = -tk[k2, "diff"], p = tk[k2, "p adj"]))
  NULL
}

ord  <- pair_order(F1, F2)
post <- do.call(rbind, lapply(ord, function(pr) {
  v <- get_pair(pr[1], pr[2])
  if (is.null(v)) return(NULL)
  n1 <- length(cells[[pr[1]]]); n2 <- length(cells[[pr[2]]])
  data.frame(
    `Post hoc`  = paste(paste0(prefix, "-", sub("^([^_]+)_(.+)$", "\\2_\\1", pr[1])),
                        "\u2212",
                        paste0(prefix, "-", sub("^([^_]+)_(.+)$", "\\2_\\1", pr[2]))),
    Difference  = round(v$d, 4),
    t           = round(v$d / sqrt(mse * (1 / n1 + 1 / n2)), 3),
    `p adj`     = signif(v$p, 4),
    Signif      = star(v$p),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}))
rownames(post) <- NULL

## raw block: one column per replicate, grouped by cell
n_max <- max(lengths(cells))
raw   <- as.data.frame(
  do.call(cbind, lapply(cells, function(v) c(v, rep(NA, n_max - length(v))))),
  check.names = FALSE
)
raw_header <- unlist(lapply(names(cells), function(nm)
  c(paste0(prefix, "-", sub("^([^_]+)_(.+)$", "\\2_\\1", nm)),
    rep("", n_max - 1))))
raw_wide <- as.data.frame(matrix(unlist(lapply(cells, function(v)
  c(v, rep(NA, n_max - length(v))))), nrow = 1), check.names = FALSE)

## fixed layout: header / values / blank / blank / statistics
r_stats <- 5                      # "statistics"
r_main  <- r_stats + 2            # Source header row
r_post  <- r_main + nrow(main) + 2  # Post hoc header row

addWorksheet(wb, sheet)
writeData(wb, sheet, t(raw_header), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, sheet, raw_wide,      startRow = 2, startCol = 1, colNames = FALSE)

writeData(wb, sheet, "statistics", startRow = r_stats, startCol = 1)
writeData(wb, sheet, t(c("Test", TEST_NAME)),
          startRow = r_stats + 1, startCol = 1, colNames = FALSE)
writeData(wb, sheet, main, startRow = r_main, startCol = 1)
writeData(wb, sheet, post, startRow = r_post, startCol = 1)

addStyle(wb, sheet, arial, rows = 1:30, cols = 1:15, gridExpand = TRUE)

print(main)

saveWorkbook(wb, OUTFILE, overwrite = TRUE)
cat("written to", OUTFILE, "\n")
