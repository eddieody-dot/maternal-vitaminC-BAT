## 01_clustering.R -------------------------------------------------------------
## snRNA-seq of neonatal interscapular brown adipose tissue.
## Reads CellRanger output for four samples, merges and integrates them,
## clusters, and exports cluster markers and between-group DEGs.
##
## Run this first; 02_apc_trajectory.R continues from the object saved here.
## -----------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------

for (pkg in c("Seurat", "SeuratObject", "patchwork", "ggplot2",
              "dplyr", "tidyr", "openxlsx")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("installing ", pkg, " ...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)


## ---- settings ---------------------------------------------------------------
## DATA_DIR holds one CellRanger output folder per sample, each containing
## filtered_feature_bc_matrix/.

DATA_DIR <- "cellranger_count"
OUT_RDS  <- "Merge.rds"

## sample name -> CellRanger run folder
SAMPLES <- c(
  VcF1_P0  = "XCSC25-RA06",
  dVcF1_P0 = "XCSC25-RA07",
  VcF1_P7  = "XCSC25-RA05",
  dVcF1_P7 = "XCSC25-RA04"
)

## order used for all split.by panels
SAMPLE_LEVELS <- c("VcF1_P0", "dVcF1_P0", "VcF1_P7", "dVcF1_P7")

## QC thresholds
MIN_CELLS    <- 3
MIN_FEATURES <- 200
MAX_FEATURES <- 6000
MAX_PERCENT_MT <- 10

## clustering
## These are the values used for the published clustering: neighbours are built
## on 10 components, UMAP and t-SNE on 15. They are kept separate deliberately
## so that the published result can be reproduced exactly.
DIMS_NEIGHBOURS <- 10
DIMS_EMBEDDING  <- 15
RESOLUTION      <- 0.3


## ---- read and merge ---------------------------------------------------------

objects <- lapply(names(SAMPLES), function(s) {
  mat <- Read10X(data.dir = file.path(DATA_DIR, SAMPLES[[s]],
                                      "filtered_feature_bc_matrix"))
  obj <- CreateSeuratObject(counts = mat, project = s,
                            min.cells = MIN_CELLS, min.features = MIN_FEATURES)
  obj$sample <- s
  obj
})
names(objects) <- names(SAMPLES)

merged <- merge(objects[[1]], y = objects[-1],
                project = "BAT", add.cell.ids = names(objects))


## ---- quality control --------------------------------------------------------

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")

VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample")

merged <- subset(merged,
                 subset = nFeature_RNA > MIN_FEATURES &
                          nFeature_RNA < MAX_FEATURES &
                          percent.mt   < MAX_PERCENT_MT)


## ---- normalisation and integration ------------------------------------------

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

merged <- IntegrateLayers(object = merged, method = CCAIntegration,
                          orig.reduction = "pca",
                          new.reduction  = "integrated.cca",
                          verbose = FALSE)
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

ElbowPlot(merged)


## ---- clustering -------------------------------------------------------------

merged <- FindNeighbors(merged, reduction = "integrated.cca",
                        dims = 1:DIMS_NEIGHBOURS)
merged <- FindClusters(merged, resolution = RESOLUTION)
merged <- RunUMAP(merged, reduction = "integrated.cca",
                  dims = 1:DIMS_EMBEDDING)
merged <- RunTSNE(merged, dims = 1:DIMS_EMBEDDING)

## The published clustering yields 11 clusters (0-10), of which cluster 0 is
## the Pdgfra+ progenitor compartment and clusters 1 and 6 are the Adipoq+
## adipocytes carried forward in 02_apc_trajectory.R.
cat("clusters found:", nlevels(merged$seurat_clusters), "\n")

## cell barcodes carry the sample prefix added by merge(); restore the column
## and fix its level order for plotting
merged$sample <- sub("_[ACGT]+-\\d+$", "", colnames(merged))
merged$sample <- factor(merged$sample, levels = SAMPLE_LEVELS)

table(merged$sample, useNA = "ifany")
table(Idents(merged), merged$sample)

ggsave("umap_clusters.pdf",
       DimPlot(merged, reduction = "umap", label = TRUE),
       width = 6, height = 6)
ggsave("umap_by_sample.pdf",
       DimPlot(merged, reduction = "umap", split.by = "sample", label = TRUE),
       width = 14, height = 5)
ggsave("tsne_by_sample.pdf",
       DimPlot(merged, reduction = "tsne", split.by = "sample", label = TRUE),
       width = 14, height = 5)


## ---- cluster composition ----------------------------------------------------

composition <- merged@meta.data %>%
  group_by(sample, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

write.xlsx(composition, "cluster_composition.xlsx")

ggsave("cluster_composition.pdf",
       ggplot(composition, aes(x = sample, y = freq, fill = seurat_clusters)) +
         geom_bar(stat = "identity", position = "fill") +
         geom_text(aes(label = scales::percent(freq, accuracy = 1)),
                   position = position_fill(vjust = 0.5),
                   size = 3, colour = "white") +
         scale_y_continuous(labels = scales::percent) +
         labs(x = "Sample", y = "Cell fraction", fill = "Cluster") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 6, height = 6)


## ---- cluster markers --------------------------------------------------------

markers <- FindAllMarkers(merged, only.pos = TRUE)
write.xlsx(markers, "cluster_markers.xlsx")

top10 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.25) %>%
  slice_head(n = 10) %>%
  ungroup()

merged <- ScaleData(merged, features = unique(top10$gene))

ggsave("cluster_marker_heatmap.pdf",
       DoHeatmap(merged, features = unique(top10$gene), label = FALSE) +
         theme(axis.text.y  = element_text(size = 3),
               axis.text.x  = element_blank(),
               axis.ticks.x = element_blank()),
       width = 12, height = 8)


## ---- between-group DEGs, cluster by cluster ---------------------------------
## Each contrast is written to one workbook, one worksheet per cluster.

CONTRASTS <- list(
  c("dVcF1_P0", "VcF1_P0"),      # effect of deficiency at P0
  c("dVcF1_P7", "VcF1_P7"),      # effect of deficiency at P7
  c("VcF1_P7",  "VcF1_P0"),      # postnatal change, replete
  c("dVcF1_P7", "dVcF1_P0")      # postnatal change, deficient
)

for (contrast in CONTRASTS) {

  wb <- createWorkbook()

  for (cl in levels(Idents(merged))) {
    message("contrast ", contrast[1], " vs ", contrast[2], ", cluster ", cl)

    deg <- try(FindMarkers(merged,
                           ident.1 = contrast[1], ident.2 = contrast[2],
                           subset.ident = cl, group.by = "sample"),
               silent = TRUE)
    if (inherits(deg, "try-error")) next      # too few cells in this cluster

    deg$gene <- rownames(deg)
    sheet <- gsub("[^A-Za-z0-9]", "_", cl)
    addWorksheet(wb, sheet)
    writeData(wb, sheet, deg)
  }

  saveWorkbook(wb,
               sprintf("DEGs_%s_vs_%s.xlsx", contrast[1], contrast[2]),
               overwrite = TRUE)
}


## ---- save -------------------------------------------------------------------

saveRDS(merged, file = OUT_RDS)
cat("saved", OUT_RDS, "\n")
