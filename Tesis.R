############################################################
# 1) Librerías
############################################################
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(patchwork)
library(tidyverse)
library(TCseq)
library(reshape2)

setwd("C:/Users/alejo/OneDrive/Escritorio/Doctorado/RNAseq")

############################################################
# 2) Funciones auxiliares
############################################################

# Procesar estadio a partir de archivos de conteo

procesar_estadio <- function(ruta, patron, nombres = NULL) {
  archivos <- list.files(ruta, pattern = patron, full.names = TRUE)
  tablas <- lapply(archivos, function(x) read.table(x, header = TRUE, comment.char = "#"))
  
  # Filtrar solo archivos válidos
  ncol_tablas <- sapply(tablas, ncol)
  tablas_filtradas <- tablas[ncol_tablas == 7]
  archivos_filtrados <- archivos[ncol_tablas == 7]
  
  if (is.null(nombres)) {
    nombres <- gsub("_counts.txt", "", basename(archivos_filtrados))
  }
  
  meta <- tablas_filtradas[[1]][, 1:6]
  
  counts <- mapply(function(tabla, nombre) {
    col <- tabla[, 7, drop = FALSE]
    colnames(col) <- nombre
    col
  }, tablas_filtradas, nombres, SIMPLIFY = FALSE)
  
  datos <- cbind(meta, do.call(cbind, counts))
  
  # Calcular RPK y TPM
  for (nombre in nombres) {
    rpk <- datos[[nombre]] / (datos$Length / 1000)
    datos[[paste0("rpk_", nombre)]] <- rpk
    datos[[paste0("tpm_", nombre)]] <- rpk / sum(rpk, na.rm = TRUE) * 1e6
  }
  
  return(datos)
}

# Parsear campos del GFF
parse_uniprot <- function(x) {
  parts <- unlist(strsplit(x, "=UniRef100_"))
  parent_data <- parts[2]
  parts <- unlist(strsplit(parent_data, ";"))
  parts[1]
}

parse_ID_data <- function(x) {
  parts <- unlist(strsplit(x, "ID="))
  parent_data <- parts[2]
  parts <- unlist(strsplit(parent_data, ";"))
  parts[1]
}

############################################################
# 3) Procesamiento de datos por estadio
############################################################
epimastigotes <- procesar_estadio(
  ruta = "conteos/epis",
  patron = "epis.*_counts.txt",
  nombres = c("epis1", "epis2", "epis3")
)

tripomastigotes <- procesar_estadio(
  ruta = "conteos/tripos",
  patron = "tripos.*_counts.txt"
)

amastigotes_intraC <- procesar_estadio(
  ruta = "conteos",
  patron = "^AI[1-4]_counts.txt",
  nombres = c("AI1", "AI2", "AI3", "AI4")
)

amastigotes_axenicos <- procesar_estadio(
  ruta = "conteos/amas_ax",
  patron = ".*_counts.txt"
)

############################################################
# 4) Anotación con GFF
############################################################
gff <- read.delim("TcDm28cT2T_manualCurated.gff", header = FALSE)
gff2 <- gff
gff2$uniref <- sapply(gff$V9, parse_uniprot)
gff2$ID <- sapply(gff$V9, parse_ID_data)
gff2$ID_clean <- sub("\\.t[0-9]+$", "", gff2$ID)

# Merge con cada estadio
epimastigotes <- merge(epimastigotes, gff2[, c("ID_clean", "uniref")],
                       by.x = "Geneid", by.y = "ID_clean", all.x = TRUE)
tripomastigotes <- merge(tripomastigotes, gff2[, c("ID_clean", "uniref")],
                         by.x = "Geneid", by.y = "ID_clean", all.x = TRUE)
amastigotes_intraC <- merge(amastigotes_intraC, gff2[, c("ID", "uniref")],
                            by.x = "Geneid", by.y = "ID", all.x = TRUE)
amastigotes_axenicos <- merge(amastigotes_axenicos, gff2[, c("ID", "uniref")],
                              by.x = "Geneid", by.y = "ID", all.x = TRUE)

############################################################
# 5) Construcción de matriz de conteos
############################################################
conteos_list <- list(
  epis1 = epimastigotes[, 7],
  epis2 = epimastigotes[, 8],
  epis3 = epimastigotes[, 9],
  tripos1 = tripomastigotes[, 7],
  tripos2 = tripomastigotes[, 8],
  tripos3 = tripomastigotes[, 9],
  tripos6 = tripomastigotes[, 10],
  AI1 = amastigotes_intraC[, 7],
  AI2 = amastigotes_intraC[, 8],
  AI3 = amastigotes_intraC[, 9],
  AI4 = amastigotes_intraC[, 10]
)

orden_ax <- c(
  paste0("A", 1:3, "D0"),
  paste0("A", 1:3, "D1"),
  paste0("A", 1:3, "D5"),
  paste0("A", 1:3, "D10")
)

for (n in orden_ax) {
  conteos_list[[n]] <- amastigotes_axenicos[, n]
}

matriz_conteos <- round(do.call(cbind, conteos_list))
rownames(matriz_conteos) <- epimastigotes$Geneid

condicion <- c(
  rep("epimastigote", 3),
  rep("trypomastigote", 4),
  rep("amastigote", 4),
  rep("axenico_D0", 3),
  rep("axenico_D1", 3),
  rep("axenico_D5", 3),
  rep("axenico_D10", 3)
)

colData <- data.frame(
  row.names = colnames(matriz_conteos),
  condition = factor(condicion, levels = unique(condicion))
)

dds <- DESeqDataSetFromMatrix(matriz_conteos, colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = TRUE)

# --- Matriz de conteos para amastigotes axénicos
counts_axenicos <- amastigotes_axenicos[, grepl("^A[1-3]D(0|1|5|10)$", colnames(amastigotes_axenicos))]
rownames(counts_axenicos) <- amastigotes_axenicos$Geneid

# Vector numérico de tiempo (días)
tiempo_num <- c(rep(0,3), rep(1,3), rep(5,3), rep(10,3))

# ColData con tiempo como numérico (para tendencia continua)
colData_ax <- data.frame(time = tiempo_num, row.names = colnames(counts_axenicos))

# DESeqDataSet con tiempo continuo
dds_axenicos <- DESeqDataSetFromMatrix(
  countData = counts_axenicos,
  colData = colData_ax,
  design = ~ time
)

# Filtrar y correr DESeq
dds_axenicos <- dds_axenicos[rowSums(counts(dds_axenicos)) > 10, ]
dds_axenicos <- DESeq(dds_axenicos)

# Transformación VST
vsd_axenicos <- vst(dds_axenicos, blind = TRUE)

# --- Matriz de conteos combinada por estadio puro
counts_estadio <- cbind(
  epimastigotes[, grep("^(epis)[0-9]+$", colnames(epimastigotes))],
  tripomastigotes[, grep("^(tripos)[0-9]+$", colnames(tripomastigotes))],
  amastigotes_intraC[, grep("^(AI)[0-9]+$", colnames(amastigotes_intraC))]
)
rownames(counts_estadio) <- epimastigotes$Geneid

# --- Condiciones ---
condiciones_estadio <- data.frame(
  row.names = colnames(counts_estadio),
  condition = factor(c(
    rep("epimastigote", ncol(epimastigotes[, grep("^(epis)[0-9]+$", colnames(epimastigotes))])),
    rep("trypomastigote", ncol(tripomastigotes[, grep("^(tripos)[0-9]+$", colnames(tripomastigotes))])),
    rep("amastigote", ncol(amastigotes_intraC[, grep("^(AI)[0-9]+$", colnames(amastigotes_intraC))]))
  ))
)

# --- Crear DESeqDataSet para heatmaps ---
dds_estadio <- DESeqDataSetFromMatrix(
  countData = round(counts_estadio),
  colData = condiciones_estadio,
  design = ~ condition
)

# --- Filtrar genes de baja cuenta ---
dds_estadio <- dds_estadio[rowSums(counts(dds_estadio)) > 10, ]
dds_estadio <- DESeq(dds_estadio)
# --- Transformación VST ---
vsd_estadio <- vst(dds_estadio, blind = TRUE)
mat_vsd_estadio <- assay(vsd_estadio)

############################################################
# 6) PCA + Clustering
############################################################
# --- PCA ---
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_coords <- pcaData[, c("PC1", "PC2")]

# --- Método del codo ---
wss <- sapply(1:10, function(k) {
  kmeans(pca_coords, centers = k, nstart = 10)$tot.withinss
})
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Número de clusters k",
     ylab = "Suma total de cuadrados intra-cluster",
     main = "Método del codo para elegir k")

abline(v = 5, lty = 2)
text(5, max(wss)*0.95, "k=5", pos = 4)

set.seed(123)
k <- 5
kmeans_result <- kmeans(pca_coords, centers = k)
pcaData$cluster <- as.factor(kmeans_result$cluster)

ggplot(pcaData, aes(PC1, PC2, color = condition, shape = cluster)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 30) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(3, 15, 19, 17, 18)) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  theme_minimal(base_size = 14) +
  ggtitle(paste("K-means clustering (k =", k, ") sobre PCA"))

# --- PCA solo axénicos ---
# Crear versión con time categórico
vsd_ax_pca <- vsd_axenicos
colData(vsd_ax_pca)$time <- factor(rep(c("D0", "D1", "D5", "D10"), each = 3),
                                   levels = c("D0", "D1", "D5", "D10"))

pcaData_axenicos <- plotPCA(vsd_ax_pca, intgroup = "time", returnData = TRUE)
percentVar_axenicos <- round(100 * attr(pcaData_axenicos, "percentVar"))

ggplot(pcaData_axenicos, aes(PC1, PC2, color = time)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_axenicos[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_axenicos[2], "% variance")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

############################################################
# 7) Heatmaps
############################################################

# 7.1) Heatmap con top 75 genes más variables
var_genes <- apply(assay(vsd), 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:75]
mat_heatmap <- assay(vsd)[top_genes, ]

annotation_col <- data.frame(Condition = colData$condition)
rownames(annotation_col) <- rownames(colData)

col_cond <- brewer.pal(length(unique(colData$condition)), "Set2")
names(col_cond) <- unique(colData$condition)

pheatmap(
  mat_heatmap, scale = "row",
  annotation_col = annotation_col,
  annotation_colors = list(Condition = col_cond),
  clustering_method = "ward.D2",
  show_rownames = TRUE, fontsize = 12,
  main = "Heatmap de expresión (top 75 genes variables)"
)

# 7.2) Heatmap con todos los genes
pheatmap(
  assay(vsd), scale = "none",
  annotation_col = annotation_col,
  annotation_colors = list(Condition = col_cond),
  clustering_method = "ward.D2",
  show_rownames = FALSE, fontsize = 10,
  main = "Heatmap de expresión (todos los genes)"
)

# 7.4) Heatmap para todos los genes por condición
condiciones <- unique(colData$condition)
listas_muestras <- split(rownames(colData), colData$condition)
submatrices <- lapply(listas_muestras, function(samples) {
  assay(vsd)[, samples, drop = FALSE]
})

# Estandarizar por fila (z-score)
z_scores <- lapply(submatrices, function(mat) t(scale(t(mat))))

# Calcular promedio real por condición
z_promedios <- lapply(z_scores, function(mat) rowMeans(mat))
z_promedios <- do.call(cbind, z_promedios)
colnames(z_promedios) <- names(z_scores)

# Colores y breaks
color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
breaks <- seq(-2.5, 2.5, length.out = 101)

heatmaps <- list()

for (cond in condiciones) {
  cat("Procesando:", cond, "\n")
  
  mat_rep <- z_scores[[cond]]
  mat_avg <- z_promedios[, cond, drop = FALSE]
  
  # Heatmap de réplicas
  if (!is.null(mat_rep) && ncol(mat_rep) > 0 && nrow(mat_rep) > 0) {
    heat_rep <- tryCatch(
      pheatmap(mat_rep,
               cluster_rows = FALSE, cluster_cols = TRUE,
               show_rownames = FALSE, main = paste(cond, "(réplicas)"),
               color = color_palette, breaks = breaks, silent = TRUE)[[4]],
      error = function(e) {
        warning(paste("Error en", cond, "(réplicas):", e$message))
        NULL
      }
    )
    heatmaps[[paste0(cond, "_rep")]] <- heat_rep
  }
  
  # Heatmap de promedio (cada fila distinta)
  if (!is.null(mat_avg) && nrow(mat_avg) > 0) {
    heat_avg <- tryCatch(
      pheatmap(mat_avg,
               cluster_rows = FALSE, cluster_cols = FALSE,
               show_rownames = FALSE, main = paste(cond, "(promedio)"),
               color = color_palette, breaks = breaks, silent = TRUE)[[4]],
      error = function(e) {
        warning(paste("Error en", cond, "(promedio):", e$message))
        NULL
      }
    )
    heatmaps[[paste0(cond, "_avg")]] <- heat_avg
  }
}

# Filtrar heatmaps válidos
heatmaps_validos <- Filter(Negate(is.null), heatmaps)

# Mostrar todos los heatmaps
grid.arrange(grobs = heatmaps_validos, ncol = length(heatmaps_validos))

############################################################
# 8) revisar e integrar
############################################################
# Comparaciones por estadio
res_epi_vs_otros <- results(dds_estadio, contrast = c("condition", "epimastigote", "trypomastigote"))
res_epi_vs_amast <- results(dds_estadio, contrast = c("condition", "epimastigote", "amastigote"))
res_tripo_vs_epi <- results(dds_estadio, contrast = c("condition", "trypomastigote", "epimastigote"))
res_tripo_vs_amas <- results(dds_estadio, contrast = c("condition", "trypomastigote", "amastigote"))
res_amas_vs_epi <- results(dds_estadio, contrast = c("condition", "amastigote", "epimastigote"))
res_amas_vs_tripo <- results(dds_estadio, contrast = c("condition", "amastigote", "trypomastigote"))

# Genes top por estadio
genes_epi <- intersect(rownames(subset(res_epi_vs_otros, padj < 0.05 & log2FoldChange > 1)),
                       rownames(subset(res_epi_vs_amast, padj < 0.05 & log2FoldChange > 1)))
genes_tripo <- intersect(rownames(subset(res_tripo_vs_epi, padj < 0.05 & log2FoldChange > 1)),
                         rownames(subset(res_tripo_vs_amas, padj < 0.05 & log2FoldChange > 1)))
genes_amas <- intersect(rownames(subset(res_amas_vs_epi, padj < 0.05 & log2FoldChange > 1)),
                        rownames(subset(res_amas_vs_tripo, padj < 0.05 & log2FoldChange > 1)))

top_epi <- head(rownames(res_epi_vs_otros[genes_epi, ][order(res_epi_vs_otros[genes_epi, ]$log2FoldChange, decreasing = TRUE), ]), 1000)
top_tripo <- head(rownames(res_tripo_vs_epi[genes_tripo, ][order(res_tripo_vs_epi[genes_tripo, ]$log2FoldChange, decreasing = TRUE), ]), 1000)
top_amas <- head(rownames(res_amas_vs_epi[genes_amas, ][order(res_amas_vs_epi[genes_amas, ]$log2FoldChange, decreasing = TRUE), ]), 1000)

top_genes <- unique(c(top_epi, top_tripo, top_amas))

# Subset y anotación
mat_heatmap <- mat_vsd_estadio[top_genes, ]
mat_heatmap <- mat_heatmap[complete.cases(mat_heatmap), ]
annotation_col <- data.frame(condition = condiciones_estadio$condition)
rownames(annotation_col) <- rownames(condiciones_estadio)

# Paleta y heatmap
palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Crear heatmap y guardar el orden real
#png("Genes_de_mayor_expresión_en_cada_estadio.png", width = 1400, height = 3000, res = 150)
p <- pheatmap(mat_heatmap,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              annotation_col = annotation_col,
              color = palette,
              scale = "row",
              fontsize_row = 6,
              fontsize = 12,
              show_rownames = F,
              main = "Genes de mayor expresión en cada estadio"
)
#dev.off()

# --- NUEVO: heatmap con todos los estadios, bien anotados ---
# Usar orden de genes del primer heatmap
orden_genes <- rownames(mat_heatmap)[p$tree_row$order]

# Filtrar genes que existan en vsd
orden_genes <- orden_genes[orden_genes %in% rownames(assay(vsd))]
if (length(orden_genes) == 0) stop("Ningún gen de orden_genes está en assay(vsd)")

# --- Definir orden correcto de muestras ---
# Estadios puros (de counts_estadio)
muestras_epis <- grep("^epis", colnames(counts_estadio), value = TRUE)
muestras_tripos <- grep("^tripos", colnames(counts_estadio), value = TRUE)
muestras_amas <- grep("^AI", colnames(counts_estadio), value = TRUE)

# Amastigotes axénicos: extraer y ordenar por día y réplica
counts_ax <- amastigotes_axenicos[, grepl("^A[1-3]D(0|1|5|10)$", colnames(amastigotes_axenicos))]
rownames(counts_ax) <- epimastigotes$Geneid

# Extraer día y réplica de nombres originales: A1D0, A2D1, etc.
dias <- gsub("^A[0-9]+(D[0-9]+)$", "\\1", colnames(counts_ax))
reps <- gsub("^A([0-9]+)D[0-9]+$", "\\1", colnames(counts_ax))

# Ordenar por día y réplica
orden_df <- data.frame(
  colname = colnames(counts_ax),
  dia = factor(dias, levels = c("D0", "D1", "D5", "D10")),
  rep = as.integer(reps)
)
orden_df <- orden_df[order(orden_df$dia, orden_df$rep), ]
counts_ax <- counts_ax[, orden_df$colname]

# Orden final: estadios puros + axénicos (con nombres originales: A1D0, A2D0, etc.)
muestras_finales <- c(
  muestras_epis,
  muestras_tripos,
  muestras_amas,
  colnames(counts_ax)  # nombres como "A1D0", que SÍ están en vsd
)

# Filtrar solo muestras que existan en vsd
muestras_finales <- muestras_finales[muestras_finales %in% colnames(assay(vsd))]
if (length(muestras_finales) == 0) stop("Ninguna muestra está en assay(vsd)")

# Extraer matriz VST
mat_vsd_sub <- assay(vsd)[orden_genes, muestras_finales, drop = FALSE]
mat_vsd_sub <- mat_vsd_sub[complete.cases(mat_vsd_sub), , drop = FALSE]

# --- Anotación: etiquetas bonitas sin cambiar nombres ---
condition <- c(
  rep("epimastigote", length(muestras_epis)),
  rep("trypomastigote", length(muestras_tripos)),
  rep("amastigote", length(muestras_amas)),
  rep("axenico_D0", sum(grepl("D0$", colnames(counts_ax)))),
  rep("axenico_D1", sum(grepl("D1$", colnames(counts_ax)))),
  rep("axenico_D5", sum(grepl("D5$", colnames(counts_ax)))),
  rep("axenico_D10", sum(grepl("D10$", colnames(counts_ax))))
)
names(condition) <- muestras_finales

# Convertir a factor con orden deseado
annotation_col <- data.frame(
  condition = factor(condition, levels = c(
    "epimastigote", "trypomastigote", "amastigote",
    "axenico_D0", "axenico_D1", "axenico_D5", "axenico_D10"
  )),
  row.names = names(condition)
)

# Paleta para heatmap
palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Colores para anotación
niveles_condicion <- c(
  "epimastigote", "trypomastigote", "amastigote",
  "axenico_D0", "axenico_D1", "axenico_D5", "axenico_D10"
)
colores_condicion <- RColorBrewer::brewer.pal(7, "Set3")
names(colores_condicion) <- niveles_condicion

# Heatmap final
#png("Expresión_de_genes_marcadores_por_estadio_en_amastigotes_axénicos.png", width = 1400, height = 3000, res = 150)
pheatmap(
  mat_vsd_sub,
  cluster_rows = F,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = list(condition = colores_condicion),
  color = palette,
  scale = "row",
  fontsize_row = 6,
  fontsize = 12,
  show_rownames = F,
  main = "Expresión de genes marcadores por estadio (incluye axénicos)"
)
#dev.off()

# --- Heatmap 1: Top 30 genes marcadores (solo estadios puros) ---
# Top 30 por estadio
top30_epi <- head(top_epi, 30)
top30_tripo <- head(top_tripo, 30)
top30_amas <- head(top_amas, 30)
genes_top30 <- unique(c(top30_epi, top30_tripo, top30_amas))

# Extraer del vsd_estadio (solo estadios puros)
mat_puro <- mat_vsd_estadio[genes_top30, , drop = FALSE]
mat_puro <- mat_puro[complete.cases(mat_puro), ]

# Anotación
annotation_puro <- data.frame(condition = condiciones_estadio$condition)
rownames(annotation_puro) <- colnames(mat_puro)

# Paleta
pal <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

#png("heatmap_top30_puros.png", width = 1000, height = 1200, res = 150)
pheatmap(
  mat_puro,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_puro,
  color = pal,
  scale = "row",
  fontsize_row = 8,
  fontsize_col = 10,
  show_rownames = TRUE,
  main = "Top 30 genes marcadores (estadios puros)"
)
#dev.off()

# --- Heatmap 2: Top 30 con amastigotes axénicos ---
# Usar nombres originales para muestras axénicas
counts_ax <- amastigotes_axenicos[, grepl("^A[1-3]D(0|1|5|10)$", colnames(amastigotes_axenicos))]
rownames(counts_ax) <- epimastigotes$Geneid

# Extraer día y réplica
dias <- gsub("^A[0-9]+(D[0-9]+)$", "\\1", colnames(counts_ax))
reps <- gsub("^A([0-9]+)D[0-9]+$", "\\1", colnames(counts_ax))

# Ordenar por día y réplica
orden_df <- data.frame(colname = colnames(counts_ax), dia = factor(dias, levels = c("D0","D1","D5","D10")), rep = as.integer(reps))
orden_df <- orden_df[order(orden_df$dia, orden_df$rep), ]
counts_ax <- counts_ax[, orden_df$colname]

# Orden final de muestras
muestras_finales <- c(colnames(counts_estadio), colnames(counts_ax))
muestras_finales <- muestras_finales[muestras_finales %in% colnames(assay(vsd))]

# Extraer del vsd completo
mat_full <- assay(vsd)[genes_top30, muestras_finales, drop = FALSE]
mat_full <- mat_full[complete.cases(mat_full), , drop = FALSE]

# Anotación
condition <- c(
  rep("epimastigote", ncol(counts_estadio[, grep("^epis", colnames(counts_estadio))]) ),
  rep("trypomastigote", ncol(counts_estadio[, grep("^tripos", colnames(counts_estadio))]) ),
  rep("amastigote", ncol(counts_estadio[, grep("^AI", colnames(counts_estadio))]) ),
  rep("axenico_D0", sum(grepl("D0$", colnames(counts_ax)))),
  rep("axenico_D1", sum(grepl("D1$", colnames(counts_ax)))),
  rep("axenico_D5", sum(grepl("D5$", colnames(counts_ax)))),
  rep("axenico_D10", sum(grepl("D10$", colnames(counts_ax))))
)
names(condition) <- colnames(mat_full)

annotation_full <- data.frame(
  condition = factor(condition, levels = c("epimastigote", "trypomastigote", "amastigote",
                                           "axenico_D0", "axenico_D1", "axenico_D5", "axenico_D10"))
)
rownames(annotation_full) <- colnames(mat_full)

# Colores anotación
colores <- RColorBrewer::brewer.pal(7, "Set3")
names(colores) <- levels(annotation_full$condition)

# Heatmap final
#png("heatmap_top30_con_axenicos.png", width = 1400, height = 1200, res = 150)
pheatmap(
  mat_full,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_full,
  annotation_colors = list(condition = colores),
  color = pal,
  scale = "row",
  fontsize_row = 8,
  fontsize_col = 10,
  show_rownames = TRUE,
  main = "Top 30 genes marcadores (con axénicos)"
)
#dev.off()

############################################################
# 9) Análisis temporal en amastigotes axénicos
############################################################

# Resultados del efecto del tiempo
res_time <- results(dds_axenicos, name = "time")
summary(res_time)  # Muestra cuántos genes up/down

# Mostrar los mejores candidatos aunque no sean significativos
sig_genes <- subset(as.data.frame(res_time), padj < 0.05)

# Si no hay significativos, tomá los de menor pvalor
if (nrow(sig_genes) == 0) {
  message("No hay genes con padj < 0.05. Mostrando el gen con menor pvalor.")
  res_sorted <- res_time[order(res_time$pvalue), ]
  gen <- rownames(res_sorted)[1]
  message(paste("Gen con menor pvalor:", gen, 
                "| log2FoldChange =", round(res_sorted$log2FoldChange[1], 3),
                "| pvalue =", format(res_sorted$pvalue[1], digits = 2)))
} else {
  gen <- rownames(sig_genes)[1]
  message(paste("Gen significativo encontrado:", gen))
}

# Graficar siempre, sin importar significancia
vsd_tmp <- vst(dds_axenicos, blind = FALSE)
expr_mat <- assay(vsd_tmp)

df <- data.frame(expression = expr_mat[gen, ], time = tiempo_num)

ggplot(df, aes(x = time, y = expression)) +
  geom_point(size = 3, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("Progresión temporal de", gen),
    subtitle = ifelse(nrow(sig_genes) > 0, "Gen significativo (padj < 0.05)", "Mejor candidato (menor pvalor)"),
    x = "Tiempo (días)",
    y = "Expresión (VST)"
  ) +
  scale_x_continuous(breaks = c(0, 1, 5, 10))

############################################################
# 10) Expresión de genes reguladores del cAMP
############################################################

# Genes de interés
genes_objetivo <- c("g3430", "g8087", "g8086", "g5111", "g5087")
genes_objetivo <- genes_objetivo[genes_objetivo %in% rownames(vsd)]

# Diccionario de nombres
nombres_genes <- c(
  g3430 = "PDEA",
  g8087 = "PDEB1",
  g8086 = "PDEB2",
  g5111 = "PDEC2",
  g5087 = "PDED"
)

# Expresión en amastigotes intracelulares (AI1-AI4)
amastigotes_intracelulares <- assay(vsd)[genes_objetivo, c("AI1", "AI2", "AI3", "AI4"), drop = FALSE]
amastigotes_df <- reshape2::melt(amastigotes_intracelulares)
colnames(amastigotes_df) <- c("Gene", "Sample", "VST_expression")

# Promedio y error estándar
amastigotes_prom <- amastigotes_df %>%
  group_by(Gene) %>%
  summarise(
    mean_expr = mean(VST_expression),
    se_expr = sd(VST_expression) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Gene = recode(Gene, !!!nombres_genes))

# Todas las muestras en orden
orden_muestras <- c(
  "tripos1", "tripos2", "tripos3", "tripos6",
  paste0("A", 1:3, "D0"), paste0("A", 1:3, "D1"),
  paste0("A", 1:3, "D5"), paste0("A", 1:3, "D10"),
  "epis1", "epis2", "epis3"
)
orden_muestras <- orden_muestras[orden_muestras %in% colnames(vsd)]

# Extraer expresión
matriz_vst <- assay(vsd)[genes_objetivo, orden_muestras, drop = FALSE]
df_largo <- reshape2::melt(matriz_vst)
colnames(df_largo) <- c("Gene", "Sample", "VST_expression")

# Clasificar por grupo
df_largo <- df_largo %>%
  mutate(
    Group = case_when(
      Sample %in% c("tripos1", "tripos2", "tripos3", "tripos6") ~ "Tripomastigote",
      grepl("^A[1-3]D0$", Sample) ~ "Axenico_D0",
      grepl("^A[1-3]D1$", Sample) ~ "Axenico_D1",
      grepl("^A[1-3]D5$", Sample) ~ "Axenico_D5",
      grepl("^A[1-3]D10$", Sample) ~ "Axenico_D10",
      Sample %in% c("epis1", "epis2", "epis3") ~ "Epimastigote",
      TRUE ~ NA_character_
    ),
    Gene = recode(Gene, !!!nombres_genes)
  )

# Factor con orden deseado
niveles_grupo <- c("Tripomastigote", "Axenico_D0", "Axenico_D1", "Axenico_D5", "Axenico_D10", "Epimastigote")
df_largo$Group <- factor(df_largo$Group, levels = niveles_grupo)

# Estadísticas por grupo
estadisticas <- df_largo %>%
  group_by(Gene, Group) %>%
  summarise(
    mean_expr = mean(VST_expression),
    se_expr = sd(VST_expression) / sqrt(n()),
    .groups = "drop"
  )

# Gráfico 1: Barras (solo amastigotes intracelulares)
p1 <- ggplot(amastigotes_prom, aes(x = Gene, y = mean_expr, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                width = 0.2, size = 0.7) +
  theme_bw() +
  labs(subtitle = "Amastigotes Intracelulares", x = NULL, y = "Expresión (VST)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(8.5, 12.25))

# Gráfico 2: Líneas (todos los estadios)
p2 <- ggplot(estadisticas, aes(x = Group, y = mean_expr, color = Gene, group = Gene)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                width = 0.2, size = 0.7) +
  theme_bw() +
  labs(subtitle = "Variación durante el proceso de Diferenciación ",
       x = "Condición",
       y = "Expresión (VST)",
       color = "Gen") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")+
  coord_cartesian(ylim = c(8.5, 12.25))

# Combinar gráficos
p1 + p2 + plot_layout(widths = c(1, 3))+
  plot_annotation("Expresión de PDEs")


############################################################
# 10.1) ANOVA y Tukey HSD para cada gen PDE
############################################################

# Lista para guardar resultados
resultados_tukey <- list()

for (g in genes_objetivo) {
  # Renombrar gen
  gen_nombre <- recode(g, !!!nombres_genes)
  
  # Datos del gen
  df_gene <- df_largo %>% filter(Gene == gen_nombre)
  
  # ANOVA
  aov_model <- aov(VST_expression ~ Group, data = df_gene)
  
  # Tukey HSD
  tukey <- TukeyHSD(aov_model)
  
  # Convertir a tabla
  tukey_tab <- as.data.frame(tukey$Group)
  tukey_tab$Comparison <- rownames(tukey_tab)
  tukey_tab$Gene <- gen_nombre
  tukey_tab <- tukey_tab %>% select(Gene, Comparison, diff, lwr, upr, `p adj`)
  colnames(tukey_tab)[colnames(tukey_tab) == "p adj"] <- "p.adj"
  
  resultados_tukey[[gen_nombre]] <- tukey_tab
}

# Combinar todos los resultados
tabla_tukey_pde <- bind_rows(resultados_tukey)

# Mostrar primeros resultados
head(tabla_tukey_pde)

# Opcional: guardar en CSV
# write.csv(tabla_tukey_pde, "TukeyHSD_PDEs_por_grupo.csv", row.names = FALSE)

############################################################
# 11) Hipótesis
############################################################
# 11.1) Niveles de cAMP por estadio (AC/PDE balance)
# Definir genes PDE
pde_genes <- c("g3430", "g8087", "g8086", "g5111", "g5087")

# Filtrar candidatos AC del GFF
ac_candidates <- gff2 %>%
  filter(grepl("adenylate_cyclase", V9, ignore.case = TRUE)) %>%
  pull(ID_clean)

# Genes AC presentes en los datos
ac_genes <- intersect(ac_candidates,
                      c(epimastigotes$Geneid,
                        tripomastigotes$Geneid,
                        amastigotes_intraC$Geneid,
                        amastigotes_axenicos$Geneid))

amastigotes_axenicos$Geneid <- sub("\\.t1$", "", amastigotes_axenicos$Geneid)
amastigotes_intraC$Geneid <- sub("\\.t1$", "", amastigotes_intraC$Geneid)

# Función para calcular índice log2(ΣAC / ΣPDE)
calcular_indice <- function(df) {
  tpm_cols <- grep("^tpm_", colnames(df), value = TRUE)
  if (length(tpm_cols) == 0) return(data.frame(muestra = character(0), indice_log2 = numeric(0)))
  
  ac_sum  <- colSums(df[df$Geneid %in% ac_genes,  tpm_cols, drop = FALSE], na.rm = TRUE)
  pde_sum <- colSums(df[df$Geneid %in% pde_genes, tpm_cols, drop = FALSE], na.rm = TRUE)
  
  indice <- log2((ac_sum + 1e-6) / (pde_sum + 1e-6))
  data.frame(muestra = tpm_cols, indice_log2 = as.numeric(indice))
}

# --- Calcular para cada estadio ---
indice_epi <- calcular_indice(epimastigotes) %>% mutate(estadio = "epimastigote")
indice_tripos <- calcular_indice(tripomastigotes) %>% mutate(estadio = "tripomastigote")
indice_amai <- calcular_indice(amastigotes_intraC) %>% mutate(estadio = "amastigote_intraC")
indice_amas_ax <- calcular_indice(amastigotes_axenicos)

# Asignar estadio por día (solo el día, no A1D0)
dias <- c("D0", "D1", "D5", "D10")
indice_amas_ax$estadio <- paste0("amastigote_ax_", dias)

# Combinar todo
df_indices <- bind_rows(indice_epi, indice_tripos, indice_amai, indice_amas_ax)

# Orden biológico
df_indices$estadio <- factor(df_indices$estadio, levels = c(
  "amastigote_intraC",
  "tripomastigote",
  "amastigote_ax_D0",
  "amastigote_ax_D1",
  "amastigote_ax_D5",
  "amastigote_ax_D10",
  "epimastigote"
))

# Colores
colores <- c(
  "amastigote_intraC" = "#FF7D03",
  "tripomastigote" = "#FF0303",
  "amastigote_ax_D0" = "#A4A607",
  "amastigote_ax_D1" = "#028357",
  "amastigote_ax_D5" = "grey",
  "amastigote_ax_D10" = "#490291",
  "epimastigote" = "#FF03EA"
)

# Gráfico final
# ggplot(df_indices, aes(x = estadio, y = indice_log2, color = estadio)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
#   scale_color_manual(values = colores) +
#   labs(x = "Estadio", y = "log2(ΣAC / ΣPDE)", color = "Estadio") +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggplot(df_indices, aes(x = estadio, y = indice_log2, fill = estadio)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = estadio)) +
  scale_fill_manual(values = colores) +
  scale_color_manual(values = colores) +
  labs(#x = "Estadio", 
    y = "log2(ΣAC / ΣPDE)", 
    color = "Estadio", fill = "Estadio") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
        axis.title.y = element_text(size = 15),
        legend.position = "none")


#11.2.1) Expresión de genes AC
# Filtrar genes AC que estén en vsd
genes_objetivo <- ac_genes[ac_genes %in% rownames(vsd)]

# Expresión en amastigotes intracelulares (AI1-AI4)
amastigotes_intracelulares <- assay(vsd)[genes_objetivo, c("AI1", "AI2", "AI3", "AI4"), drop = FALSE]
amastigotes_df <- reshape2::melt(amastigotes_intracelulares)
colnames(amastigotes_df) <- c("Gene", "Sample", "VST_expression")

amastigotes_prom <- amastigotes_df %>%
  group_by(Gene) %>%
  summarise(
    mean_expr = mean(VST_expression),
    se_expr = sd(VST_expression) / sqrt(n()),
    .groups = "drop"
  )

# Todas las muestras en orden
orden_muestras <- c(
  "tripos1", "tripos2", "tripos3", "tripos6",
  paste0("A", 1:3, "D0"), paste0("A", 1:3, "D1"),
  paste0("A", 1:3, "D5"), paste0("A", 1:3, "D10"),
  "epis1", "epis2", "epis3"
)
orden_muestras <- orden_muestras[orden_muestras %in% colnames(vsd)]

# Extraer expresión
matriz_vst <- assay(vsd)[genes_objetivo, orden_muestras, drop = FALSE]
df_largo <- reshape2::melt(matriz_vst)
colnames(df_largo) <- c("Gene", "Sample", "VST_expression")

# Clasificar por grupo
df_largo <- df_largo %>%
  mutate(
    Group = case_when(
      grepl("^tripos", Sample) ~ "Tripomastigote",
      grepl("^A[0-9]+D0$", Sample) ~ "Axenico_D0",
      grepl("^A[0-9]+D1$", Sample) ~ "Axenico_D1",
      grepl("^A[0-9]+D5$", Sample) ~ "Axenico_D5",
      grepl("^A[0-9]+D10$", Sample) ~ "Axenico_D10",
      grepl("^epis", Sample) ~ "Epimastigote",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Group = factor(Group, levels = c("Tripomastigote", "Axenico_D0", "Axenico_D1", "Axenico_D5", "Axenico_D10", "Epimastigote")))

# Estadísticas por grupo
estadisticas <- df_largo %>%
  group_by(Gene, Group) %>%
  summarise(
    mean_expr = mean(VST_expression),
    se_expr = sd(VST_expression) / sqrt(n()),
    .groups = "drop"
  )

# Gráfico 1: Barras (solo amastigotes intracelulares)
p1 <- ggplot(amastigotes_prom, aes(x = Gene, y = mean_expr, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                width = 0.2, size = 0.7) +
  theme_bw() +
  labs(x = NULL, y = "Expresión (VST)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico 2: Líneas (todos los estadios)
p2 <- ggplot(estadisticas, aes(x = Group, y = mean_expr, color = Gene, group = Gene)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                width = 0.2, size = 0.7) +
  theme_bw() +
  labs(title = "Expresión de genes AC", x = "Condición", y = "Expresión (VST)", color = "Gen") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

# Combinar gráficos
p1 + p2 + plot_layout(widths = c(1, 3))

# 11.2.2) ANOVA y Tukey por gen AC 

resultados_tukey_ac <- list()

for (g in genes_objetivo) {
  df_gene <- df_largo %>% filter(Gene == g)
  aov_model <- aov(VST_expression ~ Group, data = df_gene)
  tukey <- TukeyHSD(aov_model)
  
  tukey_tab <- as.data.frame(tukey$Group)
  tukey_tab$Comparison <- rownames(tukey_tab)
  tukey_tab$Gene <- g
  resultados_tukey_ac[[g]] <- tukey_tab  # ✅ Ahora sí: resultados_tukey_ac
}

# Renombrar y crear tabla
tabla_diferencias_ac <- bind_rows(resultados_tukey_ac) %>%
  rename(p.adj = `p adj`) %>%
  select(Gene, Comparison, diff, lwr, upr, p.adj)

head(tabla_diferencias_ac)

# 11.3) Contrivucion relativa de AC y PDE en niveles de cAMP
# Función para calcular contribución relativa con días para axénicos
calcular_contribucion_relativa <- function(df, genes, tipo, estadio_nombre, dias_map = NULL) {
  tpm_cols <- grep("^tpm_", colnames(df), value = TRUE)
  
  df_sel <- df %>%
    filter(Geneid %in% genes) %>%
    select(Geneid, all_of(tpm_cols)) %>%
    pivot_longer(cols = all_of(tpm_cols), names_to = "muestra", values_to = "TPM") %>%
    group_by(muestra) %>%
    mutate(frac_contrib = TPM / sum(TPM, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      tipo = tipo,
      estadio = if(!is.null(dias_map)) dias_map[muestra] else estadio_nombre
    )
  
  return(df_sel)
}

# Definir días de los axénicos según columnas tpm
ax_cols <- grep("^tpm_A", colnames(amastigotes_axenicos), value = TRUE)
dias_map <- gsub("tpm_.*(D[0-9]+)$", "amastigote_ax_\\1", ax_cols)
names(dias_map) <- ax_cols

# Calcular contribuciones
df_contribucion <- bind_rows(
  # Epimastigotes
  calcular_contribucion_relativa(epimastigotes, ac_genes, "AC", "epimastigote"),
  calcular_contribucion_relativa(epimastigotes, pde_genes, "PDE", "epimastigote"),
  # Tripomastigotes
  calcular_contribucion_relativa(tripomastigotes, ac_genes, "AC", "tripomastigote"),
  calcular_contribucion_relativa(tripomastigotes, pde_genes, "PDE", "tripomastigote"),
  # Amastigotes intracelulares
  calcular_contribucion_relativa(amastigotes_intraC, ac_genes, "AC", "amastigote_intraC"),
  calcular_contribucion_relativa(amastigotes_intraC, pde_genes, "PDE", "amastigote_intraC"),
  # Amastigotes axénicos con días
  calcular_contribucion_relativa(amastigotes_axenicos, ac_genes, "AC", "amastigote_ax", dias_map = dias_map),
  calcular_contribucion_relativa(amastigotes_axenicos, pde_genes, "PDE", "amastigote_ax", dias_map = dias_map)
)

# Orden de factores para eje x
df_contribucion$estadio <- factor(df_contribucion$estadio, levels = c(
  "amastigote_intraC",
  "tripomastigote",
  "amastigote_ax_D0",
  "amastigote_ax_D1",
  "amastigote_ax_D5",
  "amastigote_ax_D10",
  "epimastigote"
))

# Definir paleta de colores para los genes
genes_totales <- unique(df_contribucion$Geneid)
n_genes <- length(genes_totales)
paleta_genes <- RColorBrewer::brewer.pal(min(n_genes, 12), "Set3")
if(n_genes > 12) paleta_genes <- colorRampPalette(paleta_genes)(n_genes)
names(paleta_genes) <- genes_totales

# Para cada estadio, calcular la contribución media de cada gen
orden_genes <- df_contribucion %>%
  group_by(estadio, Geneid) %>%
  summarise(mean_frac = mean(frac_contrib, na.rm = TRUE), .groups="drop") %>%
  arrange(estadio, -mean_frac) %>%
  group_by(estadio) %>%
  mutate(Geneid_ord = fct_reorder(Geneid, mean_frac, .desc = TRUE)) %>%
  ungroup() %>%
  select(estadio, Geneid, Geneid_ord)

# Unir con df_contribucion
df_contribucion_ord <- df_contribucion %>%
  left_join(orden_genes, by=c("estadio","Geneid"))

# # Graficar con genes ordenados
# ggplot(df_contribucion_ord, aes(x=estadio, y=frac_contrib, fill=Geneid_ord)) +
#   geom_bar(stat="identity", position="fill", color="black") +
#   facet_wrap(~tipo, scales="free_y") +
#   scale_fill_manual(values = setNames(paleta_genes[levels(orden_genes$Geneid_ord)], levels(orden_genes$Geneid_ord))) +
#   labs(x="Estadio", y="Contribución relativa (%)", fill="Gen") +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) 

# Promediar contribución por gen y estadio
df_promedio <- df_contribucion_ord %>%
  group_by(estadio, tipo, Geneid) %>%
  summarise(mean_frac = mean(frac_contrib), .groups = "drop") %>%
  ungroup()

# Gráfico con una barra por gen
ggplot(df_promedio, aes(x = estadio, y = mean_frac, fill = Geneid)) +
  geom_col(position = "fill", color = "black", size = 0.5) +  # Borde negro alrededor de cada barra
  facet_wrap(~tipo, scales = "free_y") +
  scale_fill_manual(values = setNames(paleta_genes[levels(orden_genes$Geneid_ord)], levels(orden_genes$Geneid_ord))) +
  labs(x = "Estadio", y = "Contribución relativa (%)", fill = "Gen") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#####
amas_size <- colSums(amastigotes_intraC[,c("AI1", 
                                           #"AI2", 
                                           "AI3"
                                           #"AI4"
                                           )])
epis_size <- colSums(epimastigotes[,c("epis1", "epis2", "epis3")])
tripos_size <- colSums(tripomastigotes[,c("tripos1", 
                                          "tripos2", 
                                          "tripos3", 
                                          "tripos6"
                                          )])

library_sizes <- data.frame(
  size = c(epis_size, tripos_size, amas_size),
  estadio = c(
    rep("Epimastigotes", length(epis_size)),
    rep("Tripomastigotes", length(tripos_size)),
    rep("Amastigotes", length(amas_size))
  )
)
tapply(library_sizes$size, library_sizes$estadio, mean)
aggregate(size ~ estadio, library_sizes, mean)
test<-lm(size ~ estadio, data = library_sizes)
pairwise.wilcox.test(library_sizes$size, library_sizes$estadio,
                     p.adjust.method = "BH")
ggplot(library_sizes, aes(x = estadio, y = size)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Tamaño de librería por estadio",
       y = "Library size (counts totales)",
       x = "")





############################################################
# 4) Unificación de Conteos y Heatmap de Correlación
############################################################

# 1. Función auxiliar para limpiar el Geneid y extraer SOLO los conteos crudos
get_counts <- function(df) {
  df %>% 
    mutate(Geneid = gsub("\\..*", "", Geneid)) %>% 
    select(
      -c(Chr, Start, End, Strand, Length), 
      -starts_with("rpk_"), 
      -starts_with("tpm_")
    )
}

# 2. Unificar todas las tablas en una sola matriz de conteos
lista_conteos <- list(
  get_counts(epimastigotes),
  get_counts(tripomastigotes),
  get_counts(amastigotes_intraC),
  get_counts(amastigotes_axenicos)
)

counts_merged <- lista_conteos %>% 
  reduce(inner_join, by = "Geneid") %>%
  column_to_rownames("Geneid")

# 3. Crear los metadatos (colData) para DESeq2
coldata <- data.frame(row.names = colnames(counts_merged))

coldata$Condicion <- case_when(
  grepl("epis", rownames(coldata)) ~ "Epimastigote",
  grepl("tripos", rownames(coldata)) ~ "Tripomastigote",
  grepl("AI", rownames(coldata)) ~ "Amastigote_IntraC",
  grepl("D10", rownames(coldata)) ~ "Amas_Ax_D10",
  grepl("D0", rownames(coldata)) ~ "Amas_Ax_D0",
  grepl("D1", rownames(coldata)) ~ "Amas_Ax_D1",
  grepl("D5", rownames(coldata)) ~ "Amas_Ax_D5",
  TRUE ~ "Desconocido"
)

coldata$Condicion <- as.factor(coldata$Condicion)

# Validar que matriz y metadatos tengan el mismo orden
stopifnot(all(rownames(coldata) == colnames(counts_merged)))

# 4. Objeto DESeqDataSet y Transformación
dds <- DESeqDataSetFromMatrix(countData = counts_merged,
                              colData = coldata,
                              design = ~ Condicion)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

vsd <- vst(dds, blind = TRUE)

# 5. Cálculo de la matriz de correlación (Pearson)
matriz_correlacion <- cor(assay(vsd), method = "pearson")

# 6. Graficar el Heatmap
colores_heatmap <- colorRampPalette(brewer.pal(9, "Blues"))(255)

df_anotacion <- as.data.frame(colData(dds)[, "Condicion", drop = FALSE])

pheatmap(matriz_correlacion,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colores_heatmap,
         annotation_col = df_anotacion,
         annotation_row = df_anotacion,
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "Heatmap de Correlación de Muestras")

install.packages("WGCNA")
library(WGCNA)
