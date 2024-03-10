#' Find marker genes based on Seurat
#'
#' @param counts a matrix indicating the raw count expression of cells (gene x cell).
#' @param labels a vector indicating the corresponding cell type labels.
#' @param n_markers an integer indicating the number of marker genes to select for each cell type. The default is \code{NULL} and all marker genes will be output.
#'
#' @return a list containing marker genes for each cell type.
#' @export
#' @import Seurat
find_markers <- function(counts, labels, n_markers = NULL) {
  dat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 0)
  dat <- NormalizeData(dat)
  dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(dat)
  dat <- ScaleData(dat, features = all.genes)
  Idents(dat) <- labels
  markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  unique_labels <- unique(labels)
  markers_list <- lapply(unique_labels, function(x) {
    temp <- markers[markers$cluster == x, ]
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp
  })
  names(markers_list) <- unique_labels
  if (!is.null(n_markers)) {
    marker_genes <- lapply(markers_list, function(x) {
      temp <- x$gene
      n_markers <- min(n_markers, length(temp))
      temp[1:n_markers]
    })
  } else {
    marker_genes <- lapply(markers_list, function(x) {
      x$gene
    })
  }
  return(marker_genes)
}



#' Max-min sampling
#'
#' @param data a matrix of samples (sample x feature).
#' @param n_points an integer indicating the number of samples to take.
#'
#' @return a vector of the indices of the samples that have been sampled.
#' @export
max_min_sampling <- function(data, n_points = 5) {
  n_data <- nrow(data)
  data_scale <- scale(data[, colSums(data) != 0], center = TRUE, scale = TRUE)
  data_dist <- cosine_distance(data_scale)
  all_points <- which.max(rowSums(data_scale^2))
  if (n_points > 1) {
    for (i in 2:n_points) {
      all_points <- c(all_points, which.max(apply(as.matrix(data_dist[, all_points]), 1, min)))
    }
  }
  return(all_points)
}



#' Calculate the cosine distance
#'
#' @param data a matrix of samples (sample x feature).
#'
#' @return a matrix of distances between samples (sample x sample).
#' @export
cosine_distance <- function(data) {
  data <- as.matrix(data)
  temp <- data / sqrt(rowSums(data^2))
  res <- 1 - temp %*% t(temp)
  diag(res) <- 0
  return(res)
}



#' Calculate the cosine similarity
#'
#' @param data_1 a matrix of the first set of samples (sample x feature).
#' @param data_2 a matrix of the second set of samples (sample x feature).
#'
#' @return a matrix of similarity between two sets of samples.
#' @export
cosine_similarity <- function(data_1, data_2) {
  data_1 <- as.matrix(data_1)
  data_2 <- as.matrix(data_2)
  data <- rbind(data_1, data_2)
  data_scale <- scale(data[, colSums(data) != 0], center = TRUE, scale = TRUE)
  data_1 <- data_scale[1:nrow(data_1), ]
  data_2 <- data_scale[(nrow(data_1) + 1):nrow(data_scale), ]
  data_1 <- data_1 / sqrt(rowSums(data_1^2))
  data_2 <- data_2 / sqrt(rowSums(data_2^2))
  res <- data_1 %*% t(data_2)
  return(res)
}
