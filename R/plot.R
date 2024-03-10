#' Visualization of archetypes
#'
#' @param Z a matrix indicating the normalized expression of cells (cell x gene).
#' @param phi a matrix indicating the expression of archetypes (archetype x gene).
#' @param f_name a character string naming the pdf file of the figure. The default is "archetypes".
#' @param save_dir a character string indicating the directory where the pdf file is saved. The default is the current directory.
#' @param width a value indicating the width of the figure. The default is 7.
#' @param height a value indicating the height of the figure. The default is 7.
#'
#' @return a pdf file containing the figure.
#' @export
#' @importFrom stats prcomp predict
#' @importFrom grDevices pdf dev.off
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @import umap
plot_archetypes <- function(Z, phi, f_name = "archetypes", save_dir = ".", width = 7, height = 7) {
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  gene_intersect <- intersect(colnames(Z), colnames(phi))
  Z <- Z[, gene_intersect]
  phi <- phi[, gene_intersect]

  idx_gene <- which(colSums(Z) != 0)
  pca_res <- prcomp(Z[, idx_gene], center = TRUE, scale. = TRUE)
  ndims <- min(30, dim(Z)[1])
  Z_pca <- as.data.frame(pca_res$x[, 1:ndims])
  umap_settings <- umap.defaults
  umap_settings$n_neighbors <- min(15, dim(Z)[1])
  umap_res <- umap(as.matrix(Z_pca), config = umap_settings)
  Z_umap <- as.data.frame(umap_res$layout)
  colnames(Z_umap) <- c("dim1", "dim2")
  phi_pca <- as.data.frame(predict(pca_res, phi[, idx_gene])[, 1:ndims])
  phi_umap <- as.data.frame(predict(umap_res, phi_pca))
  phi_umap <- cbind(phi_umap, 1:nrow(phi))
  colnames(phi_umap) <- c("dim1", "dim2", "name")

  pdf(paste(save_dir, paste(f_name, "pdf", sep = "."), sep = "/"), width = width, height = height)
  fig <- ggplot() +
    geom_point(aes(x = dim1, y = dim2), color = "#70B2DE", shape = 16, size = 1, data = Z_umap) +
    geom_point(aes(x = dim1, y = dim2), color = "#F0525F", shape = 17, size = 2, data = phi_umap) +
    geom_text_repel(aes(x = dim1, y = dim2, label = name), data = phi_umap, fontface = "bold", size = 2) +
    labs(title = f_name, x = "UMAP_1", y = "UMAP_2") +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )
  print(fig)
  dev.off()
}
