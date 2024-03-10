#' Archetypal analysis on scRNA-seq reference
#'
#' @param sc_counts a matrix indicating the raw count expression of cells (cell x gene).
#' @param sc_labels a vector indicating the corresponding cell type labels.
#' @param n_archetypes_vec either an integer or a vector indicating the number of archetypes for each cell type. If it is an integer, the number of archetypes is the same for all cell types and is determined by this integer. If specifying a different number of archetypes for different cell types, set \code{n_archetypes_vec} to a vector with a length equal to the number of cell types and name each element according to the respective cell type. The default is 10.
#' @param tol a value indicating the convergence threshold. The default is 1e-8.
#' @param maxiter an integer indicating the maximum number of iterations. The default is 2000.
#' @param do_initialize a logical value indicating whether to perform initialization. The default is \code{TRUE}.
#' @param do_plot a logical value indicating whether to visualize the archetypes using UMAP (Uniform Manifold Approximation and Projection). The default is \code{TRUE}.
#' @param n_hvgs an integer indicating the number of highly variable genes used for archetypal analysis. The default is 2000.
#' @param n_markers an integer indicating the number of marker genes used for archetypal analysis for each cell type. The default is 20.
#' @param n_sample_cells an integer indicating the maximum number of cells contained in each cell type. Cell types with a cell count exceeding \code{n_sample_cells} will be downsampled to \code{n_sample_cells}. The default is 500.
#' @param do_parallel a logical value indicating whether to use parallel computation. The default is \code{TRUE}.
#' @param n_cores an integer indicating the number of cores used for parallel computation. The default is 20.
#' @param save_res a logical value indicating whether to save the results. The default is \code{TRUE}.
#' @param save_dir a character string indicating the directory where the results are saved. The default is the current directory.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{archetypes_all_list}}{a list containing the archetypes with all genes for each cell type.}
#'   \item{\code{archetypes_list}}{a list containing the archetypes with selected genes (highly variable genes and marker genes) for each cell type.}
#'   \item{\code{res_list}}{a list containing all results of archetypal analysis for each cell type.}
#' }
#' @export
#' @import parallel
#' @import foreach
#' @import doParallel
sc_train <- function(sc_counts, sc_labels, n_archetypes_vec = 10, tol = 1e-8, maxiter = 2000, do_initialize = TRUE, do_plot = TRUE, n_hvgs = 2000, n_markers = 20, n_sample_cells = 500, do_parallel = TRUE, n_cores = 20, save_res = TRUE, save_dir = ".") {

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  sc_data <- sc_preprocess(sc_counts, sc_labels, n_hvgs = n_hvgs, n_markers = n_markers, n_sample_cells = n_sample_cells)
  X_list <- sc_data$X_list
  Z_list <- sc_data$Z_list
  Z_all_list <- sc_data$Z_all_list
  N_list <- sc_data$N_list
  unique_labels <- sc_data$unique_labels

  if (length(n_archetypes_vec) == 1) {
    n_archetypes_vec <- rep(n_archetypes_vec, length(unique_labels))
    names(n_archetypes_vec) <- unique_labels
  }

  if (do_parallel) {
    n_cores <- min(n_cores, length(unique_labels), detectCores() - 1)

    cl <- makeCluster(n_cores, outfile = "")
    registerDoParallel(cl)

    res_list <- foreach(i = unique_labels, .export = c("sc_model", "sc_initialize", "max_min_sampling", "cosine_distance", "cosine_similarity", "plot_archetypes")) %dopar% {
      cat(i, "\t", "start", "\n")
      res_i <- sc_model(X = X_list[[i]], Z = Z_list[[i]], N = N_list[[i]], n_archetypes = as.numeric(n_archetypes_vec[i]), tol = tol, maxiter = maxiter, do_initialize = do_initialize, do_plot = do_plot, f_name = paste("archetypes_umap", i, sep = "_"), save_dir = paste(save_dir, "archetypes_umap", sep = "/"))
      rownames(res_i$phi) <- paste(i, seq(dim(res_i$phi)[1]), sep = "_")
      colnames(res_i$alpha) <- paste(i, seq(dim(res_i$alpha)[2]), sep = "_")
      rownames(res_i$W) <- paste(i, seq(dim(res_i$W)[1]), sep = "_")
      res_i$phi_all <- res_i$W %*% Z_all_list[[i]]
      cat(i, "\t", "end", "\n")
      res_i
    }

    stopImplicitCluster()
    stopCluster(cl)

    names(res_list) <- unique_labels
    archetypes_list <- lapply(res_list, function(x) {
      x$phi
    })
    archetypes_all_list <- lapply(res_list, function(x) {
      x$phi_all
    })

  } else {

    archetypes_list <- list()
    archetypes_all_list <- list()
    res_list <- list()

    for (i in unique_labels) {
      cat(i, "\t", "start", "\n")
      res_i <- sc_model(X = X_list[[i]], Z = Z_list[[i]], N = N_list[[i]], n_archetypes = as.numeric(n_archetypes_vec[i]), tol = tol, maxiter = maxiter, do_initialize = do_initialize, do_plot = do_plot, f_name = paste("archetypes_umap", i, sep = "_"), save_dir = paste(save_dir, "archetypes_umap", sep = "/"))
      rownames(res_i$phi) <- paste(i, seq(dim(res_i$phi)[1]), sep = "_")
      colnames(res_i$alpha) <- paste(i, seq(dim(res_i$alpha)[2]), sep = "_")
      rownames(res_i$W) <- paste(i, seq(dim(res_i$W)[1]), sep = "_")
      res_i$phi_all <- res_i$W %*% Z_all_list[[i]]
      archetypes_list[[i]] <- res_i$phi
      archetypes_all_list[[i]] <- res_i$phi_all
      res_list[[i]] <- res_i
      cat(i, "\t", "end", "\n")
    }
  }

  sc_results <- list(archetypes_all_list = archetypes_all_list, archetypes_list = archetypes_list, res_list = res_list)

  if (save_res) {
    saveRDS(sc_results, file = paste(save_dir, "sc_results.rds", sep = "/"))
  }

  return(sc_results)
}



#' Deconvolution of spatial transcriptomics data
#'
#' @param st_counts a matrix indicating the raw count expression of spots (spot x gene).
#' @param sc_results a list indicating the output of the \code{\link{sc_train}} function.
#' @param n_hvgs an integer indicating the number of highly variable genes selected from the intersection of genes in the scRNA-seq reference and the spatial transcriptomics data. The default is 5000.
#' @param sigma a value indicating the parameter that determine the standard deviation of the normal distribution for the platform effect. The default is 0.3.
#' @param tol a value indicating the convergence threshold. The default is 1e-8.
#' @param maxiter an integer indicating the maximum number of iterations. The default is 2000.
#' @param save_res a logical value indicating whether to save the results. The default is \code{TRUE}.
#' @param save_dir a character string indicating the directory where the results are saved. The default is the current directory.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{proportion}}{a matrix of the cell type proportions for spots (spot x cell type).}
#'   \item{\code{beta}}{a matrix of the cell type abundance for spots (spot x cell type).}
#'   \item{\code{alpha}}{a matrix of the combination coefficients for spots (spot x archetype).}
#'   \item{\code{omega}}{a matrix of the auxiliary matrix (spot x archetype).}
#'   \item{\code{gamma}}{a vector of the platform effect for genes.}
#'   \item{\code{mu}}{a list containing the cell-type-specific gene expression for spots.}
#'   \item{\code{Q}}{a one-hot matrix indicating the cell type corresponding to each archetype (cell type x archetype).}
#'   \item{\code{unique_labels}}{a vector of unique cell types included in the scRNA-seq reference.}
#'   \item{\code{sigma}}{a value of the parameter that determine the standard deviation of the normal distribution for the platform effect.}
#'   \item{\code{phi}}{a matrix of the expression of archetypes for all cell types.}
#'   \item{\code{phi_all}}{a matrix of the expression of archetypes with all genes for all cell types.}
#' }
#' @export
#' @importFrom utils write.csv
st_train <- function(st_counts, sc_results, n_hvgs = 5000, sigma = 0.3, tol = 1e-8, maxiter = 2000, save_res = TRUE, save_dir = ".") {

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  archetypes_list <- sc_results$archetypes_list
  archetypes_all_list <- sc_results$archetypes_all_list

  keep_genes_sc <- intersect(Reduce(intersect, lapply(archetypes_list, function(x) {
    colnames(x)
  })), colnames(st_counts))
  keep_genes_st <- intersect(Reduce(intersect, lapply(archetypes_all_list, function(x) {
    colnames(x)
  })), colnames(st_counts))
  st_data <- st_preprocess(st_counts, n_hvgs = n_hvgs, keep_genes = keep_genes_st)
  keep_genes_st <- colnames(st_data$Y)
  keep_genes <- intersect(keep_genes_st, keep_genes_sc)
  Y <- st_data$Y[, keep_genes]
  N <- st_data$N

  unique_labels <- names(archetypes_all_list)
  archetypes_list <- lapply(archetypes_all_list, function(x) {
    x[, keep_genes]
  })
  phi <- Reduce(rbind, archetypes_list)
  phi_all <- Reduce(rbind, archetypes_all_list)
  Q <- (as.matrix(table(rep(unique_labels, times = sapply(archetypes_list, function(x) {
    dim(x)[1]
  })), rownames(phi))) > 0) * 1
  Q <- Q[unique_labels, rownames(phi)]

  st_results <- st_model(Y = Y, N = N, phi = phi, Q = Q, sigma = sigma, tol = tol, maxiter = maxiter)
  st_results$mu <- lapply(unique_labels, function(x) {
    st_results$alpha %*% diag(Q[x, ]) %*% phi_all
  })
  names(st_results$mu) <- unique_labels

  st_results$Q <- Q
  st_results$unique_labels <- unique_labels
  st_results$sigma <- sigma
  st_results$phi <- phi
  st_results$phi_all <- phi_all

  if (save_res) {
    saveRDS(st_results, file = paste(save_dir, "st_results.rds", sep = "/"))
    write.csv(as.matrix(st_results$proportion), file = paste(save_dir, "PANDA_prop.csv", sep = "/"))
    write.csv(as.matrix(st_results$alpha), file = paste(save_dir, "PANDA_alpha.csv", sep = "/"))
    for (i in names(st_results$mu)) {
      write.csv(as.matrix(st_results$mu[[i]]), file = paste(save_dir, paste(c("PANDA", i, "mu.csv"), collapse = "_"), sep = "/"))
    }
  }

  return(st_results)
}



#' Find archetypes for one cell type
#'
#' @param X a matrix indicating the raw count expression of cells (cell x gene).
#' @param Z a matrix indicating the normalized expression of cells (cell x gene).
#' @param N a vector indicating the library size of each cell.
#' @param n_archetypes an integer indicating the number of archetypes.
#' @param tol a value indicating the convergence threshold. The default is 1e-8.
#' @param maxiter an integer indicating the maximum number of iterations. The default is 2000.
#' @param do_initialize a logical value indicating whether to perform initialization. The default is \code{TRUE}.
#' @param do_plot a logical value indicating whether to visualize the archetypes using UMAP (Uniform Manifold Approximation and Projection). The default is \code{TRUE}.
#' @param f_name a character string naming the pdf file of the UMAP. The default is "archetypes".
#' @param save_dir a character string indicating the directory where the pdf file is saved. The default is the current directory.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{phi}}{a matrix of the expression of archetypes (archetype x gene).}
#'   \item{\code{alpha}}{a matrix of combination coefficients for cells (cell x archetype).}
#'   \item{\code{W}}{a matrix of combination coefficients for archetypes (archetype x cell).}
#' }
#' @export
#' @importFrom stats var runif
sc_model <- function(X, Z, N, n_archetypes, tol = 1e-8, maxiter = 2000, do_initialize = TRUE, do_plot = TRUE, f_name = "archetypes", save_dir = ".") {

  eps <- 1e-12

  n_cells <- dim(X)[1]
  n_genes <- dim(X)[2]
  n_archetypes <- min(n_archetypes, n_cells)

  tau <- var(c(X)) * 20 / n_genes

  if (do_initialize) {
    res_init <- sc_initialize(Z, n_archetypes)
    alpha <- res_init$alpha
    W <- res_init$W
  } else {
    alpha <- matrix(runif(n_cells * n_archetypes), nrow = n_cells, ncol = n_archetypes)
    W <- matrix(runif(n_archetypes * n_cells), nrow = n_archetypes, ncol = n_cells)
    alpha <- alpha / rowSums(alpha)
    W <- W / rowSums(W)
    res_init <- NULL
  }


  loss_func <- function(X, Z, N, alpha, W, tau) {
    n_cells <- dim(X)[1]
    n_genes <- dim(X)[2]
    n_archetypes <- dim(alpha)[2]
    temp <- alpha %*% W %*% Z
    loss_1 <- (1 / (n_cells * n_genes)) * (-sum(diag(t(X) %*% log((N %*% t(rep(1, n_genes))) * temp + eps))) +
      t(N) %*% temp %*% rep(1, n_genes))
    loss_2 <- (1 / n_cells) * (-t(rep(1, n_cells)) %*% log(alpha %*% rep(1, n_archetypes) + eps) + t(rep(1, n_cells)) %*% alpha %*% rep(1, n_archetypes)) +
      (1 / n_archetypes) * (-t(rep(1, n_archetypes)) %*% log(W %*% rep(1, n_cells) + eps) + t(rep(1, n_archetypes)) %*% W %*% rep(1, n_cells))
    loss <- loss_1 + tau * loss_2
    return(list(loss = as.numeric(loss), loss_1 = as.numeric(loss_1), loss_2 = as.numeric(loss_2)))
  }

  for (iter in 1:maxiter) {
    temp <- W %*% Z
    alpha <- alpha * (1 / ((1 / (n_cells * n_genes)) * N %*% t(rep(1, n_genes)) %*% t(temp) + tau / n_cells + eps)) *
      ((1 / (n_cells * n_genes)) * (X / (alpha %*% temp + eps)) %*% t(temp) + (tau / n_cells) * (1 / (alpha %*% matrix(1, n_archetypes, n_archetypes) + eps)))
    W <- W * (1 / ((1 / (n_cells * n_genes)) * t(alpha) %*% N %*% t(rep(1, n_genes)) %*% t(Z) + tau / n_archetypes + eps)) *
      ((1 / (n_cells * n_genes)) * t(alpha) %*% (X / (alpha %*% temp + eps)) %*% t(Z) + (tau / n_archetypes) * (1 / (W %*% matrix(1, n_cells, n_cells) + eps)))
    if (iter != 1) {
      loss_res <- loss_func(X, Z, N, alpha, W, tau)
      loss <- loss_res$loss
      # cat(iter, "\t", "loss = ", loss, "\t", "loss_old = ", loss_old, "\t", "loss_old - loss =", loss_old - loss, "\t", "abs(loss - loss_old) / abs(loss_old) = ", abs(loss - loss_old) / abs(loss_old), "\t", "loss_1 = ", loss_res$loss_1, "\t", "loss_2 = ", loss_res$loss_2, "\n")
      if (abs(loss - loss_old) / abs(loss_old) < tol) {
        break
      }
      loss_old <- loss
    } else {
      loss_old <- loss_func(X, Z, N, alpha, W, tau)$loss
    }
  }

  phi <- W %*% Z

  if (do_plot) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    plot_archetypes(Z, phi, f_name = f_name, save_dir = save_dir, width = 4, height = 4)
  }

  return(list(phi = phi, alpha = alpha, W = W))
}



#' Initialize the combination coefficients for cells and archetypes
#'
#' @param Z a matrix indicating the normalized expression of cells (cell x gene).
#' @param n_archetypes an integer indicating the number of archetypes.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{alpha}}{a matrix of the initial combination coefficients for cells (cell x archetype).}
#'   \item{\code{W}}{a matrix of the initial combination coefficients for archetypes (archetype x cell).}
#' }
#' @export
sc_initialize <- function(Z, n_archetypes) {

  eps <- 1e-12

  n_cells <- dim(Z)[1]

  idx_sample <- max_min_sampling(Z, n_points = n_archetypes)

  W <- cosine_similarity(Z[idx_sample, ], Z)
  W[W < eps] <- eps
  W <- W / rowSums(W)

  alpha <- cosine_similarity(Z, W %*% Z)
  alpha[alpha < eps] <- eps
  alpha <- alpha / rowSums(alpha)

  return(list(alpha = alpha, W = W))
}



#' The core function for deconvolution
#'
#' @param Y a matrix indicating the raw count expression of spots (spot x gene).
#' @param N a vector indicating the total counts of each spot.
#' @param phi a matrix indicating the expression of archetypes (archetype x gene)
#' @param Q a one-hot matrix indicating the cell type corresponding to each archetype (cell type x archetype).
#' @param sigma a value indicating the parameter that determine the standard deviation of the normal distribution for the platform effect. The default is 0.3.
#' @param tol a value indicating the convergence threshold. The default is 1e-8.
#' @param maxiter an integer indicating the maximum number of iterations. The default is 2000.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{proportion}}{a matrix of the cell type proportions for spots (spot x cell type).}
#'   \item{\code{beta}}{a matrix of the cell type abundance for spots (spot x cell type).}
#'   \item{\code{alpha}}{a matrix of the combination coefficients for spots (spot x archetype).}
#'   \item{\code{omega}}{a matrix of the auxiliary matrix (spot x archetype).}
#'   \item{\code{gamma}}{a vector of the platform effect for genes.}
#' }
#' @export
st_model <- function(Y, N, phi, Q, sigma = 0.3, tol = 1e-8, maxiter = 2000) {

  eps <- 1e-12

  n_spots <- dim(Y)[1]
  n_genes <- dim(Y)[2]
  n_archetypes <- dim(phi)[1]

  sigma <- sigma / sqrt(n_spots)

  loss_func <- function(Y, N, phi, omega, gamma, sigma) {
    n_spots <- dim(Y)[1]
    n_genes <- dim(Y)[2]
    temp <- omega %*% phi
    loss <- -sum(diag(t(Y) %*% log((N %*% t(rep(1, n_genes))) * temp * (rep(1, n_spots) %*% t(exp(gamma))) + eps))) +
      t(N) %*% temp %*% exp(gamma) +
      (1 / (2 * (sigma^2))) * (t(rep(1, n_genes)) %*% (gamma^2))
    return(as.numeric(loss))
  }

  gamma <- rep(0, n_genes)
  omega <- cosine_similarity(Y / N, phi)
  omega[omega < eps] <- eps
  omega <- omega / rowSums(omega)

  for (iter in 1:maxiter) {
    omega <- omega * (1 / (N %*% t(phi %*% exp(gamma)) + eps)) * ((Y / (omega %*% phi + eps)) %*% t(phi))
    if (iter != 1) {
      loss <- loss_func(Y, N, phi, omega, gamma, sigma)
      # cat(iter, "\t", "loss = ", loss, "\t", "loss_old = ", loss_old, "\t", "loss_old - loss =", loss_old - loss, "\t", "abs(loss - loss_old) / abs(loss_old) =", abs(loss - loss_old) / abs(loss_old), "\n")
      if (abs(loss - loss_old) / abs(loss_old) < tol) {
        break
      }
      loss_old <- loss
    } else {
      loss_old <- loss_func(Y, N, phi, omega, gamma, sigma)
    }
  }

  initial_diff <- (t(Y) %*% rep(1, n_spots)) / (t(phi) %*% t(omega) %*% N + eps)
  keep_genes <- which((initial_diff < 2.5) & (initial_diff > 0.4))
  Y <- Y[, keep_genes]
  phi <- phi[, keep_genes]
  gamma <- rep(0, length(keep_genes))

  for (iter in 1:maxiter) {
    omega <- omega * (1 / (N %*% t(phi %*% exp(gamma)) + eps)) * ((Y / (omega %*% phi + eps)) %*% t(phi))
    temp <- exp(gamma) * (t(phi) %*% t(omega) %*% N)
    gamma <- gamma - (-t(Y) %*% rep(1, n_spots) + temp + gamma / (sigma^2)) / (temp + 1 / (sigma^2))
    if (iter != 1) {
      loss <- loss_func(Y, N, phi, omega, gamma, sigma)
      # cat(iter, "\t", "loss = ", loss, "\t", "loss_old = ", loss_old, "\t", "loss_old - loss =", loss_old - loss, "\t", "abs(loss - loss_old) / abs(loss_old) =", abs(loss - loss_old) / abs(loss_old), "\n")
      if (abs(loss - loss_old) / abs(loss_old) < tol) {
        break
      }
      loss_old <- loss
    } else {
      loss_old <- loss_func(Y, N, phi, omega, gamma, sigma)
    }
  }

  beta <- omega %*% t(Q)
  alpha <- omega / (beta %*% Q + eps)

  proportion <- beta / rowSums(beta)

  rownames(proportion) <- rownames(Y)
  rownames(beta) <- rownames(Y)
  rownames(alpha) <- rownames(Y)
  rownames(omega) <- rownames(Y)

  return(list(proportion = proportion, beta = beta, alpha = alpha, omega = omega, gamma = gamma))
}



#' Preprocess the scRNA-seq reference
#'
#' @param sc_counts a matrix indicating the raw count expression of cells (cell x gene).
#' @param sc_labels a vector indicating the corresponding cell type labels.
#' @param n_hvgs an integer indicating the number of highly variable genes to select. The default is 2000.
#' @param n_markers an integer indicating the number of marker genes to select for each cell type. The default is 20.
#' @param n_sample_cells an integer indicating the maximum number of cells contained in each cell type. Cell types with a cell count exceeding \code{n_sample_cells} will be downsampled to \code{n_sample_cells}. The default is 500.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{X_list}}{a list containing the preprocessed count matrix for each cell type.}
#'   \item{\code{Z_list}}{a list containing the preprocessed normalized expression matrix for each cell type.}
#'   \item{\code{Z_all_list}}{a list containing the preprocessed normalized expression matrix with all genes for each cell type.}
#'   \item{\code{N_list}}{a list containing the library size vector for each cell type.}
#'   \item{\code{unique_labels}}{a vector of unique cell types included in the scRNA-seq reference.}
#' }
#' @export
#' @import Seurat
sc_preprocess <- function(sc_counts, sc_labels, n_hvgs = 2000, n_markers = 20, n_sample_cells = 500) {

  sc_counts <- as.matrix(sc_counts)

  N <- rowSums(sc_counts)
  idx_keep <- which(N >= 200)
  N <- N[idx_keep]
  sc_counts <- sc_counts[idx_keep, ]
  sc_labels <- sc_labels[idx_keep]

  Z <- sc_counts / N
  X <- sc_counts
  Z_all <- Z

  n_hvgs <- min(n_hvgs, dim(sc_counts)[2])
  dat <- CreateSeuratObject(counts = t(sc_counts), min.cells = 3, min.features = 200)
  dat <- NormalizeData(dat)
  dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = n_hvgs)
  hvgs <- VariableFeatures(dat)

  markers <- find_markers(t(sc_counts), sc_labels, n_markers = n_markers)
  markers <- Reduce("union", markers)

  keep_genes <- union(hvgs, markers)

  X <- X[, keep_genes]
  Z <- Z[, keep_genes]

  unique_labels <- unique(sc_labels)
  X_list <- list()
  Z_list <- list()
  Z_all_list <- list()
  N_list <- list()

  for (i in unique_labels) {
    idx <- which(sc_labels == i)
    if (length(idx) > n_sample_cells) {
      idx_keep <- max_min_sampling(Z[idx, ], n_points = n_sample_cells)
      idx <- idx[idx_keep]
    }
    X_list[[i]] <- X[idx, ]
    Z_list[[i]] <- Z[idx, ]
    Z_all_list[[i]] <- Z_all[idx, ]
    N_list[[i]] <- N[idx]
  }

  return(list(X_list = X_list, Z_list = Z_list, Z_all_list = Z_all_list, N_list = N_list, unique_labels = unique_labels))
}



#' Preprocess the spatial transcriptomics data
#'
#' @param st_counts a matrix indicating the raw count expression of spots (spot x gene).
#' @param n_hvgs an integer indicating the number of highly variable genes to select. The default is 5000.
#' @param keep_genes a vector of genes from which the highly variable genes will be selected. If \code{NULL}, the highly variable genes will be selected from all genes in the \code{st_counts}.
#'
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{Y}}{a matrix of the preprocessed count matrix (spot x gene).}
#'   \item{\code{N}}{a vector of total counts for spot.}
#' }
#' @export
#' @import Seurat
st_preprocess <- function(st_counts, n_hvgs = 5000, keep_genes = NULL) {

  st_counts <- as.matrix(st_counts)

  N <- rowSums(st_counts)
  Y <- st_counts

  if (is.null(keep_genes)) {
    keep_genes <- colnames(st_counts)
  }

  dat <- CreateSeuratObject(counts = t(st_counts[, keep_genes]), min.cells = 0, min.features = 0)
  dat <- NormalizeData(dat)
  dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = n_hvgs)
  hvgs <- VariableFeatures(dat)

  Y <- Y[, hvgs]

  return(list(Y = Y, N = N))
}
