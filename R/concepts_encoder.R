# SAE concept encoder: weight handling, an R reference implementation, and the
# server-side Earth Engine deployment. The SAE is trained once offline
# (dev/02_train_sae.R) and ships as data; users only ever run the encoder as
# image math on Earth Engine.

# alphaearth_rescale() lives in gee_logic.R; every SAE consumer uses it.

#' Band names for the m concept bands
#' @noRd
concept_band_names <- function(m) {
  width <- max(3L, nchar(as.character(m)))
  sprintf(paste0("C%0", width, "d"), seq_len(m))
}

#' Construct and validate an SAE weights object
#'
#' Both matrices are stored latent-major (m x 64): row j of `W_enc` reads
#' latent j and row j of `W_dec` is that latent's decoder direction d_j.
#' Encoding is a = TopK(ReLU(W_enc (x - b_pre) + b_enc)); reconstruction is
#' x_hat = t(W_dec) a + b_pre.
#'
#' @param W_enc Encoder matrix, m x 64.
#' @param b_enc Encoder bias, length m.
#' @param W_dec Decoder matrix, m x 64, rows unit-norm by training.
#' @param b_pre Pre-encoder bias subtracted from the input, length 64.
#' @param k TopK sparsity parameter.
#' @param variant "topk" or "batchtopk".
#' @param metadata Free-form training metadata list.
#' @noRd
new_sae <- function(W_enc, b_enc, W_dec, b_pre, k,
                    variant = "topk", metadata = list()) {
  sae <- list(W_enc = as.matrix(W_enc), b_enc = as.numeric(b_enc),
              W_dec = as.matrix(W_dec), b_pre = as.numeric(b_pre),
              m = nrow(as.matrix(W_enc)), k = as.integer(k),
              variant = variant, metadata = metadata)
  validate_sae(sae)
}

#' @noRd
validate_sae <- function(sae) {
  need <- c("W_enc", "b_enc", "W_dec", "b_pre", "m", "k")
  missing <- setdiff(need, names(sae))
  if (length(missing) > 0)
    stop("SAE weights object is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  # Strip dimnames inherited from training data: reticulate turns a named
  # vector's as.list() into a Python dict, which ee.Array rejects.
  sae$W_enc <- unname(as.matrix(sae$W_enc)); sae$W_dec <- unname(as.matrix(sae$W_dec))
  sae$b_enc <- unname(as.numeric(sae$b_enc)); sae$b_pre <- unname(as.numeric(sae$b_pre))
  m <- sae$m
  if (!identical(dim(sae$W_enc), c(m, 64L)) && !identical(dim(sae$W_enc), as.integer(c(m, 64))))
    stop("W_enc must be m x 64; got ", paste(dim(sae$W_enc), collapse = " x "), call. = FALSE)
  if (!identical(dim(sae$W_dec), dim(sae$W_enc)))
    stop("W_dec must be m x 64; got ", paste(dim(sae$W_dec), collapse = " x "), call. = FALSE)
  if (length(sae$b_enc) != m) stop("b_enc must have length m = ", m, call. = FALSE)
  if (length(sae$b_pre) != 64L) stop("b_pre must have length 64.", call. = FALSE)
  if (sae$k < 1L || sae$k > m) stop("k must be in [1, m].", call. = FALSE)
  vals <- c(sae$W_enc, sae$b_enc, sae$W_dec, sae$b_pre)
  if (!all(is.finite(vals))) stop("SAE weights contain non-finite values.", call. = FALSE)
  sae
}

#' Load the shipped SAE weights
#'
#' @param path Path to an .rds written by dev/02_train_sae.R; NULL for the
#'   weights installed with the package.
#' @noRd
load_sae_weights <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "sae_weights_v1.rds", package = "AlphaSDM")
    if (!nzchar(path)) {
      stop("No SAE weights are installed with this copy of AlphaSDM. ",
           "The concept layer is experimental: train weights with dev/02_train_sae.R ",
           "and place them in inst/extdata/sae_weights_v1.rds, or pass `sae = `.",
           call. = FALSE)
    }
  }
  validate_sae(readRDS(path))
}

#' R reference encoder
#'
#' The ground truth the Earth Engine deployment and the torch training script
#' are both tested against. TopK keeps values >= the k-th largest, so ties at
#' the threshold all survive; the GEE helper has the same semantics.
#'
#' @param X Numeric matrix, n x 64, in the shared value convention.
#' @param sae SAE weights object.
#' @return n x m activation matrix, columns named as the concept bands.
#' @noRd
sae_encode <- function(X, sae) {
  X <- alphaearth_rescale(as.matrix(X))
  A <- sweep(X, 2L, sae$b_pre) %*% t(sae$W_enc)
  A <- sweep(A, 2L, sae$b_enc, `+`)
  A[A < 0] <- 0
  if (sae$k < sae$m) {
    for (i in seq_len(nrow(A))) {
      th <- sort(A[i, ], decreasing = TRUE)[sae$k]
      A[i, A[i, ] < th] <- 0
    }
  }
  colnames(A) <- concept_band_names(sae$m)
  A
}

#' R reference reconstruction
#' @param A n x m activation matrix from [sae_encode()].
#' @return n x 64 reconstructed embeddings.
#' @noRd
sae_reconstruct <- function(A, sae) {
  sweep(as.matrix(A) %*% sae$W_dec, 2L, sae$b_pre, `+`)
}

#' R reference ablation of one concept
#'
#' x - a_j d_j: the subtraction form, so reconstruction error never enters the
#' delta. This is the operation the GEE ablation path applies as image math.
#' @noRd
sae_ablate <- function(X, sae, j) {
  A <- sae_encode(X, sae)
  as.matrix(X) - A[, j, drop = FALSE] %*% sae$W_dec[j, , drop = FALSE]
}

#' Build an ee.Array from an R matrix
#' @noRd
ee_array_from_matrix <- function(M) {
  ee <- reticulate::import("ee")
  M <- unname(as.matrix(M))
  ee$Array(lapply(seq_len(nrow(M)), function(i) as.list(M[i, ])))
}

#' Zero every activation below the k-th largest, per pixel
#'
#' `acts` is a 1-D array image of length m. Sorting a copy ascending and
#' reading element m - k gives the k-th largest as a scalar image, which
#' broadcasts across the array in the comparison. Ties at the threshold all
#' survive, matching [sae_encode()].
#' @noRd
ee_topk_zero <- function(acts, m, k) {
  if (k >= m) return(acts)
  thresh <- acts$arraySort()$arrayGet(list(as.integer(m - k)))
  acts$multiply(acts$gte(thresh))
}

#' Concept activations as a server-side Earth Engine image
#'
#' Applies the SAE encoder to a 64-band embedding image as array-image math:
#' matrix multiply, bias, ReLU, then per-pixel TopK. Lazy; nothing is computed
#' until the result is sampled or exported.
#'
#' @param image 64-band embedding image (bands A00..A63).
#' @param sae SAE weights object.
#' @param select_concepts Optional concept band names to return. The full-m
#'   matrix multiply still runs, because the TopK threshold depends on every
#'   latent; selecting only trims what later stages carry.
#' @return m-band (or selected) image of concept activations.
#' @noRd
ee_concept_activations <- function(image, sae, select_concepts = NULL) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)
  bands <- concept_band_names(sae$m)

  centered <- alphaearth_rescale(image$select(emb_cols))$
    subtract(ee$Image$constant(as.list(unname(sae$b_pre)))$rename(emb_cols))
  arr <- centered$toArray()$toArray(1L)                          # 64 x 1

  acts <- ee$Image(ee_array_from_matrix(sae$W_enc))$
    matrixMultiply(arr)$
    add(ee$Image(ee$Array(lapply(as.list(unname(sae$b_enc)), list))))$   # m x 1
    max(ee$Image$constant(0))$
    arrayProject(list(0L))                                       # length-m vector

  acts <- ee_topk_zero(acts, sae$m, sae$k)
  out <- acts$arrayFlatten(list(as.list(bands)))
  if (!is.null(select_concepts)) out <- out$select(as.list(select_concepts))
  out
}
