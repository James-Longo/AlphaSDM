# Phase 2: train the SAE offline and freeze the shipped weights. Run once, by
# hand, after downloading the Phase 1 sample. Needs the 'torch' package
# (install.packages("torch"); torch::install_torch()).
#
# Trains TopK and BatchTopK SAEs over the (m, k) grid, plus the mandatory
# NMF and k-means dictionary baselines on the same sample. Model selection is
# NOT made here: every candidate's diagnostics are written out, Phase 4 scores
# label coherence for each, and only then is one config frozen as
# inst/extdata/sae_weights_v1.rds (subject to dead latents < 20%,
# tie-broken by reconstruction R^2).

library(AlphaSDM)
library(torch)

# ---- Config ------------------------------------------------------------------
SAMPLE_DIR  <- "~/alphasdm_sae_sample_v1"   # downloaded Phase 1 CSVs
OUT_DIR     <- "dev/artifacts/sae_candidates"
M_GRID      <- c(512L, 1024L, 2048L)
K_GRID      <- c(8L, 16L, 32L)
VAL_FRAC    <- 0.05                          # held-out fraction for diagnostics
BATCH_SIZE  <- 4096L
EPOCHS      <- 20L
LR          <- 1e-3
SEED        <- 42L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED); torch_manual_seed(SEED)
emb_cols <- sprintf("A%02d", 0:63)

# ---- Load the sample ---------------------------------------------------------
files <- list.files(path.expand(SAMPLE_DIR), pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
sample_df <- do.call(rbind, lapply(files, function(f)
  utils::read.csv(f)[, c(emb_cols, "wc_class", "longitude", "latitude")]))
sample_df <- sample_df[stats::complete.cases(sample_df[, emb_cols]), ]
message(sprintf("%d pixels loaded.", nrow(sample_df)))
message("Per-class counts (sanity checkpoint):")
print(table(sample_df$wc_class))

# Phase 1 exported raw values (see its metadata file); the shared convention
# is applied here too so a future non-identity rescale cannot be missed.
X <- AlphaSDM:::alphaearth_rescale(as.matrix(sample_df[, emb_cols]))

n_val   <- ceiling(VAL_FRAC * nrow(X))
val_idx <- sample.int(nrow(X), n_val)
X_train <- torch_tensor(X[-val_idx, ], dtype = torch_float())
X_val   <- torch_tensor(X[val_idx, ],  dtype = torch_float())

# ---- TopK / BatchTopK SAE ----------------------------------------------------
sae_module <- nn_module(
  initialize = function(m, k, batch_topk = FALSE) {
    self$m <- m; self$k <- k; self$batch_topk <- batch_topk
    self$b_pre <- nn_parameter(torch_zeros(64))
    self$enc   <- nn_linear(64, m, bias = TRUE)
    self$dec   <- nn_linear(m, 64, bias = FALSE)     # untied decoder
    self$normalize_decoder()
  },
  normalize_decoder = function() {
    with_no_grad({
      w <- self$dec$weight                            # 64 x m; columns are d_j
      self$dec$weight$copy_(w / w$norm(dim = 1, keepdim = TRUE)$clamp(min = 1e-8))
    })
  },
  encode = function(x) {
    a <- nnf_relu(self$enc(x - self$b_pre))
    if (self$batch_topk) {
      # BatchTopK: one threshold over the whole batch at k * batch_size kept
      # activations, which lets busy pixels use more latents than quiet ones.
      n_keep <- self$k * a$size(1)
      flat   <- a$flatten()
      thresh <- flat$topk(n_keep)[[1]]$min()
      a * (a >= thresh)
    } else {
      thresh <- a$topk(self$k, dim = 2)[[1]][, self$k]$unsqueeze(2)
      a * (a >= thresh)
    }
  },
  forward = function(x) self$dec(self$encode(x)) + self$b_pre
)

train_sae <- function(m, k, batch_topk = FALSE) {
  model <- sae_module(m, k, batch_topk)
  opt   <- optim_adam(model$parameters, lr = LR)
  n     <- X_train$size(1)
  losses <- numeric(0)
  for (epoch in seq_len(EPOCHS)) {
    perm <- torch_randperm(n) + 1L
    epoch_loss <- 0; nb <- 0
    for (s in seq(1L, n, by = BATCH_SIZE)) {
      idx  <- perm[s:min(s + BATCH_SIZE - 1L, n)]
      xb   <- X_train[idx, ]
      opt$zero_grad()
      loss <- nnf_mse_loss(model(xb), xb)     # TopK carries sparsity; no L1
      loss$backward()
      opt$step()
      model$normalize_decoder()
      epoch_loss <- epoch_loss + loss$item(); nb <- nb + 1
    }
    losses <- c(losses, epoch_loss / nb)
    message(sprintf("  m=%d k=%d %s epoch %d/%d loss %.6f", m, k,
                    if (batch_topk) "batchtopk" else "topk", epoch, EPOCHS,
                    losses[epoch]))
  }
  list(model = model, losses = losses)
}

diagnostics <- function(model) {
  with_no_grad({
    A_val <- model$encode(X_val)
    xhat  <- model$dec(A_val) + model$b_pre
    ss_res <- (X_val - xhat)$pow(2)$sum()$item()
    ss_tot <- (X_val - X_val$mean(dim = 1))$pow(2)$sum()$item()
    act_freq <- as.numeric((A_val > 0)$to(dtype = torch_float())$mean(dim = 1))
    W_dec <- t(as.matrix(model$dec$weight))            # m x 64, unit rows
    cos_sim <- W_dec %*% t(W_dec)
    list(r2 = 1 - ss_res / ss_tot,
         dead_frac = mean(act_freq == 0),
         act_freq = act_freq,
         # feature-splitting check: strongest off-diagonal decoder cosine
         max_offdiag_cos = max(abs(cos_sim[upper.tri(cos_sim)])))
  })
}

as_sae_object <- function(model, variant, meta) {
  AlphaSDM:::new_sae(
    W_enc = as.matrix(model$enc$weight),               # m x 64
    b_enc = as.numeric(model$enc$bias),
    W_dec = t(as.matrix(model$dec$weight)),            # stored latent-major
    b_pre = as.numeric(model$b_pre),
    k = model$k, variant = variant, metadata = meta)
}

# ---- Train the grid ----------------------------------------------------------
diag_rows <- list()
for (m in M_GRID) for (k in K_GRID) for (bt in c(FALSE, TRUE)) {
  variant <- if (bt) "batchtopk" else "topk"
  fit  <- train_sae(m, k, bt)
  d    <- diagnostics(fit$model)
  meta <- list(n_train = X_train$size(1), n_val = n_val, year = 2023,
               stratification = "WorldCover class x continent", seed = SEED,
               loss_curve = fit$losses, dead_frac = d$dead_frac, r2 = d$r2)
  sae  <- as_sae_object(fit$model, variant, meta)
  tag  <- sprintf("%s_m%d_k%d", variant, m, k)
  saveRDS(sae, file.path(OUT_DIR, paste0("sae_", tag, ".rds")))
  saveRDS(d$act_freq, file.path(OUT_DIR, paste0("actfreq_", tag, ".rds")))
  diag_rows[[tag]] <- data.frame(variant = variant, m = m, k = k, r2 = d$r2,
                                 dead_frac = d$dead_frac,
                                 max_offdiag_cos = d$max_offdiag_cos)
}

# ---- Baselines (the null hypotheses the SAE has to beat) ---------------------
X_r <- as.matrix(X_train)

# NMF needs non-negative input; embeddings span [-1, 1], so factorise the
# standard split X = [max(X,0), max(-X,0)] and fold the parts back together.
nmf_dictionary <- function(X, m, iters = 200L) {
  Xs <- cbind(pmax(X, 0), pmax(-X, 0))
  set.seed(SEED)
  W <- matrix(runif(nrow(Xs) * m), nrow(Xs), m)
  H <- matrix(runif(m * ncol(Xs)), m, ncol(Xs))
  for (i in seq_len(iters)) {
    H <- H * (t(W) %*% Xs) / (t(W) %*% W %*% H + 1e-12)
    W <- W * (Xs %*% t(H)) / (W %*% H %*% t(H) + 1e-12)
  }
  D <- H[, 1:64] - H[, 65:128]                         # m x 64 directions
  D / pmax(sqrt(rowSums(D^2)), 1e-12)
}

for (m in M_GRID) {
  D <- nmf_dictionary(X_r, m)
  saveRDS(D, file.path(OUT_DIR, sprintf("baseline_nmf_m%d.rds", m)))
  km <- stats::kmeans(X_r, centers = m, iter.max = 100L, nstart = 1L)
  Dk <- km$centers / pmax(sqrt(rowSums(km$centers^2)), 1e-12)
  saveRDS(Dk, file.path(OUT_DIR, sprintf("baseline_kmeans_m%d.rds", m)))
}

diag_df <- do.call(rbind, diag_rows)
utils::write.csv(diag_df, file.path(OUT_DIR, "diagnostics.csv"), row.names = FALSE)
print(diag_df)
message("Checkpoint: review diagnostics.csv, run dev/03_autolabel.R over the ",
        "candidates, then freeze the winner as inst/extdata/sae_weights_v1.rds.")
