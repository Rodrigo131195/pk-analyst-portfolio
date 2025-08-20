# ===================================================
# GRU Autoencoder → Embeddings → k-means → Clusters
# ===================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(cluster); library(factoextra); library(torch); library(glue)
})

# ---------------- Config (edit as needed) ----------------
set.seed(42); torch_manual_seed(42)
input_path   <- "data/public/tbl500_plus_demo.csv"   # cols: ID, ARIA, adalimumab, Days
time_max     <- 196
time_grid    <- c(0, 14, 28, 56, 84, 168)            # days
latent_dim   <- 4
learning_rate<- 0.010
epochs       <- 2000
max_clusters <- 10
use_gpu      <- FALSE                                 # set TRUE if CUDA available
out_dir      <- "derived/gru_autoencoder"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, w = 7, h = 5, dpi = 300) {
  pngf <- file.path(out_dir, paste0(filename, ".png"))
  pdff <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(pngf, plot, width = w, height = h, dpi = dpi)
  ggsave(pdff, plot, width = w, height = h, device = "pdf")
  message("Saved: ", basename(pngf), " & ", basename(pdff))
}

# ---------------- Load & prepare ----------------
stopifnot(file.exists(input_path))
tbl500 <- readr::read_csv(input_path, show_col_types = FALSE)

df <- tbl500 %>%
  select(ID, ARIA, adalimumab, Days) %>%
  rename(ADA = ARIA, CONC = adalimumab, TIME = Days) %>%
  filter(TIME <= time_max) %>%
  arrange(ID, TIME) %>%
  group_by(ID) %>%
  filter(sum(!is.na(ADA)) >= 2) %>%
  ungroup()

# Interpolation helpers (robust to NAs; rule=2 clamps ends)
interp_to_grid <- function(x_time, x_val, grid) {
  ok <- is.finite(x_time) & is.finite(x_val)
  x_time <- x_time[ok]; x_val <- x_val[ok]
  if (length(x_time) < 1) return(rep(NA_real_, length(grid)))
  if (length(unique(x_time)) == 1) {
    rep(x_val[1], length(grid))
  } else {
    approx(x_time, x_val, xout = grid, rule = 2)$y
  }
}

# Build patient sequences on the grid: [patient, time, feature(ADA,CONC)]
ids <- unique(df$ID)
seq_list <- lapply(ids, function(id) {
  sub <- df %>% filter(ID == id)
  ada  <- interp_to_grid(sub$TIME, sub$ADA,  time_grid)
  conc <- interp_to_grid(sub$TIME, sub$CONC, time_grid)
  cbind(ADA = ada, CONC = conc)
})

sequence_array <- array(unlist(seq_list), dim = c(length(ids), length(time_grid), 2))

# Standardize features across all patients/time (per feature)
feat_means <- apply(sequence_array, 3, function(m) mean(as.vector(m), na.rm = TRUE))
feat_sds   <- apply(sequence_array, 3, function(m) sd(as.vector(m),   na.rm = TRUE))
for (j in 1:2) {
  sequence_array[,,j] <- (sequence_array[,,j] - feat_means[j]) / ifelse(feat_sds[j] == 0, 1, feat_sds[j])
}

# Replace any remaining NA from interpolation with 0 (post-standardization)
sequence_array[is.na(sequence_array)] <- 0

# Torch tensors
device <- if (use_gpu && cuda_is_available()) torch_device("cuda") else torch_device("cpu")
X <- torch_tensor(sequence_array, dtype = torch_float(), device = device)

n_patients <- dim(X)[1]
time_steps <- dim(X)[2]
n_features <- dim(X)[3]

# ---------------- Model ----------------
gru_autoencoder <- nn_module(
  "GRUAutoencoder",
  initialize = function() {
    self$enc <- nn_gru(input_size = n_features, hidden_size = latent_dim, batch_first = TRUE)
    self$dec <- nn_gru(input_size = n_features, hidden_size = latent_dim, batch_first = TRUE)
    self$out <- nn_linear(latent_dim, n_features)
  },
  forward = function(x) {
    enc_out <- self$enc(x)     # list(output, h_n)
    h <- enc_out_
    