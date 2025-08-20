# =============================================================
# ADA CLUSTERING - FEATURE SELECTION + CLUSTER OPTIMIZATION
# Extended feature grid + silhouette and elbow methods
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(cluster); library(e1071); library(factoextra); library(rstatix)
  library(purrr); library(zoo); library(glue)
})

# ----------------------#
# Config (edit paths)
# ----------------------#
input_path <- "data/public/tbl500_plus_demo.csv"   # cols: ID, Days, Week.nr, ARIA, adalimumab
time_max_days <- 196
k_range <- 2:10
ada_high_30 <- 30
ada_high_100 <- 100
out_dir <- "derived/clustering_feature_sel"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, w = 7, h = 5, dpi = 300) {
  pngf <- file.path(out_dir, paste0(filename, ".png"))
  pdff <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(pngf, plot, width = w, height = h, dpi = dpi)
  ggsave(pdff, plot, width = w, height = h, device = "pdf")
  message("Saved: ", basename(pngf), " & ", basename(pdff))
}

normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# ----------------------#
# Load and Prepare Data
# ----------------------#
stopifnot(file.exists(input_path))
df <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  select(ID, Days, Week.nr, ARIA, adalimumab) %>%
  rename(TIME = Days, WEEK = Week.nr, ADA = ARIA, CONC = adalimumab) %>%
  filter(TIME <= time_max_days)

# keep subjects with >= 2 ADA measures
valid_ids <- df %>%
  filter(!is.na(ADA)) %>%
  count(ID, name = "n") %>%
  filter(n >= 2) %>%
  pull(ID)

df <- df %>% filter(ID %in% valid_ids)

# ----------------------#
# Extended Feature Extraction
# ----------------------#
compute_auc <- function(x, y) {
  if (length(x) > 1 && length(y) > 1) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
  } else NA_real_
}

safe_baseline <- function(x_time, x_val) {
  o <- order(x_time)
  x_val <- x_val[o]
  x_val[which(!is.na(x_val))[1]] %||% NA_real_
}

feature_df <- df %>%
  arrange(ID, TIME) %>%
  group_by(ID) %>%
  summarise(
    ada_baseline = safe_baseline(TIME, ADA),
    ada_max      = suppressWarnings(max(ADA, na.rm = TRUE)),
    ada_min      = suppressWarnings(min(ADA, na.rm = TRUE)),
    ada_mean     = mean(ADA, na.rm = TRUE),
    ada_var      = var(ADA, na.rm = TRUE),
    ada_skew     = e1071::skewness(ADA, na.rm = TRUE, type = 2),
    ada_kurt     = e1071::kurtosis(ADA, na.rm = TRUE, type = 2),
    ada_auc      = compute_auc(TIME, ADA),
    ada_range    = suppressWarnings(max(ADA, na.rm = TRUE) - min(ADA, na.rm = TRUE)),
    pos_ada_count   = sum(ADA > ada_high_30, na.rm = TRUE),
    prop_ada_high   = mean(ADA > ada_high_100, na.rm = TRUE),
    time_max_ada    = TIME[which.max(ADA)][1],
    duration_ada_high = if (sum(ADA > ada_high_30, na.rm = TRUE) >= 2) {
      max(TIME[ADA > ada_high_30], na.rm = TRUE) - min(TIME[ADA > ada_high_30], na.rm = TRUE)
    } else NA_real_,
    ada_slope    = if (length(na.omit(ADA)) > 1 && var(TIME, na.rm = TRUE) > 0) {
      as.numeric(coef(lm(ADA ~ TIME))[2])
    } else NA_real_,
    conc_mean    = mean(CONC, na.rm = TRUE),
    conc_min     = suppressWarnings(min(CONC, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  # Drop rows with all-NA or impossible stats
  filter(!if_any(everything(), ~ is.nan(.x))) %>%
  tidyr::drop_na()

stopifnot(nrow(feature_df) > 1)

# ----------------------#
# Distance once (Gower)
# ----------------------#
gower_dist <- cluster::daisy(feature_df %>% select(-ID), metric = "gower")

# ----------------------#
# Grid Search Across k
# ----------------------#
results <- vector("list", length(k_range)); names(results) <- as.character(k_range)
elbow_vals <- numeric(length(k_range)); names(elbow_vals) <- as.character(k_range)
sil_vals <- numeric(length(k_range)); names(sil_vals) <- as.character(k_range)

set.seed(42)
for (k in k_range) {
  pam_fit <- cluster::pam(gower_dist, diss = TRUE, k = k)
  sil_vals[as.character(k)] <- pam_fit$silinfo$avg.width
  elbow_vals[as.character(k)] <- pam_fit$objective[1]  # total dissimilarity (elbow)
  
  cluster_assignments <- pam_fit$clustering
  
  feat_long <- feature_df %>%
    mutate(cluster = factor(cluster_assignments)) %>%
    pivot_longer(-c(ID, cluster), names_to = "feature", values_to = "value")
  
  # Kruskal–Wallis per feature
  feat_stats <- feat_long %>%
    group_by(feature) %>%
    group_split() %>%
    map_df(function(df_feat) {
      p_val <- tryCatch(
        rstatix::kruskal_test(df_feat, value ~ cluster)$p,
        error = function(e) NA_real_
      )
      tibble(
        feature = unique(df_feat$feature),
        p = p_val,
        k = k,
        sil = pam_fit$silinfo$avg.width
      )
    })
  
  results[[as.character(k)]] <- feat_stats
}

# ----------------------#
# Combine and Summarize
# ----------------------#
all_results <- bind_rows(results)

summary_k <- all_results %>%
  group_by(k) %>%
  summarise(
    avg_sil = unique(sil)[1],
    n_sig   = sum(p < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    elbow = as.numeric(elbow_vals[as.character(k)]),
    sil_norm   = normalize01(avg_sil),
    elbow_norm = normalize01(elbow),
    nsig_norm  = normalize01(n_sig)
  )

# ----------------------#
# Plot Silhouette + Elbow + #Sig
# ----------------------#
p_combo <- ggplot(summary_k, aes(x = k)) +
  geom_line(aes(y = sil_norm), size = 1.2) +
  geom_point(aes(y = sil_norm), size = 2) +
  geom_line(aes(y = elbow_norm), linetype = "dashed") +
  geom_point(aes(y = elbow_norm), shape = 1) +
  geom_line(aes(y = nsig_norm), color = "firebrick") +
  geom_point(aes(y = nsig_norm), color = "firebrick", size = 2) +
  labs(
    title = "Model selection signals vs k",
    subtitle = "Silhouette (solid), Elbow (dashed), Significant features (red)",
    y = "Normalized score", x = "Number of clusters (k)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_both(p_combo, "k_selection_signals", w = 7.5, h = 5)

# ----------------------#
# Best k and Top Features
# ----------------------#
best_k <- summary_k$k[which.max(summary_k$avg_sil)]
best_features <- all_results %>%
  filter(k == best_k, !is.na(p)) %>%
  arrange(p)

message("Best k (by silhouette): ", best_k)
utils::write.csv(best_features, file.path(out_dir, glue("best_features_k{best_k}.csv")), row.names = FALSE)
utils::write.csv(summary_k,   file.path(out_dir, "clustering_summary.csv"), row.names = FALSE)

# Refit chosen model to export clusters and a couple of visuals
final_pam <- cluster::pam(gower_dist, diss = TRUE, k = best_k)
assignments <- tibble(ID = feature_df$ID, cluster = factor(final_pam$clustering))
readr::write_csv(assignments, file.path(out_dir, glue("clusters_k{best_k}.csv")))
message("Wrote cluster labels: clusters_k", best_k, ".csv")

# Silhouette plot for best_k
p_sil_full <- factoextra::fviz_silhouette(final_pam) +
  theme_minimal(base_size = 12) + labs(title = glue("Silhouette plot (k = {best_k})"))
save_both(p_sil_full, glue("silhouette_full_k{best_k}"), w = 7.5, h = 5.5)

# MDS projection
mds <- cmdscale(gower_dist, k = 2)
mds_df <- tibble(Dim1 = mds[,1], Dim2 = mds[,2]) %>%
  bind_cols(assignments %>% select(cluster))

p_mds <- ggplot(mds_df, aes(Dim1, Dim2, color = cluster)) +
  geom_point(size = 2.7, alpha = 0.85) +
  labs(title = glue("MDS projection of Gower distance (k = {best_k})"),
       x = "Dim 1", y = "Dim 2", color = "Cluster") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
save_both(p_mds, glue("mds_gower_k{best_k}"), w = 7, h = 5)

message("All done. Outputs in: ", normalizePath(out_dir, mustWork = FALSE))
