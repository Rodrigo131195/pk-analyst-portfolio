# ================================================================
# ADA Clustering (Gower + PAM): features, silhouette, MDS, boxplots
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(cluster); library(e1071); library(factoextra)
  library(glue)
})

input_path <- "data/public/tbl500_plus_demo.csv"  # columns: ID, ARIA, adalimumab, Days, Week.nr
time_max_days <- 196
ada_high_threshold <- 30
valid_weeks <- c(-2, 0, 4, 16, 28, 52)            # for weekly boxplots
k_range <- 2:10
out_dir <- "derived/clustering"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, w = 7, h = 5, dpi = 300) {
  pngf <- file.path(out_dir, paste0(filename, ".png"))
  pdff <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(pngf, plot, width = w, height = h, dpi = dpi)
  ggsave(pdff, plot, width = w, height = h, device = "pdf")
  message("Saved: ", basename(pngf), " & ", basename(pdff))
}

theme_base <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

# ---------------------- Load & prepare -----------------------
stopifnot(file.exists(input_path))
df <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  select(ID, ARIA, adalimumab, Days, Week.nr) %>%
  filter(Days <= time_max_days) %>%
  rename(ADA = ARIA, CONC = adalimumab, TIME = Days, WEEK = Week.nr)

# keep subjects with >= 2 ADA measures
valid_ids <- df %>%
  filter(!is.na(ADA)) %>%
  count(ID, name = "n") %>%
  filter(n >= 2) %>%
  pull(ID)

df <- df %>% filter(ID %in% valid_ids)

# ------------------- Feature engineering ---------------------
compute_auc <- function(x, y) {
  if (length(x) > 1 && length(y) > 1) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
  } else NA_real_
}

compute_slope <- function(x, y) {
  if (length(x) > 1 && length(y) > 1) {
    lm_fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
    if (!is.null(lm_fit)) as.numeric(coef(lm_fit)[2]) else NA_real_
  } else NA_real_
}

features <- df %>%
  group_by(ID) %>%
  summarise(
    ada_mean  = mean(ADA, na.rm = TRUE),
    ada_max   = max(ADA, na.rm = TRUE),
    ada_var   = var(ADA, na.rm = TRUE),
    ada_range = max(ADA, na.rm = TRUE) - min(ADA, na.rm = TRUE),
    ada_auc   = compute_auc(TIME, ADA),
    ada_slope = compute_slope(TIME, ADA),
    ada_skew  = e1071::skewness(ADA, na.rm = TRUE, type = 2),
    prop_ada_high = mean(ADA > ada_high_threshold, na.rm = TRUE),
    conc_mean = mean(CONC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Keep rows with complete features (simplifies distance/clustering)
  filter(!if_any(-ID, is.na))

stopifnot(nrow(features) > 1)

# -------------------- Gower + PAM clustering -----------------
gower_dist <- cluster::daisy(features %>% select(-ID), metric = "gower")

sil_tbl <- tibble(k = k_range, sil = NA_real_)
pam_models <- vector("list", length(k_range))
names(pam_models) <- as.character(k_range)

set.seed(42)
for (k in k_range) {
  fit <- cluster::pam(gower_dist, diss = TRUE, k = k)
  sil_tbl$sil[sil_tbl$k == k] <- fit$silinfo$avg.width
  pam_models[[as.character(k)]] <- fit
}

best_k <- sil_tbl$k[which.max(sil_tbl$sil)]
message("Best k (by average silhouette width): ", best_k)
final_pam <- pam_models[[as.character(best_k)]]
features$cluster <- factor(final_pam$clustering)

# Silhouette vs k plot
p_sil <- ggplot(sil_tbl, aes(k, sil)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
  labs(title = "Average silhouette width by k", x = "k (clusters)", y = "Avg silhouette width") +
  theme_base
save_both(p_sil, "silhouette_by_k", w = 6.5, h = 4.5)

# Full silhouette plot for best_k
p_sil_full <- factoextra::fviz_silhouette(final_pam) + theme_base +
  labs(title = glue("Silhouette plot (k = {best_k})"))
save_both(p_sil_full, "silhouette_full", w = 7.5, h = 5.5)

# ---------------------- MDS projection -----------------------
mds <- cmdscale(gower_dist, k = 2)
mds_df <- tibble(Dim1 = mds[,1], Dim2 = mds[,2], cluster = features$cluster)

p_mds <- ggplot(mds_df, aes(Dim1, Dim2, color = cluster)) +
  geom_point(size = 2.7, alpha = 0.8) +
  labs(title = "MDS projection of Gower distance", x = "Dim 1", y = "Dim 2", color = "Cluster") +
  theme_base
save_both(p_mds, "mds_gower", w = 6.5, h = 5)

# -------------------- Feature boxplots -----------------------
long_df <- features %>%
  pivot_longer(cols = -c(ID, cluster), names_to = "feature", values_to = "value") %>%
  mutate(feature = recode(feature,
                          ada_mean = "ADA_mean",
                          ada_max = "ADA_max",
                          ada_var = "ADA_var",
                          ada_range = "ADA_range",
                          ada_auc = "AUC",
                          ada_slope = "Slope",
                          ada_skew = "Skew",
                          prop_ada_high = "PctADA>30",
                          conc_mean = "CONC_mean"
  ))

p_feat <- ggplot(long_df, aes(cluster, value, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  labs(title = "Feature distributions by cluster", x = "Cluster", y = "Value") +
  theme_base +
  theme(strip.text = element_text(size = 9),
        panel.spacing = unit(1, "lines"))
save_both(p_feat, "feature_boxplots", w = 8, h = 7)

# -------------- ADA & CONC by week (boxplots) ----------------
df_clustered <- df %>%
  semi_join(features, by = "ID") %>%
  left_join(features %>% select(ID, cluster), by = "ID") %>%
  filter(!is.na(WEEK), WEEK %in% valid_weeks)

p_ada_wk <- ggplot(df_clustered %>% filter(!is.na(ADA)),
                   aes(x = factor(WEEK, levels = valid_weeks), y = ADA)) +
  geom_boxplot(outlier.shape = NA, fill = "firebrick", alpha = 0.6) +
  facet_wrap(~ cluster) +
  labs(title = "ADA by week and cluster", x = "Week number", y = "ADA titer") +
  theme_base
save_both(p_ada_wk, "ada_by_week", w = 8, h = 5)

p_conc_wk <- ggplot(df_clustered %>% filter(!is.na(CONC)),
                    aes(x = factor(WEEK, levels = valid_weeks), y = CONC)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ cluster) +
  labs(title = "Adalimumab concentration by week and cluster", x = "Week number", y = "Concentration") +
  theme_base
save_both(p_conc_wk, "conc_by_week", w = 8, h = 5)

# --------------------- Export labels -------------------------
lab_path <- file.path(out_dir, "ada_clusters_final.csv")
readr::write_csv(features %>% select(ID, cluster), lab_path)
message("Wrote clusters: ", lab_path)
