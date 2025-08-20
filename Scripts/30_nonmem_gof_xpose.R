# ====================================================
# NONMEM GOF plots (DV~PRED, DV~IPRED) + ETA summaries
# ====================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(rlang); library(glue)
})

# Optional: use xpose4 if you have real NONMEM tables locally
has_xpose <- requireNamespace("xpose4", quietly = TRUE)
if (!has_xpose) message("xpose4 not installed; CSV mode will be used unless you change it.")

# ------------------ Config (edit as needed) ------------------
mode    <- if (has_xpose) "csv" else "csv"   # "csv" or "xpose"
runno   <- "02-035"                           # used for xpose mode and output dir naming
csv_path <- "data/public/nonmem_sdtab_demo.csv" # expected cols: ID,TAD,CMTX,DV,PRED,IPRED,ETA1,ETA2,ETA3 (optional)
run_dir <- file.path("runs", runno)           # where NONMEM tables would live (xpose mode)
out_dir <- file.path("derived/nonmem_gof", runno)
use_log_axes <- FALSE                         # set TRUE for log10 axes

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, w = 7, h = 5, dpi = 300) {
  pngf <- file.path(out_dir, paste0(filename, ".png"))
  pdff <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(pngf, plot, width = w, height = h, dpi = dpi)
  ggsave(pdff, plot, width = w, height = h, device = "pdf")
  message("Saved: ", basename(pngf), " & ", basename(pdff))
}

# ------------------ Load data ------------------
if (mode == "xpose") {
  if (!has_xpose) stop("xpose4 not available. Install it or switch to CSV mode.")
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  stopifnot(dir.exists(run_dir))
  setwd(run_dir)
  run_df <- xpose4::read.nm.tables(runno = runno)
} else {
  stopifnot(file.exists(csv_path))
  run_df <- readr::read_csv(csv_path, show_col_types = FALSE)
}

# ------------------ Sanity checks & filters ------------------
required <- c("ID","TAD","CMTX","DV","PRED","IPRED")
miss <- setdiff(required, names(run_df))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

run_df <- run_df %>%
  filter(TAD != 0, CMTX == 2) %>%            # central compartment rows, post-dose
  filter(DV > 0, PRED > 0, IPRED > 0) %>%    # log-safe if needed
  mutate(across(c(DV, PRED, IPRED), as.numeric))

# ------------------ Plots: DV vs PRED / DV vs IPRED ------------------
base_theme <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top")

maybe_log <- function(p) {
  if (use_log_axes) {
    p + scale_x_log10() + scale_y_log10()
  } else p
}

p1 <- ggplot(run_df, aes(x = DV, y = PRED)) +
  geom_point(alpha = 0.25, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Goodness-of-Fit: DV vs PRED",
       x = "Observation (DV)", y = "Population Prediction (PRED)") +
  base_theme
p1 <- maybe_log(p1)
save_both(p1, "gof_dv_pred")

p2 <- ggplot(run_df, aes(x = DV, y = IPRED)) +
  geom_point(alpha = 0.25, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Goodness-of-Fit: DV vs IPRED",
       x = "Observation (DV)", y = "Individual Prediction (IPRED)") +
  base_theme
p2 <- maybe_log(p2)
save_both(p2, "gof_dv_ipred")

# ------------------ ETA summaries (if available) ------------------
eta_cols <- grep("^ETA\\d+$", names(run_df), value = TRUE)
if (length(eta_cols)) {
  # Map ETA1/ETA2/ETA3 → CL/Vc/KA if those exist; otherwise keep generic
  name_map <- c("ETA1" = "CL", "ETA2" = "Vc", "ETA3" = "KA")
  eta_long <- run_df %>%
    select(ID, all_of(eta_cols)) %>%
    distinct() %>%
    pivot_longer(-ID, names_to = "ETA_name", values_to = "ETA") %>%
    mutate(IIV = factor(dplyr::recode(ETA_name, !!!name_map, .default = ETA_name)))
  
  p3 <- ggplot(eta_long, aes(IIV, ETA)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(shape = 21, width = 0.12, alpha = 0.5, size = 1.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(title = "ETA distribution by parameter", x = "IIV parameter", y = "ETA") +
    base_theme
  save_both(p3, "eta_boxplots")
  
  # Focused Vc panel only if present
  if ("ETA2" %in% eta_cols) {
    vc_df <- eta_long %>% filter(IIV %in% c("Vc","ETA2"))
    p_vc <- ggplot(vc_df, a
                   