# ================================
# Initial Data Inspection — individual figures (Panels B–F)
# ================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(glue)
  # if you prefer project-rooted paths, you can also use: library(here)
})

# ---------- Config (EDIT THESE IF NEEDED) ----------
# Point to a public/synthetic file in your repo, not your internship drive.
# Supports .xlsx or .csv with columns: ID/MapNummerOld, Days, ARIA, adalimumab
input_path <- "data/public/tbl500_demo.xlsx"   # e.g., "data/public/tbl500_demo.xlsx"
input_sheet <- 1                                # ignored for CSV
out_dir <- "derived/figures/initial_inspection"
ada_thresh_report <- 70                         # ADA-positive threshold for panel F
every_n_for_x <- 4                              # label every N weeks on x-axis

# ---------- I/O helpers ----------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, w = 6, h = 4, dpi = 300) {
  fn_png <- file.path(out_dir, paste0(filename, ".png"))
  fn_pdf <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(fn_png, plot, width = w, height = h, dpi = dpi)
  ggsave(fn_pdf, plot, width = w, height = h, device = "pdf")
  message("Saved: ", basename(fn_png), " & ", basename(fn_pdf))
}

short_num <- scales::label_number(scale_cut = scales::cut_short_scale())

every_n_weeks <- function(x, n = 4) {
  ix <- suppressWarnings(as.integer(as.character(x)))
  x[!is.na(ix) & (ix %% n == 0)]
}

# ---------- Read data (xlsx or csv) ----------
read_tbl500 <- function(path, sheet = 1) {
  stopifnot(file.exists(path))
  ext <- tolower(tools::file_ext(path))
  dat <- switch(ext,
                "xlsx" = readxl::read_excel(path, sheet = sheet),
                "csv"  = readr::read_csv(path, show_col_types = FALSE),
                stop("Unsupported file extension: ", ext)
  )
  # standardize ID name if needed
  if ("MapNummerOld" %in% names(dat) && !"ID" %in% names(dat)) {
    dat <- dat %>% rename(ID = MapNummerOld)
  }
  req <- c("ID","Days","ARIA","adalimumab")
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  dat
}

tbl500 <- read_tbl500(input_path, input_sheet) %>%
  mutate(Week = floor(Days / 7))

# ---------- ADA status per patient ----------
pat_ada_status <- tbl500 %>%
  group_by(ID) %>%
  summarise(
    any_ada_measured = any(!is.na(ARIA)),
    any_ada_pos      = any(ARIA >= ada_thresh_report, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ADA_STATUS = case_when(
      !any_ada_measured               ~ "Not measured",
      any_ada_pos                     ~ "ADA+",
      any_ada_measured & !any_ada_pos ~ "ADA-"
    )
  ) %>% select(ID, ADA_STATUS)

df <- tbl500 %>%
  left_join(pat_ada_status, by = "ID") %>%
  mutate(ADA_STATUS = factor(ADA_STATUS, levels = c("ADA+","ADA-","Not measured")))

# Consistent colors
group_cols  <- c("ADA+" = "#D62728", "ADA-" = "#1F77B4", "Not measured" = "grey60")
group_fills <- c("ADA+" = alpha("#D62728", 0.35),
                 "ADA-" = alpha("#1F77B4", 0.35),
                 "Not measured" = alpha("grey60", 0.35))

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# =========================================
# Panel B — ADA over time (BOXPLOTS), log y
# =========================================
b_dat <- df %>% filter(!is.na(ARIA), !is.na(ADA_STATUS), ARIA > 0) %>% mutate(WeekF = factor(Week))

pB <- ggplot(b_dat, aes(x = WeekF, y = ARIA, fill = ADA_STATUS)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"),
               outlier.alpha = 0.25, outlier.size = 0.7) +
  scale_fill_manual(values = group_fills, drop = FALSE) +
  scale_y_log10(labels = short_num) +
  scale_x_discrete(breaks = function(x) every_n_weeks(x, n = every_n_for_x)) +
  labs(title = "ADA over time by ADA status",
       x = "Week", y = "ADA (AU/mL, log scale)", fill = "ADA status") +
  base_theme

save_both(pB, "fig1B_ada_time_box", w = 9, h = 5)

# =========================================
# Panel C — Adalimumab over time (BOXPLOTS)
# =========================================
c_dat <- df %>% filter(!is.na(adalimumab), !is.na(ADA_STATUS)) %>% mutate(WeekF = factor(Week))

pC <- ggplot(c_dat, aes(x = WeekF, y = adalimumab, fill = ADA_STATUS)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"),
               outlier.alpha = 0.25, outlier.size = 0.7) +
  scale_fill_manual(values = group_fills, drop = FALSE) +
  scale_x_discrete(breaks = function(x) every_n_weeks(x, n = every_n_for_x)) +
  labs(title = "Adalimumab over time by ADA status",
       x = "Week", y = "Adalimumab (mg/L)", fill = "ADA status") +
  base_theme

save_both(pC, "fig1C_drug_time_box", w = 9, h = 5)

# =========================================
# Panel D — Filtered ADA vs Adalimumab, log–log, with density + rho
# =========================================
d_dat <- df %>%
  filter(!is.na(ARIA), !is.na(adalimumab), !is.na(ADA_STATUS),
         ARIA >= 30, adalimumab >= 0.1,
         ADA_STATUS %in% c("ADA+","ADA-")) %>%
  mutate(logA = log10(adalimumab), logADA = log10(ARIA))

# Spearman rho per group (on log10 variables)
rho_df <- dplyr::summarise(group_by(d_dat, ADA_STATUS),
                           n = dplyr::n(),
                           rho = suppressWarnings(cor(logA, logADA, method = "spearman")),
                           .groups = "drop")

# per-facet annotation positions (use group-wise mins/maxes)
ann_pos <- d_dat %>%
  group_by(ADA_STATUS) %>%
  summarise(x_lab = min(adalimumab, na.rm = TRUE) * 1.2,
            y_lab = max(ARIA, na.rm = TRUE) / 1.2,
            .groups = "drop")

ann_df <- left_join(rho_df, ann_pos, by = "ADA_STATUS") %>%
  mutate(label = glue("Spearman \u03C1 = {sprintf('%.2f', rho)}\nN = {n}"))

pD <- ggplot(d_dat, aes(x = adalimumab, y = ARIA)) +
  geom_point(alpha = 0.25, size = 1) +
  stat_density_2d(aes(color = after_stat(level)), bins = 6, show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8, linewidth = 0.8, color = "black") +
  scale_x_log10(labels = short_num) +
  scale_y_log10(labels = short_num) +
  facet_wrap(~ ADA_STATUS, nrow = 1, scales = "fixed") +
  labs(title = "ADA vs adalimumab (filtered: ADA ≥ 30 AU/mL, drug ≥ 0.1 mg/L)",
       x = "Adalimumab (mg/L, log scale)", y = "ADA (AU/mL, log scale)") +
  base_theme +
  geom_text(data = ann_df, aes(x = x_lab, y = y_lab, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5)

save_both(pD, "fig1D_ada_vs_drug_filtered", w = 9, h = 4.5)

# =========================================
# Panel E — Sampling coverage heatmap
# =========================================
covr <- df %>%
  select(ID, Week, ARIA, adalimumab) %>%
  pivot_longer(c("ARIA","adalimumab"), names_to = "measure", values_to = "value") %>%
  mutate(measured = !is.na(value)) %>%
  group_by(measure, Week) %>%
  summarise(n = sum(measured), .groups = "drop") %>%
  mutate(measure = recode(measure,
                          "ARIA" = "ADA measured",
                          "adalimumab" = "Adalimumab measured"))

pE <- ggplot(covr, aes(x = Week, y = measure, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "grey90", high = "black", limits = c(0, max(covr$n))) +
  labs(title = "Sampling coverage", x = "Week", y = NULL, fill = "N samples") +
  base_theme +
  theme(legend.position = "right")

save_both(pE, "fig1E_sampling_coverage", w = 7, h = 3.8)

# =========================================
# Panel F — ECDF of time to first ADA-positive (≥70 AU/mL)
# =========================================
first_pos <- df %>%
  group_by(ID) %>%
  summarise(
    any_measured  = any(!is.na(ARIA)),
    any_pos       = any(ARIA >= ada_thresh_report, na.rm = TRUE),
    first_pos_week = suppressWarnings(min(Week[ARIA >= ada_thresh_report], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(any_measured, any_pos, is.finite(first_pos_week))

f_median <- median(first_pos$first_pos_week, na.rm = TRUE)
f_q1     <- quantile(first_pos$first_pos_week, 0.25, na.rm = TRUE)
f_q3     <- quantile(first_pos$first_pos_week, 0.75, na.rm = TRUE)

pF <- ggplot(first_pos, aes(x = first_pos_week)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = c(f_q1, f_median, f_q3),
             linetype = c("dotted","dashed","dotted")) +
  annotate("text", x = f_median, y = 0.05,
           label = paste0("Median = ", round(f_median, 1), " wk"),
           vjust = 0, hjust = -0.1, size = 3.5) +
  labs(title = "Time to first ADA-positive (≥70 AU/mL)",
       x = "Week of first ADA-positive", y = "Cumulative fraction of patients") +
  base_theme

save_both(pF, "fig1F_time_to_first_ada_ecdf", w = 6, h = 4)

# Preview in IDE if desired
# print(pB); print(pC); print(pD); print(pE); print(pF)
message("All figures written to: ", normalizePath(out_dir, mustWork = FALSE))
