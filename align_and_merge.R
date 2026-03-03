#!/usr/bin/env Rscript
# align_and_merge.R
# Purpose: Align 11B calibration runs to their corresponding AllCal runs for a
# specific detector/crystal, merge unique peaks, fit a linear calibration to the
# combined peak list, and characterise residual structure via a CV-selected
# smoothing spline.

## ------------------------- Configuration ---------------------------------
# Path to file containing fit funcitons and helper utilities for energy calibration
energy_fit_path <- "/home/tylermk/TUNL/Data/NRF/70Ge/energy_calibration/energy_fit.r"
# Directory to save output plots and calibration parameters
outdir <- "/home/tylermk/TUNL/Data/NRF/70Ge/energy_calibration/calibrations"
# Directory containing cubix_workspaces with peak fits for all runs (will search recursively)
cubix_files <- "/home/tylermk/TUNL/Data/NRF/70Ge/cubix_workspaces"

# Minimum common peaks to perform per-run alignment
min_matches <- 2
# Digits used when matching energies across files (≈eV precision)
energy_match_digits <- 3
# Smoothing spline grid resolution
spline_grid_n <- 400
# Maximum number of spline knots (basis dimension)
spline_max_knots <- 15
# GAM smoothness selection method: "REML", "GCV.Cp", "ML", etc.
spline_method <- "REML"
# Fixed smoothing parameter (NULL = auto-select via spline_method)
spline_sp <- 5
# Whether to save linear calibration plots
plot_linear_cal <- FALSE

# All detector/crystal combinations to process
all_detectors <- c("C1", "C3", "C5", "C7", "B1", "B2", "B3", "B5")
all_crystals  <- 1:4

## ------------------------- Dependencies & loader -------------------------
ensure_packages <- function() {
  suppressPackageStartupMessages({
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 required")
    library(ggplot2)
    if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package mgcv required")
    library(mgcv)
  })
}

load_data <- function() {
  if (!file.exists(energy_fit_path)) stop("Missing energy_fit.r at: ", energy_fit_path)
  source(energy_fit_path)
  ensure_packages()
  collect_calibrations(root = cubix_files, recursive = TRUE, verbose = FALSE)
}

## ------------------------- Utility helpers ------------------------------
energy_key <- function(energy, digits = energy_match_digits) {
  round(as.numeric(energy), digits)
}

as_df <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

dedupe_peaks <- function(df, digits = energy_match_digits) {
  df <- as_df(df)
  if (!all(c("channel", "energy") %in% names(df))) return(df[0, , drop = FALSE])
  df <- df[!is.na(df$energy) & !is.na(df$channel), , drop = FALSE]
  if (nrow(df) == 0) return(df)
  df$key <- energy_key(df$energy, digits)
  keep <- !duplicated(df$key)
  out <- df[keep, , drop = FALSE]
  out$key <- NULL
  out
}

assign_channel_error <- function(df) {
  df <- as_df(df)
  if (!"channel" %in% names(df) || nrow(df) == 0) return(df)
  errs <- ifelse(df$channel <= 8500, 0.05,
                 ifelse(df$channel <= 22000, 0.1,
                        ifelse(df$channel <= 30000, 0.7, 1.8)))
  df$xerr <- errs
  df
}

select_common_peaks <- function(run_pts, allcal_pts, digits = energy_match_digits) {
  run_df <- dedupe_peaks(run_pts, digits)
  all_df <- dedupe_peaks(allcal_pts, digits)
  if (nrow(run_df) == 0 || nrow(all_df) == 0) return(data.frame())
  run_df$energy_key <- energy_key(run_df$energy, digits)
  all_df$energy_key <- energy_key(all_df$energy, digits)
  common <- merge(run_df, all_df, by = "energy_key", suffixes = c(".run", ".all"))
  common$energy_key <- NULL
  common
}

align_columns <- function(df, cols) {
  out <- as_df(df)
  if (length(cols) == 0) return(out)
  missing <- setdiff(cols, names(out))
  if (length(missing) > 0) {
    for (m in missing) out[[m]] <- NA
  }
  out <- out[, cols, drop = FALSE]
  out
}

append_unique_peaks <- function(existing, additions, seen_keys, digits = energy_match_digits) {
  existing <- as_df(existing)
  additions <- as_df(additions)
  if (!"energy" %in% names(additions)) {
    return(list(data = existing, seen = seen_keys, added = 0))
  }
  additions <- additions[!is.na(additions$energy), , drop = FALSE]
  if (nrow(additions) == 0) return(list(data = existing, seen = seen_keys, added = 0))
  additions$energy_key <- energy_key(additions$energy, digits)
  keep <- !(additions$energy_key %in% seen_keys)
  additions <- additions[keep, , drop = FALSE]
  if (nrow(additions) == 0) return(list(data = existing, seen = seen_keys, added = 0))

  additions$energy_key <- NULL
  all_cols <- union(names(existing), names(additions))
  existing <- align_columns(existing, all_cols)
  additions <- align_columns(additions, all_cols)
  combined <- rbind(existing, additions)
  list(data = combined, seen = c(seen_keys, energy_key(additions$energy, digits)), added = nrow(additions))
}

fit_alignment <- function(common) {
  sigma_y <- if ("xerr.all" %in% names(common)) pmax(common$xerr.all, .Machine$double.eps) else NA_real_
  sigma_y2 <- ifelse(is.na(sigma_y), 0, sigma_y^2)

  sigma_x <- if ("xerr.run" %in% names(common)) pmax(common$xerr.run, .Machine$double.eps) else NA_real_
  sigma_x2 <- ifelse(is.na(sigma_x), 0, sigma_x^2)

  if ("xerr" %in% names(common) && all(is.na(sigma_x))) {
    sigma_x2 <- pmax(common$xerr, .Machine$double.eps)^2
  }
  if ("xerr" %in% names(common) && all(is.na(sigma_y))) {
    sigma_y2 <- pmax(common$xerr, .Machine$double.eps)^2
  }

  sigma_tot <- sigma_y2 + sigma_x2
  if (all(sigma_tot <= 0 | is.na(sigma_tot))) {
    lm(channel.all ~ channel.run, data = common)
  } else {
    weights <- 1 / pmax(sigma_tot, .Machine$double.eps^2)
    lm(channel.all ~ channel.run, data = common, weights = weights)
  }
}

apply_correction <- function(pts_run, fit) {
  pts_df <- as_df(pts_run)
  pred <- predict(fit, newdata = data.frame(channel.run = pts_df$channel), se.fit = TRUE)
  slope <- unname(if (!is.null(coef(fit)["channel.run"])) coef(fit)["channel.run"] else coef(fit)[2])
  meas_sigma <- if ("xerr" %in% names(pts_df)) pts_df$xerr else NA_real_
  var_meas <- ifelse(is.na(meas_sigma) | is.na(slope), 0, (slope^2) * (meas_sigma^2))
  var_fit <- pred$se.fit^2
  total_var <- var_meas + var_fit
  pts_df$channel_corrected <- pred$fit
  pts_df$channel_corrected_err <- sqrt(pmax(total_var, 0))
  pts_df$channel_corrected_err[!is.finite(pts_df$channel_corrected_err)] <- NA_real_
  pts_df
}

fit_global_cal <- function(df_combined) {
  w <- NULL
  if ("yerr" %in% names(df_combined) && any(!is.na(df_combined$yerr))) {
    w <- 1 / (pmax(df_combined$yerr, .Machine$double.eps)^2)
  }
  if (is.null(w)) {
    lm(energy ~ channel_corrected, data = df_combined)
  } else {
    lm(energy ~ channel_corrected, data = df_combined, weights = w)
  }
}

fit_spline <- function(df_combined, grid_n = spline_grid_n) {
  n_total <- nrow(df_combined)
  idx <- which(!is.na(df_combined$pred) & !is.na(df_combined$resid))
  if (length(idx) < 3) {
    empty_vec <- rep(NA_real_, n_total)
    return(list(df_smooth = NULL, spline_model = NULL, sp_data_y = empty_vec, sp_data_se = empty_vec, resid_sigma = empty_vec))
  }

  err_vec <- if ("energy_err_total" %in% names(df_combined)) df_combined$energy_err_total else if ("yerr" %in% names(df_combined)) df_combined$yerr else NA_real_
  weight_vals <- ifelse(is.na(err_vec[idx]), NA_real_, 1 / (pmax(err_vec[idx], .Machine$double.eps)^2))
  if (all(is.na(weight_vals))) weight_vals <- rep(1, length(idx))

  tmp <- data.frame(
    x = df_combined$pred[idx],
    y = df_combined$resid[idx],
    w = ifelse(is.na(weight_vals), 1, weight_vals)
  )
  tmp <- tmp[is.finite(tmp$x) & is.finite(tmp$y) & is.finite(tmp$w), , drop = FALSE]
  if (nrow(tmp) < 3) {
    empty_vec <- rep(NA_real_, n_total)
    return(list(df_smooth = NULL, spline_model = NULL, sp_data_y = empty_vec, sp_data_se = empty_vec, resid_sigma = empty_vec))
  }

  unique_x <- sort(unique(tmp$x))
  if (length(unique_x) < 3) {
    empty_vec <- rep(NA_real_, n_total)
    return(list(df_smooth = NULL, spline_model = NULL, sp_data_y = empty_vec, sp_data_se = empty_vec, resid_sigma = empty_vec))
  }

  basis_k <- min(length(unique_x) - 1, spline_max_knots)
  resid_sd <- sd(tmp$y, na.rm = TRUE)
  if (!is.finite(resid_sd) || resid_sd <= 0) {
    resid_sd <- mean(abs(tmp$y), na.rm = TRUE)
  }
  if (!is.finite(resid_sd) || resid_sd <= 0) resid_sd <- NA_real_

  sp_arg <- if (!is.null(spline_sp)) list(sp = spline_sp) else list()

  spline_model <- tryCatch({
    do.call(mgcv::gam, c(list(formula = y ~ s(x, k = basis_k, bs = "cr"), data = tmp, weights = w, method = spline_method), sp_arg))
  }, error = function(e) {
    warning("GAM fit failed with weights (", e$message, "), retrying without weights")
    tryCatch(do.call(mgcv::gam, c(list(formula = y ~ s(x, k = basis_k, bs = "cr"), data = tmp, method = spline_method), sp_arg)), error = function(e2) {
      warning("GAM fit failed: ", e2$message)
      NULL
    })
  })
  if (is.null(spline_model)) {
    empty_vec <- rep(NA_real_, n_total)
    return(list(df_smooth = NULL, spline_model = NULL, sp_data_y = empty_vec, sp_data_se = empty_vec, resid_sigma = empty_vec))
  }

  xs <- seq(min(tmp$x), max(tmp$x), length.out = grid_n)
  grid_df <- data.frame(x = xs)
  sp <- tryCatch(predict(spline_model, newdata = grid_df, se.fit = TRUE, type = "response"), error = function(e) {
    warning("GAM prediction on grid failed: ", e$message)
    NULL
  })

  df_smooth <- NULL
  if (!is.null(sp)) {
    sevec <- sp$se.fit
    if (is.null(sevec) || all(is.na(sevec))) {
      sevec <- rep(resid_sd, length(sp$fit))
    }
    if (any(sevec < 0, na.rm = TRUE)) sevec[sevec < 0] <- NA_real_
    df_smooth <- data.frame(x = xs, y = sp$fit, se = sevec)
  }

  sp_data_y <- rep(NA_real_, nrow(df_combined))
  sp_data_se <- rep(NA_real_, nrow(df_combined))
  if (length(idx) > 0) {
    obs_df <- data.frame(x = df_combined$pred[idx])
    sp_obs <- tryCatch(predict(spline_model, newdata = obs_df, se.fit = TRUE, type = "response"), error = function(e) {
      warning("GAM prediction on observations failed: ", e$message)
      NULL
    })
    if (!is.null(sp_obs)) {
      sp_data_y[idx] <- sp_obs$fit
      se_pred <- sp_obs$se.fit
      if (is.null(se_pred) || all(is.na(se_pred))) se_pred <- rep(resid_sd, length(sp_obs$fit))
      if (any(se_pred < 0, na.rm = TRUE)) se_pred[se_pred < 0] <- NA_real_
      sp_data_se[idx] <- se_pred
    }
  }

  meas_vars <- pmax(err_vec^2, 0)
  resid_sigma <- sqrt(pmax(meas_vars + sp_data_se^2, 0))

  list(df_smooth = df_smooth, spline_model = spline_model, sp_data_y = sp_data_y, sp_data_se = sp_data_se, resid_sigma = resid_sigma)
}

## ------------------------- Plotting helpers -----------------------------
plot_per_run <- function(common, rt, outdir, cfg_detector, cfg_crystal) {
  if (nrow(common) == 0) return(invisible(NULL))
  run_err_df <- if ("xerr.run" %in% names(common)) common[!is.na(common$xerr.run) & is.finite(common$xerr.run), , drop = FALSE] else common[0, , drop = FALSE]
  all_err_df <- if ("xerr.all" %in% names(common)) common[!is.na(common$xerr.all) & is.finite(common$xerr.all), , drop = FALSE] else common[0, , drop = FALSE]
  p <- ggplot(common, aes(x = channel.run, y = channel.all)) +
    geom_segment(
      data = run_err_df,
      aes(x = channel.run - xerr.run, xend = channel.run + xerr.run, y = channel.all, yend = channel.all),
      inherit.aes = FALSE,
      colour = "grey60",
      alpha = 0.6
    ) +
    geom_segment(
      data = all_err_df,
      aes(x = channel.run, xend = channel.run, y = channel.all - xerr.all, yend = channel.all + xerr.all),
      inherit.aes = FALSE,
      colour = "grey60",
      alpha = 0.6
    ) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("Align", rt, "-> AllCal", sprintf("(%sE%d)", cfg_detector, cfg_crystal)),
      x = "channel (11B run)",
      y = "channel (AllCal)"
    ) +
    theme_minimal()
  safe_rt <- gsub("[^A-Za-z0-9_-]", "_", rt)
  #out_fn <- file.path(outdir, paste0("align_", safe_rt, "_to_AllCal_", cfg_detector, "E", cfg_crystal, ".png"))
  #suppressMessages(suppressWarnings(ggsave(filename = out_fn, plot = p, width = 6.5, height = 4.5, dpi = 150)))
  invisible(NULL)
}

plot_combined <- function(df_combined, spline_df, knots_df, outdir, cfg_detector, cfg_crystal) {
  df_plot <- as_df(df_combined)
  df_plot$yerr_plot <- if ("energy_err_total" %in% names(df_plot)) df_plot$energy_err_total else if ("yerr" %in% names(df_plot)) df_plot$yerr else NA_real_

  df_plot <- df_plot[order(df_plot$channel_corrected, df_plot$energy), ]
  df_plot$run_display <- if ("run_type" %in% names(df_plot)) {
    ifelse(is.na(df_plot$run_type) | df_plot$run_type == "", "AllCal", df_plot$run_type)
  } else {
    rep("AllCal", nrow(df_plot))
  }

  p_linear <- ggplot(df_plot, aes(x = channel_corrected, y = energy, colour = run_display)) +
    geom_point(alpha = 0.7) +
    geom_line(data = df_plot, aes(x = channel_corrected, y = pred), colour = "steelblue", linewidth = 0.9, inherit.aes = FALSE) +
    labs(
      title = paste("Linear energy fit", sprintf("(%sE%d)", cfg_detector, cfg_crystal)),
      x = "Channel (AllCal-aligned)",
      y = "Energy (keV)"
    ) +
    theme_minimal() +
    scale_colour_discrete(name = "Run type")
  if (plot_linear_cal) {
    suppressMessages(suppressWarnings(ggsave(file.path(outdir, paste0(cfg_detector, "E", cfg_crystal, ".lin_cal", ".png")),
                                             p_linear, width = 6.5, height = 4.5, dpi = 150)))
  }

  p_resid <- ggplot(df_plot, aes(x = pred, y = resid)) +
    geom_errorbar(aes(ymin = resid - ifelse(is.na(yerr_plot), 0, yerr_plot), ymax = resid + ifelse(is.na(yerr_plot), 0, yerr_plot)), width = 0, colour = "grey70", alpha = 0.5) +
    geom_errorbar(aes(ymin = resid - ifelse(is.na(resid_sigma), 0, resid_sigma), ymax = resid + ifelse(is.na(resid_sigma), 0, resid_sigma)), width = 0, colour = "navy", alpha = 0.25) +
    geom_point(alpha = 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    labs(
      title = paste("Residuals of linear calibration", sprintf("(%sE%d)", cfg_detector, cfg_crystal)),
      x = "Predicted energy (keV)",
      y = "Residual (keV)"
    ) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    theme_minimal()

  if (!is.null(spline_df)) {
    if (all(is.na(spline_df$se))) {
      p_resid <- p_resid + geom_line(data = spline_df, aes(x = x, y = y), colour = "steelblue", inherit.aes = FALSE)
    } else {
      spline_df_extended <- merge(spline_df, unique(data.frame(x = df_plot$pred, resid_sigma = df_plot$resid_sigma)), by = "x", all.x = TRUE, sort = FALSE)
      if (!"resid_sigma" %in% names(spline_df_extended)) {
        spline_df_extended$resid_sigma <- NA_real_
      }
      if (all(is.na(spline_df_extended$resid_sigma))) {
        spline_df_extended$resid_sigma <- mean(df_plot$resid_sigma, na.rm = TRUE)
      } else {
        missing_idx <- which(is.na(spline_df_extended$resid_sigma))
        if (length(missing_idx)) {
          spline_df_extended$resid_sigma[missing_idx] <- mean(df_plot$resid_sigma, na.rm = TRUE)
        }
      }
      band_specs <- data.frame(
        mult = c(3, 2, 1),
        fill = c("#0b3c68", "#1f6db2", "#74a9d8"),
        alpha = c(0.12, 0.18, 0.28)
      )
      for (i in seq_len(nrow(band_specs))) {
        spec <- band_specs[i, ]
        band_df <- spline_df_extended
        total_se <- sqrt(pmax(band_df$se^2 + band_df$resid_sigma^2, 0))
        band_df$ymin <- band_df$y - spec$mult * total_se
        band_df$ymax <- band_df$y + spec$mult * total_se
        band_df <- band_df[is.finite(band_df$ymin) & is.finite(band_df$ymax), , drop = FALSE]
        if (nrow(band_df) == 0) next
        p_resid <- p_resid +
          geom_ribbon(data = band_df,
                      aes(x = x, ymin = ymin, ymax = ymax),
                      fill = spec$fill,
                      alpha = spec$alpha,
                      inherit.aes = FALSE)
      }
      p_resid <- p_resid +
        geom_line(data = spline_df, aes(x = x, y = y), colour = "steelblue", inherit.aes = FALSE)
      if (!is.null(knots_df) && nrow(knots_df) > 0) {
        p_resid <- p_resid +
          geom_point(data = knots_df, aes(x = x_energy, y = y_residual),
                     colour = "red", size = 3, inherit.aes = FALSE)
      }
    }
  }

  if ("spline_pred" %in% names(df_plot)) {
    df_plot$resid_minus_spline <- df_plot$resid - df_plot$spline_pred
    df_plot$combo_se <- sqrt(pmax(ifelse(is.na(df_plot$resid_sigma), 0, df_plot$resid_sigma)^2, 0))
    p_resid_spline <- ggplot(df_plot, aes(x = pred, y = resid_minus_spline)) +
      geom_point(alpha = 0.65) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
      labs(
        title = paste("Residuals after spline subtraction", sprintf("(%sE%d)", cfg_detector, cfg_crystal)),
        x = "Predicted energy (keV)",
        y = "Residual - spline (keV)"
      ) +
      coord_cartesian(ylim = c(-2.5, 2.5)) +
      theme_minimal()
    if (any(!is.na(df_plot$combo_se))) {
      p_resid_spline <- p_resid_spline +
        geom_errorbar(aes(ymin = resid_minus_spline - 1.96 * combo_se, ymax = resid_minus_spline + 1.96 * combo_se),
                      width = 0, colour = "grey70", alpha = 0.5)
    }
    suppressMessages(suppressWarnings(ggsave(file.path(outdir, paste0(cfg_detector, "E", cfg_crystal, ".spline_res", ".png")),
                                             p_resid_spline, width = 6.5, height = 4.5, dpi = 150)))
  }

  suppressMessages(suppressWarnings(ggsave(file.path(outdir, paste0(cfg_detector, "E", cfg_crystal, ".lin_res", ".png")),
                                           p_resid, width = 6.5, height = 4.5, dpi = 150)))
  invisible(NULL)
}

save_spline_knots <- function(spline_model, final_fit, outdir, cfg_detector, cfg_crystal) {
  if (is.null(spline_model)) {
    warning("Cannot save spline knots: spline_model is NULL")
    return(invisible(NULL))
  }

  smooth_object <- spline_model$smooth[[1]]
  if (is.null(smooth_object)) {
    warning("Cannot extract smooth object from GAM model")
    return(invisible(NULL))
  }

  knot_x <- smooth_object$xp
  if (is.null(knot_x) || length(knot_x) == 0) {
    warning("No knot positions available in smooth object")
    return(invisible(NULL))
  }

  knot_y <- predict(spline_model, newdata = data.frame(x = knot_x), type = "response")

  knots_df <- data.frame(
    x_energy = round(knot_x, 3),
    y_residual = round(knot_y, 3)
  )

  out_file <- file.path(outdir, paste0(cfg_detector, "E", cfg_crystal, ".cal_params", ".txt"))

  # Open file for writing
  file_conn <- file(out_file, "w")
  # Write linear fit parameters as header
  if (!is.null(final_fit)) {
    coeffs <- coef(final_fit)
    intercept <- coeffs[1]
    slope <- coeffs[2]
    writeLines(paste("# Linear fit parameters"), file_conn)
    writeLines(paste(intercept, slope), file_conn)
    writeLines(paste("# Spline knots: ", nrow(knots_df), sep = ""), file_conn)
  }

  close(file_conn)

  # Append knots data (without column names to avoid duplication)
  write.table(knots_df, file = out_file, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t", append = TRUE)
  message("Spline knots saved to: ", out_file)
  invisible(NULL)
}

save_models <- function(alignments, final_fit, spline_model, outdir, cfg_detector, cfg_crystal) {
  saveRDS(
    list(alignments = alignments, final_fit = final_fit, spline_model = spline_model),
    file = file.path(outdir, paste0(cfg_detector, "E", cfg_crystal, ".models", ".rds"))
  )
}

## ------------------------- Main workflow --------------------------------
run_one <- function(cfg_detector, cfg_crystal, res_all) {

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  message(sprintf("\n========== Processing %sE%d ==========", cfg_detector, cfg_crystal))

  allcal_points <- subset(res_all$points, run_type == "AllCal" & detector == cfg_detector & crystal == cfg_crystal)
  if (!nrow(allcal_points)) {
    message(sprintf("No AllCal peaks found for %sE%d — skipping", cfg_detector, cfg_crystal))
    return(invisible(NULL))
  }
  allcal_points <- dedupe_peaks(allcal_points)
  allcal_points <- assign_channel_error(allcal_points)
  allcal_points$channel_corrected <- allcal_points$channel
  allcal_points$channel_corrected_err <- if ("xerr" %in% names(allcal_points)) allcal_points$xerr else NA_real_

  seen_keys <- energy_key(allcal_points$energy)
  combined <- allcal_points

  runs_11B <- subset(res_all$summary, grepl("11B", run_type) & detector == cfg_detector & crystal == cfg_crystal)
  runs_11B <- unique(runs_11B$run_type)
  runs_11B <- runs_11B[order(runs_11B)]

  alignments <- list()

  for (rt in runs_11B) {
    pts <- subset(res_all$points, run_type == rt & detector == cfg_detector & crystal == cfg_crystal)
    pts <- dedupe_peaks(pts)
    pts <- assign_channel_error(pts)
    if (nrow(pts) < 2) {
      message(sprintf("Skipping %s: insufficient peaks", rt))
      next
    }

    common <- select_common_peaks(pts, allcal_points)
    if (nrow(common) < min_matches) {
      message(sprintf("Skipping %s: only %d matching peaks (need >= %d)", rt, nrow(common), min_matches))
      next
    }

    fit <- fit_alignment(common)
    alignments[[rt]] <- list(model = fit, matches = common)

    corrected <- apply_correction(pts, fit)
    corrected$run_type <- rt

    res_append <- append_unique_peaks(combined, corrected, seen_keys)
    combined <- res_append$data
    seen_keys <- res_append$seen

    if (res_append$added > 0) {
      message(sprintf("Added %d unique peaks from %s", res_append$added, rt))
    } else {
      message(sprintf("No new unique peaks contributed by %s", rt))
    }

    plot_per_run(common, rt, outdir, cfg_detector, cfg_crystal)
  }

  combined <- as_df(combined)
  if (!"channel_corrected" %in% names(combined)) combined$channel_corrected <- combined$channel

  cat("\nAllCal peaks (deduplicated) for ", cfg_detector, "E", cfg_crystal, ":\n", sep = "")
  print(allcal_points[order(allcal_points$energy), c("energy", "channel", "channel_corrected")])

  cat("\nCombined peak list prior to linear fit:\n")
  combined_print <- combined
  combined_print <- combined_print[order(combined_print$energy, combined_print$run_type %||% ""), ]
  print(combined_print[, intersect(c("run_type", "energy", "channel", "channel_corrected"), names(combined_print))])

  final_fit <- fit_global_cal(combined)
  summary_final <- summary(final_fit)
  cat("Final fit coefficients:\n")
  print(coef(final_fit))
  cat("\nR-squared:", summary_final$r.squared, "\n")

  pred_final <- predict(final_fit, newdata = combined, se.fit = TRUE)
  combined$pred <- pred_final$fit
  combined$resid <- combined$energy - combined$pred

  slope_final <- unname(if (!is.null(coef(final_fit)["channel_corrected"])) coef(final_fit)["channel_corrected"] else coef(final_fit)[2])
  channel_var <- if ("channel_corrected_err" %in% names(combined)) pmax(combined$channel_corrected_err^2, 0) else 0
  meas_energy_var <- if ("yerr" %in% names(combined)) pmax(ifelse(is.na(combined$yerr), 0, combined$yerr)^2, 0) else 0
  fit_var <- pmax(pred_final$se.fit^2, 0)
  slope_contrib <- ifelse(is.na(slope_final), 0, (slope_final^2) * channel_var)
  combined$energy_err_total <- sqrt(pmax(meas_energy_var + fit_var + slope_contrib, 0))
  combined$energy_err_total[!is.finite(combined$energy_err_total)] <- NA_real_

  spline_res <- fit_spline(combined)
  if (!is.null(spline_res$sp_data_y)) {
    combined$spline_pred <- spline_res$sp_data_y
    combined$spline_se <- spline_res$sp_data_se
    combined$resid_sigma <- spline_res$resid_sigma
  }

  knots_df <- NULL
  if (!is.null(spline_res$spline_model)) {
    smooth_object <- spline_res$spline_model$smooth[[1]]
    if (!is.null(smooth_object)) {
      knot_x <- smooth_object$xp
      if (!is.null(knot_x) && length(knot_x) > 0) {
        knot_y <- predict(spline_res$spline_model, newdata = data.frame(x = knot_x), type = "response")
        knots_df <- data.frame(x_energy = knot_x, y_residual = knot_y)
      }
    }
  }

  plot_combined(combined, spline_res$df_smooth, knots_df, outdir, cfg_detector, cfg_crystal)

  save_spline_knots(spline_res$spline_model, final_fit, outdir, cfg_detector, cfg_crystal)
  save_models(alignments, final_fit, spline_res$spline_model, outdir, cfg_detector, cfg_crystal)
  invisible(NULL)
}

main <- function() {
  res_all <- load_data()
  for (det in all_detectors) {
    for (crys in all_crystals) {
      tryCatch(
        run_one(cfg_detector = det, cfg_crystal = crys, res_all = res_all),
        error = function(e) message(sprintf("Error processing %sE%d: %s", det, crys, e$message))
      )
    }
  }
  message("\nAll detector/crystal combinations complete.")
}

# Run main when invoked non-interactively (Rscript)
main()
