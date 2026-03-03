## Calibration collection utilities
# Scans a cubix_workspaces-style tree for 70Ge_*_calibration.(dat|func|res)
# and returns a list with tidy data.frames: points (from .dat), fits (from .func),
# residuals (from .res) and a summary table of runs and available files.

parse_filename_meta <- function(filepath) {
  bn <- basename(filepath)
  # Expect names like: 70Ge_T_DEN_calibration.t  where DEN may include crystal (e.g. B1E1)
  m <- regexec("^70Ge_([^_]+)_([^_]+)_calibration\\.(dat|func|res)$", bn, perl = TRUE)
  res <- regmatches(bn, m)[[1]]
  if (length(res) == 0) {
    return(list(run_type = NA_character_, detector = NA_character_, crystal = NA_integer_, file_type = NA_character_))
  }
  run_type <- res[2]
  den <- res[3]
  file_type <- res[4]
  # detector may be letters+digits before an optional "E" crystal number
  crystal <- NA_integer_
  if (grepl("E[0-9]+$", den)) {
    crystal <- as.integer(sub("^.*E", "", den))
    detector <- sub("E[0-9]+$", "", den)
  } else {
    detector <- den
  }
  list(run_type = run_type, detector = detector, crystal = crystal, file_type = file_type)
}

safe_read_table <- function(path, col.names = NULL) {
  tryCatch({
    df <- read.table(path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
    if (!is.null(col.names) && ncol(df) == length(col.names))
      names(df) <- col.names
    df
  }, error = function(e) {
    # try a loose numeric extraction
    lines <- readLines(path, warn = FALSE)
    nums <- lapply(lines, function(L) as.numeric(unlist(regmatches(L, gregexpr("-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?", L, perl = TRUE)))))
    nums <- Filter(function(x) length(x) > 0, nums)
    if (length(nums) == 0) return(data.frame())
    maxlen <- max(sapply(nums, length))
    mat <- do.call(rbind, lapply(nums, function(row) {
      length(row) <- maxlen; row
    }))
    df2 <- as.data.frame(mat, stringsAsFactors = FALSE)
    if (!is.null(col.names) && ncol(df2) == length(col.names))
      names(df2) <- col.names
    df2
  })
}

parse_dat <- function(path) {
  # Read numeric table ignoring comment lines; take first two numeric columns as channel and energy
  df <- safe_read_table(path)
  if (nrow(df) == 0) return(df)
  # keep first two columns (channel, energy) and preserve optional error columns (xerr, yerr)
  if (ncol(df) >= 2) {
    if (ncol(df) >= 4) {
      df2 <- df[, 1:4, drop = FALSE]
      names(df2) <- c("channel", "energy", "xerr", "yerr")
    } else if (ncol(df) == 3) {
      df2 <- df[, 1:3, drop = FALSE]
      names(df2) <- c("channel", "energy", "xerr")
    } else {
      df2 <- df[, 1:2, drop = FALSE]
      names(df2) <- c("channel", "energy")
    }
    # coerce to numeric
    for (nm in names(df2)) df2[[nm]] <- as.numeric(df2[[nm]])
    return(df2)
  }
  # if single column, assume energy only
  if (ncol(df) == 1) {
    names(df) <- "energy"
    df$energy <- as.numeric(df$energy)
  }
  df
}

parse_func <- function(path) {
  # .func files may contain named parameters like a0, a1. Ignore comment lines starting with #.
  lines <- readLines(path, warn = FALSE)
  # try to find named a0/a1
  a0_line <- grep("^\\s*a0\\b", lines, ignore.case = TRUE, value = TRUE)
  a1_line <- grep("^\\s*a1\\b", lines, ignore.case = TRUE, value = TRUE)
  params <- list()
  if (length(a0_line) > 0) {
    val <- as.numeric(unlist(regmatches(a0_line[1], gregexpr("-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?", a0_line[1], perl = TRUE))))
    if (length(val) > 0) params$a0 <- val[1]
  }
  if (length(a1_line) > 0) {
    val <- as.numeric(unlist(regmatches(a1_line[1], gregexpr("-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?", a1_line[1], perl = TRUE))))
    if (length(val) > 0) params$a1 <- val[1]
  }
  # fallback: if no named params found, extract all numeric tokens
  if (length(params) == 0) {
    nums <- as.numeric(unlist(regmatches(paste(lines, collapse = " "), gregexpr("-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?", paste(lines, collapse = " "), perl = TRUE))))
    if (length(nums) == 0) return(list(params = numeric(0), raw = lines))
    params <- as.list(nums)
    names(params) <- paste0("param", seq_along(params))
  }
  list(params = params, raw = lines)
}

parse_res <- function(path) {
  # Residual files contain X value (energy) and Y value (residual) plus optional error columns. # nolint: line_length_linter.
  df <- safe_read_table(path)
  if (nrow(df) == 0) return(df)
  if (ncol(df) >= 2) {
    # map: energy, residual, optional xerr, yerr
    cols <- c("energy", "residual")
    if (ncol(df) >= 3) cols <- c(cols, "xerr")
    if (ncol(df) >= 4) cols <- c(cols, "yerr")
    df2 <- df[, 1:length(cols), drop = FALSE]
    names(df2) <- cols
    # coerce numeric
    for (nm in names(df2)) df2[[nm]] <- as.numeric(df2[[nm]])
    return(df2)
  }
  df
}

collect_calibrations <- function(root = "cubix_workspaces", recursive = TRUE, sample_n = NULL, verbose = TRUE) {
  # list all files then filter by basename to avoid pattern mismatches when full paths
  all_files <- list.files(path = root, recursive = recursive, full.names = TRUE)
  files <- all_files[grepl("_calibration\\.(dat|func|res)$", basename(all_files), perl = TRUE)]
  if (length(files) == 0) {
    if (verbose) message("No calibration files found under ", root)
    return(list(points = data.frame(), fits = data.frame(), residuals = data.frame(), summary = data.frame()))
  }
  if (!is.null(sample_n)) files <- head(files, sample_n)

  points_list <- list()
  fits_list <- list()
  res_list <- list()
  summary_rows <- list()

  for (f in files) {
    meta <- parse_filename_meta(f)
    meta$path <- f
    # record availability
    summary_rows[[length(summary_rows) + 1]] <- data.frame(run_type = meta$run_type, detector = meta$detector, crystal = meta$crystal, file_type = meta$file_type, path = meta$path, stringsAsFactors = FALSE)

    if (meta$file_type == "dat") {
      df <- tryCatch(parse_dat(f), error = function(e) {
        warning("Failed parsing dat: ", f, " : ", e$message); data.frame()
      })
      if (nrow(df) > 0) {
        df$run_type <- meta$run_type; df$detector <- meta$detector; df$crystal <- meta$crystal; df$path <- f
        points_list[[length(points_list) + 1]] <- df
      }
    } else if (meta$file_type == "func") {
      ff <- tryCatch(parse_func(f), error = function(e) {
        warning("Failed parsing func: ", f, " : ", e$message); list(params = numeric(0), raw = character(0))
      })
      fits_list[[length(fits_list) + 1]] <- data.frame(run_type = meta$run_type, detector = meta$detector, crystal = meta$crystal, path = f, params = list(ff$params), raw = paste(ff$raw, collapse = "\\n"), stringsAsFactors = FALSE)
    } else if (meta$file_type == "res") {
      df <- tryCatch(parse_res(f), error = function(e) {
        warning("Failed parsing res: ", f, " : ", e$message); data.frame()
      })
      if (nrow(df) > 0) {
        df$run_type <- meta$run_type; df$detector <- meta$detector; df$crystal <- meta$crystal; df$path <- f
        res_list[[length(res_list) + 1]] <- df
      }
    }
  }

  points <- if (length(points_list) > 0) do.call(rbind, lapply(points_list, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  })) else data.frame()
  fits <- if (length(fits_list) > 0) do.call(rbind, fits_list) else data.frame()
  residuals <- if (length(res_list) > 0) do.call(rbind, lapply(res_list, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  })) else data.frame()
  summary <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()

  # tidy up types
  if ("crystal" %in% names(summary)) summary$crystal <- as.integer(summary$crystal)

  list(points = points, fits = fits, residuals = residuals, summary = summary)
}

# Convenience wrapper for interactive quick inspection
collect_calibrations_example <- function() {
  root <- "cubix_workspaces"
  if (!dir.exists(root)) root <- "."
  collect_calibrations(root = root, sample_n = 20, verbose = TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## End of file