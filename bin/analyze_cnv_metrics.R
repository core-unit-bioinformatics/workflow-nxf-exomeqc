#!/usr/bin/env Rscript

# Canonical CNV metrics expected from DRAGEN cnv_metrics CSV output.
DEFAULT_CNV_METRICS <- c(
  "sample",
  "NORMAL_SEX",
  "NORMAL_SEX-relative",
  "TUMOR_SEX",
  "TUMOR_SEX-relative",
  "Number of filtered records (total)-relative",
  "OutlierBafFraction",
  "PMAD",
  "Coverage MAD",
  "Median Bin Count",
  "Post-Normalization Bin Count Sigma",
  "Number of segments",
  "Number of passing amplifications",
  "Number of passing amplifications-relative",
  "Number of passing deletions",
  "Number of passing deletions-relative"
)

DEFAULT_VCF_METRICS_RAW <- c(
  "ModelSource",
  "EstimatedTumorPurity",
  "DiploidCoverage",
  "OverallPloidy",
  "HomozygosityIndex",
  "AlternativeModelDedup",
  "AlternativeModelDup"
)

DEFAULT_VCF_METRICS <- c(
  "ModelSource",
  "EstimatedTumorPurity",
  "DiploidCoverage",
  "OverallPloidy",
  "HomozygosityIndex",
  "AlternativeModelDedupPurity",
  "AlternativeModelDedupDiploidCoverage",
  "AlternativeModelDupPurity",
  "AlternativeModelDupDiploidCoverage"
)

# Suggested defaults for DRAGEN CNV QC. Keep configurable for project calibration.
DEFAULT_THRESHOLDS <- list(
  pmad_max = NULL,
  coverage_mad_max = NULL,
  outlier_baf_fraction_max = NULL,
  diploid_coverage_min = NULL,
  purity_min = NULL,
  purity_max = NULL,
  ploidy_min = NULL,
  ploidy_max = NULL,
  post_norm_sigma_max = NULL,
  num_segments_max = NULL,
  filtered_records_relative_max = NULL,
  homozygosity_index_max = NULL,
  total_excluded_bp_max = NULL
)

stopf <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

load_required_packages <- local({
  loaded <- FALSE
  function() {
    if (loaded) {
      return(invisible(TRUE))
    }
    pkgs <- c("karyoploteR", "ggplot2", "GenomicRanges")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stopf("Missing required R packages: %s", paste(missing_pkgs, collapse = ", "))
    }

    suppressPackageStartupMessages(library(karyoploteR))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(GenomicRanges))
    loaded <<- TRUE
    invisible(TRUE)
  }
})

msgf <- function(...) {
  message(sprintf(...))
}

# Parse CLI options in the form --key=value or --key value.
parse_cli_args <- function(args) {
  positional <- character()
  options <- list()
  i <- 1

  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      positional <- c(positional, token)
      i <- i + 1
      next
    }

    keyval <- sub("^--", "", token)
    if (grepl("=", keyval, fixed = TRUE)) {
      split_idx <- regexpr("=", keyval, fixed = TRUE)
      key <- substr(keyval, 1, split_idx - 1)
      value <- substr(keyval, split_idx + 1, nchar(keyval))
      options[[key]] <- value
      i <- i + 1
      next
    }

    next_is_value <- i < length(args) && !startsWith(args[[i + 1]], "--")
    if (next_is_value) {
      options[[keyval]] <- args[[i + 1]]
      i <- i + 2
    } else {
      options[[keyval]] <- "TRUE"
      i <- i + 1
    }
  }

  list(positional = positional, options = options)
}

is_true_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

is_missing_value <- function(x) {
  if (length(x) == 0 || is.null(x)) {
    return(TRUE)
  }
  sx <- trimws(as.character(x))
  length(sx) == 0 || is.na(sx) || tolower(sx) %in% c("na", "nan", ".", "")
}

parse_alternative_model_pair <- function(x) {
  out <- c(Purity = NA_real_, DiploidCoverage = NA_real_)

  if (is_missing_value(x)) {
    return(out)
  }

  sx <- trimws(as.character(x))
  parts <- strsplit(sx, ",", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(out)
  }

  purity <- to_numeric_safe(trimws(parts[[1]]))
  diploid_cov <- to_numeric_safe(trimws(parts[[2]]))

  if (!is.na(purity)) {
    out[["Purity"]] <- purity
  }
  if (!is.na(diploid_cov)) {
    out[["DiploidCoverage"]] <- diploid_cov
  }
  out
}

ensure_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      df[[col]] <- NA
    }
  }
  df[, required_cols, drop = FALSE]
}

# Strictly assert that required columns are present.
require_columns <- function(df, required_cols, context) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stopf(
      "%s is missing required columns: %s",
      context,
      paste(missing_cols, collapse = ", ")
    )
  }
  df[, required_cols, drop = FALSE]
}

parse_thresholds <- function(options) {
  thresholds <- DEFAULT_THRESHOLDS
  for (name in names(DEFAULT_THRESHOLDS)) {
    if (!is.null(options[[name]])) {
      parsed <- to_numeric_safe(options[[name]])
      if (is.na(parsed)) {
        stopf("Invalid numeric value for --%s: %s", name, options[[name]])
      }
      thresholds[[name]] <- parsed
    }
  }

  if (!is.null(thresholds$purity_min) && (thresholds$purity_min < 0 || thresholds$purity_min > 1)) {
    stopf("Invalid purity_min=%s. Must be in [0, 1].", thresholds$purity_min)
  }
  if (!is.null(thresholds$purity_max) && (thresholds$purity_max < 0 || thresholds$purity_max > 1)) {
    stopf("Invalid purity_max=%s. Must be in [0, 1].", thresholds$purity_max)
  }
  if (!is.null(thresholds$purity_min) && !is.null(thresholds$purity_max) && thresholds$purity_min > thresholds$purity_max) {
    stopf("Invalid purity range: purity_min=%s purity_max=%s", thresholds$purity_min, thresholds$purity_max)
  }

  if (!is.null(thresholds$ploidy_min) && thresholds$ploidy_min <= 0) {
    stopf("Invalid ploidy_min=%s. Must be > 0.", thresholds$ploidy_min)
  }
  if (!is.null(thresholds$ploidy_max) && thresholds$ploidy_max <= 0) {
    stopf("Invalid ploidy_max=%s. Must be > 0.", thresholds$ploidy_max)
  }
  if (!is.null(thresholds$ploidy_min) && !is.null(thresholds$ploidy_max) && thresholds$ploidy_min > thresholds$ploidy_max) {
    stopf("Invalid ploidy range: ploidy_min=%s ploidy_max=%s", thresholds$ploidy_min, thresholds$ploidy_max)
  }

  if (!is.null(thresholds$total_excluded_bp_max) && thresholds$total_excluded_bp_max < 0) {
    stopf("Invalid total_excluded_bp_max=%s. Must be >= 0.", thresholds$total_excluded_bp_max)
  }

  thresholds
}

write_plain_csv <- function(df, outfile) {
  if ("sample" %in% names(df)) {
    ordered_cols <- c("sample", setdiff(names(df), "sample"))
    df <- df[, ordered_cols, drop = FALSE]
  }
  utils::write.csv(df, file = outfile, row.names = FALSE, quote = TRUE)
}

build_long_purity_model_df <- function(pur_cov_model_list) {
  long_rows <- lapply(names(pur_cov_model_list), function(sample_name) {
    pur_cov_model <- pur_cov_model_list[[sample_name]]
    max_purity_df <- aggregate(logl ~ purity, pur_cov_model, function(x) max(x, na.rm = TRUE))
    max_purity_df <- max_purity_df[order(max_purity_df$purity), , drop = FALSE]

    sample_max <- max(max_purity_df$logl, na.rm = TRUE)
    data.frame(
      sample = sample_name,
      purity = max_purity_df$purity,
      value = max_purity_df$logl - sample_max,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, long_rows)
  out[order(out$sample, out$purity), c("sample", "purity", "value"), drop = FALSE]
}

escape_json_string <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  gsub('"', '\\\\"', x)
}

format_json_number <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 12)
}

build_purity_model_lineplot_json <- function(pur_cov_model_long_df) {
  sample_names <- unique(pur_cov_model_long_df$sample)
  sample_entries <- lapply(sample_names, function(sample_name) {
    sample_df <- pur_cov_model_long_df[pur_cov_model_long_df$sample == sample_name, c("purity", "value"), drop = FALSE]
    sample_df <- sample_df[is.finite(sample_df$purity) & is.finite(sample_df$value), , drop = FALSE]

    pair_strings <- apply(sample_df, 1, function(row) {
      paste0("[", format_json_number(as.numeric(row[["purity"]])), ",", format_json_number(as.numeric(row[["value"]])), "]")
    })

    paste0('"', escape_json_string(sample_name), '":[', paste(pair_strings, collapse = ","), "]")
  })

  paste0(
    "{\n",
    "  \"id\": \"exomeqc_purcov_maxmodel\",\n",
    "  \"plot_type\": \"linegraph\",\n",
    "  \"data\": {", paste(sample_entries, collapse = ","), "}\n",
    "}\n"
  )
}

write_combined_purity_model_outputs <- function(pur_cov_model_list, plotfile_prefix) {
  if (length(pur_cov_model_list) == 0) {
    return(invisible(NULL))
  }

  purity_long_df <- build_long_purity_model_df(pur_cov_model_list)
  
  #skip saving as CSV file.
  #csv_outfile <- paste0(plotfile_prefix, "purcov_maxmodel.csv")
  #write_plain_csv(purity_long_df, csv_outfile)

  json_outfile <- paste0(plotfile_prefix, "purcov_maxmodel_lineplot.json")
  writeLines(build_purity_model_lineplot_json(purity_long_df), con = json_outfile)

  #invisible(list(csv = csv_outfile, json = json_outfile))
  invisible(list(json = json_outfile))
}

usage <- function() {
  paste(
    "Usage:",
    "  analyze_cnv_metrics.R --dragen=<dragen_output_dir> --prefix=<plotfile_prefix> [--outfile=<outfile_csv>] [--threshold=value ...] [--run-self-tests]",
    "",
    "Threshold options (all optional; omitted thresholds are not QC-checked):",
    "  --pmad_max=<number>",
    "  --coverage_mad_max=<number>",
    "  --outlier_baf_fraction_max=<number>",
    "  --diploid_coverage_min=<number>",
    "  --purity_min=<number>",
    "  --purity_max=<number>",
    "  --ploidy_min=<number>",
    "  --ploidy_max=<number>",
    "  --post_norm_sigma_max=<number>",
    "  --num_segments_max=<number>",
    "  --filtered_records_relative_max=<number>",
    "  --homozygosity_index_max=<number>",
    "  --total_excluded_bp_max=<number>",
    "  --plot-purcov-maxmodel=true|false",
    "  --plot-purcov-allmodel=true|false",
    "",
    "Notes:",
    "  - Required arguments are named-only: --dragen, --prefix.",
    "  - Optional named argument: --outfile.",
    "  - The per-sample max-model purity PNG is enabled by default.",
    "  - Positional arguments are not accepted.",
    "  - EstimatedTumorPurity=NA is always marked FAIL (hard-coded model confidence check).",
    "  - Optional self-tests can be run with --run-self-tests.",
    sep = "\n"
  )
}

# Discover sample-specific files and name entries by sample basename.
get_files <- function(dragen_output_dir, pattern, fail_if_empty = TRUE) {
  files_to_load <- list.files(
    path = dragen_output_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files_to_load) == 0) {
    if (fail_if_empty) {
      stopf("No files found in %s for pattern: %s", dragen_output_dir, pattern)
    }
    return(setNames(character(), character()))
  }

  sample_names <- sub(pattern, "", basename(files_to_load), perl = TRUE)
  if (anyDuplicated(sample_names) > 0) {
    duplicated_names <- unique(sample_names[duplicated(sample_names)])
    stopf(
      "Duplicate sample names after pattern stripping (%s): %s",
      pattern,
      paste(duplicated_names, collapse = ", ")
    )
  }
  names(files_to_load) <- sample_names
  msgf("Discovered %d file(s) for pattern %s", length(files_to_load), pattern)
  msgf("Sample/file map: %s", paste(sprintf("%s=%s", names(files_to_load), files_to_load), collapse = " | "))
  files_to_load
}

# Convert a list of per-sample DRAGEN metric tables into one harmonized row-wise table.
format_metrics_file <- function(list_of_dfs,
                                subset_col = "V1",
                                subset_name = NULL,
                                metric_name_col = "V3",
                                value_col = "V4",
                                perc_col = "V5",
                                remove_samplename_pattern = NULL) {
  vec_list <- lapply(names(list_of_dfs), function(df_name) {
    df <- list_of_dfs[[df_name]]

    if (!is.null(subset_name)) {
      df <- df[df[, subset_col] == subset_name, , drop = FALSE]
    }
    if (nrow(df) == 0) {
      stopf("No metrics left after filtering for sample: %s", df_name)
    }

    metric_names <- as.character(df[, metric_name_col])
    if (anyDuplicated(metric_names) > 0) {
      duplicated_metrics <- unique(metric_names[duplicated(metric_names)])
      stopf("Duplicate metric names for sample %s: %s", df_name, paste(duplicated_metrics, collapse = ", "))
    }

    val_row <- as.list(df[, value_col])
    perc_row <- as.list(df[, perc_col])
    names(val_row) <- metric_names
    names(perc_row) <- paste0(metric_names, "-relative")

    sample_name <- if (!is.null(remove_samplename_pattern)) {
      gsub(remove_samplename_pattern, "", df_name)
    } else {
      df_name
    }

    one_row <- c(list(sample = sample_name), val_row, perc_row)
    as.data.frame(one_row, check.names = FALSE, stringsAsFactors = FALSE)
  })

  all_colnames <- unique(unlist(lapply(vec_list, names)))
  vec_list_fix <- lapply(vec_list, function(x) {
    for (missing_col in setdiff(all_colnames, names(x))) {
      x[[missing_col]] <- NA
    }
    x[, all_colnames, drop = FALSE]
  })

  df_sum <- do.call(rbind, vec_list_fix)
  row.names(df_sum) <- df_sum$sample
  df_sum
}

get_cnv_metrics <- function(dragen_output_dir, cnv_metrics_names) {
  cnv_metricfiles <- get_files(dragen_output_dir, "\\.cnv_metrics\\.csv$", fail_if_empty = TRUE)
  msgf("Loading CNV metric tables...")

  cnv_metrics <- lapply(cnv_metricfiles, function(x) {
    read.csv(x, header = FALSE, row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  })
  msgf("Loaded CNV metric table dimensions: %s", paste(sprintf("%s=%dx%d", names(cnv_metrics), sapply(cnv_metrics, nrow), sapply(cnv_metrics, ncol)), collapse = " | "))

  for (i in names(cnv_metrics)) {
    cnv_metrics[[i]][cnv_metrics[[i]][, "V1"] == "SEX GENOTYPER" & grepl("N", cnv_metrics[[i]][, "V3"]), "V3"] <- "NORMAL_SEX"
    cnv_metrics[[i]][cnv_metrics[[i]][, "V1"] == "SEX GENOTYPER" & grepl("T", cnv_metrics[[i]][, "V3"]), "V3"] <- "TUMOR_SEX"
  }

  cnv_metrics_df <- format_metrics_file(
    list_of_dfs = cnv_metrics,
    remove_samplename_pattern = "\\.cnv_metrics$"
  )
  msgf("Formatted CNV metrics shape: %d x %d", nrow(cnv_metrics_df), ncol(cnv_metrics_df))

  require_columns(cnv_metrics_df, cnv_metrics_names, "CNV metrics")
}

read_vcf_header <- function(vcf_filepath) {
  con <- if (endsWith(vcf_filepath, ".gz")) gzfile(vcf_filepath, "rt") else file(vcf_filepath, "rt")
  on.exit(close(con), add = TRUE)

  vcf_header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || !startsWith(line, "#")) {
      break
    }
    vcf_header <- c(vcf_header, line)
  }
  vcf_header
}

extract_metric_from_vcfheader <- function(vcf_header, metricname, samplename, allow_missing = FALSE) {
  metric_pattern <- paste0("^##", gsub("([][{}()+*^$\\\\.|?])", "\\\\\\1", metricname), "=")
  metricstring <- grep(metric_pattern, vcf_header, value = TRUE, perl = TRUE)

  if (length(metricstring) == 0) {
    if (isTRUE(allow_missing)) {
      return(data.frame(sample = samplename, metric = metricname, value = NA_character_, stringsAsFactors = FALSE))
    }
    stopf("VCF header metric %s is missing in sample %s", metricname, samplename)
  }
  if (length(metricstring) > 1) {
    stopf("Multiple VCF header entries for %s in sample %s", metricname, samplename)
  }

  metric_value <- sub(metric_pattern, "", metricstring[[1]], perl = TRUE)
  data.frame(sample = samplename, metric = metricname, value = metric_value, stringsAsFactors = FALSE)
}

get_vcf_metrics <- function(dragen_output_dir, vcf_metrics_names) {
  cnv_vcffiles <- get_files(dragen_output_dir, "\\.cnv\\.vcf\\.gz$", fail_if_empty = TRUE)
  msgf("Loading VCF headers...")
  cnv_vcf_header <- lapply(cnv_vcffiles, read_vcf_header)
  msgf("Loaded VCF header line counts: %s", paste(sprintf("%s=%d", names(cnv_vcf_header), sapply(cnv_vcf_header, length)), collapse = " | "))

  vcf_metrics_df_list <- lapply(names(cnv_vcf_header), function(sample_name) {
    metrics_df_long <- do.call(rbind, lapply(DEFAULT_VCF_METRICS_RAW, function(vcf_metric_name) {
      extract_metric_from_vcfheader(
        vcf_header = cnv_vcf_header[[sample_name]],
        metricname = vcf_metric_name,
        samplename = sample_name,
        allow_missing = TRUE
      )
    }))

    metrics_df <- reshape(
      metrics_df_long[, c("sample", "metric", "value")],
      idvar = "sample",
      timevar = "metric",
      direction = "wide"
    )
    names(metrics_df) <- sub("^value\\.", "", names(metrics_df))

    dedup_pair <- parse_alternative_model_pair(metrics_df$AlternativeModelDedup)
    dup_pair <- parse_alternative_model_pair(metrics_df$AlternativeModelDup)

    metrics_df$AlternativeModelDedupPurity <- dedup_pair[["Purity"]]
    metrics_df$AlternativeModelDedupDiploidCoverage <- dedup_pair[["DiploidCoverage"]]
    metrics_df$AlternativeModelDupPurity <- dup_pair[["Purity"]]
    metrics_df$AlternativeModelDupDiploidCoverage <- dup_pair[["DiploidCoverage"]]

    metrics_df$AlternativeModelDedup <- NULL
    metrics_df$AlternativeModelDup <- NULL
    metrics_df
  })

  vcf_metrics_df <- do.call(rbind, vcf_metrics_df_list)
  msgf("Formatted VCF metrics shape: %d x %d", nrow(vcf_metrics_df), ncol(vcf_metrics_df))
  require_columns(vcf_metrics_df, c("sample", vcf_metrics_names), "VCF metrics")
}

# Parse excluded intervals and summarize excluded bp by reason per sample.
analyze_excluded_intervals <- function(dragen_output_dir) {
  load_required_packages()

  excluded_intervals_filepaths <- get_files(
    dragen_output_dir,
    pattern = "\\.cnv\\.excluded_intervals\\.bed\\.gz$",
    fail_if_empty = FALSE
  )
  if (length(excluded_intervals_filepaths) == 0) {
    return(data.frame(sample = character(), stringsAsFactors = FALSE))
  }

  excluded_intervals <- lapply(excluded_intervals_filepaths, function(x) {
    read.delim(x, header = FALSE, stringsAsFactors = FALSE)
  })
  msgf("Loaded excluded interval files: %s", paste(sprintf("%s=%dx%d", names(excluded_intervals), sapply(excluded_intervals, nrow), sapply(excluded_intervals, ncol)), collapse = " | "))

  bad_files <- names(excluded_intervals)[sapply(excluded_intervals, ncol) != 4]
  if (length(bad_files) > 0) {
    stopf("Excluded interval files must have exactly 4 columns. Failing files: %s", paste(bad_files, collapse = ", "))
  }

  excluded_intervals_gr <- lapply(excluded_intervals, toGRanges)
  reason_levels <- unique(unlist(sapply(excluded_intervals_gr, function(x) unique(x$V4))))

  get_total_length <- function(gr, samplename, reason_levels) {
    reason_totals <- sapply(reason_levels, function(reason) {
      sum(width(gr[gr$V4 == reason, ]))
    })
    all_values <- c(reason_totals, total_excluded_bp = sum(width(gr)), sample = samplename)
    as.data.frame(as.list(all_values), check.names = FALSE, stringsAsFactors = FALSE)
  }

  interval_lengths <- lapply(names(excluded_intervals_gr), function(sample_name) {
    get_total_length(excluded_intervals_gr[[sample_name]], sample_name, reason_levels)
  })

  out <- do.call(rbind, interval_lengths)
  msgf("Excluded interval summary shape: %d x %d", nrow(out), ncol(out))
  out
}

create_pur_cov_plot <- function(pur_cov_model, samplename, plotfile_prefix, plot_maxmodel = TRUE, plot_allmodel = FALSE) {
  load_required_packages()

  if (nrow(pur_cov_model) == 0) {
    return(invisible(NULL))
  }

  max_idx <- which.max(pur_cov_model$logl)
  titlestring <- paste0(
    "Purity:", pur_cov_model[max_idx, "purity"],
    ", Coverage:", pur_cov_model[max_idx, "coverage"]
  )

  max_pur_df <- aggregate(logl ~ purity, pur_cov_model, function(x) max(x, na.rm = TRUE))
  names(max_pur_df)[2] <- "logL"

  model_max_plot <- ggplot(max_pur_df, mapping = aes(x = purity, y = logL)) +
    geom_point() +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ggtitle(titlestring) +
    theme_bw()

  model_plot <- ggplot(data = pur_cov_model) +
    geom_point(mapping = aes(x = purity, y = logl, col = coverage)) +
    coord_cartesian() +
    theme_bw()

  if (isTRUE(plot_maxmodel)) {
    ggsave(
      filename = paste0(plotfile_prefix, "pur-cov-maxmodel-", samplename, "_mqc.png"),
      plot = model_max_plot
    )
  }
  if (isTRUE(plot_allmodel)) {
    ggsave(
      filename = paste0(plotfile_prefix, "pur-cov-allmodel-", samplename, "_mqc.png"),
      plot = model_plot
    )
  }
}

analyze_purity_cov_model <- function(dragen_output_dir, plotfile_prefix, plot_maxmodel = FALSE, plot_allmodel = FALSE) {
  pur_cov_model_files <- get_files(
    dragen_output_dir = dragen_output_dir,
    pattern = "\\.cnv\\.purity\\.coverage\\.models\\.tsv$",
    fail_if_empty = FALSE
  )
  if (length(pur_cov_model_files) == 0) {
    msgf("No purity/coverage model files found. Skipping purity/coverage plots.")
    return(invisible(NULL))
  }

  pur_cov_model_list <- lapply(pur_cov_model_files, function(pur_cov_model_file) {
    df <- read.delim(pur_cov_model_file, stringsAsFactors = FALSE)
    if (ncol(df) < 3) {
      stopf("Purity/coverage model file has <3 columns: %s", pur_cov_model_file)
    }
    df <- df[, 1:3]
    names(df) <- c("purity", "coverage", "logl")

    if (
      any(df$purity < 0 | df$purity > 1, na.rm = TRUE) ||
      any(df$coverage < 0, na.rm = TRUE)
    ) {
      stopf("Invalid purity/coverage model values in file: %s", pur_cov_model_file)
    }
    df
  })
  msgf("Loaded purity/coverage models: %s", paste(sprintf("%s=%dx%d", names(pur_cov_model_list), sapply(pur_cov_model_list, nrow), sapply(pur_cov_model_list, ncol)), collapse = " | "))

  write_combined_purity_model_outputs(pur_cov_model_list, plotfile_prefix)

  invisible(lapply(names(pur_cov_model_list), function(samplename) {
    create_pur_cov_plot(
      pur_cov_model = pur_cov_model_list[[samplename]],
      samplename = samplename,
      plotfile_prefix = plotfile_prefix,
      plot_maxmodel = plot_maxmodel,
      plot_allmodel = plot_allmodel
    )
  }))
}

create_karyoplot <- function(df_gr, ymin = NULL, ymax = NULL, outfile) {
  load_required_packages()

  if (length(df_gr) == 0 || all(is.na(df_gr$value))) {
    return(invisible(NULL))
  }

  png(filename = outfile, width = 1000, height = 500)
  kp <- plotKaryotype(genome = "hg38", plot.type = 4)
  if (!is.null(ymin) || !is.null(ymax)) {
    kpPoints(kp, data = df_gr, y = df_gr$value, ymin = ymin, ymax = ymax, cex = 0.2)
    kpAxis(kp, ymin = ymin, ymax = ymax, numticks = 10)
  } else {
    minv <- min(df_gr$value, na.rm = TRUE)
    maxv <- max(df_gr$value, na.rm = TRUE)
    kpPoints(kp, data = df_gr, y = df_gr$value, ymin = minv, ymax = maxv, cex = 0.2)
    kpAxis(kp, ymin = minv, ymax = maxv, numticks = 10)
  }
  dev.off()
}

convert_to_GRanges <- function(df) {
  load_required_packages()

  required <- c("contig", "start", "stop", "value")
  if (!all(required %in% names(df))) {
    stopf("Data frame is missing required columns for GRanges conversion: %s", paste(required, collapse = ","))
  }

  makeGRangesFromDataFrame(
    df,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "stop",
    keep.extra.columns = TRUE
  )
}

process_karyotypeplots <- function(dragen_output_dir, plotfile_prefix) {
  tn_files <- get_files(dragen_output_dir, "\\.tn\\.tsv\\.gz$", fail_if_empty = FALSE)
  if (length(tn_files) > 0) {
    names(tn_files) <- paste0(names(tn_files), "-tn")
  }

  baf_files <- get_files(dragen_output_dir, "\\.tumor\\.baf\\.bedgraph\\.gz$", fail_if_empty = FALSE)
  if (length(baf_files) > 0) {
    names(baf_files) <- paste0(names(baf_files), "-tumor.baf")
  }

  count_dfs <- lapply(tn_files, function(x) {
    df <- read.delim(x, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    if (ncol(df) < 5) {
      stopf("TN file has <5 columns: %s", x)
    }
    names(df)[5] <- "value"
    df
  })

  baf_dfs <- lapply(baf_files, function(x) {
    read.delim(
      x,
      col.names = c("contig", "start", "stop", "value"),
      header = FALSE,
      stringsAsFactors = FALSE
    )
  })

  if (length(count_dfs) > 0) {
    msgf("Loaded TN files: %s", paste(sprintf("%s=%dx%d", names(count_dfs), sapply(count_dfs, nrow), sapply(count_dfs, ncol)), collapse = " | "))
  }
  if (length(baf_dfs) > 0) {
    msgf("Loaded BAF files: %s", paste(sprintf("%s=%dx%d", names(baf_dfs), sapply(baf_dfs, nrow), sapply(baf_dfs, ncol)), collapse = " | "))
  }

  counts_grs <- lapply(count_dfs, convert_to_GRanges)
  baf_grs <- lapply(baf_dfs, convert_to_GRanges)

  invisible(lapply(names(counts_grs), function(x) {
    create_karyoplot(
      df_gr = counts_grs[[x]],
      outfile = paste0(plotfile_prefix, "karyoplot-", x, "_mqc.png")
    )
  }))

  invisible(lapply(names(baf_grs), function(x) {
    create_karyoplot(
      df_gr = baf_grs[[x]],
      outfile = paste0(plotfile_prefix, "karyoplot-", x, "_mqc.png"),
      ymin = 0,
      ymax = 1
    )
  }))
}

evaluate_qc_status <- function(all_metrics_df, thresholds) {
  get_num <- function(row, key) {
    if (!(key %in% names(row))) {
      return(NA_real_)
    }
    to_numeric_safe(row[[key]])
  }

  qc_status <- apply(all_metrics_df, 1, function(row) {
    fail_reasons <- character()

    model_source <- if ("ModelSource" %in% names(row)) trimws(as.character(row[["ModelSource"]])) else NA_character_
    estimated_purity_raw <- if ("EstimatedTumorPurity" %in% names(row)) row[["EstimatedTumorPurity"]] else NA_character_
    estimated_purity <- to_numeric_safe(estimated_purity_raw)

    # DRAGEN docs: EstimatedTumorPurity=NA means no confident model.
    if (is_missing_value(estimated_purity_raw) || is.na(estimated_purity)) {
      fail_reasons <- c(fail_reasons, "NoModelFound:EstimatedTumorPurity=NA")
    }
    if (is_missing_value(model_source)) {
      fail_reasons <- c(fail_reasons, "ModelSource:MISSING")
    }

    check_max <- function(metric_name, threshold_value, label) {
      if (is.null(threshold_value) || is.na(threshold_value)) {
        return(invisible(NULL))
      }
      val <- get_num(row, metric_name)
      if (is.na(val)) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":MISSING"))
      } else if (val > threshold_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), ">", threshold_value))
      }
    }
    check_min <- function(metric_name, threshold_value, label) {
      if (is.null(threshold_value) || is.na(threshold_value)) {
        return(invisible(NULL))
      }
      val <- get_num(row, metric_name)
      if (is.na(val)) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":MISSING"))
      } else if (val < threshold_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), "<", threshold_value))
      }
    }
    check_range <- function(metric_name, min_value, max_value, label) {
      if (is.null(min_value) && is.null(max_value)) {
        return(invisible(NULL))
      }
      val <- get_num(row, metric_name)
      if (is.na(val)) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":MISSING"))
      } else if (!is.null(min_value) && val < min_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), "<", min_value))
      } else if (!is.null(max_value) && val > max_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), ">", max_value))
      }
    }

    check_range_purity <- function(metric_name, min_value, max_value, label) {
      if (is.null(min_value) && is.null(max_value)) {
        return(invisible(NULL))
      }
      val <- get_num(row, metric_name)
      if (is.na(val)) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":MISSING"))
      } else if (!is.null(min_value) && val < min_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), "<", min_value))
      } else if (!is.null(max_value) && val > max_value) {
        fail_reasons <<- c(fail_reasons, paste0(label, ":", signif(val, 4), ">", max_value))
      }
    }

    check_max("PMAD", thresholds$pmad_max, "PMAD")
    check_max("Coverage MAD", thresholds$coverage_mad_max, "CoverageMAD")
    check_max("OutlierBafFraction", thresholds$outlier_baf_fraction_max, "OutlierBafFraction")
    check_min("DiploidCoverage", thresholds$diploid_coverage_min, "DiploidCoverage")
    check_range_purity("EstimatedTumorPurity", thresholds$purity_min, thresholds$purity_max, "EstimatedTumorPurity")
    check_range("OverallPloidy", thresholds$ploidy_min, thresholds$ploidy_max, "OverallPloidy")
    check_max("Post-Normalization Bin Count Sigma", thresholds$post_norm_sigma_max, "PostNormBinCountSigma")
    check_max("Number of segments", thresholds$num_segments_max, "NumSegments")
    check_max("Number of filtered records (total)-relative", thresholds$filtered_records_relative_max, "FilteredRecordsRelative")
    check_max("HomozygosityIndex", thresholds$homozygosity_index_max, "HomozygosityIndex")
    check_max("total_excluded_bp", thresholds$total_excluded_bp_max, "total_excluded_bp")

    if (length(fail_reasons) == 0) {
      c(QC_STATUS = "PASS", QC_FAIL_REASONS = "none")
    } else {
      c(QC_STATUS = "FAIL", QC_FAIL_REASONS = paste(unique(fail_reasons), collapse = ";"))
    }
  })

  qc_status_df <- as.data.frame(t(qc_status), stringsAsFactors = FALSE)
  cbind(qc_status_df, all_metrics_df)
}

# Merge all metric sources by sample and keep all samples for QC reporting.
merge_metrics <- function(vcf_metrics_df, cnv_metrics_df, excl_intervals_lengths) {
  merged <- merge(vcf_metrics_df, cnv_metrics_df, by = "sample", all = TRUE, sort = TRUE)
  if (nrow(excl_intervals_lengths) > 0) {
    merged <- merge(merged, excl_intervals_lengths, by = "sample", all = TRUE, sort = TRUE)
  }
  merged
}

run_self_tests <- function() {
  msgf("Running self-tests...")

  parsed <- parse_cli_args(c("--dragen=in", "--outfile=out", "--prefix=pref", "--pmad_max=0.2", "--run-self-tests"))
  stopifnot(length(parsed$positional) == 0)
  stopifnot(parsed$options$dragen == "in")
  stopifnot(parsed$options$outfile == "out")
  stopifnot(parsed$options$prefix == "pref")
  stopifnot(parsed$options$pmad_max == "0.2")
  stopifnot(parsed$options$`run-self-tests` == "TRUE")

  thr <- parse_thresholds(list(pmad_max = "0.2", ploidy_min = "1.2", ploidy_max = "6"))
  stopifnot(isTRUE(all.equal(thr$pmad_max, 0.2)))
  stopifnot(isTRUE(all.equal(thr$ploidy_min, 1.2)))

  hdr <- c("##fileformat=VCFv4.2", "##EstimatedTumorPurity=0.42", "#CHROM")
  m <- extract_metric_from_vcfheader(hdr, "EstimatedTumorPurity", "S1")
  stopifnot(m$value == "0.42")
  missing_metric_error <- tryCatch(
    {
      extract_metric_from_vcfheader(hdr, "DiploidCoverage", "S1")
      FALSE
    },
    error = function(e) TRUE
  )
  stopifnot(missing_metric_error)

  dup_hdr <- c("##PMAD=0.1", "##PMAD=0.2", "#CHROM")
  duplicate_metric_error <- tryCatch(
    {
      extract_metric_from_vcfheader(dup_hdr, "PMAD", "S1")
      FALSE
    },
    error = function(e) TRUE
  )
  stopifnot(duplicate_metric_error)

  mock <- data.frame(
    sample = c("S1", "S2", "S3"),
    PMAD = c("0.1", "0.9", "0.1"),
    `Coverage MAD` = c("0.1", "0.1", "0.1"),
    OutlierBafFraction = c("0.02", "0.02", "0.02"),
    DiploidCoverage = c("40", "40", "40"),
    EstimatedTumorPurity = c("0.5", "0.5", "NA"),
    OverallPloidy = c("2", "2", "2"),
    `Post-Normalization Bin Count Sigma` = c("0.2", "0.2", "0.2"),
    `Number of segments` = c("100", "100", "100"),
    `Number of filtered records (total)-relative` = c("0.2", "0.2", "0.2"),
    HomozygosityIndex = c("0.8", "0.8", "0.8"),
    ModelSource = c("DEPTH+BAF", "DEPTH+BAF", "DEPTH+BAF"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  qc <- evaluate_qc_status(mock, parse_thresholds(list(pmad_max = "0.2")))
  stopifnot(qc$QC_STATUS[qc$sample == "S1"] == "PASS")
  stopifnot(qc$QC_STATUS[qc$sample == "S2"] == "FAIL")
  stopifnot(qc$QC_STATUS[qc$sample == "S3"] == "FAIL")
  stopifnot(grepl("NoModelFound:EstimatedTumorPurity=NA", qc$QC_FAIL_REASONS[qc$sample == "S3"]))

  purity_models <- list(
    S1 = data.frame(purity = c(0.1, 0.2, 0.3), coverage = c(40, 40, 40), logl = c(-10, -5, -8)),
    S2 = data.frame(purity = c(0.1, 0.2), coverage = c(30, 30), logl = c(-2, -1))
  )
  purity_long <- build_long_purity_model_df(purity_models)
  stopifnot(identical(names(purity_long), c("sample", "purity", "value")))
  s1_vals <- purity_long$value[purity_long$sample == "S1"]
  s2_vals <- purity_long$value[purity_long$sample == "S2"]
  stopifnot(max(s1_vals) == 0)
  stopifnot(max(s2_vals) == 0)

  msgf("Self-tests completed successfully.")
}

main <- function() {
  parsed <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  options <- parsed$options

  if (!is.null(options$help) || !is.null(options$h)) {
    cat(usage(), "\n")
    quit(save = "no", status = 0)
  }

  if (!is.null(options$`run-self-tests`) && is_true_flag(options$`run-self-tests`)) {
    run_self_tests()
    quit(save = "no", status = 0)
  }

  if (length(parsed$positional) > 0) {
    stopf("Positional arguments are not supported. Use named args only.\n%s", usage())
  }
  required_named <- c("dragen", "prefix")
  missing_named <- required_named[!vapply(required_named, function(x) !is.null(options[[x]]), logical(1))]
  if (length(missing_named) > 0) {
    stopf("Missing required named argument(s): %s\n%s", paste(missing_named, collapse = ", "), usage())
  }

  dragen_output_dir <- options$dragen
  plotfile_prefix <- options$prefix
  outfile_csv <- if (!is.null(options$outfile) && nzchar(options$outfile)) {
    options$outfile
  } else {
    paste0(plotfile_prefix, "cnv-metrics.csv")
  }
  plot_purcov_maxmodel <- if (!is.null(options[["plot-purcov-maxmodel"]])) {
    is_true_flag(options[["plot-purcov-maxmodel"]])
  } else {
    FALSE
  }
  plot_purcov_allmodel <- !is.null(options[["plot-purcov-allmodel"]]) && is_true_flag(options[["plot-purcov-allmodel"]])

  if (!dir.exists(dragen_output_dir)) {
    stopf("DRAGEN output directory does not exist: %s", dragen_output_dir)
  }

  msgf("Input arguments: dragen=%s prefix=%s", dragen_output_dir, plotfile_prefix)
  thresholds <- parse_thresholds(options)
  enabled_thresholds <- thresholds[!vapply(thresholds, is.null, logical(1))]
  if (length(enabled_thresholds) == 0) {
    msgf("Threshold configuration: none (only hard-coded model confidence check is active)")
  } else {
    msgf("Threshold configuration: %s", paste(sprintf("%s=%s", names(enabled_thresholds), unlist(enabled_thresholds)), collapse = " | "))
  }

  msgf("Get CNV metrics files...")
  cnv_metrics_df <- get_cnv_metrics(
    dragen_output_dir = dragen_output_dir,
    cnv_metrics_names = DEFAULT_CNV_METRICS
  )
  msgf("CNV metrics files loaded.")

  msgf("Load VCF files and extract metrics...")
  vcf_metrics_df <- get_vcf_metrics(
    dragen_output_dir = dragen_output_dir,
    vcf_metrics_names = DEFAULT_VCF_METRICS
  )
  msgf("VCF file metrics loaded.")

  msgf("Create plots for purity/coverage model...")
  analyze_purity_cov_model(
    dragen_output_dir = dragen_output_dir,
    plotfile_prefix = plotfile_prefix,
    plot_maxmodel = plot_purcov_maxmodel,
    plot_allmodel = plot_purcov_allmodel
  )

  msgf("Summarize excluded intervals...")
  excl_intervals_lengths <- analyze_excluded_intervals(dragen_output_dir = dragen_output_dir)

  msgf("Create karyoploteR plots...")
  process_karyotypeplots(dragen_output_dir, plotfile_prefix)

  all_metrics_df <- merge_metrics(vcf_metrics_df, cnv_metrics_df, excl_intervals_lengths)
  msgf("Merged metrics shape: %d x %d", nrow(all_metrics_df), ncol(all_metrics_df))
  all_metrics_df <- evaluate_qc_status(all_metrics_df, thresholds)
  msgf("QC summary: %s", paste(sprintf("%s=%d", names(table(all_metrics_df$QC_STATUS)), as.integer(table(all_metrics_df$QC_STATUS))), collapse = " | "))

  msgf("Write all metrics in output file...")
  write_plain_csv(all_metrics_df, outfile = outfile_csv)
}

main()
