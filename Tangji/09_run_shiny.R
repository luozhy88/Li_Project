#!/usr/bin/env Rscript
# ============================================================================
# 09_run_shiny.R
# ThyGlycoPortal 启动脚本
# Usage: Rscript 09_run_shiny.R [port]
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
port <- ifelse(length(args) >= 1, as.integer(args[1]), 3838)

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required. Please install it first.")
}

# Check if database exists
if (!file.exists("output/thyroid_glyco_db.sqlite")) {
  stop("Database not found at output/thyroid_glyco_db.sqlite. Please run 05_create_database.py and 06_import_glycomics_data.py first.")
}

cat("============================================================\n")
cat("  ThyGlycoPortal - Thyroid Cancer Glycosylation Portal\n")
cat("============================================================\n")
cat("Starting Shiny app on port:", port, "\n")
cat("Open your browser at: http://localhost:", port, "\n", sep = "")
cat("Press Ctrl+C to stop\n")
cat("============================================================\n\n")

# Launch the app
shiny::runApp("app.R", host = "0.0.0.0", port = port, launch.browser = FALSE)
