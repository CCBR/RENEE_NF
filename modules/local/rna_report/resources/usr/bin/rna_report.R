#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-m",
  "--rmarkdown",
  type = "character",
  required = TRUE,
  help = "Required input file: parameterized rNA R Markdown file"
)

parser$add_argument(
  "-r",
  "--raw_counts",
  type = "character",
  required = TRUE,
  help = "Required input file: raw counts matrix"
)

parser$add_argument(
  "-t",
  "--tin_counts",
  type = "character",
  required = TRUE,
  help = "Required input file: TIN matrix"
)

parser$add_argument(
  "-q",
  "--qc_table",
  type = "character",
  required = TRUE,
  help = "Required input file: QC metadata table"
)

parser$add_argument(
  "-f",
  "--output_filename",
  type = "character",
  required = FALSE,
  default = "rNA.html",
  help = "Optional output HTML filename"
)

parser$add_argument(
  "-a",
  "--annotate",
  action = "store_true",
  default = FALSE,
  help = "Display sample names in the complex heatmap"
)

args <- parser$parse_args()

rmarkdown::render(
  args$rmarkdown,
  output_file = args$output_filename,
  output_dir = getwd(),
  params = list(
    raw = args$raw_counts,
    tin = args$tin_counts,
    qc = args$qc_table,
    annot = args$annotate
  )
)
