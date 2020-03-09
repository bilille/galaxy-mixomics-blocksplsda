#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Get parameters ##
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the circleCor function')

parser$add_argument('--input_rdata', dest='input_rdata', required=TRUE, help="Input RData file")
parser$add_argument('--blocks_vec', dest='blocks_vec', required=TRUE, help="Blocks vector")
parser$add_argument('--responses_var', dest='responses_var', required=TRUE, help="Responses variables")
parser$add_argument('--x_min', dest='x_min', type='double', required=TRUE, help="X min")
parser$add_argument('--x_max', dest='x_max', type='double', required=TRUE, help="X max")
parser$add_argument('--y_min', dest='y_min', type='double', required=TRUE, help="Y min")
parser$add_argument('--y_max', dest='y_max', type='double', required=TRUE, help="Y max")
parser$add_argument('--output_var', dest='output_var', required=TRUE, help="Output variables file")
parser$add_argument('--output_pdf', dest='output_pdf', required=TRUE, help="Output PDF file")

args <- parser$parse_args()

## Print parameters
print("Input RData:")
print(args$input_rdata)
print("Blocks vector:")
print(args$blocks_vec)
print("Response variables:")
print(args$responses_var)
print("X min:")
print(args$x_min)
print("X max:")
print(args$x_max)
print("Y min:")
print(args$y_min)
print("Y max:")
print(args$y_max)
print("Output variables file:")
print(args$output_var)
print("Output PDF file:")
print(args$output_pdf)

## Loading libraries
suppressPackageStartupMessages(require(ellipse))
suppressPackageStartupMessages(require(grDevices))
suppressPackageStartupMessages(require(RColorBrewer))
# suppressPackageStartupMessages(require(mixOmics))


# R script call
source_local <- function(fname)
{
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
## Loading local functions
source_local("Integration_block_splsda_fonc.R")

load(args$input_rdata)

blocks_vector = strsplit(args$blocks_vec, ",")[[1]]
response_variables = strsplit(args$responses_var, ",")[[1]]


print("liste_vec_indice_blockSelect:")
print(liste_vec_indice_blockSelect)
print("liste_vec_blocks:")
print(liste_vec_blocks)
print("Mat cor comp1:")
print(mat_cor_comp1)
print("blocks_vector:")
print(blocks_vector)


# pdf(args$output_pdf, width=12, height=9)
pdf(args$output_pdf)

mar = c(5.1, 4.1, 4.1, 9.1)
par(mar = mar)

varSelect = circleCor(liste_dataframe_cor_comp_var_global = liste_dataframe_cor_comp_var_global,
                      liste_vec_indice_blockSelect = liste_vec_indice_blockSelect,
                      mat_cor_comp1 = mat_cor_comp1,
                      mat_cor_comp2 = mat_cor_comp2,
                      vec_blocks = blocks_vector,
                      nomsVarReponses = response_variables,
                      min.X = args$x_min,
                      max.X = args$x_max,
                      min.Y = args$y_min,
                      max.Y = args$y_max,
                      cutoff = 0.85,
                      rad.in = 0.5,
                      cex = 0.7,
                      cex_legend = 0.8,
                      pos = c(1.2, 0),
                      pch = 20,
                      inset = c(-0.25, 0))

dev.off()

write(varSelect, file=args$output_var, ncolumns=1)