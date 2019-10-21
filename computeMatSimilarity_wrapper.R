#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Get parameters ##
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the computeMatSimilarity function')

parser$add_argument('--input_rdata', dest='input_rdata', required=TRUE, help="Input RData file")
parser$add_argument('--output_rdata', dest='output_rdata', required=TRUE, help="Output RData file")

args <- parser$parse_args()

## Print parameters
print("Input RData:")
print(args$input_rdata)
print("Output RData:")
print(args$output_rdata)

## Loading libraries
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

liste_res_matSimilarity_group = compute_matSimilarity(liste_dataframe_cor_comp_var_global)

liste_matSimilarity_group = liste_res_matSimilarity_group$liste_matSimilarity_group
comp = liste_res_matSimilarity_group$comp

save(liste_matSimilarity_group,
     comp,
     file = args$output_rdata)