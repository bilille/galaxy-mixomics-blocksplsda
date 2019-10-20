#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Get parameters ##
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the networkVar function')

parser$add_argument('--mat_similarity_rdata', dest='mat_similarity_rdata', required=TRUE, help="matSimilarity RData file")
parser$add_argument('--block_splsda_rdata', dest='block_splsda_rdata', required=TRUE, help="block.splsda RData file")
parser$add_argument('--var_list_file', dest='var_list_file', required=TRUE, help="Variables list file")
parser$add_argument('--blocks_vec', dest='blocks_vec', required=TRUE, help="Blocks vector")
parser$add_argument('--responses_var', dest='responses_var', required=TRUE, help="Responses variables")
parser$add_argument('--output_graph', dest='output_graph', required=TRUE, help="Output graphml")

args <- parser$parse_args()

## Print parameters
print("matSimilarity RData file:")
print(args$mat_similarity_rdata)
print("block.splsda RData file:")
print(args$block_splsda_rdata)
print("Variables list file:")
print(args$var_list_file)
print("Blocks vector:")
print(args$blocks_vec)
print("Response variables:")
print(args$responses_var)
print("Output graphml:")
print(args$output_graph)

## Loading libraries
suppressPackageStartupMessages(require(mixOmics))

# R script call
source_local <- function(fname)
{
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
## Loading local functions
source_local("Integration_block_splsda_fonc.R")

load(args$mat_similarity_rdata)
load(args$block_splsda_rdata)

blocks_vector = strsplit(args$blocks_vec, ",")
response_variables = strsplit(args$responses_var, ",")

network = networkVariableSelect(liste_matSimilarity_group = liste_matSimilarity_group,
                                comp = comp,
                                res_block_splsda = mixomics_result,
                                cutoff_comp = 0.8,
                                vec_varBlock = blocks_vector,
                                vec_varRep = response_variables)

write.graph(network$gR, file = args$output_graph, format = "graphml")