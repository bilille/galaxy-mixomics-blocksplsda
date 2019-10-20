#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Get parameters ##
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Compute the matCorEtBlockSelect and addVariablesReponses functions')

parser$add_argument('--input_rdata', dest='input_rdata', required=TRUE, help="Input RData file")
parser$add_argument('--cutoff_comp', dest='cutoff_comp', type='double', required=TRUE, help="")
parser$add_argument('--mat_block_Y_file', dest='mat_block_Y_file', required=TRUE, help="Matrix block Y filepath")
parser$add_argument('--output_rdata', dest='output_rdata', required=TRUE, help="Output RData file")
parser$add_argument('--output_blocks_comb', dest='output_blocks_comb', required=TRUE, help="Output blocks combinations file")

args <- parser$parse_args()

## Print parameters
print("Input RData:")
print(args$input_rdata)
print("Cutoff comp:")
print(args$cutoff_comp)
print("Mat Block Y file:")
print(args$mat_block_Y_file)
print("Output RData:")
print(args$output_rdata)
print("Output Blocks combinations:")
print(args$output_blocks_comb)

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

load(args$input_rdata)

liste_matCorEtBlockSelect = matCorEtBlockSelect(res_block_splsda = mixomics_result,
                                                cutoff_comp = args$cutoff_comp,
                                                comp = 1:2)

mat_cor_comp1 = liste_matCorEtBlockSelect$mat_cor_comp1
mat_cor_comp2 = liste_matCorEtBlockSelect$mat_cor_comp2
dataframe_cor_comp_var_global = liste_matCorEtBlockSelect$dataframe_cor_comp_var_global
liste_vec_indice_blockSelect = liste_matCorEtBlockSelect$liste_vec_indice_blockSelect
liste_vec_blocks = liste_matCorEtBlockSelect$liste_vec_blocks

print("Mat cor comp1:")
print(mat_cor_comp1)
print("Mat cor comp2:")
print(mat_cor_comp2)
print("dataframe_cor_comp_var_global:")
print(dataframe_cor_comp_var_global)
print("liste_vec_indice_blockSelect:")
print(liste_vec_indice_blockSelect)
print("liste_vec_blocks:")
print(liste_vec_blocks)


lapply(liste_vec_blocks, write, file=args$output_blocks_comb,  append=TRUE, ncolumns=100, sep=",")


print("Reading Mat Block Y")
mat_block_Y = read.table(args$mat_block_Y_file, header=TRUE, row.names=1)
print(mat_block_Y)


liste_dataframe_cor_comp_var_global = addVariablesReponses(res_block_splsda = mixomics_result,
                                                           dataframe_cor_comp_var_global = dataframe_cor_comp_var_global,
                                                           liste_vec_indice_blockSelect = liste_vec_indice_blockSelect,
                                                           mat_block_Y = mat_block_Y)


save(mixomics_result,
     liste_matCorEtBlockSelect,
     mat_cor_comp1,
     mat_cor_comp2,
     dataframe_cor_comp_var_global,
     liste_vec_indice_blockSelect,
     liste_vec_blocks,
     liste_dataframe_cor_comp_var_global,
     file = args$output_rdata)
