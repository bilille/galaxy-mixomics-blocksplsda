#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Main Function ##

suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the mixOmics block.splsda function')

parser$add_argument('--block', dest='blocks_list', nargs=2, action="append", required=TRUE, help="Block file")
parser$add_argument('--samples', dest='samples_file', required=TRUE, help="Samples description file")
parser$add_argument('--ncomp', dest='ncomp', type='integer', required=TRUE, help="Number of components to include in the model")
parser$add_argument('--correlation', dest='correlation', action="store_true", help="Add correlation between all blocks")
parser$add_argument('--scheme', dest='scheme', required=TRUE, help="Scheme")
parser$add_argument('--mode', dest='mode', required=TRUE, help="Mode")
parser$add_argument('--maxiter', dest='maxiter', type='integer', required=TRUE, help="Maximum number of iterations")
parser$add_argument('--outrdata', dest='output_rdata', required=TRUE, help="Output Rdata file")

args <- parser$parse_args()

##
print("Blocks:")
print(args$blocks_list)
print("Sample description file:")
print(args$samples_file)
print("Number of components:")
print(args$ncomp)
print("Compute correlation between all blocks:")
print(args$correlation)
print("Scheme:")
print(args$scheme)
print("Mode:")
print(args$mode)
print("Max nb of iterations:")
print(args$maxiter)
print("Output Rdata file:")
print(args$output_rdata)

# loading libraries
suppressPackageStartupMessages(require(mixOmics))

list_X <- c()

for(i in 1:nrow(args$blocks_list))
{
    list_X[[args$blocks_list[i,1]]] <- read.table(args$blocks_list[i,2], sep='\t', header=TRUE, row.names=1)
}

# print(list_X)

samples_descriptions <- read.table(args$samples_file, sep='\t', header=TRUE, row.names=1)

print("Samples descriptions matrix:")
print(samples_descriptions)

Y <- factor(samples_descriptions[[1]])

print("Y factor matrix:")
print(Y)

block_nb <- nrow(args$blocks_list)

design <- matrix(0, nrow = block_nb, ncol = block_nb)

if(args$correlation)
{
    design <- matrix(1, nrow = block_nb, ncol = block_nb)
    diag(design) <- 0
}

print("Design matrix:")
print(design)

res_block_splsda <- block.splsda(X = list_X,
                                 Y = Y,
                                 ncomp = args$ncomp,
                                 design = design,
                                 scheme = args$scheme,
                                 mode = args$mode,
                                 scale = TRUE,
                                 init = "svd",
                                 tol = 1e-06,
                                 max.iter = args$maxiter,
                                 near.zero.var = FALSE,
                                 all.outputs = TRUE)

print("Block.splsda object:")
print(res_block_splsda)

save(res_block_splsda, file=args$output_rdata)