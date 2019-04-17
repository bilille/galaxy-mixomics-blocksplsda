#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Main Function ##

library(argparse)

parser <- ArgumentParser(description='Run the mixOmics block.splsda function')

parser$add_argument('--block', dest='blocks_list', nargs=2, action="append", required=TRUE, help="Block file")
parser$add_argument('--samples', dest='samples_file', required=TRUE, help="Samples description file")
parser$add_argument('--ncomp', dest='ncomp', type='integer', required=TRUE, help="Samples description file")

args <- parser$parse_args()

##
print(args$blocks_list)
print(args$samples_file)
print(args$ncomp)

# loading libraries
require(mixOmics)

list_X <- c()

summary(args$blocks_list)

print(args$blocks_list[3,2])

for(i in 1:nrow(args$blocks_list))
{
    print(i)
    list_X[[args$blocks_list[i,1]]] <- read.table(args$blocks_list[i,2], sep='\t', header=TRUE, row.names=1)
}

print(list_X)

samples_descriptions <- read.table(args$samples_file, sep='\t', header=TRUE, row.names=1)

print(samples_descriptions)

Y <- factor(samples_descriptions[[1]])

print(Y)

block_nb <- nrow(args$blocks_list)

design <- matrix(0, nrow = block_nb, ncol = block_nb)

print(design)

res_block_splsda <- block.splsda(X = list_X,
                                 Y = Y,
                                 ncomp = args$ncomp,
                                 design = design)

print(res_block_splsda)























