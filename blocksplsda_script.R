#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Main Function ##

library(argparse)

parser <- ArgumentParser(description='Run the mixOmics block.splsda function')

parser$add_argument('--block', dest='blocks_list', action="append", required=TRUE, help="Block file")
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

for(i in 1:length(args$blocks_list))
{
    list_X[[i]] <- read.table(args$blocks_list[i], sep='\t')
}

print(list_X)


























