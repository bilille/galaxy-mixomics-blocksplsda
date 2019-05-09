#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Main Function ##

suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the mixOmics plotIndiv function')

parser$add_argument('--input_rdata', dest='input_rdata', required=TRUE, help="Input RData file")
parser$add_argument('--legend', dest='legend', action="store_true", help="Display the legend")
parser$add_argument('--ellipse', dest='ellipse', action="store_true", help="Plot ellipse plots")
parser$add_argument('--output_pdf', dest='output_pdf', required=TRUE, help="Output PDF file")

args <- parser$parse_args()

##
print("Input RData:")
print(args$input_rdata)
print("Plot legend:")
print(args$legend)
print("Plot ellipse plots:")
print(args$ellipse)
print("Output PDF file:")
print(args$output_pdf)

# loading libraries
suppressPackageStartupMessages(require(mixOmics))

load(args$input_rdata)

pdf(args$output_pdf)

plotIndiv(mixomics_result,
          legend = args$legend,
          ellipse = args$ellipse)

for(k in 1:(length(mixomics_result$names[[3]])-1))
{
    name_block = mixomics_result$names[[3]][k]

    plotIndiv(mixomics_result,
              blocks = k,
              legend = args$legend,
              ellipse = args$ellipse)
}

dev.off()