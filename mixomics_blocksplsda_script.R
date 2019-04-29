#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Main Function ##

suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description='Run the mixOmics block.splsda function')

parser$add_argument('--block', dest='blocks_list', nargs=4, action="append", required=TRUE, help="Block name + nb variables to select + data matrix file + variables metadata file")
parser$add_argument('--sample_metadata_in', dest='sample_metadata_in', required=TRUE, help="Samples metadata file")
parser$add_argument('--sample_description_col', dest='sample_description_col', type='integer', required=TRUE, help="Sample description column number")
parser$add_argument('--ncomp', dest='ncomp', type='integer', required=TRUE, help="Number of components to include in the model")
parser$add_argument('--correlation', dest='correlation', action="store_true", help="Add correlation between all blocks")
parser$add_argument('--scheme', dest='scheme', required=TRUE, help="Scheme")
parser$add_argument('--mode', dest='mode', required=TRUE, help="Mode")
parser$add_argument('--maxiter', dest='maxiter', type='integer', required=TRUE, help="Maximum number of iterations")
parser$add_argument('--scale', dest='scale', action="store_true", help="Each block is standardized to zero means and unit variances")
parser$add_argument('--init', dest='init', required=TRUE, help="Init (svd or svd.single)")
parser$add_argument('--tol', dest='tol', type='double', required=TRUE, help="Convergence stopping value")
parser$add_argument('--nearzerovar', dest='nearzerovar', action="store_true", help="Should be set in particular for data with many zero values")
parser$add_argument('--rdata_out', dest='rdata_out', required=TRUE, help="Output Rdata file")
parser$add_argument('--sample_metadata_out', dest='sample_metadata_out', required=TRUE, help="Output sample metadata file")
parser$add_argument('--variable_metadata_outdir', dest='variable_metadata_outdir', required=TRUE, help="Output variable metadata directory")

args <- parser$parse_args()

##
print("Blocks:")
print(args$blocks_list)
print("Sample metadata file:")
print(args$sample_metadata_in)
print("Sample description column number:")
print(args$sample_description_col)
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
print("Scale:")
print(args$scale)
print("Tol:")
print(args$tol)
print("near.zero.var:")
print(args$nearzerovar)
print("Output Rdata file:")
print(args$rdata_out)
print("Output sample metadata file:")
print(args$sample_metadata_out)
print("Output variable metadata directory:")
print(args$variable_metadata_outdir)

# loading libraries
suppressPackageStartupMessages(require(mixOmics))

list_X <- c()
keepX <- c()

for(i in 1:nrow(args$blocks_list))
{
    block_name <- args$blocks_list[i,1]
    block_keep <- strtoi(args$blocks_list[i,2])
    block_data_matrix <- args$blocks_list[i,3]
    # block_meta_var <- args$blocks_list[i,4]
    list_X[[block_name]] <- t(read.table(block_data_matrix, sep='\t', header=TRUE, row.names=1)) # transpose the matrix
    nb_variables = ncol(list_X[[block_name]])
    if(block_keep > 0)
    {
        keepX[[block_name]] <- rep(block_keep, args$ncomp)
    }
    else
    {
        keepX[[block_name]] <- rep(nb_variables, args$ncomp)
    }
    print(sprintf("Block %s contains %d variables and %d will be selected", block_name, nb_variables, block_keep))
}

# print(list_X)

sample_metadata <- read.table(args$sample_metadata_in, sep='\t', header=TRUE, row.names=1)

print("Sample metadata matrix:")
print(head(sample_metadata))

description_column <- ncol(sample_metadata)
if(args$sample_description_col > 0)
{
    description_column <- args$sample_description_col
}

Y <- factor(sample_metadata[[description_column]])

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
                                 keepX = keepX,
                                 design = design,
                                 scheme = args$scheme,
                                 mode = args$mode,
                                 scale = args$scale,
                                 init = args$init,
                                 tol = args$tol,
                                 max.iter = args$maxiter,
                                 near.zero.var = args$nearzerovar,
                                 all.outputs = TRUE)

print("Block.splsda object:")
print(res_block_splsda)

save(res_block_splsda, file=args$rdata_out)

# print("Block.splsda variates:")
# print(res_block_splsda$variates)

for(bname in names(res_block_splsda$variates))
{
    # print(bname)
    # print(res_block_splsda$variates[[bname]])
    colnames(res_block_splsda$variates[[bname]]) <- paste("block.splsda", bname, gsub(" ", "_", colnames(res_block_splsda$variates[[bname]])), sep = "_")
    # print(res_block_splsda$variates[[bname]])
    sample_metadata <- cbind2(sample_metadata, res_block_splsda$variates[[bname]])
}

# print(sample_metadata)

write.table(sample_metadata, file = args$sample_metadata_out, quote = TRUE, sep = "\t", row.names = TRUE, col.names = NA)

# print("Block.splsda loadings:")
# print(res_block_splsda$loadings)

for(i in 1:nrow(args$blocks_list))
{
    block_name <- args$blocks_list[i,1]
    # block_keep <- strtoi(args$blocks_list[i,2])
    # block_data_matrix <- args$blocks_list[i,3]
    block_meta_var <- args$blocks_list[i,4]

    meta_variable <- res_block_splsda$loadings[[block_name]]

    # print(head(meta_variable))

    if(block_meta_var != "None")
    {
        input_meta_variable <- read.table(block_meta_var, sep='\t', header=TRUE, row.names=1)
        # print(head(input_meta_variable))
        meta_variable <- cbind2(input_meta_variable, meta_variable)
    }

    print(head(meta_variable))

    block_meta_var_output_filename <- paste("mixomics_blocksplsda_block_", block_name, "_variable_metadata.tsv", sep="")

    write.table(meta_variable, file = paste(args$variable_metadata_outdir,block_meta_var_output_filename, sep='/'), quote = TRUE, sep = "\t", row.names = TRUE, col.names = NA)
}