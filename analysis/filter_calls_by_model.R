#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("Usage: filter_calls_by_model model_file input.vcf output.vcf rpath [threshold]")
} else {
    model.file <- args[1] 
    in.vcf <- args[2]    
    out.vcf <- args[3]    
    rpath <- args[4]
    if (length(args) > 4) {
        thresh <- as.numeric(args[5])
    } else {
        thresh <- 0.71
    }
}

load(model.file)

if (!exists("model")) {
    stop("Model file must contain a model named model")
}
if (!exists("predict_function")) {
    stop("Model file must contain a prediction function named predict_function")
}
if (!exists("formula")) {
    stop("Model file must contain a formula named formula")
}

source(paste0(rpath, 'estimate_accuracies.R'))
source(paste0(rpath, 'apply_model.R'))

filterVcfByModel(in.vcf, out.vcf, model, formula, predict_function, thresh)
