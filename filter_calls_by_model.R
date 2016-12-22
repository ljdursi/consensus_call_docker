#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("Usage: filter_calls_by_model model_file input.vcf output.vcf [threshold]")
} else {
    model.file <- args[1] 
    in.vcf <- args[2]    
    out.vcf <- args[3]    
    if (length(args) > 3) {
        thresh <- as.numeric(args[4])
    } else {
        thresh <- 0.66
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

source('analysis/estimate_accuracies.R')
source('analysis/apply_model.R')

filterVcfByModel(in.vcf, out.vcf, model, formula, predict_function, thresh)
