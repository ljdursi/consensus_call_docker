require(VariantAnnotation)
require(glmnet)

vcfToDataFrame <- function(vcf) {
  callers <- sapply(vcfInfo(vcf)$Callers, function(x) paste(x, collapse=" "))
  
  broad_snowman <- rep(0, length(callers))
  smufin <- rep(0, length(callers))
  sanger <- rep(0, length(callers))
  dkfz <- rep(0, length(callers))

  broad_snowman[grep("broad", callers)] <- 1
  smufin[grep("smufin", callers)] <- 1
  sanger[grep("sanger", callers)] <- 1
  dkfz[grep("dkfz", callers)] <- 1

  df <- data.frame(broad_snowman, smufin, sanger, dkfz)
  df$broad_snowman <- as.integer(broad_snowman)
  df$smufin <- as.integer(smufin)
  df$sanger <- as.integer(sanger)
  df$dkfz <- as.integer(dkfz)
  
  ncore <- df$smufin + df$sanger + df$dkfz 
  df$union      <- ifelse(ncore >= 1, 1, 0)
  df$intersect2 <- ifelse(ncore >= 2, 1, 0)
  df$intersect3 <- ifelse(ncore >= 3, 1, 0)
  
  df$wgs_nvardepth <- vcfInfo(vcf)$NormalVarDepth
  df$wgs_tvardepth <- vcfInfo(vcf)$TumorVarDepth
  df$wgs_ndepth <- as.integer(vcfInfo(vcf)$NormalTotalDepth)
  df$wgs_tdepth <- as.integer(vcfInfo(vcf)$TumorTotalDepth)
  df$wgs_nvaf <- vcfInfo(vcf)$NormalVAF
  df$wgs_tvaf <- vcfInfo(vcf)$TumorVAF
  
  df$dbsnp <- as.numeric(ifelse(is.na(vcfInfo(vcf)$dbsnp), 0, 1))
  df$cosmic <- as.numeric(ifelse(is.na(vcfInfo(vcf)$cosmic), 0, 1))
  df$repeat_masker <- as.numeric(ifelse(is.na(vcfInfo(vcf)$repeat_masker), 0, 1))
  df$thousand_genomes <- as.numeric(ifelse(is.na(vcfInfo(vcf)[["1000genomes_ID"]]), 0, 1))
  
  df$repeat_count <- as.numeric(vcfInfo(vcf)$RepeatRefCount)
  
  reflen <- sapply(ref(vcf), length)
  altlen <- sapply(alt(vcf), function(x) length(x[[1]]))
  df$varlen <- as.integer(altlen-reflen)
  
  foo <- as.data.frame(rowRanges(vcf))
  df$chrom <- as.character(foo$seqnames)
  df$pos <- foo$start
  df$ref <- foo$REF
  df$alt <- sapply(alt(vcf), function(x) as.character(x[[1]]))
  return(df)
}

filterVcfByModel <- function(infilename, outfilename, model, formula, predict_function, threshold) {
  vcf_records <- readVcf(infilename, 'hg19')
  vcf_df <- vcfToDataFrame(vcf_records)
  
  all.prediction <- rep("NOCALL", nrow(vcf_df))
  prediction.prob <- predict_function(model, formula, vcf_df)
  prediction.prob[is.na(prediction.prob)] <- 0.
  threshold <- as.double(threshold)
  all.prediction <- ifelse(prediction.prob >= threshold, ".", "LOWSUPPORT")
  tstdf <- data.frame(prediction.prob, all.prediction) 
  if (class(all.prediction) == "matrix") {
    all.prediction <- all.prediction[,1]
  }

  filt(vcf_records) <- all.prediction
  
  rownames(vcf_records) <- rep(".", nrow(vcf_records))
  vcfInfo(vcf_records)$model_score <- prediction.prob
  writeVcf(vcf_records, outfilename, index=FALSE)
}
