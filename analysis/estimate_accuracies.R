dummy.lines <- function(template, n, vafprobs=c(1,0,0,0,0)) {
  
  tmp <- template[1:n,]
  tmp$chrom <- "chr0"
  tmp$val_tvaf <- -1
  tmp$val_nvaf <- -1
  tmp$val_tdepth <- -1
  tmp$val_ndepth <- -1
  tmp$wgs_tvaf <- -1
  tmp$wgs_nvaf <- -1
  tmp$binned_wgs_tvaf <- sample(levels(template$binned_wgs_tvaf), n, replace=TRUE, prob=vafprobs)

  tmp$repeat_count <- -1
  
  tmp$concordance <- 0
  tmp$adiscan <- 0
  tmp$broad_mutect <- 0
  tmp$dkfz <- 0
  tmp$lohcomplete <- 0
  tmp$mda_hgsc_gatk_muse <- 0
  tmp$oicr_bl <- 0
  tmp$oicr_sga <- 0
  tmp$sanger <- 0
  tmp$smufin <- 0
  tmp$wustl <- 0
  tmp$common_sample <- TRUE
  tmp$union <- 0
  tmp$intersect2 <- 0
  tmp$intersect3 <- 0
  
  tmp
}

made.up.validation.calls <- function(template, passall=TRUE, privateaccuracies=NA, npersample=50, callername="new_mutect", vafprobs=c(1,0,0,0,0)) {
  sample_levels <- levels(template$sample)
  result <- data.frame()
  
  for (isample in seq(length(levels(template$sample)))) {
    sample <- sample_levels[isample]
    tmp <- dummy.lines(template, npersample, vafprobs)
    
    tmp$sample <- sample
    tmp[[callername]] <- 1
    
    if (!is.na(privateaccuracies)) {
      frac <- privateaccuracies[sample,]$precision
      tmp$validate_true <- runif(npersample) > frac
      tmp$status <- ifelse(tmp$validate_true, "PASS", "GERMLINE")
    } else if (passall) {
      tmp$status <- "PASS"
      tmp$validate_true <- TRUE
    } else {
      tmp$status <- "NOTSEEN"
      tmp$validate_true <- FALSE
    }
    
    result <- rbind(result, tmp)
  }
  result
}

estimates <- function(snvs, snv_calls,
                      filename='~/Desktop/otherbroad/new_mutect.csv',
                      newcaller='new_mutect',
                      originalcaller='broad_mutect',
                      shortcaller='mutect') {
  
  privatecaller <- paste0(shortcaller, '_private')
  newcaller <- paste0('new_',shortcaller)
  
  snvs[[privatecaller]] <- snvs[[originalcaller]] & snvs$concordance == 1
  snv_calls[[privatecaller]] <- snv_calls[[originalcaller]] & snv_calls$concordance == 1
  
  old_mutect_private_acc_by_sample <- corrected.accuracies(snvs, snv_calls, privatecaller, combine=c('concordance'))
  results <- corrected.accuracies(snvs, snv_calls, originalcaller, combine=c('concordance'))
  
  results$caller <- paste0(shortcaller,"_original")
  results<- melt(results, id=c('sample','caller'))
  
  dkfz <- corrected.accuracies(snvs, snv_calls, 'dkfz', combine=c('concordance'))
  sanger <- corrected.accuracies(snvs, snv_calls, 'sanger', combine=c('concordance'))
  broad <- corrected.accuracies(snvs, snv_calls, 'broad_mutect', combine=c('concordance'))
  others <- rbind(dkfz, sanger, broad)
  others <- melt(others, id=c('sample','caller'))
  results <- rbind(results, others)
  
  new_mutect <- read.csv(filename)
  new_mutect <- new_mutect[new_mutect$sample %in% levels(snvs$sample),]
  new_mutect$ref <- as.character(new_mutect$ref)
  new_mutect$alt <- as.character(new_mutect$alt)
  
  # validation data; we're not going to change these, except add a 'new_mutect' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_mutect, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0

  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_mutect, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=TRUE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  new_snv_calls$concordance[is.na(new_snv_calls$concordance)] <- 0
  
  # accuracies assuming all novel calls are wrong
  new_snvs_with_novel_false <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=FALSE, callername=newcaller))
  accuracies <- corrected.accuracies(new_snvs_with_novel_false, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_bad")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)
  
  # accuracies assuming all novel calls are right
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=TRUE, callername=newcaller))
  accuracies <- corrected.accuracies(new_snvs_with_novel_true, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_good")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)
  
  # accuracies extrapolating the per-sample private call accuracy
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, callername=newcaller, privateaccuracies=old_mutect_private_acc_by_sample))
  accuracies <- corrected.accuracies(new_snvs_with_novel_true, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_interpolated")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)

  # sort samples by number of intersect2 calls:
  sorted_samples <- names(sort(tapply(snv_calls$intersect2, snv_calls$sample, FUN=sum)))
  results$sample <- factor(results$sample, levels=sorted_samples)
  
  results$base_caller <- as.character(results$caller)
  results$base_caller[grep(shortcaller, results$base_caller)] <- shortcaller
  results$base_caller <- as.factor(results$base_caller)
  results$estimated <- FALSE
  results$estimated[results$base_caller == shortcaller] <- TRUE
  results$estimated[results$caller == paste0(shortcaller,"_original")] <- FALSE
  results
}

estimates.by.vaf <- function(snvs, snv_calls, snv_callers,
                             filename='~/Desktop/otherbroad/new_mutect.csv',
                             originalcaller='broad_mutect',
                             shortcaller='mutect',
                             vafprobs=c(1,0,0,0,0)) {
  
  privatecaller <- paste0(shortcaller, '_private')
  newcaller <- paste0('new_',shortcaller)
  
  snvs[[privatecaller]] <- snvs[[originalcaller]] & snvs$concordance == 1
  snv_calls[[privatecaller]] <- snv_calls[[originalcaller]] & snv_calls$concordance == 1
  old_caller_private_acc_by_sample <- corrected.accuracies(snvs, snv_calls, privatecaller, combine=c('concordance'))
  
  results <- corrected.accuracies.by.caller.by.vaf(snvs, snv_calls, snv_callers)

  new_calls <- read.csv(filename)
  new_calls <- new_calls[new_calls$sample %in% levels(snvs$sample),]
  new_calls$ref <- as.character(new_calls$ref)
  new_calls$alt <- as.character(new_calls$alt)
  
  # validation data; we're not going to change these, except add a 'new_mutect' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0
  
  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=TRUE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  new_snv_calls$concordance[is.na(new_snv_calls$concordance)] <- 0
  new_snv_calls$binned_wgs_tvaf[is.na(new_snv_calls$binned_wgs_tvaf)] <- sample(levels(snv_calls$binned_wgs_tvaf),
                                                                                sum(is.na(new_snv_calls$binned_wgs_tvaf)),
                                                                                replace=TRUE,
                                                                                prob=vafprobs)
  
  # accuracies assuming all novel calls are wrong
  new_snvs_with_novel_false <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=FALSE, callername=newcaller, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_false, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller,"_novel_calls_bad")
  results <- rbind(results, accuracies)
  
  # accuracies assuming all novel calls are right
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=TRUE, callername=newcaller, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_true, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller, "_novel_calls_good")
  results <- rbind(results, accuracies)
  
  # accuracies extrapolating the per-sample private call accuracy
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, callername=newcaller, privateaccuracies=old_caller_private_acc_by_sample, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_true, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller, "_novel_calls_iterpolated")
  results <- rbind(results, accuracies)
  
  results <- melt(results, id=c('VAF','caller'))
  results$base_caller <- as.character(results$caller)
  results$base_caller[grep(shortcaller, results$base_caller)] <- shortcaller
  results$estimated <- FALSE
  results$estimated[results$base_caller == shortcaller] <- TRUE
  results$estimated[results$caller == originalcaller] <- FALSE
  results
}

muse.estimates <- function(snvs, snv_calls, filename='~/Desktop/othermuse/new_muse.csv') {
  estimates(snvs, snv_calls, filename=filename, newcaller='new_muse',
            originalcaller='mda_hgsc_gatk_muse', shortcaller='muse')
}

muse.vaf.estimates <- function(snvs, snv_calls, snv_callers, 
                               filename='~/Desktop/newmuse/new_muse.csv', 
                               origcaller='mda_hgsc_gatk_muse', 
                               shortcaller='muse', 
                               vafprobs=c(0.5,0.5,0,0,0)) {
  estimates.by.vaf(snvs, snv_calls, snv_callers, filename=filename, originalcaller=origcaller, shortcaller=shortcaller, vafprobs=vafprobs)
}

add.newcaller.as.feature <- function(snvs, snv_calls, newcaller='new_muse', filename='~/Desktop/newmuse/new_muse.csv') {
  
  new_calls <- read.csv(filename)
  new_calls <- new_calls[new_calls$sample %in% levels(snvs$sample),]
  new_calls$ref <- as.character(new_calls$ref)
  new_calls$alt <- as.character(new_calls$alt)
  
  # validation data; we're not going to change these, except add a 'newcaller' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0
  
  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  
  return(list(new_snvs=new_snvs, new_snv_calls=new_snv_calls))
}

library(reshape2)
#library(ggplot2)

boxplots <- function(snvs, snv_calls, callers=c('broad_mutect','dkfz','sanger'), samples=levels(snvs$sample), title="Distribution of Per-Sample Accuracies") {
  result.list <- lapply(callers, function(x) corrected.accuracies(snvs, snv_calls, x, combine=c('concordance')))
  results <- do.call(rbind, result.list)
  results <- results[results$sample %in% samples,]
  results <- melt(results,id=c('caller','sample'))
  results = results[results$value != 0,]
  
  print(filter(results, value < 0.25))
  print(summary(filter(results, caller=="stackedLogisticRegression", variable=="sensitivity")))
  print(summary(filter(results, caller=="stackedLogisticRegression", variable=="precision")))
  print(summary(filter(results, caller=="stackedLogisticRegression", variable=="f1")))
  
  ggplot(results) + geom_boxplot(aes(x=caller,y=value,color=caller)) + 
    facet_grid(variable ~ .) + ggtitle(title) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=20)) + ylab('Accuracy')
}

require(dplyr)

calc_weights <- function (validated_calls, all_calls) {
  # calculate the weights for each validated call
   
  validated_counts <- validated_calls %>% dplyr::select(concordance, sample) %>% group_by(concordance, sample) %>% summarise(nvalidated=n())
  called_counts <- all_calls %>% dplyr::select(concordance, sample) %>% group_by(concordance, sample) %>% summarise(ncalled=n())
  allcounts <- merge(validated_counts, called_counts, all=TRUE)
  allcounts[is.na(allcounts)] <- 0
  allcounts$weight <- allcounts$ncalled/(allcounts$nvalidated+1.)
  allcounts <- dplyr::select(allcounts, concordance, sample, weight)
    
  validated_calls <- merge(validated_calls, allcounts, by=c("concordance","sample"), all.x=TRUE)

  return(validated_calls)
}

applymodel <- function(snvs, snv_calls, learn_model, predict_model, formula, model_name, 
                       samples=levels(droplevels(snvs$sample)), seed=NA, n_at_a_time=13) {
  if (!is.na(seed)) set.seed(seed)
  
  shuffled_samples <- sample(samples)
  sample_bin <- as.integer((1:length(samples))/n_at_a_time)

  copy_snvs <- snvs
  copy_snvs$row <- 1:nrow(snvs)
  copy_snvs[[model_name]] <- 0
  copy_snvs$sample_bin <- 0
  
  copy_snv_calls <- snv_calls
  copy_snv_calls$row <- 1:nrow(snv_calls)
  copy_snv_calls[[model_name]] <- 0
  copy_snvs$sample_bin <- 0
  
  if (!"weight" %in% names(copy_snvs)) {
    copy_snvs <- calc_weights(copy_snvs, snv_calls)
  }
  
  for (i in 1:length(shuffled_samples)) {
    copy_snvs$sample_bin[copy_snvs$sample == shuffled_samples[i]] <-sample_bin[i]
    copy_snv_calls$sample_bin[copy_snv_calls$sample == shuffled_samples[i]] <-sample_bin[i]
  }
  
  models <- list()
  for (i in 0:max(sample_bin)) {
    print(paste(i+1,"of",max(sample_bin)+1))
    train_snvs <- filter(copy_snvs, sample_bin != i)
    train_snv_calls <- filter(copy_snv_calls, sample_bin != i)
    test_snvs <- filter(copy_snvs, sample_bin == i)
    test_snv_calls <- filter(copy_snv_calls, sample_bin == i)
    
    # learn on the train and apply to test
    models[[i+1]] <- learn_model(formula, train_snvs)    
    test_snv_calls_predictions <- predict_model(models[[i+1]], formula, test_snv_calls)
    test_snvs_predictions <- predict_model(models[[i+1]], formula, test_snvs)
    
    copy_snvs[[model_name]][test_snvs$row] <- test_snvs_predictions
    copy_snv_calls[[model_name]][test_snv_calls$row] <- test_snv_calls_predictions
  }
  copy_snv_calls[[model_name]][copy_snv_calls$union == 0] <- 0
  copy_snvs[[model_name]][copy_snvs$union == 0] <- 0
  
#  if ("broad_snowman" %in% names(copy_snvs)) {
#    copy_snv_calls[[model_name]][copy_snv_calls$intersect2 == 0 & copy_snv_calls$broad_snowman == 1] <- 0
#    copy_snvs[[model_name]][copy_snvs$intersect2 == 0 & copy_snvs$broad_snowman == 1] <- 0
#  }
  
  return(list(snv_calls=copy_snv_calls, snvs=copy_snvs, models=models))
}

require(glmnet)
require(party)
require(kernlab)
require(randomForest)

valid.rows <- function(data) {
  return(complete.cases(data))
}

only.required.columns <- function(formula, data, include.answer=TRUE, answer="validate_true") {
  variables <- unique(sort(unlist(strsplit(gsub("[+:|()*]"," ",as.character(formula)[[3]])," "))))
  variables <- variables[-1]
  if (include.answer) {
    variables <- c(answer, variables)
  }
  subsetdata <- data[variables]
#  print(summary(subsetdata))
  return(subsetdata)
}

svmlearn <- function(formula, data) { 
  subsetdata <- only.required.columns(formula, data)
  myrows <- valid.rows(subsetdata)
  ksvm(formula, subsetdata[myrows,], type="C-svc", prob.model=TRUE) 
}

svmpredict <- function(model, formula, data) { 
  subsetdata <- only.required.columns(formula, data, FALSE)
  myrows <- valid.rows(subsetdata)
  
  results <- rep(0,nrow(data))
  prediction <- predict(model, data, type="probabilities")[,2]
  
  results[myrows] <- prediction
  return(results)
}

randomforestlearn <- function(formula, data) { 
  subsetdata <- only.required.columns(formula, data)
  subsetdata$weight <- data$weight
  myrows <- valid.rows(subsetdata)
  randomForest(formula, subsetdata[myrows,], ntree=250)
}

randomforestpredict <- function(model, formula, data) { 
  subsetdata <- only.required.columns(formula, data, FALSE)
  myrows <- valid.rows(subsetdata)
  
  results <- rep(0,nrow(data))
  prediction <- predict(model, subsetdata[myrows,], type="prob")[,2]
  
  results[myrows] <- prediction
  return(results)
}

glmlearn <- function(formula, data) { 
  data$intweights = as.integer(10*data$weight)
  glm(formula, data, weights=intweights, family="binomial") 
}

glmpredict <- function(model, formula, data) { predict(model, data, type="response") }

ctreelearn <- function(formula, data) { ctree(formula, newdata=data, weights=data$weight) }
ctreepredict <- function(model, formula, data) { Predict(model, newdata=data) }

glmnetlearn <- function(formula, data) {
  subsetdata <- only.required.columns(formula, data)
  subsetdata$weight <- data$weight
  myrows <- valid.rows(subsetdata)
  X <- model.matrix(formula, data=subsetdata[myrows,])
  y <- subsetdata$validate_true[myrows]
 
  cv.glmnet(X, y, weights=subsetdata$weight[myrows], family="binomial", foldid=as.integer(droplevels(data$sample[myrows])))
}

glmnetpredict <- function(model, formula, data) {
  subsetdata <- only.required.columns(formula, data, FALSE)
  subsetdata$validate_true <- FALSE
  myrows <- valid.rows(subsetdata)
  
  results <- rep(0,nrow(data))
  X <- model.matrix(formula, data=subsetdata[myrows,])
  prediction <- predict.cv.glmnet(model, X, type="response", s="lambda.1se")[,1]

  results[myrows] <- prediction
  return(results)
}

binarize <- function(data, thresh=0.5) {
    ifelse(data > thresh, 1, 0)
}

binarize_both <- function(l, probs, caller, thresh=0.5) {
  l$snvs[[caller]] <- binarize(l$snvs[[probs]], thresh)
  l$snv_calls[[caller]] <- binarize(l$snv_calls[[probs]], thresh)
  return(l)
}

simple_indel_models <- function(validated.calls, all.calls, seed=1, threshold=0.5) {
  
  print("Logistic Regression")
  l <- applymodel(validated.calls, all.calls, glmlearn, glmpredict, 
 #                as.formula("validate_true ~ repeat_count +  wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + broad_snowman*dkfz*sanger"), 
                  as.formula("validate_true ~ repeat_count +  wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + dkfz*sanger*smufin*broad_snowman"), 
                 "logistic_regression", seed=1)
  
  print("Decision Tree")
  l <- applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
#                 as.formula("validate_true ~  repeat_count + wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes+ broad_snowman*dkfz*sanger"), 
                  as.formula("validate_true ~ repeat_count +  wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + dkfz*sanger*smufin*broad_snowman"), 
                  "decision_tree", seed=1)
  
  print("SVM")
  l <- applymodel(l$snvs, l$snv_calls, svmlearn, svmpredict, 
#                  as.formula("validate_true ~  repeat_count + wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes+ broad_snowman*dkfz*sanger"), 
                  as.formula("validate_true ~ repeat_count +  wgs_tvaf + wgs_tvardepth + wgs_nvaf + wgs_nvardepth + varlen + cosmic + dbsnp + thousand_genomes + dkfz*sanger*smufin*broad_snowman"), 
                  "svm_prob", seed=1)
  
  print("Stacked Logistic Regression")
  l <- applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
#                  as.formula("validate_true ~ (broad_snowman + dkfz + sanger)*(wgs_nvaf + wgs_tvaf  + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count)"), 
                   as.formula("validate_true ~ (dkfz + sanger + smufin)*(wgs_tvaf + wgs_tvardepth + wgs_nvaf + wgs_nvardepth + varlen + cosmic + dbsnp + thousand_genomes + repeat_masker + repeat_count + broad_snowman)"),
                   "stacked_logistic_regression", seed=1) 
  
  print("RandomForest")
  l <- applymodel(l$snvs, l$snv_calls, randomforestlearn, randomforestpredict, 
                  as.formula("as.factor(validate_true) ~  repeat_count + wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic*dbsnp + thousand_genomes+ dkfz*sanger*broad_snowman*smufin"), 
#                  as.formula("as.factor(validate_true) ~  repeat_count + wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic*dbsnp + thousand_genomes+ broad_snowman*dkfz*sanger"), 
                  "random_forest", seed=1) 

  l$snvs$logisticRegression <- binarize(l$snvs$logistic_regression, thresh=threshold)
  l$snv_calls$logisticRegression <- binarize(l$snv_calls$logistic_regression, thresh=threshold)
  
  l$snvs$decisionTree <- binarize(l$snvs$decision_tree, thresh=threshold)
  l$snv_calls$decisionTree <- binarize(l$snv_calls$decision_tree, thresh=threshold)
  
  l$snv_calls$two_plus <- ifelse(l$snv_calls$sanger + l$snv_calls$broad_mutect + l$snv_calls$dkfz + l$snv_calls$smufin >= 2, 1, 0)
  l$snvs$two_plus <- ifelse(l$snvs$sanger + l$snvs$broad_mutect + l$snvs$dkfz + l$snvs$smufin >= 2, 1, 0)
  
  l$snvs$stackedLogisticRegression <- binarize(l$snvs$stacked_logistic_regression, thresh=threshold)
  l$snv_calls$stackedLogisticRegression <- binarize(l$snv_calls$stacked_logistic_regression, thresh=threshold)
  
  l$snvs$svm <- binarize(l$snvs$svm_prob, thresh=threshold)
  l$snv_calls$svm <- binarize(l$snv_calls$svm_prob, thresh=threshold)
  
  l$snvs$randomForest <- binarize(l$snvs$random_forest, thresh=threshold)
  l$snv_calls$randomForest <- binarize(l$snv_calls$random_forest, thresh=threshold)
  return(l)
}

vary_ncallers_tree <- function(validated.calls, all.calls, seed=1, threshold=0.5) {
  
  print("Two")
  l <- applymodel(validated.calls, all.calls, ctree, ctreepredict,  
                  as.formula("validate_true ~ (dkfz*sanger + wgs_tvaf*wgs_tvardepth*wgs_tdepth + wgs_nvaf*wgs_nvardepth*wgs_ndepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_masker + repeat_count)"), 
                  "two_callers", seed=1) 
  l <- binarize_both(l, "two_callers", "twoCallers", threshold)
  
  print("All four")
  l <- applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                  as.formula("validate_true ~ (dkfz*sanger*broad_snowman*smufin + wgs_tvaf*wgs_tvardepth*wgs_tdepth + wgs_nvaf*wgs_nvardepth*wgs_ndepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_masker + repeat_count)"), 
                  "four_callers", seed=1) 
  l <- binarize_both(l, "four_callers", "fourCallers", threshold)
  
  print("two + smufin")
  l <- applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                  as.formula("validate_true ~ (dkfz*sanger*smufin + wgs_tvaf*wgs_tvardepth*wgs_tdepth + wgs_nvaf*wgs_nvardepth*wgs_ndepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_masker + repeat_count)"), 
                  "two_plus_smufin", seed=1) 
  l <- binarize_both(l, "two_plus_smufin", 'twoPlusSmufin', threshold)

  print("two + snowman")
  l <- applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                  as.formula("validate_true ~ (dkfz*sanger*broad_snowman + wgs_tvaf*wgs_tvardepth*wgs_tdepth + wgs_nvaf*wgs_nvardepth*wgs_ndepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_masker + repeat_count)"), 
                  "two_plus_snowman_caller", seed=1) 
  l <- binarize_both(l, "two_plus_snowman_caller", 'twoPlusSnowmanCaller', threshold)
  
  l$snv_calls$two_plus <- ifelse(l$snv_calls$sanger + l$snv_calls$broad_snowman + l$snv_calls$dkfz + l$snv_calls$smufin >= 2, 1, 0)
  l$snvs$two_plus <- ifelse(l$snvs$sanger + l$snvs$broad_mutect + l$snvs$dkfz + l$snvs$smufin >= 2, 1, 0)

  return(l)
}

vary_ncallers <- function(validated.calls, all.calls, seed=1, threshold=0.5) {
  
  print("Two")
  l <- applymodel(validated.calls, all.calls, glmnetlearn, glmnetpredict, 
                  as.formula("validate_true ~ (dkfz + sanger)*(wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count)"), 
                  "two_callers", seed=1) 
  l <- binarize_both(l, "two_callers", "twoCallers", threshold)
  
  print("All four")
  l <- applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                  as.formula("validate_true ~ (dkfz + sanger + smufin)*(wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count + broad_snowman)"), 
                  "four_callers", seed=1) 
  l <- binarize_both(l, "four_callers", "fourCallers", threshold)
  
  print("two + smufin")
  l <- applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                  as.formula("validate_true ~ (dkfz + sanger + smufin)*(wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count)"), 
                  "two_plus_smufin", seed=1) 
  l <- binarize_both(l, "two_plus_smufin", 'twoPlusSmufin', threshold)
  
  print("two + snowman")
  l <- applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                  as.formula("validate_true ~ (dkfz + sanger)*(wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count + broad_snowman)"), 
                  "two_plus_snowman_feature", seed=1) 
  l <- binarize_both(l, "two_plus_snowman_feature", 'twoPlusSnowmanFeature', threshold)
  
  print("two + snowman caller")
  l <- applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                  #                  as.formula("validate_true ~ (broad_snowman + dkfz + sanger)*(wgs_nvaf + wgs_tvaf  + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count)"), 
                  as.formula("validate_true ~ (dkfz + sanger + broad_snowman)*(wgs_tvaf*wgs_tvardepth + wgs_nvaf*wgs_nvardepth + varlen + cosmic:dbsnp + thousand_genomes + repeat_count)"), 
                  "two_plus_snowman_caller", seed=1) 
  l <- binarize_both(l, "two_plus_snowman_caller", 'twoPlusSnowmanCaller', threshold)
  
  
  l$snv_calls$two_plus <- ifelse(l$snv_calls$sanger + l$snv_calls$broad_snowman + l$snv_calls$dkfz + l$snv_calls$smufin >= 2, 1, 0)
  l$snvs$two_plus <- ifelse(l$snvs$sanger + l$snvs$broad_mutect + l$snvs$dkfz + l$snvs$smufin >= 2, 1, 0)
  
  return(l)
}

simple_snv_models <- function(validated.calls, all.calls, seed=1, threshold=0.5) {
  print("Logistic Regression")
  l <-  applymodel(validated.calls, all.calls, glmlearn, glmpredict, 
                   as.formula("validate_true ~ union + intersect2 + intersect3 + wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                   "logistic_regression", seed=1)
  
  print("Decision Tree")
  l <-  applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                   as.formula("validate_true ~ union + intersect2 + intersect3 +  wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                   "decision_tree", seed=1)
  
  print("Stacked Logistic Regression")
  l <-  applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                   as.formula("validate_true ~ (broad_mutect + dkfz + sanger)*(wgs_nvaf + wgs_tvaf + gencode + cosmic + dbsnp + muse_feature)"), 
                   "stacked_logistic_regression", seed=1) 
  
  print("SVM")
  l <- applymodel(l$snvs, l$snv_calls, svmlearn, svmpredict, 
                  as.formula("validate_true ~ union + intersect2 + intersect3 +  wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                  "svm_prob", seed=1)
  
  print("RandomForest")
  l <- applymodel(l$snvs, l$snv_calls, randomforestlearn, randomforestpredict, 
                  as.formula("as.factor(validate_true) ~ union + intersect2 + intersect3 +  wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_nvaf*wgs_nvardepth + wgs_tvaf*wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                  "random_forest", seed=1) 
  
  l$snvs$logisticRegression <- binarize(l$snvs$logistic_regression, thresh=threshold)
  l$snv_calls$logisticRegression <- binarize(l$snv_calls$logistic_regression, thresh=threshold)
  
  l$snvs$decisionTree <- binarize(l$snvs$decision_tree)
  l$snv_calls$decisionTree <- binarize(l$snv_calls$decision_tree)
  
  l$snvs$stackedLogisticRegression <- binarize(l$snvs$stacked_logistic_regression, thresh=threshold)
  l$snv_calls$stackedLogisticRegression <- binarize(l$snv_calls$stacked_logistic_regression, thresh=threshold)
  
  l$snvs$svm <- binarize(l$snvs$svm_prob, thresh=threshold)
  l$snv_calls$svm <- binarize(l$snv_calls$svm_prob, thresh=threshold)
  
  l$snvs$randomForest <- binarize(l$snvs$random_forest, thresh=threshold)
  l$snv_calls$randomForest <- binarize(l$snv_calls$random_forest, thresh=threshold)
  
  return(l)
}

indel_model_callers <- c('broad_snowman','dkfz', 'sanger', 'union', 'intersect2', 'two_plus', 'intersect3',
                   'logisticRegression', 'decisionTree'
                   ,'stackedLogisticRegression', 'svm', 'randomForest')

doboxplot <- function(l, callers=c('broad_mutect','dkfz', 'sanger', 'union', 'intersect2', 'intersect3', 'two_plus',
                                 'logisticRegression', 'decisionTree'
                                 ,'stackedLogisticRegression', 'svm', 'randomForest'), ...) {
  boxplots(l$snvs, l$snv_calls, callers, ...)
}

dovafplot <- function(l, variant="SNV", callers=c('union', 'intersect2', 'intersect3', 'two_plus', 'logisticRegression', 'stackedLogisticRegression', 'decisionTree', 'svm', 'randomForest')) {
  results <- corrected.accuracies.by.caller.by.vaf(l$snvs, l$snv_calls, callers) 
  results <- melt(results, id=c('caller','VAF'))
  results$feature <- FALSE
  results$feature[results$caller %in% c("logisticRegression","decisionTree","stackedLogisticRegression","svm","randomForest")] <- TRUE
  
  print(results %>% filter(variable=="f1", VAF=="[0,0.1]") %>% arrange(value))
  
  ggplot(results, aes(x=VAF,y=value)) + 
    geom_line(aes(group=caller,color=caller,linetype=feature)) + 
    facet_grid(variable~.)  +
    ggtitle(paste0("Accuracies of ", variant, " merge by VAF")) +
    ylab("Accuracy")
}

dosampleplot <- function(l, callers=c("union","intersect2","intersect3","two_plus","logisticRegression","stackedLogisticRegression","decisionTree","svm","randomForest"), variant="SNV") {
  results <- corrected.accuracies.by.caller.by.sample(l$snvs, l$snv_calls, callers)
  results <- melt(results, id=c('caller','sample'))
  results$feature <- FALSE
  results$feature[results$caller %in% c("logisticRegression","decisionTree","stackedLogisticRegression","svm","randomForest")] <- TRUE
  results <- results[results$value > 0,]
  
  samples_levels <- names(sort(table(l$snv_calls$sample)))
  results$sample <- factor(results$sample, levels=samples_levels)
  
  ggplot(results, aes(x=sample,y=value)) + 
    geom_line(aes(group=caller,color=caller,linetype=feature)) + 
    facet_grid(variable~.) +
    ggtitle(paste0("Accuracies of ", variant, " merge by sample")) +
    ylab("Accuracy") + xlab("Sample")
}

dorocplot <- function(l, title="SNV Models", xlim=c(0,1), ylim=c(0,1)) {
  rocplot(l$snvs, l$snv_calls, 
          c('broad_mutect','dkfz', 'sanger', 'union', 'intersect2', 'intersect3', 'two_plus'),  
          c('logistic_regression', 'decision_tree', 'stacked_logistic_regression', 'svm_prob', 'random_forest'),
          c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
          title, xlim=xlim, ylim=ylim)
}

require(tidyr)
frac.by.status <- function(validated.calls, all.calls, caller, statuses=c("GERMLINE","NOTSEEN","PASS")) {
  all.caller.calls <- filter(all.calls, all.calls[[caller]]==1)
  validated.caller.calls <- filter(validated.calls, validated.calls[[caller]]==1)
  validated.caller.calls$sample <- factor(as.character(validated.caller.calls$sample), levels=levels(all.caller.calls$sample))
  
  all.matrix <- as.matrix(xtabs(~concordance+sample, all.caller.calls, drop.unused.levels=FALSE))
  all.matrix[is.na(all.matrix)] <- 0
  
  all.bysample <- colSums(all.matrix)
  
  matsize <- dim(all.matrix)
  validated.matrix <- matrix(0., nrow=matsize[1], ncol=matsize[2])
  vm  <- as.matrix(xtabs(~concordance+sample, validated.caller.calls, drop.unused.levels=FALSE))
  validated.matrix[1:dim(vm)[1],] <- vm
  validated.matrix[is.na(validated.matrix)] <- 0
  
  results <- data.frame()
  
  for (stat in statuses) {
    validated.status.caller.calls <- filter(validated.caller.calls, status==stat)
    validated.status.matrix <- matrix(0., nrow=matsize[1], ncol=matsize[2])
    vm <- as.matrix(xtabs(~concordance+sample, validated.status.caller.calls, drop.unused.levels=FALSE))
    validated.status.matrix[1:dim(vm)[1],] <- vm
    
    frac <- validated.status.matrix / validated.matrix
    frac0 <- frac
    frac0[is.na(frac0)] <- 0
    frac1 <- frac
    frac1[is.na(frac1)] <- 1
    
    low.num <- colSums(frac0*all.matrix)
    low.frac <- low.num/all.bysample
    
    hi.num <- colSums(frac1*all.matrix)
    hi.frac <- hi.num/all.bysample
    
    hi.stat.df <- data.frame(sample=names(hi.num), num=hi.num, frac=hi.frac)
    lo.stat.df <- data.frame(sample=names(low.num), num=low.num, frac=low.frac)
    stat.df <- rbind(hi.stat.df, lo.stat.df)
    stat.df$status <- stat
    stat.df$sample <- as.factor(stat.df$sample)
  
    results <- rbind(results, stat.df)  
  }
  
  results$caller <- caller
  
  badsamples <- c("a34f1dba-5758-45c8-b825-1c888d6c4c13", 'ae1fd34f-6a0f-43db-8edb-c329ccf3ebae')
  results <- filter(results, !sample %in% badsamples )
  return(results[complete.cases(results),])
}

frac.by.status.by.caller <- function(validated.calls, all.calls, callers,  statuses=c("GERMLINE","NOTSEEN","PASS")) {
results <- lapply(callers, function(caller){frac.by.status(validated.calls, all.calls, caller, statuses)})
  return(do.call(rbind, results))
}

hets_and_homs <- function(validated.calls, hetfrac=0.5, hetsd=2, homfrac=1., homsd=2) {
    newdata <- validated.calls
    newdata$status <- as.character(newdata$status)
    newdata$status[newdata$status == "GERMLINE"] <- "NORMAL_EVIDENCE"
    het <- (newdata$val_nvaf >= hetfrac - hetsd*(hetfrac)/sqrt(newdata$val_ndepth)) & (newdata$val_nvaf <= hetfrac + hetsd*(hetfrac)/sqrt(newdata$val_ndepth))
    hom <- (newdata$val_nvaf >= homfrac - homsd*(homfrac)/sqrt(newdata$val_ndepth))
    newdata$status[het] <- "GERMLINE"
    newdata$status[hom] <- "GERMLINE"
    newdata$status <- factor(newdata$status)
    return(newdata)
}