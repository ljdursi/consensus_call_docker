
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
#require(party)
#require(kernlab)
#require(randomForest)

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
