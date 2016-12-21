source('analysis/precision_recall.R')
source('analysis/estimate_accuracies.R')

allfour <- as.formula("validate_true ~ (dkfz + sanger + smufin)*(wgs_tvaf + wgs_tvardepth + wgs_nvaf + wgs_nvardepth + varlen + cosmic + dbsnp + thousand_genomes + repeat_masker + repeat_count + broad_snowman)")
nosmufin<- as.formula("validate_true ~ (dkfz + sanger)*(wgs_tvaf + wgs_tvardepth + wgs_nvaf + wgs_nvardepth + varlen + cosmic + dbsnp + thousand_genomes + repeat_masker + repeat_count + broad_snowman)")
nobroad <- as.formula("validate_true ~ (dkfz + sanger + smufin)*(wgs_tvaf + wgs_tvardepth + wgs_nvaf  + wgs_nvardepth + varlen + cosmic + dbsnp + thousand_genomes + repeat_masker + repeat_count)")

train_stacked_model <- function(model.formula, indels, indel_calls) {
  copy_indels <- calc_weights(indels, indel_calls)
  model <- glmnetlearn(model.formula, copy_indels)
}

train_stacked_models <- function(indels, indel_calls, dir='../consensus/final_consensus/models/') {
  formula <- allfour
  model <- train_stacked_model(formula, indels, indel_calls)
  predict_function <- glmnetpredict
  
  save(formula, model, predict_function, file=paste0(dir,"stacked-logistic-all-four.RData"))
  
  
  formula <- nosmufin
  model <- train_stacked_model(formula, indels, indel_calls)
  predict_function <- glmnetpredict
  
  save(formula, model, predict_function, file=paste0(dir, "stacked-logistic-no-smufin.RData"))
  
  
  formula <- nobroad
  model <- train_stacked_model(formula, indels, indel_calls)
  predict_function <- glmnetpredict
  
  save(formula, model, predict_function, file=paste0(dir,"stacked-logistic-no-broad.RData"))
}
