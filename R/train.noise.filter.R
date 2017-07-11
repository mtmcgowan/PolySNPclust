#' Train noise filter
#'
#' Reads in completed training data and trains a gradient boosted model outputting prediction accuracy statistics
#'
#' @param filepath the path to a completed .csv file output from the make.training.data() function that has been scored manually
#'
#' @return A gbm model object
#'

train.noise.filter <- function(trainingpath = './training.data/training.csv', seed = NULL) {
  # Reading in the training data
  cat('\n', 'Reading in the training data and preparing for modeling', '\n')
  train.data <- read.table(trainingpath, header = T, sep = '\t', stringsAsFactors = F)

  # Converting the classifier into a factor
  train.data$class <- ifelse(train.data$class==1,'real','noise')
  train.data$class <- as.factor(train.data$class)

  # Splitting the data into a training and test set
  if (!is.null(seed)) {set.seed(seed)}
  sample <- sample(1:nrow(train.data), round(0.1 * nrow(train.data)))
  subdata.train <- train.data[-sample,]
  subdata.test <- train.data[sample,]

  # Configuring Caret to automatically control the resampling/testing of my data (training it 10 times in this case)
  objControl <- trainControl(method='cv', number=10, returnResamp='none', summaryFunction = twoClassSummary, classProbs = TRUE)

  # Setting outcome and predictors
  outcome_name <- 'class'
  predictor_names <- c('prop', 'xcoord', 'ycoord', 'var1', 'var2')

  # print update
  cat('\n', 'Training a CV model with a subset of training data', '\n')
  # Opening a redirect pipe to write stdout to a log file
  sink("gbm.training.log.txt")
  # Running the cross-validation model
  model_gbm <- train(subdata.train[,predictor_names], subdata.train[,outcome_name],
                     method = 'gbm',
                     trControl=objControl,
                     metric = 'ROC',
                     preProc = c('center', 'scale')
  )
  sink()
  # Printing CV model summary and accuracies
  summary(model_gbm)

  predictions <- predict(object=model_gbm, subdata.test[,predictor_names], type = 'raw')
  print(postResample(pred = predictions, obs = as.factor(subdata.test[,outcome_name])))

  cat('\n', 'Training the final model with all training data', '\n')
  sink("gbm.training.log.txt")
  # Training the complete model
  model_gbm <- train(train.data[,predictor_names], train.data[,outcome_name],
                     method = 'gbm',
                     trControl=objControl,
                     metric = 'ROC',
                     preProc = c('center', 'scale')
  )
  sink()

  cat('Modeling Finished!')
  return(model_gbm)
}
