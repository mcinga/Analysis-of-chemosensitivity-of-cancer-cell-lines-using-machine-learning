library(tidyverse)
library(pROC)
library(caret)

#===============================================================================
gd_09_17<-read_csv("GDSC_09_CRISPR_t5000.csv") 
unique(gd_09_17$CELL_LINE_NAME)
unique(gd_09_17$DRUG_NAME)
dim(gd_09_17)

#scaling and centering the data
gd_09_17<-gd_09_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))


#training the models
set.seed(123)
inTraining <- createDataPartition(gd_09_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_09_17[ inTraining,]
testing  <- gd_09_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

##CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)

models<-c("knn","svmRadial","gbm","rf")
results_table_09 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_09.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_09[results_table_09$models %in% i, "Precision"] <- precision
  results_table_09[results_table_09$models %in% i, "Recall"] <- recall
  results_table_09[results_table_09$models %in% i, "F1score"] <- f1
  results_table_09[results_table_09$models %in% i, "Accuracy"] <- accuracy
  results_table_09[results_table_09$models %in% i, "AUC"] <- auc
}
write_csv(results_table_09, "results_table_09_t5000_2(1).csv")

#===============================================================================
gd_08_17<-read_csv("GDSC_08_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_08_17$CELL_LINE_NAME)
unique(gd_08_17$DRUG_NAME)
dim(gd_08_17)

#scaling and centering the data
gd_08_17<-gd_08_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_08_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_08_17[ inTraining,]
testing  <- gd_08_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_08 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_08.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_08[results_table_08$models %in% i, "Precision"] <- precision
  results_table_08[results_table_08$models %in% i, "Recall"] <- recall
  results_table_08[results_table_08$models %in% i, "F1score"] <- f1
  results_table_08[results_table_08$models %in% i, "Accuracy"] <- accuracy
  results_table_08[results_table_08$models %in% i, "AUC"] <- auc
  
}
write_csv(results_table_08, "results_table_08_t5000_2.csv")


#===============================================================================
gd_07_17<-read_csv("GDSC_07_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_07_17$CELL_LINE_NAME)
unique(gd_07_17$DRUG_NAME)
dim(gd_07_17)

#scaling and centering the data
gd_07_17<-gd_07_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_07_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_07_17[ inTraining,]
testing  <- gd_07_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_07 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_07.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_07[results_table_07$models %in% i, "Precision"] <- precision
  results_table_07[results_table_07$models %in% i, "Recall"] <- recall
  results_table_07[results_table_07$models %in% i, "F1score"] <- f1
  results_table_07[results_table_07$models %in% i, "Accuracy"] <- accuracy
  results_table_07[results_table_07$models %in% i, "AUC"] <- auc
  
}
write_csv(results_table_07, "results_table_07_t5000_2.csv")
#===============================================================================
gd_06_17<-read_csv("GDSC_06_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_06_17$CELL_LINE_NAME)
unique(gd_06_17$DRUG_NAME)
dim(gd_06_17)

#scaling and centering the data
gd_06_17<-gd_06_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_06_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_06_17[ inTraining,]
testing  <- gd_06_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_06 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_06.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_06[results_table_06$models %in% i, "Precision"] <- precision
  results_table_06[results_table_06$models %in% i, "Recall"] <- recall
  results_table_06[results_table_06$models %in% i, "F1score"] <- f1
  results_table_06[results_table_06$models %in% i, "Accuracy"] <- accuracy
  results_table_06[results_table_06$models %in% i, "AUC"] <- auc
}
write_csv(results_table_06, "results_table_06_t5000_2.csv")

#===============================================================================
gd_05_17<-read_csv("GDSC_05_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_05_17$CELL_LINE_NAME)
unique(gd_05_17$DRUG_NAME)
dim(gd_05_17)

#scaling and centering the data
gd_05_17<-gd_05_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_05_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_05_17[ inTraining,]
testing  <- gd_05_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_05 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_05.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_05[results_table_05$models %in% i, "Precision"] <- precision
  results_table_05[results_table_05$models %in% i, "Recall"] <- recall
  results_table_05[results_table_05$models %in% i, "F1score"] <- f1
  results_table_05[results_table_05$models %in% i, "Accuracy"] <- accuracy
  results_table_05[results_table_05$models %in% i, "AUC"] <- auc
  
}
write_csv(results_table_05, "results_table_05_t5000_2.csv")

#==============================================================================
gd_04_17<-read_csv("GDSC_04_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_04_17$CELL_LINE_NAME)
unique(gd_04_17$DRUG_NAME)
dim(gd_04_17)

#scaling and centering the data
gd_04_17<-gd_04_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_04_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_04_17[ inTraining,]
testing  <- gd_04_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_04 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_04.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_04[results_table_04$models %in% i, "Precision"] <- precision
  results_table_04[results_table_04$models %in% i, "Recall"] <- recall
  results_table_04[results_table_04$models %in% i, "F1score"] <- f1
  results_table_04[results_table_04$models %in% i, "Accuracy"] <- accuracy
  results_table_04[results_table_04$models %in% i, "AUC"] <- auc
  
}
write_csv(results_table_04, "results_table_04_t5000_2.csv")
#=====================================================================

gd_03_17<-read_csv("GDSC_03_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_03_17$CELL_LINE_NAME)
unique(gd_03_17$DRUG_NAME)
dim(gd_03_17)

#scaling and centering the data
gd_03_17<-gd_03_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_03_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_03_17[ inTraining,]
testing  <- gd_03_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_03 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_03.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_03[results_table_03$models %in% i, "Precision"] <- precision
  results_table_03[results_table_03$models %in% i, "Recall"] <- recall
  results_table_03[results_table_03$models %in% i, "F1score"] <- f1
  results_table_03[results_table_03$models %in% i, "Accuracy"] <- accuracy
  results_table_03[results_table_03$models %in% i, "AUC"] <- auc  
}
write_csv(results_table_03, "results_table_03_t5000_2.csv")

#===============================================================================
gd_02_17<-read_csv("GDSC_02_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_02_17$CELL_LINE_NAME)
unique(gd_02_17$DRUG_NAME)
dim(gd_02_17)

#scaling and centering the data
gd_02_17<-gd_02_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_02_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_02_17[ inTraining,]
testing  <- gd_02_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_02 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_02.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_02[results_table_02$models %in% i, "Precision"] <- precision
  results_table_02[results_table_02$models %in% i, "Recall"] <- recall
  results_table_02[results_table_02$models %in% i, "F1score"] <- f1
  results_table_02[results_table_02$models %in% i, "Accuracy"] <- accuracy
  results_table_02[results_table_02$models %in% i, "AUC"] <-auc
  
}
write_csv(results_table_02, "results_table_02_t5000_2.csv")

#=============================================================================
gd_01_17<-read_csv("GDSC_01_CRISPR_t5000.csv") #this is for the full dataset
unique(gd_01_17$CELL_LINE_NAME)
unique(gd_01_17$DRUG_NAME)
dim(gd_01_17)

#scaling and centering the data
gd_01_17<-gd_01_17%>%
  mutate(across(where(is.numeric), \(x) scale(x)[,1]))
#training the models
set.seed(123)
inTraining <- createDataPartition(gd_01_17$BIOACTIVITY, p = .7, list = FALSE)
training <- gd_01_17[ inTraining,]
testing  <- gd_01_17[-inTraining,]

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10,repeats = 5, verboseIter = F, classProbs = T)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                         stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)


models<-c("knn","svmRadial","gbm","rf")
results_table_01 <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(BIOACTIVITY~., data = training, method = i,
                       trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  
  probabilities <- predict(model_train, newdata = testing, type = "prob")
  roc_obj <- roc(testing$BIOACTIVITY, probabilities[, 2])
  auc <- auc(roc_obj)
   pdf(paste0(i, "_roc_curve_01.pdf"))
  plot(roc_obj, main = paste0(i, " - ROC Curve"), print.auc = TRUE,
       auc.polygon = TRUE, grid = TRUE)
  dev.off()
  
  table_mat<-table(testing$BIOACTIVITY, predictions)
  accuracy<-sum(diag(table_mat))/sum(table_mat)
  precision<-posPredValue(predictions, testing$BIOACTIVITY)
  recall <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision* recall) / (precision + recall)
  
  # put that in the results table
  results_table_01[results_table_01$models %in% i, "Precision"] <- precision
  results_table_01[results_table_01$models %in% i, "Recall"] <- recall
  results_table_01[results_table_01$models %in% i, "F1score"] <- f1
  results_table_01[results_table_01$models %in% i, "Accuracy"] <- accuracy
  results_table_01[results_table_01$models %in% i, "AUC"] <- auc
  
}
write_csv(results_table_01, "results_table_01_t5000_2.csv")

