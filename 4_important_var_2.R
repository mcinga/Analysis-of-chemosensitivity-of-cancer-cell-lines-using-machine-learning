library(tidyverse)
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
  
}
#write_csv(results_table_09, "results_table_09_t5000_2(test).csv")
#WILL USE THE PERMUTATION FEATURE METHOD HERE
importance1 <- varImp(model_train, useModel=F, scale = FALSE, type = 1, test=testing)

importance_df <- data.frame(variable = rownames(importance1$importance),
                            importance = importance1$importance[, 1], row.names = NULL)
importance_df <- importance_df[order(-importance_df$importance),]
df_imp<-importance_df[1:300,]

write_csv(df_imp,"top_300_variables_(2).csv")