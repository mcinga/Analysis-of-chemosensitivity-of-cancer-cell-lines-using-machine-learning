library(tidyverse)
df1<-read_csv("GDSC1_fitted_dose_response_25Feb20.csv")
df2<-read_csv("GDSC2_fitted_dose_response_25Feb20.csv")


#WANT TO FOCUS ON THE TISSUE TYPE THAT IS ONLY BRCA 
#(BREAST CANCER) SO I HAVE TO FILTER.
BRCA1<-df1%>%
  filter(TCGA_DESC == "BRCA")%>%
  dplyr::select(DATASET,CELL_LINE_NAME, 
                DRUG_NAME,PATHWAY_NAME,LN_IC50, Z_SCORE)%>%
  select(-DATASET)

BRCA2<-df2%>%
  filter(TCGA_DESC == "BRCA")%>%
  dplyr::select(DATASET,CELL_LINE_NAME,
                DRUG_NAME, PATHWAY_NAME, LN_IC50,Z_SCORE)%>%
  select(-DATASET)

#MERGE DATA RETAIN ONLY THOSE THAT ARE NOT DUPLICATES IN THE DATA.
library(sqldf)
BRCA<-sqldf("select * from BRCA1 union select * from BRCA2")
unique(BRCA$CELL_LINE_NAME) #51 CELL LINES


#SEPERATING THE DRUGS ACCORDING TO SENSISTIVITY, 
#USING THE Z_SCORE
BRCA<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.5 ~ "SENSITIVE",
                                 Z_SCORE >= 0.5 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))
#I will need to explain the above and state that with more intermediate drug responses the model did not perform well at all.

#SEPERATING THE DRUGS ACCORDING TO SENSISTIVITY, 
#USING THE Z_SCORE
BRCA<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.3 ~ "SENSITIVE",
                                 Z_SCORE >= 0.3 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))


#LOAD THE GDSC DRUG INFORMATION
drug1<-read_csv("GDSC1 drugs.csv")
drug2<-read_csv("GDSC2 drugs.csv")

drug1<-drug1%>%
  select(-c(1,3,6))%>%
  rename(DRUG_NAME=drug_name,
         PATHWAY_NAME=pathway_name,TARGETS=targets)

drug2<-drug2%>%
  select(-c(1,3,6))%>%
  rename(DRUG_NAME=drug_name,
         PATHWAY_NAME=pathway_name,TARGETS=targets)

#USE INNER JOIN BECAUSE IT RETURNS ONLY ROWS FOUND IN BOTH DATASETS
drugs<-inner_join(drug1,drug2, by=c("DRUG_NAME","PATHWAY_NAME","TARGETS"))

#REMOVE ALL THOSE THAT HAVE GOT AN INTERMEDIATE RESPONSE.
gdsc2<-gdsc[gdsc$BIOACTIVITY !="INTERMEDIATE",]

#IN THIS DATASET, I WANT TO REMOVE THE DUPLICATED ROWS.
#(1) ALL THE ROWS WHERE ALL THE CONTENTS ARE THE SAME,
#FROM CELL LINE  UP TO THE BIOACTIVITY, I WANT TO RANDOMMLY PICK ONE,
#AND NOT JUST PICK THE FIRST ONE.
gdsc3<-gdsc2%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

#(2)ROWS WHERE THE CELL LINE, IS TREATED BY THE SAME DRUG,AFFECTS THE SAME PATHWAY,
#HAS THE SAME TARGET BUT THE BIOACTIVITY IS DIFFERENT.
#I WILL CHOOSE THE SECOND ONE, AND OMIT THE FIRST ROW.
library(dplyr)
gdsc4<-gdsc3 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

#REORDERING THE COLUMNS
library(datawizard)
GDSC<-data_relocate(gdsc4, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_CRISPR<-data_relocate(GDSC_CRISPR, select = "PATHWAY_NAME", before = "DRUG_NAME")

#write.csv(GDSC, "C:\\Users\\School EC\\Desktop\\MSc Stuff\\Datasets\\GDSC_CLEAN.csv",row.names = F)

#=================================================================================================================================================================
#WE ARE GOING TO READ THE ACHILLES-CRSIPR DATA.
crispr<-read_csv("CRISPR_gene_effect.csv")
sample_inf<-read_csv("sample_info.csv")
sample_inf2<-sample_inf%>%
  select(c(1,2)) #wANT TO GET THE CELL_LINES AND DEPMAP_ID. 

crispr<-crispr%>%inner_join(sample_inf2,by="DepMap_ID")
crispr<-crispr%>%
  select(cell_line_name,everything())%>%
  rename(CELL_LINE_NAME =cell_line_name) 

#crispr<-crispr%>%
crispr<-crispr%>%
  select(-DepMap_ID)

#CLEANING MY VARIABLE NAMES IN THE NEW CRISPR TABLE
#I WILL START BY RESHAPING THE DATA FROM WIDE TO LONG FORMAT.
crispr2<-crispr

crispr2<-crispr2%>%
  gather(key = "GENES", value = "EXP_VAL", -CELL_LINE_NAME)

#NOW THE NEXT STEP IS TO GSUB THE PARENTHISIS IN THE GENES COLUMN,
#THIS WILL THEN BE FOLLOWED BY SPREADING THE DATA AGAN TO WIDE FROM LONG FORMAT.
library(stringr)
crispr3<-crispr2%>%
  mutate(GENES = str_remove(GENES, "\\(.*"))
write(crispr3,"CRISPR_3.csv") #need to find the script that i used to generate this.

#LETS NOW TRY TO CHANGE THE FORMAT FROM LONG TO WIDE. 
#RUN ON THE CLUSTER. CHECK RESHAPE.R
#TAKES A LOT OF COMPUTATIONAL TIME, HAVE TO RUN IT ON THE CLUSTER.
library(datawizard)
crispr4<-data_to_wide(data = crispr3,id_cols ="CELL_LINE_NAME", values_from = "EXP_VAL", names_from = "GENES")
#ON THE CLUSTER THIS THE CODE USED.
library(tidyverse)
crispr4<-read_csv("CRISPR_3.csv")
crispr4<-crispr4%>%select(-1)

library(datawizard)
crispr5<-data_to_wide(data=crispr4,id_cols =NULL, values_from ="EXP_VAL",names_from ="GENES")
write.csv(crispr5, "CRISPR_5.csv")

#================================================================================================================================================================
#GENE DEPENDENCY DATASET IS HERE. 
#PERFORM FEATURE SELECTION ON THIS DATA
#WILL DO THAT ON MATLAB.  

CRISPR5<-read_csv("CRISPR_5.csv")
CRISPR5<-CRISPR5%>%
  select(-1)
sum(is.na(CRISPR5)) # 3391 missing values

CRISPR_6<-CRISPR5
is.na(CRISPR_6$CELL_LINE_NAME) # returns true, on index 77
##removing the index 77 which has a missing value.
CRISPR_7<-CRISPR_6%>%
  slice(-c(77))
## Remove columns and rows with more than 50% NA
CRISPR_7<-CRISPR_7[which(rowMeans(!is.na(CRISPR_7)) > 0.5),
                   which(colMeans(!is.na(CRISPR_7)) > 0.5)]
#=================================================================================================================================================================
library(caret)
library(RANN)
#IMPUTE THE MISSING VALUES. 
#IMPUTING THE MISSING VALUES
CRISP_IMP<-preProcess(as.data.frame(CRISPR_7), method = 'knnImpute', k=10)
crispr_imp<-predict(CRISP_IMP,newdata = CRISPR_7)
sum(is.na(crispr_imp)) 
is.na(crispr_imp) #RETURNS NO MISSING VALUES. 

crispr_imp_2<-crispr_imp
#write.csv(crispr_imp_2, "C:\\Users\\School EC\\Desktop\\MSc Stuff\\Datasets\\IMPUTED_CRISPR.csv",row.names = F)
#===============================================================================================================================================================
#FILTER OUT GENES THAT HAVE A HIGH CORRELATION TO EACH OTHER.
#RUN IT N THE CLUSTER BECASUE IT TAKES A LOT OF COMPUTATIONAL TIME. USED the CORRDATA script.
library(caret)
library(tidyverse)
crispr_imp_3<-read_csv("IMPUTED_CRISPR.csv")
crispr_imp_3<-crispr_imp_3%>%
  column_to_rownames(var = "CELL_LINE_NAME")
#REMOVING HIGHLY CORRELATED VARIABLES (GENES)
#USE A THRESHOLD WE WANT TO DEEM CORRELATION AS TOO HIGH. CUT OFF BEING 0.8.
nonColinearData<-crispr_imp_3[, -findCorrelation(cor(crispr_imp_3), cutoff = .8)]
#write.csv(nonColinearData, "C:\\Users\\School EC\\Desktop\\MSc Stuff\\Datasets\\COL_DATA.csv",row.names = F)

#===============================================================================================================================================================
#MERGING THE GDSC AND THE ACHILLES-CRISPR DATASETS
#PERFORM FEATURE SELECTION ON THIS DATA
#WILL DO THAT ON MATLAB. 
GDSC<-read_csv("GDSC_CLEAN.csv")
crispr_imp_4<-read_csv("COL_DATA.csv") # 70 GENES REMOVED

crispr_imp_4<-crispr_imp_4%>%
  rename(CELL_LINE_NAME = ...1)
str(crispr_imp_4)
class(crispr_imp_4)

GDSC_CRISPR<-merge(x=GDSC,y=crispr_imp_4,by="CELL_LINE_NAME")
str(GDSC_CRISPR)

library(datawizard)
GDSC_CRISPR<-data_relocate(GDSC_CRISPR, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
write.csv(GDSC_CRISPR, "C:\\Users\\School EC\\Desktop\\MSc Stuff\\Datasets\\GDSC_CRISPR2.csv",row.names = F)

#===============================================================================================================================================================
#WE ARE GOING TO BE DOING FEATURE SELECTION HERE TO REDUCE THE SIZE OF OUR DATASET.
# AND WE WILL CHOOSE THE TOP 10 000 FEATURES TO TRAIN OUR MODEL.
# WE WILL PERFORM THIS FEATURE SELECTIN USING MATLAB ON THE UCT CLUSTER.
#LOAD MATLAB ON R, IS ANOTHER OPTION USING :
#system('matlab -nodisplay -r "a=2; b=1; display(a+b); exit"').

#%Select the top 10000 predictors
#features=readtable("GDSC_CRISPR.csv");
#features(1:10,1:10);
 
#%Feature selection will use the MRMR alogorithm;
#[idx,scores]=fscmrmr(features,"BIOACTIVITY");

#%Plotting the importance scores.
#bar(sscores(idx(1:15)))
#xlabel("Predictor rankings")
#ylabel("Importance score")
#title("Top 15 features)

#numoffeatures=10 000
#data=features(:,idx(1:numoffeatures))
#data(1;10,1:7)
#data.BIOACTIVITY=[]
#DATA=horzcat(features(:,1),data); 
#writetable(DATA,"TOP_10000_FEATURES.xlsx")

#==================================================================================================================================================================
GD_TRAIN<-read_csv("TOP_10000_FEATURES.csv")

unique(GD_TRAIN$CELL_LINE_NAME)
unique(GD_TRAIN$DRUG_NAME)
sum(is.na(GD_TRAIN))

library(xgboost)
library(gbm)
library(randomForest)
library(kernlab) #for SVM
library(caret)

set.seed(123)
inTraining <- createDataPartition(DATA$BIOACTIVITY, p = .70, list = FALSE)
training <- features[ inTraining,]
testing  <- features[-inTraining,]

#===================================================================================================================================================================
#TRAINING THE MODELS

# Run algorithms using 10-fold cross validation
control <- trainControl(method="repeatedcv", number=10, repeats=10)
metric <- "Accuracy"
#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
training<- as.data.frame(unclass(training),                     
                          stringsAsFactors = TRUE)

#CHANGING THE CHARACTERS INTO FACTORS VARAIBLES
testing <- as.data.frame(unclass(testing),                     
                         stringsAsFactors = TRUE)

str(training)
#KNN
set.seed(123)
fit.knn <- train(BIOACTIVITY~., data=training, method="knn", metric=metric, trControl=control)
print(fit.knn)

predictions <- predict(fit.knn, testing)
confusionMatrix(predictions, testing$BIOACTIVITY, mode ="prec_recall")

#SVM

fit.svm <- train(BIOACTIVITY~., data=training, method="svmRadial", metric=metric, trControl=control)
print(fit.svm)

predictions <- predict(fit.svm, testing)
confusionMatrix(predictions, testing$BIOACTIVITY, mode ="prec_recall")

#RF
fit.rf <- train(BIOACTIVITY~., data=training, method="rf", metric=metric, trControl=control)
print(fit.rf)

predictions <- predict(fit.rf, testing)
confusionMatrix(predictions,testing$BIOACTIVITY, mode="prec_recall")

#GBM
set.seed(123)
fit.gbm<- train(BIOACTIVITY~., data=training, method="gbm", metric=metric, trControl=control, verbose=F)
print(fit.gbm)

predictions <- predict(fit.gbm, testing)
confusionMatrix(predictions,testing$BIOACTIVITY, mode="prec_recall")

#XGBM
set.seed(123)
fit.xgbm<- train(BIOACTIVITY~., data=training, method="xgbTree", metric=metric, trControl=control)
print(fit.xgbm)

predictions <- predict(fit.xgbm, testing)
confusionMatrix(predictions,testing$BIOACTIVITY, mode ="prec_recall")
#I might actually plot the ROC/AUC. I will think about it.

#===========================================================================================================================================
#NOW TO PERFORM REGRESSION ANANLYSIS
#SVM model
library(e1071)
library(caret)

model_reg = svm(class~., data=train)
print(model_reg)

pred = predict(model_reg, test)

x=1:length(test$class)
plot(x, test$class, pch=18, col="red")
lines(x, pred, lwd="1", col="blue")

# accuracy check 
#mse = MSE(test$medv, pred)
mae = MAE(test$class, pred)
rmse = RMSE(test$class, pred)
r2 = R2(test$class, pred, form = "traditional")

cat(" MAE:", mae,  "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)

#=========================================================================================================================================================
#RANDOM FOREST
library(caret)
rf<-train(class~.,data=train,method="rf")
print(rf)

#make predictions.
pred<-predict(rf,test)
#Accuracy check
mae = MAE(test$class, pred)
rmse = RMSE(test$class, pred)
r2 = R2(test$class, pred, form = "traditional")

cat(" MAE:", mae,  "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)
#====================================================================================================================================================
#GBM(2)
gbm<-train(class~.,data=train,verbose=F)
print(gbm)
#prediction
pred2 = predict(gbm, test)
#accuracy 
#mse = MSE(test$medv, pred2)
mae = MAE(test$class, pred2)
rmse = RMSE(test$class, pred2)
r2 = R2(test$class, pred2, form = "traditional")

cat(" MAE:", mae,  "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)
#==================================================================================================================================================
#XGBM
#==================================================================================================================================================
#KNN


