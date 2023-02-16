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

crispr<-crispr%>%
  select(-DepMap_ID)

#CLEANING MY VARIABLE NAMES IN THE NEW CRISPR TABLE
#I WILL START BY RESHAPING THE DATA FROM WIDE TO LONG FORMAT.
crispr2<-crispr
crispr2<-crispr2%>%
  gather(key = "GENES", value = "EXP_VAL", -CELL_LINE_NAME)

#NOW THE NEXT STEP IS TO GSUB THE PARENTHISIS IN THE GENES COLUMN,
library(stringr)
crispr3<-crispr2%>%
  mutate(GENES = str_remove(GENES, "\\(.*"))
#write(crispr3,"CRISPR_3.csv") #need to find the script that i used to generate this.

#LETS NOW TRY TO CHANGE THE FORMAT FROM LONG TO WIDE. 
#RUN ON THE CLUSTER. CHECK RESHAPE.R
#TAKES A LOT OF COMPUTATIONAL TIME, HAVE TO RUN IT ON THE CLUSTER.
crispr4<-data_to_wide(data = crispr3,id_cols ="CELL_LINE_NAME", values_from = "EXP_VAL", names_from = "GENES")
crispr4<-crispr4%>%select(-1)
sum(is.na(crispr4))
is.na(crispr$CELL_LINE_NAME)
crispr4<-crispr4%>%
  slice(-c(77))
## Remove columns and rows with more than 50% NA
crispr4<-crispr4[which(rowMeans(!is.na(crispr4)) > 0.5),
                   which(colMeans(!is.na(crispr4)) > 0.5)]  
#IMPUTING MISSING VALUES
library(caret)
library(RANN)
CRISP_IMP<-preProcess(as.data.frame(crispr4), method = 'knnImpute', k=10)
crispr_imp<-predict(CRISP_IMP,newdata = crispr4)
sum(is.na(crispr_imp)) 
is.na(crispr_imp) #RETURNS NO MISSING VALUES. 

#FILTER OUT GENES THAT HAVE A HIGH CORRELATION TO EACH OTHER.
#RUN IT N THE CLUSTER BECASUE IT TAKES A LOT OF COMPUTATIONAL TIME. USED the CORRDATA script.
crispr_imp<-crispr_imp%>%
  column_to_rownames(var = "CELL_LINE_NAME")
#REMOVING HIGHLY CORRELATED VARIABLES (GENES)
#Might create a correlation plot.

#USE A THRESHOLD WE WANT TO DEEM CORRELATION AS TOO HIGH. CUT OFF BEING 0.8.
non_corr<-crispr_imp[, -findCorrelation(cor(crispr_imp), cutoff = .8)]

#MERGING THE GDSC AND THE ACHILLES-CRISPR DATASETS
#PERFORM FEATURE SELECTION ON THIS DATA
#WILL DO THAT ON MATLAB. 

non_corr<-non_corr%>%
  rename(CELL_LINE_NAME = ...1)
str(non_corr)
class(non_corr)

GDSC_CRISPR<-merge(x=GDSC,y=non_corr,by="CELL_LINE_NAME")
str(GDSC_CRISPR)
#REORDER THE COLUMNS OF MERGED DATA
GDSC_CRISPR<-data_relocate(GDSC_CRISPR, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
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

set.seed(123)
inTraining <- createDataPartition(DATA$BIOACTIVITY, p = .70, list = FALSE)
training <- features[ inTraining,]
testing  <- features[-inTraining,]

#===================================================================================================================================================================
#TRAINING THE MODELS
#THESE WILL BE CLASSIFICATION MODELS

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
set.seed(123)

models<-c("knn","rf","svmRadial","gbm","xgbTree")
results_table <- data.frame(models = models, stringsAsFactors = F)

for (i in models){
  model_train <- train(class~., data = training, method = i,
                     trControl= control, metric = "Accuracy")
  assign("fit", model_train)
  predictions <- predict(model_train, newdata = testing)
  precision_ <-posPredValue(predictions, testing$BIOACTIVITY)
  recall_ <- sensitivity(predictions, testing$BIOACTIVITY)
  f1 <- (2*precision_ * recall_) / (precision_ + recall_)
  
  # put that in the results table
  results_table[results_table$models %in% i, "Precision"] <- precision_
  results_table[results_table$models %in% i, "Recall"] <- recall_
  results_table[results_table$models %in% i, "F1score"] <- f1
  }
  

