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
                                 T ~ "INTERMEDIATE"))


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

#IN THIS DATASET, I WANT TO REMOVE THE DUPLICATED ROWS.
#(1) ALL THE ROWS WHERE ALL THE CONTENTS ARE THE SAME,
#FROM CELL LINE  UP TO THE BIOACTIVITY, I WANT TO RANDOMMLY PICK ONE,
#AND NOT JUST PICK THE FIRST ONE.
gdsc2<-gdsc%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

#(2)ROWS WHERE THE CELL LINE, IS TREATED BY THE SAME DRUG,AFFECTS THE SAME PATHWAY,
#HAS THE SAME TARGET BUT THE BIOACTIVITY IS DIFFERENT.
#I WILL CHOOSE THE SECOND ONE, AND OMIT THE FIRST ROW.
library(dplyr)
gdsc3<-gdsc2 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

#write.csv(gdsc3, "C:\\Users\\School EC\\Desktop\\MSc Stuff\\Datasets\\GDSC_CLEAN.csv",row.names = F)

#=======================================================================================================================================
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
#PERFORM FEATURE SELECTION ON THIS DATA
#WILL DO THAT ON MATLAB. 

#===============================================================================================================================================================
#MERGE THE GDSC AND THE FEATURE SELECTION RESULTS SO THAT I CAN TRAIN THE MODELS.


