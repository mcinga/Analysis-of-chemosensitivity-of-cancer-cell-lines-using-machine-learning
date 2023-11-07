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

#RETAIN ONLY THOSE THAT ARE NOT DUPLICATES IN THE DATA.

library(sqldf)
BRCA<-sqldf("select * from BRCA1 union select * from BRCA2")
unique(BRCA$CELL_LINE_NAME) #51 CELL LINES

#SEPERATING THE DRUGS ACCORDING TO SENSISTIVITY, 
#USING THE Z_SCORE. WE WILL START WITH 0.5 AND SEE THE PERFORMANCE OF MODEL


BRCA_09<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.9 ~ "SENSITIVE",
                                 Z_SCORE >= 0.9 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_08<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.8 ~ "SENSITIVE",
                                 Z_SCORE >= 0.8 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_07<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.7 ~ "SENSITIVE",
                                 Z_SCORE >= 0.7 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_06<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.6 ~ "SENSITIVE",
                                 Z_SCORE >= 0.6 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_05<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.5 ~ "SENSITIVE",
                                 Z_SCORE >= 0.5 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_04<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.4 ~ "SENSITIVE",
                                 Z_SCORE >= 0.4 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_03<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.3 ~ "SENSITIVE",
                                 Z_SCORE >= 0.3 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_02<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.2 ~ "SENSITIVE",
                                 Z_SCORE >= 0.2 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

BRCA_01<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.1 ~ "SENSITIVE",
                                 Z_SCORE >= 0.1 ~ "RESISTANT",
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

#NOW TO MERGE THESE TWO DATASETS AND CREATE ONE GDSC DATASET THAT HAS DRUG INFORMATION TOO
gdsc_09<-inner_join(BRCA_09,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_08<-inner_join(BRCA_08,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_07<-inner_join(BRCA_07,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_06<-inner_join(BRCA_06,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_05<-inner_join(BRCA_05,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_04<-inner_join(BRCA_04,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_03<-inner_join(BRCA_03,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_02<-inner_join(BRCA_02,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))
gdsc_01<-inner_join(BRCA_01,drugs,by=c("DRUG_NAME","PATHWAY_NAME"))

#DELETING ROWS THAT HAVE OR SHOW AN INTERMEDIATE DRUG RESPONSE.

gdsc_09<-gdsc_09[gdsc_09$BIOACTIVITY !="INTERMEDIATE",]
gdsc_08<-gdsc_08[gdsc_08$BIOACTIVITY !="INTERMEDIATE",]
gdsc_07<-gdsc_07[gdsc_07$BIOACTIVITY !="INTERMEDIATE",]
gdsc_06<-gdsc_06[gdsc_06$BIOACTIVITY !="INTERMEDIATE",]
gdsc_05<-gdsc_05[gdsc_05$BIOACTIVITY !="INTERMEDIATE",]
gdsc_04<-gdsc_04[gdsc_04$BIOACTIVITY !="INTERMEDIATE",]
gdsc_03<-gdsc_05[gdsc_03$BIOACTIVITY !="INTERMEDIATE",]
gdsc_02<-gdsc_05[gdsc_02$BIOACTIVITY !="INTERMEDIATE",]
gdsc_01<-gdsc_05[gdsc_01$BIOACTIVITY !="INTERMEDIATE",]

# REMOVE THE DUPLICATED ROWS.
#(1) ALL THE ROWS WHERE ALL THE CONTENTS ARE THE SAME,
#FROM CELL LINE UP TO THE BIOACTIVITY, I WANT TO PICK ONE,
#AND NOT JUST PICK THE FIRST ONE.

gdsc_09<-gdsc_09%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_08<-gdsc_08%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_07<-gdsc_07%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_06<-gdsc_06%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_05<-gdsc_05%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_04<-gdsc_04%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_03<-gdsc_03%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_02<-gdsc_02%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

gdsc_01<-gdsc_01%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

#(2)ROWS WHERE THE CELL LINE, IS TREATED BY THE SAME DRUG,AFFECTS THE SAME PATHWAY,
#HAS THE SAME TARGET BUT THE BIOACTIVITY IS DIFFERENT.
#I WILL CHOOSE THE SECOND ONE, AND OMIT THE FIRST ROW.
gdsc_09<-gdsc_09 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_08<-gdsc_08 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_07<-gdsc_07 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_06<-gdsc_06 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_05<-gdsc_05 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_04<-gdsc_04 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_03<-gdsc_03 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_02<-gdsc_02 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

gdsc_01<-gdsc_01 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

#REORDERING THE COLUMNS
library(datawizard)

GDSC_09<-data_relocate(gdsc_09, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_09<-data_relocate(GDSC_09, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_08<-data_relocate(gdsc_08, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_08<-data_relocate(GDSC_08, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_07<-data_relocate(gdsc_07, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_07<-data_relocate(GDSC_07, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_06<-data_relocate(gdsc_06, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_06<-data_relocate(GDSC_06, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_05<-data_relocate(gdsc_05, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_05<-data_relocate(GDSC_05, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_04<-data_relocate(gdsc_04, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_04<-data_relocate(GDSC_04, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_03<-data_relocate(gdsc_03, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_03<-data_relocate(GDSC_03, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_02<-data_relocate(gdsc_02, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_02<-data_relocate(GDSC_02, select = "PATHWAY_NAME", before = "DRUG_NAME")

GDSC_01<-data_relocate(gdsc_01, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_01<-data_relocate(GDSC_01, select = "PATHWAY_NAME", before = "DRUG_NAME")

#Checking for class imbalance
GDSC_05<- as.data.frame(unclass(GDSC_05),                     
                         stringsAsFactors = TRUE)

pdf("CLASS_DISTRIBUTION_01_05.pdf")
par(mfrow=c(3,2))
barplot(prop.table(table(GDSC_05$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.5" )


barplot(prop.table(table(GDSC_04$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.4" )


barplot(prop.table(table(GDSC_03$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.3" )

barplot(prop.table(table(GDSC_02$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.2" )

barplot(prop.table(table(GDSC_01$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.1" )
dev.off()

pdf("CLASS_DISTRIBUTION_06_09.pdf")
par(mfrow=c(3,2))
barplot(prop.table(table(GDSC_06$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.6" )


barplot(prop.table(table(GDSC_07$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.7" )


barplot(prop.table(table(GDSC_08$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.8" )

barplot(prop.table(table(GDSC_09$BIOACTIVITY)),
        col = rainbow(2),ylim = c(0,07), main = "Class Distribution at 0.9" )
dev.off()
#===============================================================================
#WE ARE GOING TO READ THE ACHILLES-CRSIPR DATA.
crispr<-read.csv("CRISPR_gene_effect.csv")
sample_inf<-read.csv("sample_info.csv")
sample_inf<-sample_inf%>%
  select(c(1,2)) #wANT TO GET THE CELL_LINES AND DEPMAP_ID. 

crispr<-crispr%>%
  inner_join(sample_inf,by="DepMap_ID")%>%
  select(cell_line_name,everything())%>%
  rename(CELL_LINE_NAME =cell_line_name) 

crispr<-crispr%>%
  select(-DepMap_ID)

#CLEANING MY VARIABLE NAMES IN THE NEW CRISPR TABLE
#I WILL START BY RESHAPING THE DATA FROM WIDE TO LONG FORMAT.
library(tidyr)
crispr<-crispr%>%
  gather(key = "GENES", value = "EXP_VAL", -CELL_LINE_NAME)

#NOW THE NEXT STEP IS TO GSUB THE PARENTHISIS IN THE GENES COLUMN,
#THIS WILL THEN BE FOLLOWED BY SPREADING THE DATA AGAN TO WIDE FROM LONG FORMAT.
library(stringr)
crispr<-crispr%>%
  mutate(GENES = str_remove(GENES, "\\..*"))

#write(crispr3,"CRISPR_3.csv") #need to find the script that i used to generate this.
#LETS NOW TRY TO CHANGE THE FORMAT FROM LONG TO WIDE. 
#RUN ON THE CLUSTER. CHECK RESHAPE.R
#TAKES A LOT OF COMPUTATIONAL TIME, HAVE TO RUN IT ON THE CLUSTER.
library(datawizard)
crispr<-data_to_wide(data = crispr,id_cols =NULL, values_from = "EXP_VAL", names_from = "GENES")
#write_csv(crispr,"crispr_wide_3.csv")

#========================================================================================================================================
#PERFORM FEATURE SELECTION ON THIS DATA
#WILL PERFORM THIS USING MATLAB TERMINAL ON THE CLUSTER.
#========================================================================================================================================
#GENE DEPENDENCY DATASET IS HERE. 
sum(is.na(crispr)) 
is.na(crispr$CELL_LINE_NAME) 
##removing the index 77 which has a missing value.
crispr<-crispr%>%
  slice(-c(38,77,587))

## REMOVE COLUMNS AND ROWS WITH MORE THAN 50% NA's.
crispr<-crispr[which(rowMeans(!is.na(crispr)) > 0.5),
               which(colMeans(!is.na(crispr)) > 0.5)]

#=========================================================================================
library(caret)
library(RANN)
#IMPUTE THE MISSING VALUES. 
#IMPUTING THE MISSING VALUES
imputed<-preProcess(as.data.frame(crispr), method = 'knnImpute', k=10)
crispr_imp<-predict(imputed,newdata = crispr)
sum(is.na(crispr_imp)) 
is.na(crispr_imp) #RETURNS NO MISSING VALUES. 

#write_csv(crispr_imp,"crispr_imp_3.csv")

#========================================================================================================================================
#FILTER OUT GENES THAT HAVE A HIGH CORRELATION TO EACH OTHER.
#RUN IT N THE CLUSTER BECASUE IT TAKES A LOT OF COMPUTATIONAL TIME.
crispr_imp<-crispr_imp%>%
  column_to_rownames(var = "CELL_LINE_NAME")

#REMOVING HIGHLY CORRELATED VARIABLES (GENES)
#USE A THRESHOLD WE WANT TO DEEM CORRELATION AS TOO HIGH. CUT OFF BEING 0.7.
crispr_imp2<-cor(crispr_imp)
hc<-findCorrelation(crispr_imp2,cutoff=0.7)
hc<-sort(hc)
correlated_Data<-crispr_imp[,-c(hc)]

crispr_cell_name<-crispr%>%
  select(CELL_LINE_NAME)

newData<-data.frame(crispr_cell_name, correlated_Data)
#newData<-merge(x=crispr_cell_name,y=correlated_Data)

dim(newData)
#==========================================================================================================

GDSC_09_CRISPR<-merge.data.frame(x=GDSC_09,y=newData,by="CELL_LINE_NAME")
GDSC_08_CRISPR<-merge.data.frame(x=GDSC_08,y=newData,by="CELL_LINE_NAME")
GDSC_07_CRISPR<-merge.data.frame(x=GDSC_07,y=newData,by="CELL_LINE_NAME")
GDSC_06_CRISPR<-merge.data.frame(x=GDSC_06,y=newData,by="CELL_LINE_NAME")
GDSC_05_CRISPR<-merge.data.frame(x=GDSC_05,y=newData,by="CELL_LINE_NAME")
GDSC_04_CRISPR<-merge.data.frame(x=GDSC_04,y=newData,by="CELL_LINE_NAME")
GDSC_03_CRISPR<-merge.data.frame(x=GDSC_03,y=newData,by="CELL_LINE_NAME")
GDSC_02_CRISPR<-merge.data.frame(x=GDSC_02,y=newData,by="CELL_LINE_NAME")
GDSC_01_CRISPR<-merge.data.frame(x=GDSC_01,y=newData,by="CELL_LINE_NAME")

#REORDERING THE COLUMNS

GDSC_09_CRISPR<-data_relocate(GDSC_09_CRISPR, select= "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_08_CRISPR<-data_relocate(GDSC_08_CRISPR, select= "BIOACTIVITY", before="CELL_LINE_NAME")
GDSC_07_CRISPR<-data_relocate(GDSC_07_CRISPR, select= "BIOACTIVITY",before = "CELL_LINE_NAME") 
GDSC_06_CRISPR<-data_relocate(GDSC_06_CRISPR, select= "BIOACTIVITY", before="CELL_LINE_NAME")
GDSC_05_CRISPR<-data_relocate(GDSC_05_CRISPR, select ="BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_04_CRISPR<-data_relocate(GDSC_04_CRISPR, select= "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_03_CRISPR<-data_relocate(GDSC_03_CRISPR, select= "BIOACTIVITY", before="CELL_LINE_NAME")
GDSC_02_CRISPR<-data_relocate(GDSC_02_CRISPR, select= "BIOACTIVITY",before = "CELL_LINE_NAME") 
GDSC_01_CRISPR<-data_relocate(GDSC_01_CRISPR, select= "BIOACTIVITY", before="CELL_LINE_NAME")

#===============================================================================

write_csv(GDSC_09_CRISPR,"GDSC_09_CRISPR.csv")
write_csv(GDSC_08_CRISPR,"GDSC_08_CRISPR.csv")
write_csv(GDSC_07_CRISPR,"GDSC_07_CRISPR.csv")
write_csv(GDSC_06_CRISPR,"GDSC_06_CRISPR.csv")
write_csv(GDSC_05_CRISPR,"GDSC_05_CRISPR.csv")
write_csv(GDSC_04_CRISPR,"GDSC_04_CRISPR.csv")
write_csv(GDSC_03_CRISPR,"GDSC_03_CRISPR.csv")
write_csv(GDSC_02_CRISPR,"GDSC_02_CRISPR.csv")
write_csv(GDSC_01_CRISPR,"GDSC_01_CRISPR.csv")


#FINDING THE DIMENSIONS
unique(GDSC_05_CRISPR$CELL_LINE_NAME)
unique(GDSC_05_CRISPR$DRUG_NAME)
#===============================================================================
