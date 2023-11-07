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

gdsc_09<-gdsc_09[gdsc_09$BIOACTIVITY !="INTERMEDIATE",]

# REMOVE THE DUPLICATED ROWS.
#(1) ALL THE ROWS WHERE ALL THE CONTENTS ARE THE SAME,
#FROM CELL LINE UP TO THE BIOACTIVITY, I WANT TO PICK ONE,
#AND NOT JUST PICK THE FIRST ONE.

gdsc_09<-gdsc_09%>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,BIOACTIVITY,TARGETS)%>%
  slice_sample(n=1)

#(2)ROWS WHERE THE CELL LINE, IS TREATED BY THE SAME DRUG,AFFECTS THE SAME PATHWAY,
#HAS THE SAME TARGET BUT THE BIOACTIVITY IS DIFFERENT.
#I WILL CHOOSE THE SECOND ONE, AND OMIT THE FIRST ROW.
gdsc_09<-gdsc_09 %>%
  group_by(CELL_LINE_NAME,DRUG_NAME,PATHWAY_NAME,TARGETS) %>%
  slice(min(2, n())) %>%
  ungroup

#REORDERING THE COLUMNS
library(datawizard)

GDSC_09<-data_relocate(gdsc_09, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC_09<-data_relocate(GDSC_09, select = "PATHWAY_NAME", before = "DRUG_NAME")

#WE ARE GOING TO READ THE CCLE-GENE EXPRESSION DATA.
ccle<-read_csv("CCLE_expression.csv")
sample_inf<-read.csv("sample_info.csv")
sample_inf<-sample_inf%>%
  select(c(1,2)) #wANT TO GET THE CELL_LINES AND DEPMAP_ID.

ccle<-ccle%>%
  rename(DepMap_ID = ...1)

ccle<-ccle%>%
  inner_join(sample_inf,by="DepMap_ID")%>%
  select(cell_line_name,everything())%>%
  rename(CELL_LINE_NAME =cell_line_name)

ccle<-ccle%>%
  select(-DepMap_ID)
#removing special characters in the column names
colnames(ccle) <- sub("\\s+\\(\\d+.*", "", colnames(ccle))
  
sum(is.na(ccle)) 
is.na(ccle$CELL_LINE_NAME) 
sum(is.na(ccle$CELL_LINE_NAME))
#Removing those cell lines that have NA's in the celline column.
ccle <- ccle[!is.na(ccle$CELL_LINE_NAME),]
## REMOVE COLUMNS AND ROWS WITH MORE THAN 50% NA's.
ccle<-ccle[which(rowMeans(!is.na(ccle)) > 0.5),
               which(colMeans(!is.na(ccle)) > 0.5)]

is.na(ccle)
which(is.na(ccle))

# Check the character length of each value in the column
lengths <- nchar(ccle$CELL_LINE_NAME)

# Identify rows with values that have length 0 or only contain whitespace
empty_rows <- which(lengths == 0 | grepl("^\\s*$", ccle$CELL_LINE_NAME))
#removing those rows that had empty or whitesapces.
ccle <- ccle[-empty_rows, ]

#FILTER OUT GENES THAT HAVE A HIGH CORRELATION TO EACH OTHER.
#RUN IT N THE CLUSTER BECASUE IT TAKES A LOT OF COMPUTATIONAL TIME.
# Check for duplicate row names
duplicated_rows <- duplicated(ccle$CELL_LINE_NAME)
duplicated_names <- ccle$CELL_LINE_NAME[duplicated_rows]

#Name of cell line that is a duplicate I want to rremove.
dup_cell<- "U-251 MG"
# Remove the duplicate's first occurrence while keeping the second one
ccle2<- ccle %>%
  group_by(CELL_LINE_NAME) %>%
  mutate(dup_count = row_number()) %>%
  filter(!(CELL_LINE_NAME == dup_cell & dup_count == 1)) %>%
  ungroup() %>%
  select(-dup_count)

#REMOVING HIGHLY CORRELATED VARIABLES (GENES)
#USE A THRESHOLD WE WANT TO DEEM CORRELATION AS TOO HIGH. CUT OFF BEING 0.7.

ccle2 <- ccle %>%
  mutate(row_id = make.unique(CELL_LINE_NAME)) %>%
  column_to_rownames(var = "row_id")

ccle2<-ccle2%>%
  select(-1)

ccle3<-cor(ccle2)
hc<-findCorrelation(ccle3,cutoff=0.7)
hc<-sort(hc)
correlated_Data<-ccle2[,-c(hc)]

ccle_cell_name<-ccle%>%
  select(CELL_LINE_NAME)

newData<-data.frame(ccle_cell_name, correlated_Data)

write_csv(correlated_Data, "ccle_correlation.csv")
write_csv(newData, "new_corr_ccle_data.csv")
#===============================================================================
GDSC_09_CCLE<-merge.data.frame(x=GDSC_09,y=newData,by="CELL_LINE_NAME")
#Reordering the columns.
GDSC_09_CCCLE<-data_relocate(GDSC_09_CRISPR, select= "BIOACTIVITY", before = "CELL_LINE_NAME")
#
write_csv(GDSC_09_CCLE,"GDSC_09_CCLE.csv")

dim(newData)
