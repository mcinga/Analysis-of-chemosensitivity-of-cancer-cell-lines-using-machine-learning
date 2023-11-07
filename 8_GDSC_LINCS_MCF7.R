library(tidyverse)

MCF7<-read_csv("Level4_MCF7.csv")
names(MCF7)[1]<-"GENES"
MCF7<-t(MCF7)
new_colnames <- as.character(MCF7[1, ])
MCF7<- MCF7[-1, ]
colnames(MCF7) <- new_colnames
MCF7<-MCF7 %>%
  as.data.frame()%>%
  rownames_to_column( var ="DRUG_NAME")

#==============================================================================
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

BRCA_09<-BRCA%>%
  mutate(BIOACTIVITY = case_when(Z_SCORE <= -0.9 ~ "SENSITIVE",
                                 Z_SCORE >= 0.9 ~ "RESISTANT",
                                 T ~ "INTERMEDIATE"))%>%
  select(-c(4,5))

#=========================================================================
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

GDSC<-data_relocate(gdsc_09, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
GDSC<-data_relocate(GDSC, select = "PATHWAY_NAME", before = "DRUG_NAME")

#==============================================================================
GDSC_MCF7= GDSC %>% inner_join(MCF7, by="DRUG_NAME")
write_csv(GDSC_MCF7,"D:\\MSc Stuff\\GDSC_LINCS_MCF7.csv")





