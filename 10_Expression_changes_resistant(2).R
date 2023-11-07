library(tidyverse)
gd <- read_csv("GDSC_LINCS_MCF7.csv")
gd<-gd%>%
  filter(BIOACTIVITY == "RESISTANT")

library(datawizard)
gd_ccle <- read_csv("GDSC_09_CCLE.csv")
#removing those cell -lines that are not present in the lincs.
#
gd_ccle<-subset(gd_ccle, !(CELL_LINE_NAME %in% c("BT-549","CAL-148","CAL-85-1",
                                                 "EFM-19","EFM-192A","HCC1500",
                                                 "HCC1569","HCC2157","HDQ-P1",
                                                 "MDA-MB-361","MDA-MB-453","MDA-MB-468")))

gd_ccle<-data_relocate(gd_ccle, select = "BIOACTIVITY", before = "CELL_LINE_NAME")
gd_ccle <-gd_ccle%>%
  filter(BIOACTIVITY == "RESISTANT")
top_genes<- read_csv("new_top300_500.csv")
top_genes<-top_genes[-107,]

# get the names of the column names into a vector:
gene_names<-top_genes$variable
#check that the LINCS dataset has these genes:
missing_columns <- setdiff(gene_names, colnames(gd))
if (length(missing_columns) > 0) {
  stop(paste("The following columns are missing in GDSC_LINCS:",
             paste(missing_columns, collapse = ", ")))
}


#Extract those that are present in the GDSC_LINCS dataset.
common_genes <- intersect(gene_names, colnames(gd))
is.vector(common_genes)
# Filter LINCS dataset for common genes
GDSC_LINCS_clean <- gd[, common_genes]
#Extract these 197 genes too from the Achilles dataset.
common_genes2 <- intersect(gene_names, colnames(gd_ccle))
Achilles<-gd_ccle[,common_genes]

#===============================================================================
common_cols <- intersect(colnames(gd), colnames(gd_ccle))

lincs <- gd[,common_cols] #lincs
ccle <- gd_ccle[,common_cols] #ccle

# Identify numeric columns
num_lincs <- sapply(lincs, is.numeric) #lincs
num_ccle<-sapply(ccle,is.numeric) #ccle
# Normalize numeric columns of the lincs
lincs_norm <- lincs
str(lincs_norm)
lincs_norm[, num_lincs] <- scale(lincs_norm[, num_lincs])

ccle_norm <-ccle
str(ccle_norm)
ccle_norm[,num_ccle]<-scale(ccle_norm[,num_ccle])

# So now both datasets have the same genes on either side,
#Now I need to find the geneas that I had found to be the most important between the two.
#Extract these 197 genes too from the Achilles dataset.
#=====LINCS=====================================================================
CG1 <- intersect(gene_names, colnames(lincs_norm))
LINCS_clean <- lincs_norm[, CG1, drop = FALSE]
l<-gd%>%
  select(c("BIOACTIVITY","CELL_LINE_NAME","PATHWAY_NAME","DRUG_NAME","TARGETS"))
LINCS_CLEAN2<-cbind(l,LINCS_clean)
#======CCLE=====================================================================
CG2<- intersect(gene_names, colnames(ccle_norm))
CCLE_clean<-ccle_norm[,CG2, drop =F]

c<-gd_ccle%>%
  select(c("BIOACTIVITY","CELL_LINE_NAME","PATHWAY_NAME","DRUG_NAME","TARGETS"))

CCLE_CLEAN2<-cbind(c,CCLE_clean)

#===============================================================================
#Genes of interest
gene_list <- c( "TMEM126B","EXD3","OLAH","XPNPEP2","TARBP1","MTERF3","SAC3D1",    
                "ABCA4","STIP1","SYN1","PTPN12","HNRNPH3","RBM19","EIF3D",     
                "NRDC","TGM5","HYAL1","FIG4","MRPS30","ANXA5","MRPL19",    
                "PDLIM2","MRPS34","RETSAT","MRPL15","GADD45GIP1","MRPL28","ADAM28",    
                "SAV1","ARHGEF2","ARAF","POLR2F","TSPAN3","EDA2R","HSD17B8",   
                "RACK1","ADAM23","KCNC2","ENPP2","EXTL2", "PCDHB13","ZNF672",    
                "NDP","ASNSD1","RBP4","SYPL1","KIAA0319L","EML3","PTPRR",     
                "EIF4EBP2","MATN2","NABP2","IRF8","COX7A2","HLF","PER3",      
                "MUS81","KLF8","ABCF3","VIP","MYEF2","COPA","CNNM1",     
                "AUH","COX6C","MRPS28","KDM4C","ACSL4","CHD1L","IGSF1",     
                "OSER1","COPB1","CPT1A","MYOZ2","ABCA12","FCGRT","WASF1",     
                "EXTL3","MIF","NDUFB7","ADGRV1","ADRA1A","IER3IP1","GADD45G",   
                "KRT85","P3H1","CAP2","AQP2","COX10","STK39","SHOC2",     
                "ANKRD12","ATP9A","FGF8","SLCO2A1","LRP12","PARD6B","SYNGR1",    
                "CNNM2","PSMB1","PRDX3","EIF4EBP1","TP63","GPM6B","UQCRQ",    
                "AUNIP","AURKAIP1","HKDC1","NDUFA3","MRPL52","CTNNBIP1","NDUFB1",    
                "ZNF33B","SH3BGRL3","CLUL1","ANKH","TNFSF14","SDS","PTDSS1",    
                "NELFE","SLC25A1","TFB2M","VANGL1","MPPED2","TNIK","MORC2",     
                "GTPBP8","PTBP3","DMBT1","TRAPPC2L","NDUFA5","FPGS","ETNK1",     
                "TLE6","GPR137B","RPL8","ICAM5","CLTCL1","HGS","CTF1",      
                "COX4I1","GPR88","POU6F1","MYH15" )

expression_data <- list(
  TMEM126B = list(ccle_before = CCLE_CLEAN2$TMEM126B,lincs_after = LINCS_CLEAN2$TMEM126B),
  EXD3 =list(ccle_before = CCLE_CLEAN2$EXD3,lincs_after = LINCS_CLEAN2$EXD3),
  OLAH =list(ccle_before = CCLE_CLEAN2$OLAH, lincs_after=LINCS_CLEAN2$OLAH),
  XPNPEP2 =list(ccle_before=CCLE_CLEAN2$XPNPEP2,lincs_after=LINCS_CLEAN2$XPNPEP2),
  TARBP1 = list(ccle_before =CCLE_CLEAN2$TARBP1, lincs_after=LINCS_CLEAN2$TARBP1),
  MTERF3 =list(ccle_before=CCLE_CLEAN2$MTERF3, lincs_after=LINCS_CLEAN2$MTERF3),
  SAC3D1 =list(ccle_before=CCLE_CLEAN2$SAC3D1, lincs_after=LINCS_CLEAN2$SAC3D1),
  ABCA4=list(ccle_before=CCLE_CLEAN2$ABCA4, lincs_after=LINCS_CLEAN2$ABCA4),
  STIP1 =list(ccle_before=CCLE_CLEAN2$STIP1, lincs_after=LINCS_CLEAN2$STIP1),
  SYN1 =list(ccle_before=CCLE_CLEAN2$SYN1, lincs_after=LINCS_CLEAN2$SYN1),
  PTPN12 =list(ccle_before=CCLE_CLEAN2$PTPN12, lincs_after=LINCS_CLEAN2$PTPN12),
  HNRNPH3 =list(ccle_before=CCLE_CLEAN2$HNRNPH3, lincs_after=LINCS_CLEAN2$HNRNPH3),
  RBM19 =list(ccle_before=CCLE_CLEAN2$RBM19, lincs_after=LINCS_CLEAN2$RBM19),
  EIF3D =list(ccle_before=CCLE_CLEAN2$EIF3D, lincs_after=LINCS_CLEAN2$EIF3D),
  NRDC =list(ccle_before=CCLE_CLEAN2$NRDC, lincs_after=LINCS_CLEAN2$NRDC),
  TGM5 =list(ccle_before=CCLE_CLEAN2$TGM5, lincs_after=LINCS_CLEAN2$TGM5),
  HYAL1 =list(ccle_before=CCLE_CLEAN2$HYAL1, lincs_after=LINCS_CLEAN2$HYAL1),
  FIG4=list(ccle_before=CCLE_CLEAN2$FIG4, lincs_after=LINCS_CLEAN2$FIG4),
  MRPS30 =list(ccle_before=CCLE_CLEAN2$MRPS30, lincs_after=LINCS_CLEAN2$MRPS30),
  ANXA5 =list(ccle_before=CCLE_CLEAN2$ANXA5, lincs_after=LINCS_CLEAN2$ANXA5),
  MRPL19 =list(ccle_before=CCLE_CLEAN2$MRPL19, lincs_after=LINCS_CLEAN2$MRPL19),
  PDLIM2=list(ccle_before=CCLE_CLEAN2$PDLIM2, lincs_after=LINCS_CLEAN2$PDLIM2),
  MRPS34 =list(ccle_before=CCLE_CLEAN2$MRPS34, lincs_after=LINCS_CLEAN2$MRPS34),
  RETSAT =list(ccle_before=CCLE_CLEAN2$RETSAT, lincs_after=LINCS_CLEAN2$RETSAT),
  MRPL15=list(ccle_before=CCLE_CLEAN2$MRPL15, lincs_after=LINCS_CLEAN2$MRPL15),
  GADD45GIP1=list(ccle_before=CCLE_CLEAN2$GADD45GIP1, lincs_after=LINCS_CLEAN2$GADD45GIP1),
  MRPL28 =list(ccle_before=CCLE_CLEAN2$MRPL28, lincs_after=LINCS_CLEAN2$MRPL28),
  ADAM28= list(ccle_before=CCLE_CLEAN2$ADAM28, lincs_after=LINCS_CLEAN2$ADAM28),
  SAV1 =list(ccle_before=CCLE_CLEAN2$SAV1, lincs_after=LINCS_CLEAN2$SAV1),
  ARHGEF2=list(ccle_before=CCLE_CLEAN2$ARHGEF2, lincs_after=LINCS_CLEAN2$ARHGEF2),
  ARAF =list(ccle_before=CCLE_CLEAN2$ARAF, lincs_after=LINCS_CLEAN2$ARAF),
  POLR2F =list(ccle_before=CCLE_CLEAN2$POLR2F, lincs_after=LINCS_CLEAN2$POLR2F),
  TSPAN3=list(ccle_before=CCLE_CLEAN2$TSPAN3, lincs_after=LINCS_CLEAN2$TSPAN3),
  EDA2R =list(ccle_before=CCLE_CLEAN2$EDA2R, lincs_after=LINCS_CLEAN2$EDA2R),
  HSD17B8=list(ccle_before=CCLE_CLEAN2$HSD17B8, lincs_after=LINCS_CLEAN2$HSD17B8),   
  RACK1 =list(ccle_before=CCLE_CLEAN2$RACK1, lincs_after=LINCS_CLEAN2$RACK1),
  ADAM23 =list(ccle_before=CCLE_CLEAN2$ADAM23, lincs_after=LINCS_CLEAN2$ADAM23),
  KCNC2=list(ccle_before=CCLE_CLEAN2$KCNC2, lincs_after=LINCS_CLEAN2$KCNC2),
  ENPP2=list(ccle_before=CCLE_CLEAN2$ENPP2, lincs_after=LINCS_CLEAN2$ENPP2),
  EXTL2 =list(ccle_before=CCLE_CLEAN2$EXTL2, lincs_after=LINCS_CLEAN2$EXTL2),
  PCDHB13=list(ccle_before=CCLE_CLEAN2$PCDHB13, lincs_after=LINCS_CLEAN2$PCDHB13),
  ZNF672=list(ccle_before=CCLE_CLEAN2$ZNF672, lincs_after=LINCS_CLEAN2$ZNF672),    
  NDP=list(ccle_before=CCLE_CLEAN2$NDP, lincs_after=LINCS_CLEAN2$NDP),
  ASNSD1=list(ccle_before=CCLE_CLEAN2$ASNSD1, lincs_after=LINCS_CLEAN2$ASNSD1),
  RBP4=list(ccle_before=CCLE_CLEAN2$RBP4, lincs_after=LINCS_CLEAN2$RBP4),
  SYPL1=list(ccle_before=CCLE_CLEAN2$SYPL1, lincs_after=LINCS_CLEAN2$SYPL1),
  KIAA0319L=list(ccle_before=CCLE_CLEAN2$KIAA0319L, lincs_after=LINCS_CLEAN2$KIAA0319L),
  EML3 =list(ccle_before=CCLE_CLEAN2$EML3, lincs_after=LINCS_CLEAN2$EML3),
  PTPRR=list(ccle_before=CCLE_CLEAN2$ABCA4, lincs_after=LINCS_CLEAN2$PTPRR),     
  EIF4EBP2 =list(ccle_before=CCLE_CLEAN2$EIF4EBP2, lincs_after=LINCS_CLEAN2$EIF4EBP2),
  MATN2=list(ccle_before=CCLE_CLEAN2$MATN2, lincs_after=LINCS_CLEAN2$MATN2),
  NABP2 =list(ccle_before=CCLE_CLEAN2$NABP2, lincs_after=LINCS_CLEAN2$NABP2),
  IRF8 =list(ccle_before=CCLE_CLEAN2$IRF8, lincs_after=LINCS_CLEAN2$IRF8),
  COX7A2=list(ccle_before=CCLE_CLEAN2$COX7A2, lincs_after=LINCS_CLEAN2$COX7A2),
  HLF=list(ccle_before=CCLE_CLEAN2$HLF, lincs_after=LINCS_CLEAN2$HLF),
  PER3 =list(ccle_before=CCLE_CLEAN2$PER3, lincs_after=LINCS_CLEAN2$PER3),      
  MUS81 =list(ccle_before=CCLE_CLEAN2$MUS81, lincs_after=LINCS_CLEAN2$MUS81),
  KLF8 =list(ccle_before=CCLE_CLEAN2$KLF8, lincs_after=LINCS_CLEAN2$KLF8),
  ABCF3= list(ccle_before=CCLE_CLEAN2$ABCF3, lincs_after=LINCS_CLEAN2$ABCF3),
  VIP=list(ccle_before=CCLE_CLEAN2$VIP, lincs_after=LINCS_CLEAN2$VIP),
  MYEF2 =list(ccle_before=CCLE_CLEAN2$MYEF2, lincs_after=LINCS_CLEAN2$MYEF2),
  COPA =list(ccle_before=CCLE_CLEAN2$COPA, lincs_after=LINCS_CLEAN2$COPA),
  CNNM1 =list(ccle_before=CCLE_CLEAN2$CNNM1, lincs_after=LINCS_CLEAN2$CNNM1),     
  AUH =list(ccle_before=CCLE_CLEAN2$AUH, lincs_after=LINCS_CLEAN2$AUH),
  COX6C =list(ccle_before=CCLE_CLEAN2$COX6C, lincs_after=LINCS_CLEAN2$COX6C),
  MRPS28 =list(ccle_before=CCLE_CLEAN2$MRPS28, lincs_after=LINCS_CLEAN2$MRPS28),
  KDM4C =list(ccle_before=CCLE_CLEAN2$KDM4C, lincs_after=LINCS_CLEAN2$KDM4C),
  ACSL4 =list(ccle_before=CCLE_CLEAN2$ACSL4, lincs_after=LINCS_CLEAN2$ACSL4),
  CHD1L=list(ccle_before=CCLE_CLEAN2$CHD1L, lincs_after=LINCS_CLEAN2$CHD1L),
  IGSF1=list(ccle_before=CCLE_CLEAN2$IGSF1, lincs_after=LINCS_CLEAN2$IGSF1),     
  OSER1=list(ccle_before=CCLE_CLEAN2$OSER1, lincs_after=LINCS_CLEAN2$OSER1),
  COPB1 =list(ccle_before=CCLE_CLEAN2$COPB1, lincs_after=LINCS_CLEAN2$COPB1),
  CPT1A =list(ccle_before=CCLE_CLEAN2$CPT1A, lincs_after=LINCS_CLEAN2$CPT1A),
  MYOZ2 =list(ccle_before=CCLE_CLEAN2$MYOZ2, lincs_after=LINCS_CLEAN2$MYOZ2),
  ABCA12 =list(ccle_before=CCLE_CLEAN2$ABCA12, lincs_after=LINCS_CLEAN2$ABCA12),
  FCGRT =list(ccle_before=CCLE_CLEAN2$FCGRT, lincs_after=LINCS_CLEAN2$FCGRT),
  WASF1 =list(ccle_before=CCLE_CLEAN2$WASF1, lincs_after=LINCS_CLEAN2$WASF1),     
  EXTL3 =list(ccle_before=CCLE_CLEAN2$EXTL3, lincs_after=LINCS_CLEAN2$EXTL3),
  MIF =list(ccle_before=CCLE_CLEAN2$MIF, lincs_after=LINCS_CLEAN2$MIF),
  NDUFB7 =list(ccle_before=CCLE_CLEAN2$NDUFB7, lincs_after=LINCS_CLEAN2$NDUFB7),
  ADGRV1 =list(ccle_before=CCLE_CLEAN2$ADGRV1, lincs_after=LINCS_CLEAN2$ADGRV1),
  ADRA1A =list(ccle_before=CCLE_CLEAN2$ADRA1A, lincs_after=LINCS_CLEAN2$ADRA1A),
  IER3IP1 =list(ccle_before=CCLE_CLEAN2$IER3IP1, lincs_after=LINCS_CLEAN2$IER3IP1),
  GADD45G =list(ccle_before=CCLE_CLEAN2$GADD45G, lincs_after=LINCS_CLEAN2$GADD45G),   
  KRT85 =list(ccle_before=CCLE_CLEAN2$KRT85, lincs_after=LINCS_CLEAN2$KRT85),
  P3H1 =list(ccle_before=CCLE_CLEAN2$P3H1, lincs_after=LINCS_CLEAN2$P3H1),
  CAP2 =list(ccle_before=CCLE_CLEAN2$CAP2, lincs_after=LINCS_CLEAN2$CAP2),
  AQP2 =list(ccle_before=CCLE_CLEAN2$AQP2, lincs_after=LINCS_CLEAN2$AQP2),
  COX10 =list(ccle_before=CCLE_CLEAN2$COX10, lincs_after=LINCS_CLEAN2$COX10),
  STK39 =list(ccle_before=CCLE_CLEAN2$STK39, lincs_after=LINCS_CLEAN2$STK39),
  SHOC2 =list(ccle_before=CCLE_CLEAN2$SHOC2, lincs_after=LINCS_CLEAN2$SHOC2),     
  ANKRD12 =list(ccle_before=CCLE_CLEAN2$ANKRD12, lincs_after=LINCS_CLEAN2$ANKRD12),
  ATP9A =list(ccle_before=CCLE_CLEAN2$ATP9A, lincs_after=LINCS_CLEAN2$ATP9A),
  FGF8 =list(ccle_before=CCLE_CLEAN2$FGF8, lincs_after=LINCS_CLEAN2$FGF8),
  SLCO2A1 =list(ccle_before=CCLE_CLEAN2$SLCO2A1, lincs_after=LINCS_CLEAN2$SLCO2A1),
  LRP12=list(ccle_before=CCLE_CLEAN2$LRP12, lincs_after=LINCS_CLEAN2$LRP12),
  PARD6B =list(ccle_before=CCLE_CLEAN2$PARD6B, lincs_after=LINCS_CLEAN2$PARD6B),
  SYNGR1=list(ccle_before=CCLE_CLEAN2$SYNGR1, lincs_after=LINCS_CLEAN2$SYNGR1),    
  CNNM2 =list(ccle_before=CCLE_CLEAN2$CNNM2, lincs_after=LINCS_CLEAN2$CNNM2),
  PSMB1=list(ccle_before=CCLE_CLEAN2$PSMB1, lincs_after=LINCS_CLEAN2$PSMB1),
  PRDX3=list(ccle_before=CCLE_CLEAN2$PRDX3, lincs_after=LINCS_CLEAN2$PRDX3),
  EIF4EBP1 =list(ccle_before=CCLE_CLEAN2$EIF4EBP1, lincs_after=LINCS_CLEAN2$EIF4EBP1),
  TP63=list(ccle_before=CCLE_CLEAN2$TP63, lincs_after=LINCS_CLEAN2$TP63),
  GPM6B=list(ccle_before=CCLE_CLEAN2$GPM6B, lincs_after=LINCS_CLEAN2$GPM6B),
  UQCRQ=list(ccle_before=CCLE_CLEAN2$ABCA4, lincs_after=LINCS_CLEAN2$ABCA4),    
  AUNIP =list(ccle_before=CCLE_CLEAN2$AUNIP, lincs_after=LINCS_CLEAN2$AUNIP),
  AURKAIP1=list(ccle_before=CCLE_CLEAN2$AURKAIP1, lincs_after=LINCS_CLEAN2$AURKAIP1),
  HKDC1=list(ccle_before=CCLE_CLEAN2$HKDC1, lincs_after=LINCS_CLEAN2$HKDC1),
  NDUFA3=list(ccle_before=CCLE_CLEAN2$NDUFA3, lincs_after=LINCS_CLEAN2$NDUFA3),
  MRPL52 =list(ccle_before=CCLE_CLEAN2$MRPL52, lincs_after=LINCS_CLEAN2$MRPL52),
  CTNNBIP1=list(ccle_before=CCLE_CLEAN2$CTNNBIP1, lincs_after=LINCS_CLEAN2$CTNNBIP1),
  NDUFB1 =list(ccle_before=CCLE_CLEAN2$NDUFB1, lincs_after=LINCS_CLEAN2$NDUFB1),    
  ZNF33B=list(ccle_before=CCLE_CLEAN2$ZNF33B, lincs_after=LINCS_CLEAN2$ZNF33B),
  SH3BGRL3 =list(ccle_before=CCLE_CLEAN2$SH3BGRL3, lincs_after=LINCS_CLEAN2$SH3BGRL3),
  CLUL1=list(ccle_before=CCLE_CLEAN2$CLUL1, lincs_after=LINCS_CLEAN2$CLUL1),
  ANKH =list(ccle_before=CCLE_CLEAN2$ANKH, lincs_after=LINCS_CLEAN2$ANKH),
  TNFSF14=list(ccle_before=CCLE_CLEAN2$TNFSF14, lincs_after=LINCS_CLEAN2$TNFSF14),
  SDS=list(ccle_before=CCLE_CLEAN2$SDS, lincs_after=LINCS_CLEAN2$SDS),
  PTDSS1 =list(ccle_before=CCLE_CLEAN2$PTDSS1, lincs_after=LINCS_CLEAN2$PTDSS1),    
  NELFE =list(ccle_before=CCLE_CLEAN2$NELFE, lincs_after=LINCS_CLEAN2$NELFE),
  SLC25A1 =list(ccle_before=CCLE_CLEAN2$SLC25A1, lincs_after=LINCS_CLEAN2$SLC25A1),
  TFB2M=list(ccle_before=CCLE_CLEAN2$TFB2M, lincs_after=LINCS_CLEAN2$TFB2M),
  VANGL1 =list(ccle_before=CCLE_CLEAN2$VANGL1, lincs_after=LINCS_CLEAN2$VANGL1),
  MPPED2=list(ccle_before=CCLE_CLEAN2$MPPED2, lincs_after=LINCS_CLEAN2$MPPED2),
  TNIK=list(ccle_before=CCLE_CLEAN2$TNIK, lincs_after=LINCS_CLEAN2$TNIK),
  MORC2=list(ccle_before=CCLE_CLEAN2$MORC2, lincs_after=LINCS_CLEAN2$MORC2),     
  GTPBP8=list(ccle_before=CCLE_CLEAN2$GTPBP8, lincs_after=LINCS_CLEAN2$GTPBP8),
  PTBP3=list(ccle_before=CCLE_CLEAN2$PTBP3, lincs_after=LINCS_CLEAN2$PTBP3),
  DMBT1=list(ccle_before=CCLE_CLEAN2$DMBT1, lincs_after=LINCS_CLEAN2$DMBT1),
  TRAPPC2L=list(ccle_before=CCLE_CLEAN2$TRAPPC2L, lincs_after=LINCS_CLEAN2$TRAPPC2L),
  NDUFA5=list(ccle_before=CCLE_CLEAN2$NDUFA5, lincs_after=LINCS_CLEAN2$ABCA4),
  FPGS=list(ccle_before=CCLE_CLEAN2$FPGS, lincs_after=LINCS_CLEAN2$FPGS),
  ETNK1=list(ccle_before=CCLE_CLEAN2$ETNK1, lincs_after=LINCS_CLEAN2$ETNK1),     
  TLE6=list(ccle_before=CCLE_CLEAN2$TLE6, lincs_after=LINCS_CLEAN2$TLE6),
  GPR137B =list(ccle_before=CCLE_CLEAN2$GPR137B, lincs_after=LINCS_CLEAN2$GPR137B),
  RPL8 =list(ccle_before=CCLE_CLEAN2$ABCA4, lincs_after=LINCS_CLEAN2$ABCA4),
  ICAM5 =list(ccle_before=CCLE_CLEAN2$ICAM5, lincs_after=LINCS_CLEAN2$ICAM5),
  CLTCL1=list(ccle_before=CCLE_CLEAN2$CLTCL1, lincs_after=LINCS_CLEAN2$CLTCL1),
  HGS =list(ccle_before=CCLE_CLEAN2$HGS, lincs_after=LINCS_CLEAN2$HGS),
  CTF1=list(ccle_before=CCLE_CLEAN2$CTF1, lincs_after=LINCS_CLEAN2$CTF1),      
  COX4I1 =list(ccle_before=CCLE_CLEAN2$COX4I1, lincs_after=LINCS_CLEAN2$COX4I1),
  GPR88 =list(ccle_before=CCLE_CLEAN2$GPR88, lincs_after=LINCS_CLEAN2$GPR88),
  POU6F1=list(ccle_before=CCLE_CLEAN2$POU6F1, lincs_after=LINCS_CLEAN2$POU6F1),
  MYH15=list(ccle_before=CCLE_CLEAN2$MYH15, lincs_after=LINCS_CLEAN2$MYH15))

# Create an empty list to store the results for each gene
results_list <- list()


# Loop over the genes
for (gene in gene_list) {
  # Check if both datasets have at least one expression value
  if (length(expression_data[[gene]]$ccle_before) > 0 && length(expression_data[[gene]]$lincs_after) > 0) {
    # Subset the gene expression data for the specific gene
    ccle_before <- expression_data[[gene]]$ccle_before
    lincs_after <- expression_data[[gene]]$lincs_after
    
    # Determine the common number of samples (rows) between the datasets
    num_samples <- min(length(ccle_before), length(lincs_after))
    
    # Subset the datasets to have the same number of samples
    ccle_before <- ccle_before[1:num_samples]
    lincs_after <- lincs_after[1:num_samples]
    
    # Perform t-test
    t_result <- t.test(ccle_before, lincs_after)
    
    # Extract the relevant information
    p_value <- t_result$p.value
    mean_difference <- mean(lincs_after) - mean(ccle_before)
    
    # Store the results in a data frame
    gene_result <- data.frame(
      Gene = gene,
      Mean_Difference = mean_difference,
      P_Value = p_value,
      stringsAsFactors = FALSE
    )
    
    # Add the gene result to the list
    results_list[[gene]] <- gene_result
    
    # Create a data frame with aligned expression values
    aligned_data <- data.frame(
      Gene = rep(gene, num_samples * 2),
      Treatment = c(rep("pre_treatment", num_samples), rep("post_treatment", num_samples)),
      Expression = c(ccle_before, lincs_after),
      stringsAsFactors = FALSE
    )
    
    # Plot the box plot
    png(paste("boxplot_resistant",gene,".png",sep =""))
    boxplot(Expression ~ Treatment, data = aligned_data, main = paste("Gene:", gene), xlab = "Treatment", ylab = "Expression")
    dev.off()
  } else {
    cat("No expression values found for gene", gene, "\n")
  }
}


# Combine the results for all genes into a single data frame
results_table <- do.call(rbind, results_list)
write_csv(results_table,"D:\\Msc Stuff\\Gene_EXP_change_resistant(2).csv")
