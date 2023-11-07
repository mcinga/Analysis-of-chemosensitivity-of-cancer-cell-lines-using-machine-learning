library(clusterProfiler)
library(tidyverse)
gene_list <- c("TMEM126B", "RPL26", "EXD3", "OLAH", "XPNPEP2","TARBP1",
               "MTERF3", "SAC3D1","MS4A1","ABCA4","STIP1","SYN1","ELK1","PTPN12",
               "SNIP","HNRNPH3","RBM19",
               "EIF3D","NRDC","TGM5","TROAP","CIITA","HYAL1","FIG4",
               "MRPS30","ANXA5","MRPL19","PDLIM2","MRPS34","RETSAT","MRPL15",
               "ELAVL1","GADD45GIP1","MRPL28","GFI1B","ADAM28","SAV1","PAXBP1",
               "ARHGEF2","ARAF","PSMA2","SAE1","POLR2F","TSPAN3","EDA2R","RPS11",
               "HSD17B8","SPRR1B","RACK1","ADAM23","PDE6H","RPL30", "KCNC2","ENPP2",
               "EXTL2","PCDHB13","ZNF672","NDP", "ASNSD1","ZZEF1","PRIM1","RBM25",
               "ZNF7","RBP4","SYPL1","KIAA0319L","EML3","PTPRR","EIF4EBP2","ACAP2",
               "MATN2","NABP2","IRF8","COX7A2","HLF","RPS9","PER3","MUS81","KLF8",
               "ABCF3","VIP","SMARCA5","ATAD5","PMPCA","MYEF2","DCLRE1B","OSBPL10",
               "COPA","CNNM1","AUH","COX6C","MRPS28","SEZ6L2","KDM4C","ACSL4","CHD1L",
               "IGSF1","OSER1","COPB1","CPT1A","MYOZ2","ABCA12","INTS7","FCGRT","WASF1",
               "EXTL3","MIF","NDUFB7","ADGRV1","PDC","ADRA1A","IER3IP1","TTI1","TET3",
               "GADD45G","KRT85","P3H1","CAP2","AQP2","COX10","NRXN1","STK39","SHOC2",
               "ANKRD12","ATP9A","FGF8","SLCO2A1","LRP12","PARD6B","SYNGR1","CNNM2","PSMB1",
               "PRDX3","EIF4EBP1","TP63","NCBP1","HOXB3","GPM6B","UQCRQ","AUNIP","AURKAIP1",
               "HKDC1","NDUFA3","ZNF124","MRPL52","BTK","CTNNBIP1","PRR14L","AOC2","CCT8","NDUFB1",
               "ZNF33B","RBM15","SH3BGRL3","CLUL1","ANKH","TNFSF14","IRF4","SDS","PTDSS1","RPS20",
               "NELFE","SLC25A1","TRIM31","TFB2M","VANGL1","MPPED2","TNIK","SENP5","MORC2",
               "LIG3","GTPBP8","PTBP3","DMBT1","TRAPPC2L","NDUFA5","FPGS","ETNK1","HMGB2","TBRG4",
               "TLE6","GPR137B","SHARPIN","POLA1","RBBP4","RPL8","ICAM5","CLTCL1","PLEKHA1","DOK3",
               "HGS","CTF1","COX4I1","GPR88","POU6F1","MYH15",
               "RPUSD3", "CXorf38", "GNB4", "C4orf46", "POM121L12", "TAS2R30", "SP9", "DOCK7",
               "ZNF608", "SRFBP1", "ZNF607", "KBTBD3", "POC5", "MRPL38", "SMARCAD1", "MEDAG", 
               "KRT82", "ALYREF", "MRM2", "ODF2L", "NLRC4", "TCEAL7", "BMT2", "TPRG1", "PARS2", "VARS2",
               "DNAJC5G", "U2AF1L4", "CCDC120", "OCLN", "ZSCAN10", "NAPEPLD", "SNX31", "XAGE5", 
               "TMEM69", "TMEM241", "MRPL37", "DAPL1", "TSPOAP1", "ZSCAN20", "L3MBTL2", "MAGEE2", 
               "MOSMO", "YBEY", "MPEG1", "CCDC153", "XAGE3", "ZSWIM4", "MITD1", "MRPL53", "CDCA7L", "ATP23",
               "KRT80", "MRPL55", "ZNF28","RPF2","ATL2","NPW","SRRM4","GLOD5","FADS6","CTTNBP2NL", 
               "ZNF181", "SHISA4", "TEX37","PCDHGA6","EFCAB9","SLC8A3","DCAF8L1","HDX","CWF19L2", 
               "MRPL54", "ANKRD13B", "LCTL", "CAPZA3", "C16orf87", "ISL2", "LMOD3", "SL",
               "C22A23", "CNTNAP5",
               "ZDHHC16","ANKRD31","TLR9", "HCAR3", "MRGPRG", "METTL24", "NSG2", "GASK1B", "UFD1",
               "SANBR", "IGLL5", "MIEN1", "OLFM2", "DIPK1A", "MRPL50", "GBGT1", "ATP5MC2","TMEM54",
               "SNX32","GSDME","BTF3L4","H2BC5","SRGAP1")

# Install and load the necessary packages
library(org.Hs.eg.db)
# Convert HUGO gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, keytype = "SYMBOL", column = "ENTREZID")

# Print the mapping results
#print(entrez_ids)
ans.kegg <- enrichKEGG(gene = entrez_ids,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)

tab.kegg <- as.data.frame(ans.kegg)

#extract results.
class(ans.kegg)
head(ans.kegg@result)
head(ans.kegg@result$Description)
class(ans.kegg@result$GeneRatio)
#format results
enriched<-ans.kegg@result%>%
  separate(BgRatio, into=c("size.term", "size.category"), sep="/")%>%
  separate(GeneRatio, into=c("size.overlap.term","size.overlap.category"), sep="/")%>%
  mutate_at(c("size.overlap.term","size.overlap.category","size.term", "size.category"),as.numeric)%>%
  mutate(k.K = size.overlap.term/size.term)

#visualise  the results.
enriched2<-enriched
top_20_enriched <- head(enriched2, 20)

top_20_enriched %>%
  ggplot(aes(x = forcats::fct_reorder(Description, p.adjust), y = p.adjust)) +
  geom_bar(stat = "identity", fill = "#f68060", alpha = 0.6, width = 0.4) +
  coord_flip() +
  xlab("") +
  ylab("-log10(pvalue)") +
  theme_bw()




