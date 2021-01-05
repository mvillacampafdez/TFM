##############
### RNAseq ###
##############
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(EDASeq)
library(edgeR)
query <- GDCquery(project = "TCGA-KIRP",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)
GDCdownload(query, method = "api", files.per.chunk = 10)
kirp.exp <- GDCprepare(query)

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = "KIRP")

# get clinical data
dataClin <- GDCquery_clinic(project = "TCGA-KIRP","clinical") 

# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"NT")


## txt matrix

EXP <- matrix(kirp.exp@rowRanges$gene_id, ncol=1, nrow=length(kirp.exp@rowRanges$gene_id))
EXP <- cbind(EXP,kirp.exp@assays@data$normalized_count)
dim(EXP)
colnames(EXP) <- c('gene_name',kirp.exp$sample)
write.table(EXP, file= 'EXPR.txt', sep='\t',row.names = F,quote = F)


########################
### Immunophenoscore ###
########################

# format input expression file 
x <- read.table("EXPR.txt", header=TRUE, sep="\t", dec = ".",check.names=FALSE)
x1 <- x[!duplicated(x$gene_name),]
x1$gene_name <- gsub("-", "_", x1$gene_name)
table(is.na(x1))

# guardo per tenir-la per altres coses
write.table(x1, "EXPR_mod.txt", col.names = T, row.names = F, sep = "\t", quote= F)


### IPS
gene_expression <- read.table("EXPR_mod.txt",row.names=1, header=TRUE, sep="\t", dec = ".",check.names=FALSE) 
sample_names<-names(gene_expression) 
ipsg <- read.table("IPS_genes2.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE) #els HLAs han d'estar tambe separats per _ (no -)
#ipsg <- ipsg[ipsg$GENE%in%rownames(gene_expression),]
unique_ips_genes <- as.vector(unique(ipsg$GENE))
geneEXP <- gene_expression[ipsg$GENE%in%unique_ips_genes,]

#Iincialitzem les variables
for (i in 1){  
  B2M<-NULL
  TAP1<-NULL
  TAP2<-NULL
  HLA_A<-NULL
  HLA_B<-NULL
  HLA_C<-NULL
  HLA_DPA1<-NULL
  HLA_DPB1<-NULL
  HLA_E<-NULL
  HLA_F<-NULL
  
  MHC<-NULL
  
  PDCD1<-NULL
  CTLA4<-NULL
  LAG3<-NULL
  TIGIT<-NULL
  HAVCR2<-NULL
  CD274<-NULL
  PDCD1LG2<-NULL
  CD27<-NULL
  ICOS<-NULL
  IDO1<-NULL
  
  CP<-NULL
  
  AIM2<-NULL
  BIRC3<-NULL
  BRIP1<-NULL
  CCL20<-NULL
  CCL4<-NULL
  CCL5<-NULL
  CCNB1<-NULL
  CCR7<-NULL
  DUSP2<-NULL
  ESCO2<-NULL
  ETS1<-NULL
  EXO1<-NULL
  EXOC6<-NULL
  IARS<-NULL
  KIF11<-NULL
  KNTC1<-NULL
  NUF2<-NULL
  PRC1<-NULL
  PSAT1<-NULL
  RGS1<-NULL
  RTKN2<-NULL
  SAMSN1<-NULL
  SELL<-NULL
  TRAT1<-NULL
  
  Act_CD4<-NULL
  
  ADRM1<-NULL
  AHSA1<-NULL
  C1GALT1C1<-NULL
  CCT6B<-NULL
  CD37<-NULL
  CD3D<-NULL
  CD3E<-NULL
  CD3G<-NULL
  CD69<-NULL
  CD8A<-NULL
  CETN3<-NULL
  CSE1L<-NULL
  GEMIN6<-NULL
  GNLY<-NULL
  GPT2<-NULL
  GZMA<-NULL
  GZMH<-NULL
  GZMK<-NULL
  IL2RB<-NULL
  LCK<-NULL
  MPZL1<-NULL
  NKG7<-NULL
  PIK3IP1<-NULL
  PTRH2<-NULL
  TIMM13<-NULL
  ZAP70<-NULL
  
  Act_CD8<-NULL
  
  ATM<-NULL
  CASP3<-NULL
  CASQ1<-NULL
  CD300E<-NULL
  DARS<-NULL
  DOCK9<-NULL
  EXOSC9<-NULL
  EZH2<-NULL
  GDE1<-NULL
  IL34<-NULL
  NCOA4<-NULL
  NEFL<-NULL
  PDGFRL<-NULL
  PTGS1<-NULL
  REPS1<-NULL
  SCG2<-NULL
  SDPR<-NULL
  SIGLEC14<-NULL
  SIGLEC6<-NULL
  TAL1<-NULL
  TFEC<-NULL
  TIPIN<-NULL
  TPK1<-NULL
  UQCRB<-NULL
  USP9Y<-NULL
  WIPF1<-NULL
  ZCRB1<-NULL
  
  TEM_CD4<-NULL
  
  ACAP1<-NULL
  APOL3<-NULL
  ARHGAP10<-NULL
  ATP10D<-NULL
  C3AR1<-NULL
  CCR5<-NULL
  CD160<-NULL
  CD55<-NULL
  CFLAR<-NULL
  CMKLR1<-NULL
  DAPP1<-NULL
  FCRL6<-NULL
  FLT3LG<-NULL
  GZMM<-NULL
  HAPLN3<-NULL
  HLA_DMB<-NULL
  HLA_DPA1<-NULL
  HLA_DPB1<-NULL
  IFI16<-NULL
  LIME1<-NULL
  LTK<-NULL
  NFKBIA<-NULL
  SETD7<-NULL
  SIK1<-NULL
  TRIB2<-NULL
  
  TEM_CD8<-NULL #Valor per + Tem CD8 (del grup EC)
  EC<-NULL
  
  CCR2<-NULL
  CD14<-NULL
  CD2<-NULL
  CD86<-NULL
  CXCR4<-NULL
  FCGR2A<-NULL
  FCGR2B<-NULL
  FCGR3A<-NULL
  FERMT3<-NULL
  GPSM3<-NULL
  IL18BP<-NULL
  IL4R<-NULL
  ITGAL<-NULL
  ITGAM<-NULL
  PARVG<-NULL
  PSAP<-NULL
  PTGER2<-NULL
  PTGES2<-NULL
  S100A8<-NULL
  S100A9<-NULL
  
  MDSC<-NULL #Valor per - MDSC (del grup SC)
  
  CCL3L1<-NULL
  CD72<-NULL
  CLEC5A<-NULL
  FOXP3<-NULL
  ITGA4<-NULL
  L1CAM<-NULL
  LIPA<-NULL
  LRP1<-NULL
  LRRC42<-NULL
  MARCO<-NULL
  MMP12<-NULL
  MNDA<-NULL
  MRC1<-NULL
  MS4A6A<-NULL
  PELO<-NULL
  PLEK<-NULL
  PRSS23<-NULL
  PTGIR<-NULL
  ST8SIA4<-NULL
  STAB1<-NULL
  
  Treg<-NULL #Valor per - Treg (del grup SC)
  
  SC <- NULL
  
  AZ<-NULL
  IPS<-NULL
}  
#Recorrem la taula per a calcular les variables
for (i in 1:length(sample_names)) {
  GE<-gene_expression[[i]]
  #mGE<-mean(GE, trim=0.025)
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(ipsg$GENE),i]-mGE)/sGE
  W1<-ipsg$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    #MIG[k]<- mean(Z1[which (as.vector(ipsg$NAME)==gen)]) #Per fer els immunophenograms per class
    #WEIGHT[k]<- mean(W1[which (as.vector(ipsg$NAME)==gen)])
    #k<-k+1
    MIG[k]<- mean(Z1[which (as.vector(ipsg$GENE)==gen)]) #Per fer els immunophenograms per IPS_gene
    WEIGHT[k]<- mean(W1[which (as.vector(ipsg$GENE)==gen)])
    k<-k+1
  }
  
  WG<-MIG*WEIGHT
  
  B2M[i]<-WG[1] #els que no tenen mean son single genes
  TAP1[i]<-WG[2]
  TAP2[i]<-WG[3]
  HLA_A[i]<-WG[4]
  HLA_B[i]<-WG[5]
  HLA_C[i]<-WG[6]
  HLA_DPA1[i]<-WG[7]
  HLA_DPB1[i]<-WG[8]
  HLA_E[i]<-WG[9]
  HLA_F[i]<-WG[10] 
  
  MHC[i]<-mean(WG[1:10]) #APG
  
  PDCD1[i]<-WG[11]
  CTLA4[i]<-WG[12]
  LAG3[i]<-WG[13]
  TIGIT[i]<-WG[14]
  HAVCR2[i]<-WG[15]
  CD274[i]<-WG[16]
  PDCD1LG2[i]<-WG[17]
  CD27[i]<-WG[18]
  ICOS[i]<-WG[19]
  IDO1[i]<-WG[20]
  
  CP[i]<-mean(WG[11:20])
  
  AIM2[i]<-WG[21]
  BIRC3[i]<-WG[22]
  BRIP1[i]<-WG[23]
  CCL20[i]<-WG[24]
  CCL4[i]<-WG[25]
  CCL5[i]<-WG[26]
  CCNB1[i]<-WG[27]
  CCR7[i]<-WG[28]
  DUSP2[i]<-WG[29]
  ESCO2[i]<-WG[30]
  ETS1[i]<-WG[31]
  EXO1[i]<-WG[32]
  EXOC6[i]<-WG[33]
  IARS[i]<-WG[34]
  KIF11[i]<-WG[35]
  KNTC1[i]<-WG[36]
  NUF2[i]<-WG[37]
  PRC1[i]<-WG[38]
  PSAT1[i]<-WG[39]
  RGS1[i]<-WG[40]
  RTKN2[i]<-WG[41]
  SAMSN1[i]<-WG[42]
  SELL[i]<-WG[43]
  TRAT1[i]<-WG[44]
  
  Act_CD4[i]<-mean(WG[21:44]) #Valor per + Act CD4 (del grup EC)
  
  ADRM1[i]<-WG[45]
  AHSA1[i]<-WG[46]
  C1GALT1C1[i]<-WG[47]
  CCT6B[i]<-WG[48]
  CD37[i]<-WG[49]
  CD3D[i]<-WG[50]
  CD3E[i]<-WG[51]
  CD3G[i]<-WG[52]
  CD69[i]<-WG[53]
  CD8A[i]<-WG[54]
  CETN3[i]<-WG[55]
  CSE1L[i]<-WG[56]
  GEMIN6[i]<-WG[57]
  GNLY[i]<-WG[58]
  GPT2[i]<-WG[59]
  GZMA[i]<-WG[60]
  GZMH[i]<-WG[61]
  GZMK[i]<-WG[62]
  IL2RB[i]<-WG[63]
  LCK[i]<-WG[64]
  MPZL1[i]<-WG[65]
  NKG7[i]<-WG[66]
  PIK3IP1[i]<-WG[67]
  PTRH2[i]<-WG[68]
  TIMM13[i]<-WG[69]
  ZAP70[i]<-WG[70]
  
  Act_CD8[i]<-mean(WG[45:70]) #Valor per + Act CD8 (del grup EC)
  
  ATM[i]<-WG[71]
  CASP3[i]<-WG[72]
  CASQ1[i]<-WG[73]
  CD300E[i]<-WG[74]
  DARS[i]<-WG[75]
  DOCK9[i]<-WG[76]
  EXOSC9[i]<-WG[77]
  EZH2[i]<-WG[78]
  GDE1[i]<-WG[79]
  IL34[i]<-WG[80]
  NCOA4[i]<-WG[81]
  NEFL[i]<-WG[82]
  PDGFRL[i]<-WG[83]
  PTGS1[i]<-WG[84]
  REPS1[i]<-WG[85]
  SCG2[i]<-WG[86]
  SDPR[i]<-WG[87]
  SIGLEC14[i]<-WG[88]
  SIGLEC6[i]<-WG[89]
  TAL1[i]<-WG[90]
  TFEC[i]<-WG[91]
  TIPIN[i]<-WG[92]
  TPK1[i]<-WG[93]
  UQCRB[i]<-WG[94]
  USP9Y[i]<-WG[95]
  WIPF1[i]<-WG[96]
  ZCRB1[i]<-WG[97]
  
  TEM_CD4[i]<-mean(WG[c(71:81,83:97)]) #Valor per + Tem CD4 (del grup EC)
  
  ACAP1[i]<-WG[98]
  APOL3[i]<-WG[99]
  ARHGAP10[i]<-WG[100]
  ATP10D[i]<-WG[101]
  C3AR1[i]<-WG[102]
  CCR5[i]<-WG[103]
  CD160[i]<-WG[104]
  CD55[i]<-WG[105]
  CFLAR[i]<-WG[106]
  CMKLR1[i]<-WG[107]
  DAPP1[i]<-WG[108]
  FCRL6[i]<-WG[109]
  FLT3LG[i]<-WG[110]
  GZMM[i]<-WG[111]
  HAPLN3[i]<-WG[112]
  HLA_DMB[i]<-WG[113]
  HLA_DPA1[i]<-WG[7]
  HLA_DPB1[i]<-WG[8]
  IFI16[i]<-WG[114]
  LIME1[i]<-WG[115]
  LTK[i]<-WG[116]
  NFKBIA[i]<-WG[117]
  SETD7[i]<-WG[118]
  SIK1[i]<-WG[119]
  TRIB2[i]<-WG[120]
  
  TEM_CD8[i]<-mean(WG[c(7,8,98:114,116:120)]) #Valor per + Tem CD8 (del grup EC)
  EC[i]<-mean(Act_CD4[i],Act_CD8[i],TEM_CD4[i],TEM_CD8[i])
  
  CCR2[i]<-WG[121]
  CD14[i]<-WG[122]
  CD2[i]<-WG[123]
  CD86[i]<-WG[124]
  CXCR4[i]<-WG[125]
  FCGR2A[i]<-WG[126]
  FCGR2B[i]<-WG[127]
  FCGR3A[i]<-WG[128]
  FERMT3[i]<-WG[129]
  GPSM3[i]<-WG[130]
  IL18BP[i]<-WG[131]
  IL4R[i]<-WG[132]
  ITGAL[i]<-WG[133]
  ITGAM[i]<-WG[134]
  PARVG[i]<-WG[135]
  PSAP[i]<-WG[136]
  PTGER2[i]<-WG[137]
  PTGES2[i]<-WG[138]
  S100A8[i]<-WG[139]
  S100A9[i]<-WG[140]
  
  MDSC[i]<-mean(WG[121:140]) #Valor per - MDSC (del grup SC)
  
  CCL3L1[i]<-WG[141]
  CD72[i]<-WG[142]
  CLEC5A[i]<-WG[143]
  FOXP3[i]<-WG[144]
  ITGA4[i]<-WG[145]
  L1CAM[i]<-WG[146]
  LIPA[i]<-WG[147]
  LRP1[i]<-WG[148]
  LRRC42[i]<-WG[149]
  MARCO[i]<-WG[150]
  MMP12[i]<-WG[151]
  MNDA[i]<-WG[152]
  MRC1[i]<-WG[153]
  MS4A6A[i]<-WG[154]
  PELO[i]<-WG[155]
  PLEK[i]<-WG[156]
  PRSS23[i]<-WG[157]
  PTGIR[i]<-WG[158]
  ST8SIA4[i]<-WG[159]
  STAB1[i]<-WG[160]
  
  Treg[i]<-mean(WG[142:160]) #Valor per - Treg (del grup SC)
  SC[i] <- mean(MDSC[i],Treg[i])
  
  #AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  #IPS[i]<-ipsmap(AZ[i])
  
  
}  


#Generem la nova taula

DF <- data.frame(SAMPLE=sample_names, B2M=B2M, TAP1=TAP1, TAP2=TAP2, HLA_A=HLA_A, HLA_B=HLA_B, HLA_C=HLA_C, HLA_DPA1=HLA_DPA1, HLA_DPB1=HLA_DPB1, HLA_E=HLA_E, HLA_F=HLA_F,
                 PDCD1=PDCD1, CTLA4=CTLA4, LAG3=LAG3, TIGIT=TIGIT, HAVCR2=HAVCR2, CD274=CD274, PDCD1LG2=PDCD1LG2, CD27=CD27, ICOS=ICOS, IDO1=IDO1,
                 AIM2=AIM2, BIRC3=BIRC3, BRIP1=BRIP1, CCL20=CCL20, CCL4=CCL4, CCL5=CCL5, CCNB1=CCNB1, CCR7=CCR7, DUSP2=DUSP2, ESCO2=ESCO2, ETS1=ETS1, EXO1=EXO1, EXOC6=EXOC6, 
                 IARS=IARS, KIF11=KIF11, KNTC1=KNTC1, NUF2=NUF2, PRC1=PRC1, PSAT1=PSAT1, RGS1=RGS1, RTKN2=RTKN2, SAMSN1=SAMSN1, SELL=SELL, TRAT1=TRAT1, Act_CD4=Act_CD4,
                 ADRM1=ADRM1, AHSA1=AHSA1, C1GALT1C1=C1GALT1C1, CCT6B=CCT6B, CD37=CD37, CD3D=CD3D, CD3E=CD3E, CD3G=CD3G, CD69=CD69, CD8A=CD8A, CETN3=CETN3, CSE1L=CSE1L, 
                 GEMIN6=GEMIN6, GNLY=GNLY, GPT2=GPT2, GZMA=GZMA, GZMH=GZMH, GZMK=GZMK, IL2RB=IL2RB, LCK=LCK, MPZL1=MPZL1, NKG7=NKG7, PIK3IP1=PIK3IP1, PTRH2=PTRH2, TIMM13=TIMM13,
                 ZAP70=ZAP70, Act_CD8=Act_CD8, ATM=ATM, CASP3=CASP3, CASQ1=CASQ1, CD300E=CD300E, DARS=DARS, DOCK9=DOCK9, EXOSC9=EXOSC9, EZH2=EZH2, GDE1=GDE1, IL34=IL34, NCOA4=NCOA4, NEFL=NEFL, PDGFRL=PDGFRL, 
                 PTGS1=PTGS1, REPS1=REPS1, SCG2=SCG2, SDPR=SDPR, SIGLEC14=SIGLEC14, SIGLEC6=SIGLEC6, TAL1=TAL1, TFEC=TFEC, TIPIN=TIPIN, TPK1=TPK1, UQCRB=UQCRB, USP9Y=USP9Y, WIPF1=WIPF1, 
                 ZCRB1=ZCRB1, TEM_CD4=TEM_CD4, ACAP1=ACAP1, APOL3=APOL3, ARHGAP10=ARHGAP10, ATP10D=ATP10D, C3AR1=C3AR1, CCR5=CCR5, CD160=CD160, CD55=CD55, CFLAR=CFLAR, CMKLR1=CMKLR1, 
                 DAPP1=DAPP1, FCRL6=FCRL6,FLT3LG=FLT3LG, GZMM=GZMM, HAPLN3=HAPLN3, HLA_DMB=HLA_DMB, HLA_DPA1=HLA_DPA1, HLA_DPB1=HLA_DPB1, IFI16=IFI16, LIME1=LIME1, LTK=LTK, NFKBIA=NFKBIA, 
                 SETD7=SETD7, SIK1=SIK1, TRIB2=TRIB2, TEM_CD8=TEM_CD8,CCR2=CCR2, CD14=CD14, CD2=CD2, CD86=CD86, CXCR4=CXCR4, FCGR2A=FCGR2A, FCGR2B=FCGR2B, FCGR3A=FCGR3A, FERMT3=FERMT3, 
                 GPSM3=GPSM3, IL18BP=IL18BP, IL4R=IL4R, ITGAL=ITGAL, ITGAM=ITGAM, PARVG=PARVG, PSAP=PSAP, PTGER2=PTGER2, PTGES2=PTGES2, S100A8=S100A8, S100A9=S100A9, MDSC=MDSC, 
                 CCL3L1=CCL3L1, CD72=CD72, CLEC5A=CLEC5A, FOXP3=FOXP3, ITGA4=ITGA4, L1CAM=L1CAM, LIPA=LIPA, LRP1=LRP1, LRRC42=LRRC42, MARCO=MARCO, MMP12=MMP12, MNDA=MNDA, MRC1=MRC1, 
                 MS4A6A=MS4A6A, PELO=PELO, PLEK=PLEK, PRSS23=PRSS23, PTGIR=PTGIR, ST8SIA4=ST8SIA4, STAB1=STAB1, Treg=Treg, MHC=MHC,EC=EC,SC=SC,CP=CP
)

DF_IPS <-  data.frame(SAMPLE=sample_names, MHC=MHC, CP=CP, Act_CD4=Act_CD4, Act_CD8=Act_CD8, TEM_CD4=TEM_CD4, TEM_CD8=TEM_CD8, EC=EC, MDSC=MDSC, Treg=Treg, SC=SC) 
write.table(DF, "IPS_results.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(DF_IPS, "IPS_results_abrev.txt", row.names = F, col.names = T, sep = "\t", quote = F)


###################
### MCP Counter ###
###################

install.packages(c("devtools", "curl"))
library(devtools)
library(curl)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)

## MODIFIED FUNCTION:
MCPcounter.estimate=function(
  expression,
  featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID","ENSEMBL_ID")[1],
  probesets=read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
  genes=read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
){
  ## marker.names=c("T cells","CD8 T cells","Cytotoxic lymphocytes","NK cells","B lineage","Monocytic lineage","Myeloid dendritic cells","Neutrophils","Endothelial cells","Fibroblasts")
  
  
  if(featuresType=="affy133P2_probesets"){
    features=probesets
    markers.names = unique(features[, 2])
    features=split(features[,1],features[,2])
    features=lapply(features,intersect,x=rownames(expression))
    features=features[sapply(features,function(x)length(x)>0)]
    missing.populations=setdiff(markers.names,names(features))
    features=features[intersect(markers.names,names(features))]
    
  } else {
    markersG=genes
  }
  
  if(featuresType=="HUGO_symbols"){
    features=subset(markersG,get("HUGO symbols")%in%(expression[,1]))
    markers.names = unique(features[, "Cell population"])
    features=split(features[,"HUGO symbols"],features[,"Cell population"])
    missing.populations=setdiff(markers.names,names(features))
    features=features[intersect(markers.names,names(features))]
    
  }
  
  
  
  if(length(missing.populations)>0){
    warning(paste("Found no markers for population(s):",paste(missing.populations,collapse=", ")))
  }
  appendSignatures <- function (xp, markers) 
  {
    names <- xp[,1]
    xp <- xp[,-1]
    res = as.data.frame(do.call(cbind, lapply(markers, function(x) {
      apply(xp[expression[,1]%in%intersect(names, x),, drop = F], 2, 
            mean, na.rm = T)
    })))
    res
  }
  t(appendSignatures(expression,features))
}


expression <- read.delim("EXPR.txt", header = T, sep = "\t")  


#use TCGA Level 3 *.rsem.genes.normalized_results
#first column with "Hugo_symbol" or "gene_name", and samples in the subsequent columns

estimates <- MCPcounter.estimate(expression, featuresType = "HUGO_symbols",
                                 probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                 genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))

write.table(estimates, "mcp_estimates_tcga.txt", row.names = T, col.names = T, sep = "\t", quote = F)
estimates <- read.table('mcp_estimates_tcga.txt', sep = "\t")
#para obtener un plot de "unsupervised clustering" teniendo en cuenta los resultados:
heatmap(as.matrix(estimates),col=colorRampPalette(c("blue","white","red"))(100))
