rm(list=ls())
library("edgeR")
library("baySeq")
library(RColorBrewer)
library(gplots)
library(gdata)
library(dplyr)
library(tidyr)
library(matrixStats)
# install.packages("rstudioapi") # run this if it's your first time using it to install
library(rstudioapi) # load it
# the following line is for getting the path of your current open file
current_path <- getActiveDocumentContext()$path 
# The next line set the working directory to the relevant one:
setwd(dirname(current_path ))
# you can make sure you are in the right directory
print( getwd() )

source("getDEgenes.R")

keep=c(1,2,3,4,5,6) # keep all the samples

#log fold change threshold

#treatment groups
group <- c('plusFe','plusFe','plusFe','minusFe','minusFe','minusFe')
#logfold change threshold
LFC = 0.75
#read the files 
cm <- read.delim("../Data/GeneCountTable/Counts_B.txt",row.names = "GeneID")
B <- getDEgenes(cm,group=group,batch_factor=c(1,1,2,1,1,2),sample="B",lfc=LFC)
B$t6 <-rep(1,nrow(B))
names(B) <- c("GeneId","qval6hr","t6")


#read the files for C
cm <- read.delim("../Data/GeneCountTable/Counts_C.txt",row.names = "GeneID")
cm <-cm[,keep]
C <- getDEgenes(cm,group=group,batch_factor= c(1,1,2,1,1,2),sample="C",lfc=LFC)
C$t12 <-rep(1,nrow(C))
names(C) <- c("GeneId","qval12hr","t12")


#read file for D 
cm <- read.delim("../Data/GeneCountTable/Counts_D.txt",row.names = "GeneID")
cm <-cm[,keep]
D <- getDEgenes(cm,group=group,batch_factor=batch_factor <- c(1,1,2,1,1,2),sample="D",lfc=LFC)
D$t18 <-rep(1,nrow(D))
names(D) <- c("GeneId","qval18hr","t18")

#read file for E 
cm <- read.delim("../Data/GeneCountTable/Counts_E.txt",row.names = "GeneID")
E <- getDEgenes(cm,group=group,batch_factor=batch_factor <- c(1,1,2,1,1,2),sample="E",lfc=LFC)
E$t24 <-rep(1,nrow(E))
names(E) <- c("GeneId","qval24hr","t24")



#read file for F 
cm <- read.delim("../Data/GeneCountTable/Counts_F.txt",row.names = "GeneID")
F <- getDEgenes(cm,group=group,batch_factor=batch_factor <- c(1,1,2,1,1,2),sample="F",lfc=LFC)
F$t30 <-rep(1,nrow(F))
names(F) <- c("GeneId","qval30hr","t30")


#read file for G
group <- c('plusFe','plusFe','plusFe','minusFe','minusFe','minusFe')
cm <- read.delim("../Data/GeneCountTable/Counts_G.txt",row.names = "GeneID")
G <- getDEgenes(cm,group=group,batch_factor=batch_factor <- c(1,2,3,1,1,3),sample="G",lfc=LFC)
G$t36 <-rep(1,nrow(G))
names(G) <- c("GeneId","qval36hr","t36")
  
  
#merge all DE genes
DEGTable <-full_join(B,C)
DEGTable <-full_join(DEGTable,D)
DEGTable <-full_join(DEGTable,E)
DEGTable <-full_join(DEGTable,F)
DEGTable <-full_join(DEGTable,G)

DEGTable[is.na(DEGTable$t6),"t6"] <-0
DEGTable[is.na(DEGTable$t12),"t12"] <-0
DEGTable[is.na(DEGTable$t18),"t18"] <-0
DEGTable[is.na(DEGTable$t24),"t24"] <-0
DEGTable[is.na(DEGTable$t30),"t30"] <-0
DEGTable[is.na(DEGTable$t36),"t36"] <-0

#calculate number of timepoints a gene was found to be DE
DEGTable <-DEGTable %>%   mutate(totalDEtimpoint = t6 + t12+ t18+ t24+ t30+ t36)


DEGTable_unique <- DEGTable %>% dplyr::select (-c(t6, t12, t18, t24, t30, t36))
DEGTable_unique[is.na(DEGTable_unique)] <-""
#add gene names
alias <- read.delim(file="../Data/gene_aliases.txt", sep="\t")
alias <-alias[,c(1,2,3)]
names(alias) <- c("GeneId","Alias","description")
DEGTable <-left_join(DEGTable,alias)
DEGTable$Alias=as.character(DEGTable$Alias)
DEGTable$description=as.character(DEGTable$description)
DEGTable[is.na(DEGTable$Alias),"Alias"] <- "unknown"
DEGTable[is.na(DEGTable$description),"description"] <- "unknown"
tmp <- DEGTable %>% 
  group_by(GeneId) %>% 
  summarise(Alias = paste(Alias, collapse = ", "),description = paste(description, collapse = ", ")) 

DEGTable_unique <-left_join(DEGTable_unique,tmp)
write.csv(DEGTable_unique, file="../Data/GeneCountTable/DEGenes_glmQtest.csv", row.names = FALSE)




#Find genes of intereset or known or control genes in DEG_table
DEGTable_unique<-read.csv(file="../Data/GeneCountTable/DEGenes_glmQtest.csv")
control = read.xls('../Data/Controls.xlsx',stringsAsFactors=FALSE)
control <- data.frame(toupper(trimws(control[,1])))
colnames(control)<-'GeneId'
inner_join(control,DEGTable_unique)
counter =0

controlgenes=control$GeneId
counter=0
i=0;
isControl=rep('unknown',nrow(DEGTable_unique))

for(gene in controlgenes)
{
  #ind=which(grepl(toupper(gene),toupper(DEGTable_unique$GeneId))==TRUE)
  gene=trimws(gene) #remove whitespace before and after the AGIs
  
  ind=grep(toupper(gene),toupper(DEGTable_unique$GeneId))
  if(length(ind)>0){
    counter=counter+1
    print(gene)
    #print(DEGTable_unique$Alias[ind])
    isControl[grep(toupper(gene),DEGTable_unique$GeneId)]<-'known'
  }
  i=i+1
}
print(paste("Number of known genes found in DE analysis:", counter))





##read RPKM_DE_genes and meanRPKM DE Genes

rpkm <- read.csv("../Data/GeneCountTable/RPKM_allgenes.csv",stringsAsFactors = FALSE)
de <- read.csv("../Data/GeneCountTable/DEGenes_glmQtest.csv",stringsAsFactors = FALSE)


de<-data.frame(de[,"GeneId"])
colnames(de)<- "GeneId"

derpkm <- left_join(de,rpkm)

write.csv(derpkm,"../Data/GeneCountTable/RPKM_DEgenes.csv",row.names = FALSE)

meanRPKM <- read.csv("../Data/GeneCountTable/meanRPKM_allgenes.csv",stringsAsFactors = FALSE)

demeanrpkm <- left_join(de,meanRPKM)

#the difference signal between -Fe and +Fe
difference<-demeanrpkm[,c(2,9:14)]- demeanrpkm[,2:8]
colnames(difference)<-c('diff0hr','diff6hr','diff12hr','diff18hr','diff24hr','diff30hr','diff36hr')

#this file should contain meanRPKM and std of RPKM values at each timepoint under plusFe and minusFe
write.csv(demeanrpkm,"../Data/GeneCountTable/meanRPKM_DEgenes.csv",row.names = FALSE)





feplus <- matrix(0,nrow = nrow(derpkm), ncol = 7)
feminus <- matrix(0,nrow = nrow(derpkm), ncol = 7)

#calculate normalized time points required for visualization and clustering the gene expression patterns.

for(i in 1:nrow(demeanrpkm)){
  tmp =(demeanrpkm[i,2:8] - rowMeans(demeanrpkm[i,2:8]))/rowSds(as.matrix(demeanrpkm[i,2:8]))
  i
  feplus[i,] = as.matrix(tmp)
  tmp =(demeanrpkm[i,c(2,9:14)] - rowMeans(demeanrpkm[i,c(2,9:14)]))/rowSds(as.matrix(demeanrpkm[i,c(2,9:14)]))
  feminus[i,] = as.matrix(tmp)
  
}

normexp <- cbind(feplus,feminus)
normexp <- data.frame(normexp)
colnames(normexp)<-c('Anorm','Bnorm','Cnorm','Dnorm','Enorm','Fnorm','Gnorm','Aminusnorm','Bminusnorm','Cminusnorm','Dminusnorm','Eminusnorm','Fminusnorm','Gminusnorm')
normexp$GeneId <- demeanrpkm$GeneId
normexp <- normexp[,c(15,1:14)]
write.csv(normexp,"../Data/GeneCountTable/normexp_DEGenes.csv",row.names = FALSE)


for(i in 1:nrow(demeanrpkm)){
  tmp =demeanrpkm[i,2:8]/max(demeanrpkm[i,2:8])
  feplus[i,] = as.matrix(tmp)
  tmp =demeanrpkm[i,c(2,9:14)]/max(demeanrpkm[i,c(2,9:14)])
  feminus[i,] = as.matrix(tmp)
}

maxnormexp <- cbind(feplus,feminus)
maxnormexp <- data.frame(maxnormexp)
colnames(maxnormexp)<-c('Amaxnorm','Bmaxnorm','Cmaxnorm','Dmaxnorm','Emaxnorm','Fmaxnorm','Gmaxnorm','Aminusmaxnorm','Bminusmaxnorm','Cminusmaxnorm','Dminusmaxnorm','Eminusmaxnorm','Fminusmaxnorm','Gminusmaxnorm')
maxnormexp$GeneId <- demeanrpkm$GeneId
maxnormexp <- maxnormexp[,c(15,1:14)]
write.csv(maxnormexp,"../Data/GeneCountTable/maxnormexp_DEGenes.csv",row.names = FALSE)



for(i in 1:nrow(demeanrpkm)){
  tmp =(demeanrpkm[i,2:8]-rowMeans(demeanrpkm[i,2:8]))/(max(demeanrpkm[i,2:8]))
  feplus[i,] = as.matrix(tmp)
  tmp =(demeanrpkm[i,c(2,9:14)]-rowMeans(demeanrpkm[i,c(2,9:14)]))/(max(demeanrpkm[i,c(2,9:14)]))
  feminus[i,] = as.matrix(tmp)
}

sumnormexp <- cbind(feplus,feminus)
sumnormexp <- data.frame(sumnormexp)
colnames(sumnormexp)<-c('Ameanmaxnorm','Bmeanmaxnorm','Cmeanmaxnorm','Dmeanmaxnorm','Emeanmaxnorm','Fmeanmaxnorm','Gmeanmaxnorm','Aminusmeanmaxnorm','Bminusmeanmaxnorm','Cminusmeanmaxnorm','Dminusmeanmaxnorm','Eminusmeanmaxnorm','Fminumeanmaxnorm','Gminusmeanmaxnorm')
sumnormexp$GeneId <- demeanrpkm$GeneId
sumnormexp <- sumnormexp[,c(15,1:14)]
write.csv(sumnormexp,"../Data/GeneCountTable/meanmaxnormexp_DEGenes.csv",row.names = FALSE)

########################################################################################################################


#deTime points
de <- read.csv("../Data/GeneCountTable/DEGenes_glmQtest.csv",stringsAsFactors = FALSE)
de<-de[,1:8]
colnames(de)<-c('GeneId','isDE6hr','isDE12hr','isDE18hr','isDE24hr','isDE30hr','isDE36hr','totalDETimepoints')
for(i in 2:7){
  de[!is.na(de[,i]),i] <- 1
  de[is.na(de[,i]),i] <- 0
  }

#replace NA values with unknown
demeanrpkm$Alias[which(is.na(demeanrpkm$Alias))]<-"unknown"
demeanrpkm$description[is.na(demeanrpkm$description)]='unknown'


#merging all information in a large table
fulltable <-cbind(demeanrpkm,normexp[,2:15],maxnormexp[,2:15],sumnormexp[,2:15],difference,de[,2:8])

#adding the clustering results to data table
#these files should be generated from hmm clustering results before running these lines
hmmclusterMinusFe <-read.csv("../Data/hmm_clustering_glmq_minus_Fe_20_filtered_may14.csv")
hmmclusterplusFe <-read.csv("../Data/hmm_clustering_glmq_plus_Fe_20_filtered_may14.csv")

cluster <-cbind(hmmclusterMinusFe$cluster,hmmclusterplusFe$cluster)
colnames(cluster)<-c('hmmclusterMinusFe','hmmclusterplusFe')
fulltable <-cbind(fulltable,cluster)


hmmclusterMinusFemaxnorm <-read.csv("../Data/hmm_clustering_glmq_minus_Fe_20_filtered_June21.csv")
hmmclusterplusFemaxnorm <-read.csv("../Data/hmm_clustering_glmq_plus_Fe_20_filtered_June21.csv")

cluster <-cbind(hmmclusterMinusFemaxnorm$cluster,hmmclusterplusFemaxnorm$cluster)
colnames(cluster)<-c('hmmclusterMinusFemaxnorm','hmmclusterplusFemaxnorm')

fulltable <-cbind(fulltable,cluster)

hmmclusterminusmeanmaxnorm <-read.csv("../Data/hmm_clustering_glmq_minus_Fe_20_meanmax.csv")
hmmclusterplusmeanmaxnorm <-read.csv("../Data/hmm_clustering_glmq_plus_Fe_20_meanmax.csv")

DPGP_FeMinus<-read.table("../Data/DPGP_clustering_mean_max_FeMinus.txt",sep="\t",header = TRUE)
names(DPGP_FeMinus)<-c('DPGPcluster','GeneId')
fulltable<-left_join(fulltable,DPGP_FeMinus)
cluster <-cbind(hmmclusterminusmeanmaxnorm$cluster,hmmclusterplusmeanmaxnorm$cluster)
colnames(cluster)<-c('hmmclusterminusmeanmaxnorm','hmmclusterplusmeanmaxnorm')

fulltable <-cbind(fulltable,cluster)

fulltable$isControl<-isControl



##FIND TFs and known genes  ### change this when we get a full list of TFs in A. thaliana


tf <- read.table("../Data/Ath_TF_list.txt",sep='\t',stringsAsFactors=FALSE,header = TRUE);
tf<-tf[,c(2,3)]
colnames(tf)<-c('GeneId','TF_Family')
tf<-inner_join(tf,de)
tf<-tf[,c(1,2)]



fulltable$isTF=rep('other',nrow(fulltable))
fulltable$TFFamily=rep(' ',nrow(fulltable))

for(i in 1:nrow(tf)) {
  ind <-grep(toupper(tf[i,"GeneId"]),fulltable$GeneId)
  print(ind)
  fulltable$isTF[ind] <- 'TF'
  fulltable$TFFamily[ind] <- tf[i,"TF_Family"]
  
}


#direction activated or repressed


direction<-rowSums(as.matrix(difference))
ind1 <-which(direction>0)
ind2 <-which(direction<0)
direction[ind1]<-'Activated'
direction[ind2]<-'Repressed'
fulltable$direction <-direction

write.csv(fulltable,'../Data/GeneCountTable/DE_Genes_AT_epidermis_featurecounts_June21.csv',row.names = FALSE)
#write.xlsx(fulltable, '../Data/GeneCountTable/DE_Genes_AT_epidermis_featurecounts.xlsx', sheetName="Sheet1",  col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
# gene ontology enrichment codes 

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}


require("biomaRt")
require("topGO")
mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
head (GTOGO)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))

all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
int.genes <- DEGTable_unique$GeneId # some random genes 
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
              , allGenes = int.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO
)

results <- runTest(go.obj, algorithm = "elim", statistic = "fisher")

results.tab <- GenTable(object = go.obj, elimFisher = results)


showSigOfNodes(go.obj, score(results), firstSigNodes = 5, useInfo = 'all')









