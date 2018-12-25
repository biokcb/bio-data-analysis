#####################################################################################################
#                                     Quick Script to analyze Diff Exp                              #




####################################################################################################

#setwd to wherever .CEL Files currently located
setwd("~/Desktop/Lab Notes/HCQ Gene Expression/R Stuff/CEL Files")

library(affy)

#create affy object
Data <- ReadAffy()

#process by RMA
Data_RMA <- rma(Data)

#check out phenodata
pData(Data_RMA)

#assign a treatment to each sample
pData(Data_RMA) <- cbind(pData(Data_RMa), treatment=factor(c("vehicle", "HCQ", "vehicle", "HCQ", "HCQ", "vehicle", "HCQ", "vehicle", "HCQ", "vehicle", "HCQ", "vehicle")))

#make sure it worked
pData(Data_RMA)

#1. By Limma
library(limma)

Info <- pData(Data_RMA)

#creates design matrix for linear model
design <- model.matrix(~ -1 + factor(Info$treatment,levels=unique(Info$treatment)))
colnames(design) <- unique(Info$treatment)
design

#set up contrast matrix
cont.matrix <- makeContrasts(HvV=HCQ-vehicle, levels=design)
cont.matrix

#run limma
fit <- lmFit(Data_RMA, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#Volcano Plot - visualize the differentially expressed genes
plot(fit2$coef,-log10(fit2$p.value),xlab="logFC",ylab="-log10(p)")
title("Volcano Plot")
abline(h=3,col="red");abline(v=+1,col="red");abline(v=-1,col="red")

#Make the nice table
AllGenes<-topTable(fit2,number=43035,genelist=NULL,sort.by="none",adjust.method="BH",p.value=1)
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)
Annotation<-read.csv("Canine_2.na34.annot.csv",skip=19,header=TRUE) # or appropriate file
SmAnnot<-Annotation[,c(1,11,14,15)] # or appropriate info
AllGenes<-merge(AllGenes,SmAnnot,by.x="Probe.Set.ID",by.y="Probe.Set.ID")
AllGenes.sorted <- AllGenes[order(AllGenes$P.Value),]

#2. By SAM
library(siggenes)

#check out phenodata 
pData(Data_RMA)

#If I hadn't already performed this step, would assign samples a treatment, but since we did that already, move on to next step

#Give treatment names a numerical value
groups <- as.numeric(pData(Data_RMA)$treament == "vehicle")

#Check to make sure it worked
groups

#Run SAM
samout <- sam(Data_RMA, groups)
samout

#Decide Delta 
plot(samout, desired delta)

#working on making a nice table of sig genes pulled out with annontation. 
