##################################################################
#     A sample script for normalization of microarray data       #
#                                                                #
#                  R/COXEN Summer Camp 2014                      #
#                                                                #
#                                                                #
##################################################################


##### ---   1. Setting up directories / files  --- ######

setwd("/Users/Kristen/Documents/GustafsonLab") #set working/home directory 

############ Batch correction for ACC29 files. 

cel.dir <- "./ACC29 files/"
canine<- list.files(path=cel.dir, pattern=".CEL") #to get all the .CEL files in a list, sometimes the pattern=".cel". Check your file types
canine #check the list
length(canine) #check the number of files you're processing

caffy <- ReadAffy(filenames=paste(cel.dir,canine,sep=''), phenoData=NULL) #creates an AffyBatch object based on files from above. Can have phenoData, I typically don't & just use a separate file for extra information

slotNames(caffy) # looking at what is in the AffyBatch
protocol <- caffy@protocolData # pulling out data from one of the "slots"
slotNames(protocol) # which also contains "slots" 
date.info <- protocol@data # this is where we pull out scan date information
date.info #from this we can manually decide how many batches there are and decide on whether or not to do batch correction

#batches! first, pull off dates
dates <- as.data.frame(strsplit(date.info$ScanDate[19],' ', fixed=T)) # separating the scan time and date
dates <- as.matrix(as.character(dates[1,])) #back to columns and only keeping dates
rownames(dates) <- rownames(date.info)[19] #making sure filenames remain attached to dates
dates 

dates1 <- as.data.frame(strsplit(date.info$ScanDate[c(1:18,20:29)],'T', fixed=T)) # separating the scan time and date
dates1 <- t(dates1[1,]) #back to columns and only keeping dates
rownames(dates1) <-rownames(date.info)[c(1:18,20:29)]
dates1

dates <- matrix(c(dates1[1:18,],dates,dates1[19:28,]), nrow=29, ncol=1)
rownames(dates) <- rownames(date.info)
dates
# IF all the dates are the same, move onto RMA...

#if there is only one batch, we cannot do batch correction without any more information. If there are any batches that only have one .CEL file, we can't use combat
nums <- aggregate(dates, by=list(dates), FUN=length) #count the numbers of files that correspond to each date
nums
nums <- nums[nums[,2]!=1,] # remove dates with only one file
dates <- nums[match(dates, nums[,1]),1] #match those dates in the dates object
names(dates) <- rownames(date.info) #keep track of file info....
dates <- dates[complete.cases(dates)] #only keep non-NA values
dates
length(dates) #look at the number of files left

# next, find out how many unique dates there were
batch <- unique(dates)
batch <- cbind(batch, 1:length(batch)) #assign each unique date a batch number
batch

# next, assign those batch numbers to the list of files
batches <- batch[match(dates,batch), 2]
batches <- as.matrix(as.numeric(batches))
rownames(batches) <- names(dates)
batches

# if you're doing batch correction, remove files that won't be included due to only having one file in a batch, before doing RMA
caffy <- caffy[,rownames(batches)] #subset by the filenames

slotNames(caffy)
pheno <- caffy@phenoData #pulling out sample names
samples <- pheno@data #to double check
identical(rownames(samples), rownames(batches)) #if false, then something went wrong when subsetting. Check again.
samples

# RMA 
caffyRMA <- rma(caffy)
caffy1 <- list(caffyRMA[,names(batches[(batches[,1] %in% 1),])],caffyRMA[,names(batches[(batches[,1] %in% 2),])])
cexprRMA <- exprs(caffyRMA)

# let's look at a pc plot
pca.can <- prcomp(t(cexprRMA)) #get principal components
summary(pca.can) #look at them

#plot pc1 vs 2
plot(pca.can$x[,1:2], type="n", main="PCA Summary for ACC29 - RMA only")
pc1 <- pca.can$x[,1]
pc2 <- pca.can$x[,2]
points(pc1[batches[,1] %in% 1],pc2[batches[,1] %in% 1], pch=20, col=2)
points(pc1[batches[,1] %in% 2],pc2[batches[,1] %in% 2], pch=20, col=3)

# Using ComBat for batch corrections
library("sva")

# in addition to the batch information, we may want to include additional confounding factors, such as tissue type. 
# I have a separate file that indicates the cell type for my dataset, so I read it in
canineInfo <- read.table("ACC29_PhenoData.txt", header=TRUE, sep="\t")
rownames(canineInfo) <- canineInfo$Filename
canineInfo

#match the filenames to the ones in your affy object
canineInfo <- canineInfo[rownames(samples),]
canineInfo
identical(rownames(samples), rownames(canineInfo))
tissue <- canineInfo$type
tissue <- as.matrix(as.factor(tissue))
rownames(tissue)<-rownames(canineInfo)
tissue
# now we can apply ComBat to either RMA
ccombatRMA <- ComBat(dat=cexprRMA, batch=batches, mod=tissue)

# let's look at a pc plot
pca.can2 <- prcomp(t(ccombatRMA)) #get principal components
summary(pca.can2) #look at them

#plot pc1 vs 2, 
plot(pca.can2$x[,1:2], type="n", main="PCA Summary for ACC29 - ComBat")
pc1 <- pca.can2$x[,1]
pc2 <- pca.can2$x[,2]
points(pc1[batches[,1] %in% 1],pc2[batches[,1] %in% 1], pch=20, col=2)
points(pc1[batches[,1] %in% 2],pc2[batches[,1] %in% 2], pch=20, col=3)

# alternatively we can use InSilicoMerge
library("inSilicoMerging")

cRMAmerge <- merge(caffy1, method="COMBAT")

# let's look at a pc plot
pca.can3 <- prcomp(t(exprs(cRMAmerge))) #get principal components
summary(pca.can3) #look at them

#plot pc1 vs 2
plot(pca.can3$x[,1:2], type="n", main="PCA Summary for ACC29 - Merge/ComBat")
pc1 <- pca.can3$x[,1]
pc2 <- pca.can3$x[,2]
points(pc1[batches[,1] %in% 1],pc2[batches[,1] %in% 1], pch=20, col=2)
points(pc1[batches[,1] %in% 2],pc2[batches[,1] %in% 2], pch=20, col=3)


############

#Pretty much the same with the NCI60. 
#Some objects we used before will be overwritten, but we won't need them again so it will be fine. 
#I chose to save the .rda files/load them in because there are 170 .CEL files to process

cel.dir <- ".NCI60 HGU133 Plus2/CEL-files/" # where the .CEL files are stored, the './' means the path starts from the working directory
save.dir <- "./NCI60/" # if you want to save files in a directory different from the home directory

#note: you can also just type these into functions directly, this is just how I do it so I can change directories and whatever else at the beginning of the code instead of in each function

filenames <- list.files(path=cel.dir, pattern=".cel") #to get all the .CEL files in a list, sometimes the pattern=".cel". Check your file types
filenames #check the list
length(filenames) #check the number of files you're processing


##### ---   2. Expression Set Creation/Processing  --- ######

library(affy) #where the affy related functions are stored. Need to load this library to continue
#affy <- ReadAffy(filenames=paste(cel.dir,filenames,sep=''), phenoData=NULL) #creates an AffyBatch object based on files from above. Can have phenoData, I typically don't & just use a separate file for extra information
#save(affy, file="NCI60.Affydata.rda")

load("NCI60.Affydata.rda")

slotNames(affy) # looking at what is in the AffyBatch
protocol <- affy@protocolData # pulling out data from one of the "slots"
slotNames(protocol) # which also contains "slots" 
date.info <- protocol@data # this is where we pull out scan date information
date.info #from this we can manually decide how many batches there are and decide on whether or not to do batch correction

#batches! first, pull off dates
dates <- as.data.frame(strsplit(date.info$ScanDate,' ', fixed=T)) # separating the scan time and date
dates <- t(dates[1,]) #back to columns and only keeping dates
rownames(dates) <- rownames(date.info) #making sure filenames remain attached to dates
dates 

# IF all the dates are the same, move onto RMA...

#if there is only one batch, we cannot do batch correction without any more information. If there are any batches that only have one .CEL file, we can't use combat
nums <- aggregate(dates, by=list(dates), FUN=length) #count the numbers of files that correspond to each date
nums
nums <- nums[nums[,2]!=1,] # remove dates with only one file
dates <- nums[match(dates, nums[,1]),1] #match those dates in the dates object
names(dates) <- rownames(date.info) #keep track of file info....
dates <- dates[complete.cases(dates)] #only keep non-NA values
dates
length(dates) #look at the number of files left

# next, find out how many unique dates there were
batch <- unique(dates)
batch <- cbind(batch, 1:length(batch)) #assign each unique date a batch number
batch

# next, assign those batch numbers to the list of files
batches <- batch[match(dates,batch), 2]
batches <- as.matrix(as.numeric(batches))
rownames(batches) <- names(dates)
batches

# if you're doing batch correction, remove files that won't be included due to only having one file in a batch, before doing RMA
#affy <- affy[,rownames(batches)] #subset by the filenames

pheno <- affy1[[1]]@phenoData #pulling out sample names
samples <- pheno@data #to double check
identical(rownames(samples), names(batches[(batches[,1] %in% 1),])) #if false, then something went wrong when subsetting. Check again.
samples
batches[(batches[,1] %in% 1),]

# RMA 
#affyRMA <- rma(affy)
#save(affyRMA, file="NCI60.AffyRMA.rda")
load("NCI60.AffyRMA.rda")
affy1 <- list(affyRMA[,names(batches[(batches[,1] %in% 1),])],affyRMA[,names(batches[(batches[,1] %in% 2),])], affyRMA[,names(batches[(batches[,1] %in% 3),])], affyRMA[,names(batches[(batches[,1] %in% 4),])], affyRMA[,names(batches[(batches[,1] %in% 5),])],affyRMA[,names(batches[(batches[,1] %in% 6),])],affyRMA[,names(batches[(batches[,1] %in% 7),])])
exprRMA <- exprs(affyRMA)

# frozen RMA
library("frma")
library("hgu133plus2frmavecs") #or load whichever frma vectors that correspond to your platform. NO canine frma vectors exist at this time.
#affyFRMA <- frma(affy)
#save(affyFRMA,file="NCI60.AffyFRMA.rda")
load("NCI60.AffyFRMA.rda")
affy2 <- list(affyFRMA[,names(batches[(batches[,1] %in% 1),])],affyFRMA[,names(batches[(batches[,1] %in% 2),])], affyFRMA[,names(batches[(batches[,1] %in% 3),])], affyFRMA[,names(batches[(batches[,1] %in% 4),])], affyFRMA[,names(batches[(batches[,1] %in% 5),])],affyFRMA[,names(batches[(batches[,1] %in% 6),])],affyFRMA[,names(batches[(batches[,1] %in% 7),])])
exprFRMA <- exprs(affyFRMA)

# Using ComBat for batch corrections
library("sva")

# in addition to the batch information, we may want to include additional confounding factors, such as tissue type. 
# I have a separate file that indicates the cell type for my dataset, so I read it in
fileInfo <- read.csv("filename.match.csv", header=TRUE)
rownames(fileInfo) <- fileInfo$filename
fileInfo

#match the filenames to the ones in your affy object
fileInfo <- fileInfo[rownames(samples),]
identical(rownames(samples), rownames(fileInfo))
tissue <- fileInfo$cell.type
tissue <- as.matrix(as.factor(tissue))
rownames(tissue)<-rownames(fileInfo)
tissue
# now we can apply ComBat to either RMA or fRMA
combatRMA <- ComBat(dat=exprRMA, batch=batches, mod=tissue)

combatFRMA <- ComBat(dat=exprFRMA, batch=batches, mod=tissue)

# alternatively we can use InSilicoMerge
library("inSilicoMerging")

RMAmerge <- merge(affy1, method="COMBAT")
FRMAmerge <- merge(affy2, method="COMBAT")

load("NCI60.RMA.COMBAT.batchrm.collapsed.rda") # this file was already combat processed by me and averaged by cell line

###############
input.dir<-"C:/Users/jsfowles/Documents/Final canine coxen"
drug <- "ptx"
ref.set <- "NCI59"
cox.set <- "ACC29"
ref.annot <- "HG-U133_Plus_2.na34.annot.csv"
ref.annot.skip <- 25
cox.annot <- "Canine_2.na34.annot.csv"
cox.annot.skip <- 19
ACCdrugTable<-sprintf("%s%s%s.%s","ACC29",drug,"rank","txt") #ACC29doxrank.txt

#loading in R packages for analysis
library(WGCNA)
library(qvalue)
source("collapseRows_04_11_11.R")
source("coxen.functions.R")

#Step 1
#loading NCI60RMA and ACC30 data
#load in ComBat processed reference set
ref.data <- combat.rma.collapse.2
cox.data <- ccombatRMA

#I need gene symbol file for both NCI60RMA and ACC29RMA files

ref_Anno <- read.csv(file=ref.annot,skip=ref.annot.skip,header=TRUE)
ref_Anno <- ref_Anno[,c(1,15)]

cox_Anno <- read.csv(file=cox.annot,skip=cox.annot.skip,header=TRUE)
cox_Anno <- cox_Anno[,c(1,15)]

#NCI60 data collapsing
ref.GenSym <- ref_Anno[,2]#this gives me the gene symbol for each row in the datset
dim(ref.data)
ref.data <- as.data.frame(ref.data)
ref.rowIDs <-rownames(ref.data)
ref.collapsed <- collapseRows(datET=ref.data,rowGroup=ref.GenSym,rowID=ref.rowIDs,method="MaxMean",connectivityBasedCollapsing=FALSE)
refgene <-ref.collapsed$datETcollapsed #this is the collapsed data, rows are labled by Gen symbol now instead of probe ID

#ACC30 data collapsing
cox.GenSym <- cox_Anno[,2]#this gives me the gene symbol for each row in the datset
dim(cox.data)
cox.data <- as.data.frame(cox.data)
cox.rowIDs <-rownames(cox.data)
cox.collapsed <- collapseRows(datET=cox.data,rowGroup=cox.GenSym,rowID=cox.rowIDs,method="MaxMean",connectivityBasedCollapsing=FALSE)
coxgene <-cox.collapsed$datETcollapsed #this is the collapsed data, rows are labled by Gen symbol now instead of probe ID


# Step 2: filtering to remove genes that are not represented on both human and canine chip

#using the %in% method for matching

ref.input <- refgene[rownames(refgene)%in%rownames(coxgene),]
cox.input <- coxgene[rownames(coxgene)%in%rownames(refgene),]

ref.input[1:20,1:2]
cox.input[1:20,1:2]

cox.length <- ncol(cox.input)

# Step 3:  run analysis for differenial gene expression identification in reference set.

##sort reference set based on drug sensitivity

ref.input<-t(ref.input)

drugdata <- read.csv(file=paste("NCI59_",drug,"_logvalues.csv",sep=''),dec=".",header=TRUE) ##where just the column of drug log GI50 values is saved as a csv file(the values are ordered according to cell line alphabetically)
drugdata <- drugdata[,1]
drugdata
ref.input <- cbind(drugdata,ref.input)
ref.input <- ref.input[order(ref.input[,1]),]
dim(ref.input)
drugdata.sorted <-ref.input[,1]
ref.input <- ref.input[,2:ncol(ref.input)]
ref.input <- t(ref.input)
ref.input[1:5,1:5]
dim(ref.input)


# selecting the size of the sensitive and resistant groups
nsen <- 12# the 12 most sensitive
nres <- (ncol(ref.input)-12)+1# the 12 most resistant

##Performing t.test.discovery function to identify DEGs
DEGs <- t.test.discovery(x=ref.input[,c(1:nsen)],y=ref.input[,c(nres:ncol(ref.input))])
DEGs <- DEGs[order(DEGs[,2]),]
DEGs <- DEGs[DEGs[,3]<0.1,]  


# extracting DEGs from ref.input and cox.input
drug_psid <- rownames(DEGs)
A <- drug_psid
A <- as.vector(A)
refIDs <- ref.input[A,]
refIDs[1:5,1:5]
dim(refIDs)

coxIDs<-cox.input[A,]
coxIDs[1:5,1:5]
dim(coxIDs)

##to standardize the NCI60RMA and carbayo89 dataset  
##make sure the columns are the probeset IDs/genesymbols##

ref.scaled <- apply(refIDs,1,scale)
ref.scaled<-t(ref.scaled)
colnames(ref.scaled)<-colnames(refIDs)
ref.scaled[1:5,1:5]

cox.scaled <- apply(coxIDs,1,scale)
cox.scaled<-t(cox.scaled)
colnames(cox.scaled)<-colnames(coxIDs)
cox.scaled[1:5,1:5]
#Step 4:  CO-expression step between reference and coxen sets.  Use the "coxen.step" function from Theodorescu hardrive


ModelGenes <-coxen.step(ref.scaled,cox.scaled,gene.final=rownames(ref.scaled))
write.csv(ModelGenes,file=paste(save.dir,"/","ModelGenes.",ref.set,".",drug,b.d.meth,".csv",sep=''))
