##################################################################
#     A sample script for normalization of microarray data       #
#                                                                #
#                  R/COXEN Summer Camp 2014                      #
#                                                                #
#                                                                #
##################################################################


##### ---   1. Setting up directories / files  --- ######

setwd("/Users/Kristen/Documents/GustafsonLab") #set working/home directory 
cel.dir <- "./CEL-files/" # where the .CEL files are stored, the './' means the path starts from the working directory
save.dir <- "./NCI60/" # if you want to save files in a directory different from the home directory

#note: you can also just type these into functions directly, this is just how I do it so I can change directories and whatever else at the beginning of the code instead of in each function

filenames <- list.files(path=cel.dir, pattern=".cel") #to get all the .CEL files in a list, sometimes the pattern=".cel". Check your file types
filenames #check the list
length(filenames) #check the number of files you're processing



##### ---   2. Expression Set Creation/Processing  --- ######

library(affy) #where the affy related functions are stored. Need to load this library to continue
affy <- ReadAffy(filenames=paste(cel.dir,filenames,sep=''), phenoData=NULL) #creates an AffyBatch object based on files from above. Can have phenoData, I typically don't & just use a separate file for extra information

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
affy <- affy[,rownames(batches)] #subset by the filenames

affy1 <- affy[,names(batches[(batches[,1] %in% 1),])]
affy2 <- affy[,names(batches[(batches[,1] %in% 2),])]
affy3 <- affy[,names(batches[(batches[,1] %in% 3),])]
affy4 <- affy[,names(batches[(batches[,1] %in% 4),])]
affy5 <- affy[,names(batches[(batches[,1] %in% 5),])]

pheno <- affy1@phenoData #pulling out sample names
samples <- pheno@data #to double check
identical(rownames(samples), rownames(batches)) #if false, then something went wrong when subsetting. Check again.
samples
batches

# RMA 
affyRMA <- rma(affy)
affyRMA <- exprs(affyRMA)

# frozen RMA
library("frma")
library("hgu133plus2frmavecs") #or load whichever frma vectors that correspond to your platform. NO canine frma vectors exist at this time.
affyFRMA <- frma(affy)
affyFRMA <- exprs(affy)

# Using ComBat for batch corrections
library("sva")

# in addition to the batch information, we may want to include additional confounding factors, such as tissue type. 
# I have a separate file that indicates the cell type for my dataset, so I read it in
fileInfo <- read.csv("./COXEN/NCI60 Batch Effects/filename.match.csv", header=TRUE)
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
combatRMA <- ComBat(dat=affyRMA, batch=batches, mod=tissue)

combatFRMA <- ComBat(dat=affyRMA, batch=batches, mod=tissue)

