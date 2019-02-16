setwd("~/Desktop/gonad_data/deseq")
par(mfrow=c(1,1))

for (i in 1:length(sampleFiles)){
  a<-read.table(sampleFiles[i])
  a[,2]<-round(a[,2])
  a <- a[!duplicated(a[,1]),]
  write.table(a, file=sampleFiles[i], quote = FALSE, row.names = FALSE, col.names = FALSE)
}


sampleFiles <- grep("*TOT2",list.files(),value=TRUE) 
sampleCondition <- c("N2","alg5","rde1","alg5rde1","N2","alg5","rde1","alg5rde1","N2","alg5","rde1","alg5rde1","N2","alg5","rde1","alg5rde1")
sampleTable <- data.frame(sampleName = sampleFiles[1:16],fileName = sampleFiles[1:16],condition = sampleCondition[1:16])
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory="./",design= ~ condition)
ddsHTSeq@colData

dds <- DESeq(ddsHTSeq)

#
res <- results(dds, contrast=c("condition","N2","alg5"))
summary(res)
plotMA(res, ylim=c(-4,4))
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
summary(resSig)

resm <-     data.frame(mean = res$baseMean,lfc = res$log2FoldChange,isDE = ifelse(is.na(res$padj), FALSE, res$padj < 0.05), row.names = res@rownames)

#
res <- results(dds, contrast=c("condition","N2","alg5rde1"))
summary(res)
plotMA(res, ylim=c(-4,4))
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
summary(resSig)
as.data.frame(resSig)

resm <-     data.frame(mean = res$baseMean,lfc = res$log2FoldChange,isDE = ifelse(is.na(res$padj), FALSE, res$padj < 0.05), row.names = res@rownames)

#
res <- results(dds, contrast=c("condition","rde1","N2"))
summary(res)
plotMA(res, ylim=c(-4,4))
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)
summary(resSig)
as.data.frame(resSig)

resm <-     data.frame(mean = res$baseMean,lfc = res$log2FoldChange,isDE = ifelse(is.na(res$padj), FALSE, res$padj < 0.05), row.names = res@rownames)

