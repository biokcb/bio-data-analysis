##### OSX microarray analysis for Ashok #####
##### A log/history of commands used for exploratory data analysis 
##### Mainly looking at batch effects, PCA, differential gene expression
 
setwd("C:/Users/kcbrown3/Desktop/")
cel.dir <- "./Ashok/"

##### read in & look at cel files #####
filenames <- list.files(path=cel.dir, pattern=".CEL")

affy <- ReadAffy(filenames=paste(cel.dir,'/',filenames,sep=''), phenoData=NULL)

slotNames(affy)
dates <- affy@protocolData
slotNames(dates)
date.info <- dates@data
date.info

hist(affy, col=2:8)
boxplot(affy, col=2:8)

##### use a batch file based on run date #####
batches <- read.csv(file=paste(cel.dir,'batches.csv',sep=''),header=TRUE)
batches
batches.rep <- batches[,4]
mod.type <- batches[,c(6,7,8)]

##### create different processed expression matrices #####
#rma
bmp2.rma <- rma(affy)
bmp2.rma <- exprs(bmp2.rma)
boxplot(bmp2.rma, col=2:8)

#frozen rma
bmp2.frma <- frma(affy)
bmp2.frma <- exprs(bmp2.frma)
boxplot(bmp2.frma, col=2:8)

# rma + combat, batches = scan date for a total of 2
bmp2.combat <- ComBat(dat=bmp2.rma, batch=batches.rep, mod=mod.type)

# frozen rma + combat, batches same as above
bmp2.fcombat <- ComBat(dat=bmp2.frma, batch=batches.rep, mod=mod.type)

##### comparison by pca #####
# rma
pca.rma <- prcomp(t(bmp2.rma))
summary(pca.rma)

# trt vs control
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: trt vs control")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# batch by day
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: batch 1 vs 2")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[batches.rep==1],pc2[batches.rep==1], pch=20, col=2)
points(pc1[batches.rep==2],pc2[batches.rep==2], pch=20, col=3)

# patient a vs b
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: patient a vs b")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: time")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

# combat rma
pca.combat <- prcomp(t(bmp2.combat))
summary(pca.combat)

# trt vs control
plot(pca.combat$x[,1:2], type="n", main="PCA Summary for RMA +COMBAT data: trt vs control")
pc1 <- pca.combat$x[,1]
pc2 <- pca.combat$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# batch by day
plot(pca.combat$x[,1:2], type="n", main="PCA Summary for RMA +COMBAT data: batch 1 vs 2")
pc1 <- pca.combat$x[,1]
pc2 <- pca.combat$x[,2]
points(pc1[batches.rep==1],pc2[batches.rep==1], pch=20, col=2)
points(pc1[batches.rep==2],pc2[batches.rep==2], pch=20, col=3)

# patient a vs b
plot(pca.combat$x[,1:2], type="n", main="PCA Summary for RMA +COMBAT data: patient a vs b")
pc1 <- pca.combat$x[,1]
pc2 <- pca.combat$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.combat$x[,1:2], type="n", main="PCA Summary for RMA +COMBAT data: time")
pc1 <- pca.combat$x[,1]
pc2 <- pca.combat$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

# frozen rma
pca.frma <- prcomp(t(bmp2.frma))
summary(pca.frma)

# trt vs control
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: trt vs control")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# batch by day
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: batch 1 vs 2")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[batches.rep==1],pc2[batches.rep==1], pch=20, col=2)
points(pc1[batches.rep==2],pc2[batches.rep==2], pch=20, col=3)

# patient a vs b
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: patient a vs b")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: time")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

# combat frma
pca.fcombat <- prcomp(t(bmp2.fcombat))
summary(pca.fcombat)

# trt vs control
plot(pca.fcombat$x[,1:2], type="n", main="PCA Summary for fRMA +COMBAT data: trt vs control")
pc1 <- pca.fcombat$x[,1]
pc2 <- pca.fcombat$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# batch by day
plot(pca.fcombat$x[,1:2], type="n", main="PCA Summary for fRMA +COMBAT data: batch 1 vs 2")
pc1 <- pca.fcombat$x[,1]
pc2 <- pca.fcombat$x[,2]
points(pc1[batches.rep==1],pc2[batches.rep==1], pch=20, col=2)
points(pc1[batches.rep==2],pc2[batches.rep==2], pch=20, col=3)

# patient a vs b
plot(pca.fcombat$x[,1:2], type="n", main="PCA Summary for fRMA +COMBAT data: patient a vs b")
pc1 <- pca.fcombat$x[,1]
pc2 <- pca.fcombat$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.fcombat$x[,1:2], type="n", main="PCA Summary for fRMA +COMBAT data: time")
pc1 <- pca.fcombat$x[,1]
pc2 <- pca.fcombat$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

##### differential gene summary #####
##### establish groups to compare #####
# trt vs control
trt.fcombat <- bmp2.fcombat[,batches[,7]=="bmp2"]
con.fcombat <- bmp2.fcombat[,batches[,7]=="control"]

trt.frma <- bmp2.frma[,batches[,7]=="bmp2"]
con.frma <- bmp2.frma[,batches[,7]=="control"]

trt.rma <- bmp2.rma[,batches[,7]=="bmp2"]
con.rma <- bmp2.rma[,batches[,7]=="control"]

trt.combat <- bmp2.combat[,batches[,7]=="bmp2"]
con.combat <- bmp2.combat[,batches[,7]=="control"]

trt.time <- batches[batches[,7]=="bmp2",6]
con.time <- batches[batches[,7]=="control",6]

# patient a vs patient b
a.fcombat <- bmp2.fcombat[,batches[,8]=="a"]
b.fcombat <- bmp2.fcombat[,batches[,8]=="b"]

a.frma <- bmp2.frma[,batches[,8]=="a"]
b.frma <- bmp2.frma[,batches[,8]=="b"]

a.rma <- bmp2.frma[,batches[,8]=="a"]
b.rma <- bmp2.rma[,batches[,8]=="b"]

a.combat <-bmp2.combat[,batches[,8]=="a"]
b.combat <-bmp2.combat[,batches[,8]=="b"]

# 24 hr trt vs 24 hr control
trt24.fcombat <- bmp2.fcombat[,c(5,13)]
con24.fcombat <- bmp2.fcombat[,c(19,23)]

trt24.frma <- bmp2.frma[,c(5,13)]
con24.frma <- bmp2.frma[,c(19,23)]

trt24.rma <- bmp2.rma[,c(5,13)]
con24.rma <- bmp2.rma[,c(19,23)]

trt24.combat <- bmp2.combat[,c(5,13)]
con24.combat <- bmp2.combat[,c(19,23)]


# 16-24hr trt vs 16-24hr control
trtend.fcombat <- bmp2.fcombat[,c(2,4,5,10,12,13)]
conend.fcombat <- bmp2.fcombat[,c(17,19,21,23)]

trtend.frma <- bmp2.frma[,c(2,4,5,10,12,13)]
conend.frma <- bmp2.frma[,c(17,19,21,23)]

trtend.rma <- bmp2.rma[,c(2,4,5,10,12,13)]
conend.rma <- bmp2.rma[,c(17,19,21,23)]

trtend.combat <- bmp2.combat[,c(2,4,5,10,12,13)]
conend.combat <- bmp2.combat[,c(17,19,21,23)]

##### sam #####
# rma
#trt vs. control
Y <- c(rep(2,12), rep(1,12))
samdata<-list(x=cbind(a.rma,b.rma), y=Y, geneid=as.character(1:nrow(a.rma)), genenames=row.names(a.rma), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

#pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table

delta.table <- samr.compute.delta.table(samr.pac,dels=seq(0.024, 0.038, 0.002))
delta.table


##### DEGs (Differentially Expressed Genes) Summary. #####
siggenes.table<-samr.compute.siggenes.table(samr.pac, 0.032, samdata, delta.table, all.genes=FALSE)

up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo

dim(up)
dim(down)

# frma
#trt vs. control
Y <- c(rep(2,12), rep(1,12))
samdata<-list(x=cbind(a.frma,b.frma), y=Y, geneid=as.character(1:nrow(a.frma)), genenames=row.names(a.frma), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

#pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table

delta.table <- samr.compute.delta.table(samr.pac,dels=seq(1.192, 1.344, 0.002))
delta.table


##### DEGs (Differentially Expressed Genes) Summary. #####
siggenes.table<-samr.compute.siggenes.table(samr.pac, 1.316, samdata, delta.table, all.genes=FALSE)

up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo

dim(up)
dim(down)

# frma +combat
#trt vs. control
Y <- c(rep(2,12), rep(1,12))
samdata<-list(x=cbind(a.fcombat,b.fcombat), y=Y, geneid=as.character(1:nrow(a.fcombat)), genenames=row.names(a.fcombat), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

#pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table

delta.table <- samr.compute.delta.table(samr.pac,dels=seq(1.498, 1.558, 0.002))
delta.table


##### DEGs (Differentially Expressed Genes) Summary. #####
siggenes.table<-samr.compute.siggenes.table(samr.pac, 0.339, samdata, delta.table, all.genes=FALSE)

up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo

dim(up)
dim(down)

# rma +combat
#trt vs. control
Y <- c(rep(2,12), rep(1,12))
samdata<-list(x=cbind(a.combat,b.combat), y=Y, geneid=as.character(1:nrow(a.combat)), genenames=row.names(a.combat), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

#pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table

delta.table <- samr.compute.delta.table(samr.pac,dels=seq(0.256, 0.269, 0.002))
delta.table


##### DEGs (Differentially Expressed Genes) Summary. #####
siggenes.table<-samr.compute.siggenes.table(samr.pac, 0.266, samdata, delta.table, all.genes=FALSE)

up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo

dim(up)
dim(down)

##### t-test #####
#rma
trt.ttest <- t.test.discovery(trt.rma,con.rma)
trt.ttest <- trt.ttest[trt.ttest[,3]<0.001,]
dim(trt.ttest)

last.ttest <- t.test.discovery(trt24.rma,con24.rma)
last.ttest <- last.ttest[last.ttest[,3]<0.001,]
dim(last.ttest)

end.ttest <- t.test.discovery(trtend.rma,conend.rma)
end.ttest <- end.ttest[end.ttest[,3]<0.001,]
dim(end.ttest)

patient.ttest <- t.test.discovery(a.rma,b.rma)
patient.ttest <- patient.ttest[patient.ttest[,3]<0.001,]
dim(patient.ttest)

#frma
trt.ttest <- t.test.discovery(trt.frma,con.frma)
trt.ttest <- trt.ttest[trt.ttest[,3]<0.001,]
dim(trt.ttest)

last.ttest <- t.test.discovery(trt24.frma,con24.frma)
last.ttest <- last.ttest[last.ttest[,3]<0.001,]
dim(last.ttest)

end.ttest <- t.test.discovery(trtend.frma,conend.frma)
end.ttest <- end.ttest[end.ttest[,3]<0.001,]
dim(end.ttest)

patient.ttest <- t.test.discovery(a.frma,b.frma)
patient.ttest <- patient.ttest[patient.ttest[,3]<0.001,]
dim(patient.ttest)

# rma+combat
trt.ttest <- t.test.discovery(trt.combat,con.combat)
trt.ttest <- trt.ttest[trt.ttest[,3]<0.001,]
dim(trt.ttest)

last.ttest <- t.test.discovery(trt24.combat,con24.combat)
last.ttest <- last.ttest[last.ttest[,3]<0.001,]
dim(last.ttest)

end.ttest <- t.test.discovery(trtend.combat,conend.combat)
end.ttest <- end.ttest[end.ttest[,3]<0.001,]
dim(end.ttest)

patient.ttest <- t.test.discovery(a.combat,b.combat)
patient.ttest <- patient.ttest[patient.ttest[,3]<0.001,]
dim(patient.ttest)

# frma+combat
trt.ttest <- t.test.discovery(trt.fcombat,con.fcombat)
trt.ttest <- trt.ttest[trt.ttest[,3]<0.001,]
dim(trt.ttest)

last.ttest <- t.test.discovery(trt24.fcombat,con24.fcombat)
last.ttest <- last.ttest[last.ttest[,3]<0.001,]
dim(last.ttest)

end.ttest <- t.test.discovery(trtend.fcombat,conend.fcombat)
end.ttest <- end.ttest[end.ttest[,3]<0.001,]
dim(end.ttest)

patient.ttest <- t.test.discovery(a.fcombat,b.fcombat)
patient.ttest <- patient.ttest[patient.ttest[,3]<0.001,]
dim(patient.ttest)

##### correlation test - trt time #####
#rma
trt.cor <- cor.test.discovery(trt.rma,trt.time)
trt.cor <- trt.cor[trt.cor[,3]<0.001,]
dim(trt.cor)

con.cor <- cor.test.discovery(con.rma,con.time)
con.cor <- con.cor[con.cor[,3]<0.001,]
dim(con.cor)

#frma
trt.cor <- cor.test.discovery(trt.frma,trt.time)
trt.cor <- trt.cor[trt.cor[,3]<0.001,]
dim(trt.cor)

con.cor <- cor.test.discovery(con.frma,con.time)
con.cor <- con.cor[con.cor[,3]<0.001,]
dim(con.cor)

#rma + combat
trt.cor <- cor.test.discovery(trt.combat,trt.time)
trt.cor <- trt.cor[trt.cor[,3]<0.001,]
dim(trt.cor)

con.cor <- cor.test.discovery(con.combat,con.time)
con.cor <- con.cor[con.cor[,3]<0.001,]
dim(con.cor)

#frma + combat
trt.cor <- cor.test.discovery(trt.fcombat,trt.time)
trt.cor <- trt.cor[trt.cor[,3]<0.001,]
dim(trt.cor)

con.cor <- cor.test.discovery(con.fcombat,con.time)
con.cor <- con.cor[con.cor[,3]<0.001,]
dim(con.cor)
#############################################################################################
##### look at cel files #####
setwd("C:/Users/kcbrown3/Desktop/")
cel.dir <- "./Ashok/bmp6/"

filenames <- list.files(path=cel.dir, pattern=".CEL")

affy <- ReadAffy(filenames=paste(cel.dir,'/',filenames,sep=''), phenoData=NULL)

slotNames(affy)
dates <- affy@protocolData
slotNames(dates)
date.info <- dates@data
date.info

hist(affy, col=2:8)
boxplot(affy, col=2:8)

##### use a batch file based on run date #####
batches <- read.csv(file=paste(cel.dir,'bmp6batches.csv',sep=''),header=TRUE)
batches
batches.rep <- batches[,4]
mod.type <- batches[,c(6,7,8)]

##### create different processed expression matrices #####

#rma
bmp6.rma <- rma(affy)
bmp6.rma <- exprs(bmp6.rma)
boxplot(bmp6.rma, col=2:8)

#frozen rma
bmp6.frma <- frma(affy)
bmp6.frma <- exprs(bmp6.frma)
boxplot(bmp6.frma, col=2:8)

##### comparison by pca #####
# rma
pca.rma <- prcomp(t(bmp6.rma))
summary(pca.rma)

# trt vs control
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: trt vs control")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# patient a vs b
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: patient a vs b")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data: time")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

# frozen rma
pca.frma <- prcomp(t(bmp6.frma))
summary(pca.frma)

# trt vs control
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: trt vs control")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[1:16],pc2[1:16], pch=20, col=2)
points(pc1[17:24],pc2[17:24], pch=20, col=3)

# patient a vs b
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: patient a vs b")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[batches[,8]=="a"],pc2[batches[,8]=="a"], pch=20, col=2)
points(pc1[batches[,8]=="b"],pc2[batches[,8]=="b"], pch=20, col=3)

# trt time
plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data: time")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[batches[,6]==1],pc2[batches[,6]==1], pch=20, col=2)
points(pc1[batches[,6]==2],pc2[batches[,6]==2], pch=20, col=3)
points(pc1[batches[,6]==4],pc2[batches[,6]==4], pch=20, col=4)
points(pc1[batches[,6]==8],pc2[batches[,6]==8], pch=20, col=5)
points(pc1[batches[,6]==12],pc2[batches[,6]==12], pch=20, col=6)
points(pc1[batches[,6]==16],pc2[batches[,6]==16], pch=20, col=7)
points(pc1[batches[,6]==20],pc2[batches[,6]==20], pch=20, col=8)
points(pc1[batches[,6]==24],pc2[batches[,6]==24], pch=20, col=9)
text(pc1,pc2, labels=batches[,6], col=1, pos=4, offset=0.1, cex=0.8)

#######################################################################################
##### look at cel files #####
setwd("C:/Users/kcbrown3/Desktop/")
cel.dir <- "./Ashok/CEL Files"

##### KNOCK DOWN OF OSTERIX #####
filenames <- list.files(path=cel.dir, pattern=".CEL")

affy <- ReadAffy(filenames=paste(cel.dir,'/',filenames,sep=''), phenoData=NULL)

slotNames(affy)
dates <- affy@protocolData
slotNames(dates)
date.info <- dates@data
date.info

hist(affy, col=2:8)
boxplot(affy, col=2:8)

#rma
KD.rma <- rma(affy)
boxplot(KD.rma,col=2:8)
KD.rma <- exprs(KD.rma)
KD.rma[1:5,1:8]

pca.rma <- prcomp(t(KD.rma))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data", xlim=c(-50,90))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc1[c(1,5)],pc2[c(1,5)], labels=c("OSX KDA", "OSX KDB"), col=1, pos=4, offset=0.2, cex=0.7)
points(pc1[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc1[2],pc2[2], labels=c("OSX KDA+BMP"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[6],pc2[6], labels=c("OSX KDB+BMP"), col=1, pos=1, offset=0.2, cex=0.7)
points(pc1[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc1[3],pc2[3], labels=c("ContA"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[7],pc2[7], labels=c( "ContB"), col=1, pos=1, offset=0.2, cex=0.7)
points(pc1[c(4,8)],pc2[c(4,8)], pch=20, col=5)
text(pc1[c(4,8)],pc2[c(4,8)], labels=c("ContA+BMP", "ContB+BMP"), col=1, pos=4, offset=0.2, cex=0.7)

#frma
KD.frma <- frma(affy)
boxplot(KD.frma,col=2:8)
KD.frma <- exprs(KD.frma)
KD.frma[1:5,1:8]

pca.frma <- prcomp(t(KD.frma))
summary(pca.frma)

plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data", xlim=c(-50,80))
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc1[c(1,5)],pc2[c(1,5)], labels=c("OSX KDA", "OSX KDB"), col=1, pos=4, offset=0.2, cex=0.7)
points(pc1[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc1[2],pc2[2], labels=c("OSX KDA+BMP"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[6],pc2[6], labels=c("OSX KDB+BMP"), col=1, pos=1, offset=0.2, cex=0.7)
points(pc1[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc1[3],pc2[3], labels=c("ContA"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[7],pc2[7], labels=c( "ContB"), col=1, pos=1, offset=0.2, cex=0.7)
points(pc1[c(4,8)],pc2[c(4,8)], pch=20, col=5)
text(pc1[c(4,8)],pc2[c(4,8)], labels=c("ContA+BMP", "ContB+BMP"), col=1, pos=4, offset=0.2, cex=0.7)

##### OVER EXPRESSION OF OSTERIX #####
filenames <- list.files(path="./Ashok", pattern="ST.CEL")

affy1 <- read.celfiles(filenames=paste("./Ashok",'/',filenames,sep=''))
affy1 <- ReadAffy(filenames=paste("./Ashok",'/',filenames,sep=''))

slotNames(affy1)
dates <- affy1@protocolData
slotNames(dates)
date.info <- dates@data
date.info

hist(affy1, col=2:8)
boxplot(affy1, col=2:8)

# rma
OE.rma1 <- rma(affy1)
boxplot(OE.rma,col=2:8)
OE.rma1<-exprs(OE.rma1)

pca.rma <- prcomp(t(OE.rma1))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data", xlim=c(-35,50))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc1[c(1,5)],pc2[c(1,5)], labels=c("Cont1", "Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc1[c(2,6)],pc2[c(2,6)], labels=c("Cont1+BMP", "Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc1[c(3,7)],pc2[c(3,7)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8)],pc2[c(4,8)], pch=20, col=5)
text(pc1[c(4,8)],pc2[c(4,8)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

plot3d(pca.rma$x[,1:3], type="n", main="PCA Summary for RMA data")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
pc3 <- pca.rma$x[,3]
spheres3d(pc1[c(1,5)],pc2[c(1,5)],pc3[c(1,5)], radius=0.3, color=c("blue"), alpha=1, shininess=20)
text3d(pc1[c(1,5)],pc2[c(1,5)],pc3[c(1,5)], c("Cont1", "Cont2"), adj=c(0.5,1))

spheres3d(pc1[c(2,6)],pc2[c(2,6)],pc3[c(2,6)], radius=0.3, color=c("red"), alpha=1, shininess=20)
text3d(pc1[c(2,6)],pc2[c(2,6)],pc3[c(2,6)], c("Cont1+BMP", "Cont2+BMP"), adj=c(0.5,1))

spheres3d(pc1[c(3,7)],pc2[c(3,7)],pc3[c(3,7)], radius=0.3, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(3,7)],pc2[c(3,7)],pc3[c(3,7)], c("OSX OE1", "OSX OE2"), adj=c(0.5,1))

spheres3d(pc1[c(4,8)],pc2[c(4,8)],pc3[c(4,8)], radius=0.3, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(4,8)],pc2[c(4,8)],pc3[c(4,8)], c("OSX OE1+BMP", "OSX OE2+BMP"), adj=c(0.5,1))

snapshot3d(paste(cel.dir,'/',"3D-OE-PCARMA1.png",sep=''), fmt="png")
     
# frma
OE.frma1 <- frma(affy1)
OE.frma1 <- exprs(OE.frma1)
OE.frma1[1:5,1:5]

pca.frma <- prcomp(t(OE.frma1))
summary(pca.frma)

plot(pca.frma$x[,c(3,2)], type="n", main="PCA Summary for fRMA data")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
pc3 <- pca.frma$x[,3]
points(pc3[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc3[c(1,5)],pc2[c(1,5)], labels=c("Cont1", "Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc3[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc3[c(2,6)],pc2[c(2,6)], labels=c("Cont1+BMP", "Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc3[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc3[c(3,7)],pc2[c(3,7)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc3[c(4,8)],pc2[c(4,8)], pch=20, col=5)
text(pc3[c(4,8)],pc2[c(4,8)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

plot3d(pca.frma$x[,1:3], type="n", main="PCA Summary for fRMA data")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
pc3 <- pca.frma$x[,3]
spheres3d(pc1[c(1,5)],pc2[c(1,5)],pc3[c(1,5)], radius=2, color=c("blue"), alpha=1, shininess=20)
text3d(pc1[c(1,5)],pc2[c(1,5)],pc3[c(1,5)], c("Cont1", "Cont2"), adj=c(0.5,1))

spheres3d(pc1[c(2,6)],pc2[c(2,6)],pc3[c(2,6)], radius=2, color=c("red"), alpha=1, shininess=20)
text3d(pc1[c(2,6)],pc2[c(2,6)],pc3[c(2,6)], c("Cont1+BMP", "Cont2+BMP"), adj=c(0.5,1))

spheres3d(pc1[c(3,7)],pc2[c(3,7)],pc3[c(3,7)], radius=2, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(3,7)],pc2[c(3,7)],pc3[c(3,7)], c("OSX OE1", "OSX OE2"), adj=c(0.5,1))

spheres3d(pc1[c(4,8)],pc2[c(4,8)],pc3[c(4,8)], radius=2, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(4,8)],pc2[c(4,8)],pc3[c(4,8)], c("OSX OE1+BMP", "OSX OE2+BMP"), adj=c(0.5,1))

snapshot3d(paste(cel.dir,'/',"3D-OE-PCAfRMA11.png",sep=''), fmt="png")

# rma w/o control 2 (col 5)
filenames <- list.files(path=cel.dir, pattern="ST.CEL")
filenames <- filenames[c(1:4,6:8)]

affy2 <- read.celfiles(filenames=paste(cel.dir,'/',filenames,sep=''))

OE.rma <- rma(affy2)
boxplot(OE.rma,col=2:8)
OE.rma<-exprs(OE.rma)
OE.rma[1:5,1:7]

pca.rma <- prcomp(t(OE.rma))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data", xlim=c(-40,50))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[1],pc2[1], pch=20, col=2)
text(pc1[1],pc2[1], labels=c("Cont1"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,5)],pc2[c(2,5)], pch=20, col=3)
text(pc1[c(2,5)],pc2[c(2,5)], labels=c("Cont1+BMP", "Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,6)],pc2[c(3,6)], pch=20, col=4)
text(pc1[c(3,6)],pc2[c(3,6)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,7)],pc2[c(4,7)], pch=20, col=5)
text(pc1[c(4,7)],pc2[c(4,7)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

plot3d(pca.rma$x[,1:3], type="n", main="PCA Summary for RMA data")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
pc3 <- pca.rma$x[,3]
spheres3d(pc1[1],pc2[1],pc3[1], radius=0.3, color=c("blue"), alpha=1, shininess=20)
text3d(pc1[1],pc2[1],pc3[1], c("Cont1"), adj=c(0.5,1))

spheres3d(pc1[c(2,5)],pc2[c(2,5)],pc3[c(2,5)], radius=0.3, color=c("red"), alpha=1, shininess=20)
text3d(pc1[c(2,5)],pc2[c(2,5)],pc3[c(2,5)], c("Cont1+BMP", "Cont2+BMP"), adj=c(0.5,1))

spheres3d(pc1[c(3,6)],pc2[c(3,6)],pc3[c(3,6)], radius=0.3, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(3,6)],pc2[c(3,6)],pc3[c(3,6)], c("OSX OE1", "OSX OE2"), adj=c(0.5,1))

spheres3d(pc1[c(4,7)],pc2[c(4,7)],pc3[c(4,7)], radius=0.3, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(4,7)],pc2[c(4,7)],pc3[c(4,7)], c("OSX OE1+BMP", "OSX OE2+BMP"), adj=c(0.5,1))

snapshot3d(paste(cel.dir,'/',"3D-OE-PCA.png",sep=''), fmt="png")

# frma
OE.frma <- frma(affy2, verbose=TRUE, summarize="median_polish", input.vecs=hugene.1.0.st.v1frmavecs)
boxplot(OE.frma,col=2:8)
OE.frma <- exprs(OE.frma)
OE.frma[1:5,1:5]

pca.frma <- prcomp(t(OE.frma))
summary(pca.frma)

plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data", xlim=c(-150,190))
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[1],pc2[1], pch=20, col=2)
text(pc1[1],pc2[1], labels=c("Cont1"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,5)],pc2[c(2,5)], pch=20, col=3)
text(pc1[c(2,5)],pc2[c(2,5)], labels=c("Cont1+BMP", "Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,6)],pc2[c(3,6)], pch=20, col=4)
text(pc1[c(3,6)],pc2[c(3,6)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,7)],pc2[c(4,7)], pch=20, col=5)
text(pc1[c(4,7)],pc2[c(4,7)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

plot3d(pca.frma$x[,1:3], type="n", main="PCA Summary for fRMA data")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
pc3 <- pca.frma$x[,3]
spheres3d(pc1[1],pc2[1],pc3[1], radius=2, color=c("blue"), alpha=1, shininess=20)
text3d(pc1[1],pc2[1],pc3[1], c("Cont1"), adj=c(0.5,1))

spheres3d(pc1[c(2,5)],pc2[c(2,5)],pc3[c(2,5)], radius=2, color=c("red"), alpha=1, shininess=20)
text3d(pc1[c(2,5)],pc2[c(2,5)],pc3[c(2,5)], c("Cont1+BMP", "Cont2+BMP"), adj=c(0.5,1))

spheres3d(pc1[c(3,6)],pc2[c(3,6)],pc3[c(3,6)], radius=2, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(3,6)],pc2[c(3,6)],pc3[c(3,6)], c("OSX OE1", "OSX OE2"), adj=c(0.5,1))

spheres3d(pc1[c(4,7)],pc2[c(4,7)],pc3[c(4,7)], radius=2, color=c("green"), alpha=1, shininess=20)
text3d(pc1[c(4,7)],pc2[c(4,7)],pc3[c(4,7)], c("OSX OE1+BMP", "OSX OE2+BMP"), adj=c(0.5,1))

snapshot3d(paste(cel.dir,'/',"3D-OE-PCA-FRMA.png",sep=''), fmt="png")

control1 <- OE.rma1[,c(1,5)]
cont.bmp1 <- OE.rma1[,c(2,6)]
OE1 <- OE.rma1[,c(3,7)]
OE.bmp1 <- OE.rma1[,c(4,8)]

contvbmp <- t.test.discovery(control1, cont.bmp1)
contvbmp <- contvbmp[contvbmp[,3]<0.001,]

contvoe <- t.test.discovery(control1, OE1)
contvoe <- contvoe[contvoe[,3]<0.001,]

##### SCAN.... frma processed data is missing the osterix row...? #####
OE.scan <- SCAN(paste("./Ashok",'/','*ST.CEL', sep=''), verbose=TRUE)
OE.scan <- exprs(OE.scan)
OE.scan[1:5,1:5]
OE.scan[rownames(OE.scan) %in% 7963664,] # osx
OE.scan[rownames(OE.scan) %in% 8120043,] #runx2

pca.scan <- prcomp(t(OE.scan))
summary(pca.scan)

plot(pca.scan$x[,1:2], type="n", main="PCA Summary for scan data", xlim=c(-15,20))
pc1 <- pca.scan$x[,1]
pc2 <- pca.scan$x[,2]
points(pc1[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc1[c(1,5)],pc2[c(1,5)], labels=c("Cont1", "Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc1[c(2,6)],pc2[c(2,6)], labels=c("Cont1+BMP", "Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc1[c(3,7)],pc2[c(3,7)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8)],pc2[c(4,8)], pch=20, col=5)
text(pc1[c(4,8)],pc2[c(4,8)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

KD.scan <- SCAN(paste(cel.dir,'/','*0.CEL', sep=''), verbose=TRUE)
KD.scan <- exprs(KD.scan)
KD.scan[1:5,1:5]
KD.scan[rownames(KD.scan) %in% "1552340_at",]

pca.scan <- prcomp(t(KD.scan))
summary(pca.scan)

plot(pca.scan$x[,1:2], type="n", main="PCA Summary for scan data, Knockdown", xlim=c(-15,25))
pc1 <- pca.scan$x[,1]
pc2 <- pca.scan$x[,2]
points(pc1[c(1,5)],pc2[c(1,5)], pch=20, col=2)
text(pc1[c(1,5)],pc2[c(1,5)], labels=c("OSX KDA", "OSX KDB"), col=1, pos=4, offset=0.2, cex=0.7)
points(pc1[c(2,6)],pc2[c(2,6)], pch=20, col=3)
text(pc1[2],pc2[2], labels=c("OSX KDA+BMP"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[6],pc2[6], labels=c("OSX KDB+BMP"), col=1, pos=1, offset=0.3, cex=0.7)
points(pc1[c(3,7)],pc2[c(3,7)], pch=20, col=4)
text(pc1[3],pc2[3], labels=c("ContA"), col=1, pos=4, offset=0.2, cex=0.7)
text(pc1[7],pc2[7], labels=c( "ContB"), col=1, pos=1, offset=0.2, cex=0.7)
points(pc1[4],pc2[4], pch=20, col=5)
text(pc1[4],pc2[4], labels=c("ContA+BMP"), col=1, pos=2, offset=0.2, cex=0.7)
points(pc1[8],pc2[8], pch=20, col=5)
text(pc1[8],pc2[8], labels=c("ContB+BMP"), col=1, pos=4, offset=0.2, cex=0.7)

##### look at these - OSX, MMP13, RUN2X MMP1 #####
##### correlation tests.... #####
# OSX, "1552340_at" && 

#A
cor.OE.rma<- cor.test.discovery(OE.rma1[-10915,c(1:4)], OE.rma1[rownames(OE.rma1) %in% 7963664,c(1:4)])
cor.OE.rma <- cor.OE.rma[cor.OE.rma[,3]<0.00001,]
dim(cor.OE.rma)

cor.OE.scan<-cor.test.discovery(OE.scan[-10915,c(1:4)], OE.scan[rownames(OE.scan) %in% 7963664,c(1:4)], cormethod="pearson")
cor.OE.scan <- cor.OE.scan[order(cor.OE.scan[,3]),]
cor.OE.scan <- cor.OE.scan[cor.OE.scan[,3]<0.01,]
dim(cor.OE.scan)

cor.KD.rma<-cor.test.discovery(KD.rma[-71,c(1:4)], KD.rma[rownames(KD.rma) %in% "1552340_at",c(1:4)])
cor.KD.rma <- cor.KD.rma[order(cor.KD.rma[,3]),]
cor.KD.rma <- cor.KD.rma[cor.KD.rma[,3]<0.0001,]
dim(cor.KD.rma)

cor.KD.frma <-cor.test.discovery(KD.frma[-71,c(1:4)], KD.frma[rownames(KD.frma) %in% "1552340_at",c(1:4)])
cor.KD.frma <- cor.KD.frma[order(cor.KD.frma[,3]),]
cor.KD.frma <- cor.KD.frma[cor.KD.frma[,3]<0.01,]
dim(cor.KD.frma)

cor.KD.scan <-cor.test.discovery(KD.scan[-71,c(1:4)], KD.scan[rownames(KD.scan) %in% "1552340_at",c(1:4)])
cor.KD.scan <- cor.KD.scan[cor.KD.scan[,3]<0.01,]
dim(cor.KD.scan)

# B
cor.OE.rmab<- cor.test.discovery(OE.rma1[-10915,c(5:8)], OE.rma1[rownames(OE.rma1) %in% 7963664,c(5:8)])
cor.OE.rmab <- cor.OE.rmab[cor.OE.rmab[,3]<0.01,]
dim(cor.OE.rmab)

cor.OE.scanb<-cor.test.discovery(OE.scan[-10915,c(5:8)], OE.scan[rownames(OE.scan) %in% 7963664,c(5:8)])
cor.OE.scanb <- cor.OE.scanb[order(cor.OE.scanb[,3]),]
cor.OE.scanb <- cor.OE.scanb[cor.OE.scanb[,3]<0.01,]
dim(cor.OE.scanb)

cor.KD.rmab<-cor.test.discovery(KD.rma[-71,c(5:8)], KD.rma[rownames(KD.rma) %in% "1552340_at",c(5:8)])
cor.KD.rmab <- cor.KD.rmab[order(cor.KD.rmab[,3]),]
cor.KD.rmab <- cor.KD.rmab[cor.KD.rmab[,3]<0.0001,]
dim(cor.KD.rmab)

cor.KD.frmab <-cor.test.discovery(KD.frma[-71,c(5:8)], KD.frma[rownames(KD.frma) %in% "1552340_at",c(5:8)])
cor.KD.frmab <- cor.KD.frmab[order(cor.KD.frmab[,3]),]
cor.KD.frmab <- cor.KD.frmab[cor.KD.frmab[,3]<0.01,]
dim(cor.KD.frmab)

cor.KD.scanb <-cor.test.discovery(KD.scan[-71,c(5:8)], KD.scan[rownames(KD.scan) %in% "1552340_at",c(5:8)])
cor.KD.scanb <- cor.KD.scanb[cor.KD.scanb[,3]<0.01,]
dim(cor.KD.scanb)

#overlap
dim(cor.OE.rmab[rownames(cor.OE.rmab) %in% rownames(cor.OE.rma),])
dim(cor.OE.scan[rownames(cor.OE.scan) %in% rownames(cor.OE.scanb),])

oe.over <-cor.OE.rmab[rownames(cor.OE.rmab) %in% rownames(cor.OE.rma),]

dim(cor.KD.rma[rownames(cor.KD.rma) %in% rownames(cor.KD.rmab),])
over <- cor.KD.rma[rownames(cor.KD.rma) %in% rownames(cor.KD.rmab),]
lap <- cor.KD.frmab[rownames(cor.KD.frmab) %in% rownames(cor.KD.frma),]
overlap <-over[rownames(over) %in% rownames(lap),]

scanover <- cor.KD.scan[rownames(cor.KD.scan) %in% rownames(cor.KD.scanb),] 
scanover[rownames(scanover) %in% c("236028_at", "207370_at"),]
lap[rownames(lap) %in% c("209875_s_at", "1568574_x_at"),]
scanover[rownames(scanover) %in% c("205959_at"),]

dim(cor.KD.frmab[rownames(cor.KD.frmab) %in% rownames(cor.KD.frma),])
dim(cor.KD.scan[rownames(cor.KD.scan) %in% rownames(cor.KD.scanb),])

cor.test(KD.rma[15406,], rank(KD.rma[71,]))
cor.test(KD.rma[15406,], KD.rma[71,], method="spearman")
rbind(KD.rma[15406,c(5:8)], KD.rma[71,c(5:8)])
rbind(rank(KD.rma[15406,c(1:4)]), rank(KD.rma[71,c(1:4)]))

plot(KD.rma[15406,], KD.rma[71,])

# t-tests....
KD.cont.rma <- t.test.discovery(KD.rma[,c(1,5)], KD.rma[,c(3,7)])
KD.cont.rma <- KD.cont.rma[KD.cont.rma[,3]<0.0001,]
dim(KD.cont.rma)

KD.cont.frma <- t.test.discovery(KD.frma[,c(1,5)], KD.frma[,c(3,7)])
KD.cont.frma <- KD.cont.frma[KD.cont.frma[,3]<0.0001,]
dim(KD.cont.frma)

KD.cont.scan <- t.test.discovery(KD.scan[,c(1,5)], KD.scan[,c(3,7)])
KD.cont.scan <- KD.cont.scan[KD.cont.scan[,3]<0.0001,]
dim(KD.cont.scan)

bmp6.rma <- t.test.discovery(KD.rma[,c(2,6)], KD.rma[,c(4,8)])
bmp6.rma <- bmp6.rma[bmp6.rma[,3]<0.0001,]
dim(bmp6.rma)

bmp6.frma <- t.test.discovery(KD.frma[,c(2,6)], KD.frma[,c(4,8)])
bmp6.frma <- bmp6.frma[bmp6.frma[,3]<0.0001,]
dim(bmp6.frma)

bmp6.scan <- t.test.discovery(KD.scan[,c(2,6)], KD.scan[,c(4,8)])
bmp6.scan <- bmp6.scan[bmp6.scan[,3]<0.0001,]
dim(bmp6.scan)

cont.rma <- t.test.discovery(KD.rma[,c(3,7)], KD.rma[,c(4,8)])
cont.rma <- cont.rma[cont.rma[,3]<0.0001,]
dim(cont.rma)

cont.frma <- t.test.discovery(KD.frma[,c(3,7)], KD.frma[,c(4,8)])
cont.frma <- cont.frma[cont.frma[,3]<0.0001,]
dim(cont.frma)

cont.scan <- t.test.discovery(KD.scan[,c(3,7)], KD.scan[,c(4,8)])
cont.scan <- cont.scan[cont.scan[,3]<0.0001,]
dim(cont.scan)

KD.bmp6.rma <- t.test.discovery(KD.rma[,c(2,6)], KD.rma[,c(3,7)])
KD.bmp6.rma <- KD.bmp6.rma[KD.bmp6.rma[,3]<0.0001,]
dim(KD.bmp6.rma)

KD.bmp6.frma <- t.test.discovery(KD.frma[,c(2,6)], KD.frma[,c(3,7)])
KD.bmp6.frma <- KD.bmp6.frma[KD.bmp6.frma[,3]<0.0001,]
dim(KD.bmp6.frma)

KD.bmp6.scan <- t.test.discovery(KD.scan[,c(2,6)], KD.scan[,c(3,7)])
KD.bmp6.scan <- KD.bmp6.scan[KD.bmp6.scan[,3]<0.0001,]
dim(KD.bmp6.scan)

# OE
OE.cont.rma <- t.test.discovery(OE.rma1[,c(1,5)], OE.rma1[,c(3,7)])
OE.cont.rma <- OE.cont.rma[OE.cont.rma[,3]<0.0001,]
dim(OE.cont.rma)

OE.cont.frma <- t.test.discovery(OE.frma1[,c(1,5)], OE.frma1[,c(3,7)])
OE.cont.frma <- OE.cont.frma[OE.cont.frma[,3]<0.0001,]
dim(OE.cont.frma)

OE.cont.scan <- t.test.discovery(OE.scan[,c(1,5)], OE.scan[,c(3,7)])
OE.cont.scan <- OE.cont.scan[OE.cont.scan[,3]<0.0001,]
dim(OE.cont.scan)

bmp6.rma <- t.test.discovery(OE.rma1[,c(2,6)], OE.rma1[,c(4,8)])
bmp6.rma <- bmp6.rma[bmp6.rma[,3]<0.0001,]
dim(bmp6.rma)

bmp6.frma <- t.test.discovery(OE.frma1[,c(2,6)], OE.frma1[,c(4,8)])
bmp6.frma <- bmp6.frma[bmp6.frma[,3]<0.0001,]
dim(bmp6.frma)

bmp6.scan <- t.test.discovery(OE.scan[,c(2,6)], OE.scan[,c(4,8)])
bmp6.scan <- bmp6.scan[bmp6.scan[,3]<0.0001,]
dim(bmp6.scan)

cont.rma <- t.test.discovery(OE.rma1[,c(1,5)], OE.rma1[,c(2,6)])
cont.rma <- cont.rma[cont.rma[,3]<0.0001,]
dim(cont.rma)

cont.frma <- t.test.discovery(OE.frma1[,c(1,5)], OE.frma1[,c(2,6)])
cont.frma <- cont.frma[cont.frma[,3]<0.0001,]
dim(cont.frma)

cont.scan <- t.test.discovery(OE.scan[,c(1,5)], OE.scan[,c(2,6)])
cont.scan <- cont.scan[cont.scan[,3]<0.0001,]
dim(cont.scan)

OE.bmp6.rma <- t.test.discovery(OE.rma1[,c(2,6)], OE.rma1[,c(3,7)])
OE.bmp6.rma <- OE.bmp6.rma[OE.bmp6.rma[,3]<0.0001,]
dim(OE.bmp6.rma)

OE.bmp6.frma <- t.test.discovery(OE.frma1[,c(2,6)], OE.frma1[,c(3,7)])
OE.bmp6.frma <- OE.bmp6.frma[OE.bmp6.frma[,3]<0.0001,]
dim(OE.bmp6.frma)

OE.bmp6.scan <- t.test.discovery(OE.scan[,c(2,6)], OE.scan[,c(3,7)])
OE.bmp6.scan <- OE.bmp6.scan[OE.bmp6.scan[,3]<0.0001,]
dim(OE.bmp6.scan)

################### KO hugene #####
filenames <- list.files(path=cel.dir, pattern=".CEL")

affy <- read.celfiles(filenames=paste(cel.dir,'/',filenames,sep=''))
affy2 <- ReadAffy(filenames=paste(cel.dir,'/',filenames,sep=''))

slotNames(affy2)
dates <- affy2@protocolData
slotNames(dates)
date.info <- dates@data
date.info

hist(affy, col=2:8)
boxplot(affy, col=2:8)

# rma
KD.rma1 <- rma(affy)
KD.rma1<-exprs(KD.rma1)
KD.rma1[1:5,1:5]

pca.rma <- prcomp(t(KD.rma1))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data")
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

# frma
KD.frma1 <- frma(affy)
KD.frma1 <- exprs(KD.frma1)
KD.frma1[1:5,1:5]

pca.frma <- prcomp(t(KD.frma1))
summary(pca.frma)

plot(pca.frma$x[,1:2], type="n", main="PCA Summary for fRMA data")
pc1 <- pca.frma$x[,1]
pc2 <- pca.frma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)


##### SCAN.... frma processed data is missing the osterix row...?
KD.scan <- SCAN(paste(cel.dir,'/','*.CEL', sep=''), verbose=TRUE)
KD.scan <- exprs(KD.scan)
KD.scan[1:5,1:5]
KD.scan[rownames(KD.scan) %in% 7963664,] # osx
KD.scan[rownames(KD.scan) %in% 8120043,] #runx2

pca.scan <- prcomp(t(KD.scan))
summary(pca.scan)

plot(pca.scan$x[,1:2], type="n", main="PCA Summary for scan data", xlim=c(-15,30))
pc1 <- pca.scan$x[,1]
pc2 <- pca.scan$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

##### RMA processed separately, then PCA together....
pca.rma <- prcomp(t(cbind(KD.rma1, OE.rma1)))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data", xlim=c(-250,450))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

points(pc1[c(13,17)],pc2[c(13,17)], pch=20, col=2)
text(pc1[c(13,17)],pc2[c(13,17)], labels=c("OE Cont1", "OE Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(14,18)],pc2[c(14,18)], pch=20, col=3)
text(pc1[c(14,18)],pc2[c(14,18)], labels=c("OE Cont1+BMP", "OE Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(15,19)],pc2[c(15,19)], pch=20, col=4)
text(pc1[c(15,19)],pc2[c(15,19)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(16,20)],pc2[c(16,20)], pch=20, col=5)
text(pc1[c(16,20)],pc2[c(16,20)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

##### SCAN then put together for PCA...
pca.scan <- prcomp(t(cbind(KD.scan,OE.scan)))
summary(pca.scan)

plot(pca.scan$x[,1:2], type="n", main="PCA Summary for scan data", xlim=c(-35,60))
pc1 <- pca.scan$x[,1]
pc2 <- pca.scan$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

points(pc1[c(13,17)],pc2[c(13,17)], pch=20, col=2)
text(pc1[c(13,17)],pc2[c(13,17)], labels=c("OE Cont1", "OE Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(14,18)],pc2[c(14,18)], pch=20, col=3)
text(pc1[c(14,18)],pc2[c(14,18)], labels=c("OE Cont1+BMP", "OE Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(15,19)],pc2[c(15,19)], pch=20, col=4)
text(pc1[c(15,19)],pc2[c(15,19)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(16,20)],pc2[c(16,20)], pch=20, col=5)
text(pc1[c(16,20)],pc2[c(16,20)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)


##### process all the data together. #####
filenames1 <- list.files(path=cel.dir, pattern=".CEL")
filenames2<- list.files(path="./Ashok", pattern="ST.CEL")

affy3 <- read.celfiles(filenames=c(paste0(cel.dir,'/',filenames1), paste0("./Ashok/", filenames2)))

# rma
rma2 <- rma(affy3)
rma2<-exprs(rma2)

# RMA processed separately, then PCA together....
pca.rma <- prcomp(t(rma2))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data", xlim=c(-100, 150))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,5,9)],pc2[c(1,5,9)], labels=c("Cont1", "Cont2", "Cont3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

points(pc1[c(13,17)],pc2[c(13,17)], pch=20, col=2)
text(pc1[c(13,17)],pc2[c(13,17)], labels=c("OE Cont1", "OE Cont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(14,18)],pc2[c(14,18)], pch=20, col=3)
text(pc1[c(14,18)],pc2[c(14,18)], labels=c("OE Cont1+BMP", "OE Cont2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(15,19)],pc2[c(15,19)], pch=20, col=4)
text(pc1[c(15,19)],pc2[c(15,19)], labels=c("OSX OE1", "OSX OE2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(16,20)],pc2[c(16,20)], pch=20, col=5)
text(pc1[c(16,20)],pc2[c(16,20)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

###### scan processing should remain the same due to the fact that it uses only 
###### information within each array to normalize instead of depending on the group

###### Combat processing
# batch info...
fileinfo <- read.csv(file="/Volumes/DanT-2TB/batcheinfo.csv", header=TRUE)
fileinfo
rownames(fileinfo) <- fileinfo[,1]

batches <- fileinfo[,2,drop=FALSE]
batches

mod <- fileinfo[,4:5,drop=FALSE]
mod <- fileinfo[,3,drop=FALSE]
mod

##### combat function using data RMA processed separately
combat.sep1 <- ComBat(dat=cbind(OE.rma1, KD.rma1) , batch=batches, mod=mod)
combat.sep[1:5,1:5]
combat.sep1 <- cbind(combat.sep1[,9:20], combat.sep1[,1:8])

combat.sep[rownames(combat.sep) %in% 7963664,] # osx
combat.sep[rownames(combat.sep) %in% 8120043,] #runx2

pca.rma <- prcomp(t(combat.sep))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data",xlim=c(-55, 35))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,9)],pc2[c(1,9)], labels=c("KDCont1", "KDCont3"), col=1, pos=3, offset=0.1, cex=0.8)
text(pc1[c(5)],pc2[c(5)], labels=c("KDCont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,6,10)],pc2[c(2,6,10)], labels=c("Cont1+BMP", "Cont2+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=4)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=5)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

points(pc1[c(13,17)],pc2[c(13,17)], pch=20, col=2)
text(pc1[c(13,17)],pc2[c(13,17)], labels=c("OE Cont1", "OE Cont2"), col=1, pos=2, offset=0.1, cex=0.8)
points(pc1[c(14,18)],pc2[c(14,18)], pch=20, col=3)
text(pc1[c(14,18)],pc2[c(14,18)], labels=c("OE Cont1+BMP", "OE Cont2+BMP"), col=1, pos=1, offset=0.1, cex=0.8)
points(pc1[c(15,19)],pc2[c(15,19)], pch=20, col=4)
text(pc1[c(15)],pc2[c(15)], labels=c("OSX OE1"), col=1, pos=4, offset=0.2, cex=0.8)
text(pc1[c(19)],pc2[c(19)], labels=c("OSX OE2"), col=1, pos=3, offset=0.2, cex=0.8)
points(pc1[c(16,20)],pc2[c(16,20)], pch=20, col=5)
text(pc1[c(16,20)],pc2[c(16,20)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=3, offset=0.2, cex=0.8)

##### combat function using data RMA processed together
combat.tog1 <- ComBat(dat=cbind(rma2[,13:20], rma2[,1:12]), batch=batches, mod=mod)
combat.tog[1:5,1:5]
combat.tog1 <- cbind(combat.tog1[,9:20], combat.tog1[,1:8])

combat.tog[rownames(combat.tog) %in% 7963664,] # osx
combat.tog[rownames(combat.tog) %in% 8120043,] #runx2

pca.rma <- prcomp(t(combat.tog))
summary(pca.rma)

plot(pca.rma$x[,1:2], type="n", main="PCA Summary for RMA data",xlim=c(-55, 35))
pc1 <- pca.rma$x[,1]
pc2 <- pca.rma$x[,2]
points(pc1[c(1,5,9)],pc2[c(1,5,9)], pch=20, col=2)
text(pc1[c(1,9)],pc2[c(1,9)], labels=c("KDCont1", "KDCont3"), col=1, pos=1, offset=0.1, cex=0.8)
text(pc1[c(5)],pc2[c(5)], labels=c("KDCont2"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(2,6,10)],pc2[c(2,6,10)], pch=20, col=3)
text(pc1[c(2,10)],pc2[c(2,10)], labels=c("Cont1+BMP", "Cont3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
text(pc1[c(6)],pc2[c(6)], labels=c("Cont2+BMP"), col=1, pos=1, offset=0.1, cex=0.8)
points(pc1[c(3,7,11)],pc2[c(3,7,11)], pch=20, col=6)
text(pc1[c(3,7,11)],pc2[c(3,7,11)], labels=c("OSX KD1", "OSX KD2", "OSX KD3"), col=1, pos=4, offset=0.1, cex=0.8)
points(pc1[c(4,8,12)],pc2[c(4,8,12)], pch=20, col=7)
text(pc1[c(4,8,12)],pc2[c(4,8,12)], labels=c("OSX KD1+BMP", "OSX KD2+BMP", "OSX KD3+BMP"), col=1, pos=4, offset=0.1, cex=0.8)

points(pc1[c(13,17)],pc2[c(13,17)], pch=20, col=2)
text(pc1[c(13,17)],pc2[c(13,17)], labels=c("OE Cont1", "OE Cont2"), col=1, pos=2, offset=0.1, cex=0.8)
points(pc1[c(14,18)],pc2[c(14,18)], pch=20, col=3)
text(pc1[c(14)],pc2[c(14)], labels=c("OE Cont1+BMP"), col=1, pos=4, offset=0.1, cex=0.8)
text(pc1[c(18)],pc2[c(18)], labels=c("OE Cont2+BMP"), col=1, pos=1, offset=0.1, cex=0.8)
points(pc1[c(15,19)],pc2[c(15,19)], pch=20, col=4)
text(pc1[c(15)],pc2[c(15)], labels=c("OSX OE1"), col=1, pos=4, offset=0.2, cex=0.8)
text(pc1[c(19)],pc2[c(19)], labels=c("OSX OE2"), col=1, pos=3, offset=0.2, cex=0.8)
points(pc1[c(16,20)],pc2[c(16,20)], pch=20, col=5)
text(pc1[c(16,20)],pc2[c(16,20)], labels=c("OSX OE1+BMP", "OSX OE2+BMP"), col=1, pos=3, offset=0.2, cex=0.8)

######
rbind(combat.tog[rownames(combat.tog) %in% 7963664,], combat.sep[rownames(combat.sep) %in% 7963664,], (combat.tog[rownames(combat.tog) %in% 7963664,] - combat.sep[rownames(combat.sep) %in% 7963664,]))
rbind(combat.tog1[rownames(combat.tog1) %in% 7963664,], combat.sep1[rownames(combat.sep1) %in% 7963664,], (combat.tog1[rownames(combat.tog1) %in% 7963664,] - combat.sep1[rownames(combat.sep1) %in% 7963664,]))

save(combat.tog, combat.sep, combat.tog1, combat.sep1, file="combatobjects.rda")
save(KD.rma1, OE.rma1, file="rmaobjects.rda")
save(affy, affy1, affy3, file="affyobjects.rda")

##### limma. #####
library(limma)
Info <- read.csv(file="/Volumes/DanT-2TB/batcheinfo.csv", header=TRUE)
design <- model.matrix(~ -1 + factor(Info$condition,levels=unique(Info$condition)))#creates design matrix for linear model
colnames(design) <- unique(Info$condition)
design

#adding bmp
cont.matrix <- makeContrasts(control-controlbmp,levels=design)#creates contrast matrix 
#kd
cont.matrix <- makeContrasts(controlbmp-kdbmp,levels=design)
# overexpression
cont.matrix <- makeContrasts(controlbmp-oebmp,levels=design)

fit <- lmFit(combat.tog, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

ModtPvals<-data.frame(PSID=rownames(fit2$p.value),modtpvals=fit2$p.value,row.names=NULL)
save(ModtPvals,file="ModtPvals.rda")

length(fit2$p.value[fit2$p.value<0.01])#find out how many genes pass p-value cutoffs
hist(fit2$p.value,breaks=50)

# Volcano Plot - visualize the differentially expressed genes
plot(fit2$coef,-log10(fit2$p.value),xlab="logFC",ylab="-log(p)")
title("Volcano Plot")
abline(h=3,col="red");abline(v=+1,col="red");abline(v=-1,col="red")#creates cutoff lines at the p<0.001 and logFC>1 marks on plot

AllGenes<-topTable(fit2,number=nrow(combat.tog),genelist=NULL,lfc=.59, sort.by="logFC",adjust.method="BH",p.value=0.05)
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)

write.csv(AllGenes, file="control-controlbmpLimma.csv")
write.csv(AllGenes, file="oebmp-controlbmpLimma.csv")
write.csv(AllGenes, file="kdbmp-controlbmpLimma.csv")

# we want to bind gene annotation onto this data table
Annotation<-read.csv("hugeneinfo.csv",header=TRUE)#these files are available at Affymetrix.com for each array typ

# Merge results, sample means and annotation
AllGenes<-merge(AllGenes,Annotation,by.x="Probe.Set.ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
AllGenes.sorted <- AllGenes[order(AllGenes$P.Value),]

# Heatmap
library(gplots)
TopList<-combat.tog[AllGenes$Probe.Set.ID,]
dim(TopList)
colnames(TopList) <- Info[match(colnames(TopList),Info[,1]),3]
colnames(TopList)
heatmap.2(TopList[,colnames(TopList) %in% c("oebmp", "controlbmp")],col="redgreen",trace="none")

## ttest
library(multtest)
CL<-(Info$condition=="control")*(Info$condition=="control")
tstats<-mt.teststat(combat.tog,CL,test="t.equalvar")
tpvals<-2*pt(abs(tstats),df=54,lower.tail=FALSE)
OrdtPvals<-data.frame(PSID=row.names(combat.tog),tpvals)
OrdtPvals<- OrdtPvals[order(OrdtPvals[,2]),]

# control v controlbmp
Infocont <- Info[Info$cond1=="control",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,5), rep(2,5))
genematrix <- combat.tog[,match(Infocont$Filename,colnames(combat.tog))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 0.592
del.min <- 0.498
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 0.592

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="controlbmp-controlSAM.csv")

# controlbmp v oebmp
Infocont <- Info[Info$condition=="controlbmp" | Info$condition=="oebmp",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,5), rep(2,2))
genematrix <- combat.tog[,match(Infocont$Filename,colnames(combat.tog))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 0.587
del.min <- 0.515
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 0.519

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="controlbmp-oebmpSAM.csv")

# controlbmp v kdbmp
Infocont <- Info[Info$condition=="controlbmp" | Info$condition=="kdbmp",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,5), rep(2,3))
genematrix <- combat.tog[,match(Infocont$Filename,colnames(combat.tog))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 0.855
del.min <- 0.799
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 0.841

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="controlbmp-kdbmpSAM.csv")

########## Now for RMA processed only files...
rm(combat.sep, combat.sep1, combat.tog1)
load("RMAobjects.rda")

####### 
# control v controlbmp
Infocont <- Info[Info$cond1=="control",]
Infocont <- Infocont[Infocont$batch=="1",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,3), rep(2,3))
genematrix <- KD.rma1[,match(Infocont$Filename,colnames(KD.rma1))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)
#pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 0.819
del.min <- 0.719
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 0.807
# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="RMA.KDcontrolbmp-controlSAM.csv")

# control v controlbmp
Infocont <- Info[Info$cond1=="control",]
Infocont <- Infocont[Infocont$batch=="2",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,2), rep(2,2))
genematrix <- OE.rma1[,match(Infocont$Filename,colnames(OE.rma1))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 0.590
del.min <- 0.488
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 0.538

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="RMA.OEcontrolbmp-controlSAM.csv")

# controlbmp v oebmp
Infocont <- Info[Info$condition=="controlbmp" | Info$condition=="oebmp",]
Infocont <- Infocont[Infocont$batch=="2",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,2), rep(2,2))
genematrix <- OE.rma1[,match(Infocont$Filename,colnames(OE.rma1))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 1.12
del.min <- 1.02
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 1.062

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="RMAcontrolbmp-oebmpSAM.csv")

#controlbmp v kdbmp
Infocont <- Info[Info$condition=="controlbmp" | Info$condition=="kdbmp",]
Infocont <- Infocont[Infocont$batch=="1",]
Infocont <- Infocont[order(Infocont$condition),]
Infocont

CL <- c(rep(1,3), rep(2,3))
genematrix <- KD.rma1[,match(Infocont$Filename,colnames(KD.rma1))]
samdata<-list(x=genematrix, y=CL, geneid=as.character(1:nrow(genematrix)), genenames=row.names(genematrix), logged2=TRUE)
samr.pac<-samr(samdata, resp.typ="Two class unpaired", nperms=500)

# pick out delta value closest to 0.1
delta.table <- samr.compute.delta.table(samr.pac)
delta.table
del.max <- 1.281
del.min <- 1.177
delta.table <- samr.compute.delta.table(samr.pac,dels=seq(del.min, del.max, 0.002))
delta.table
delta.opt <- 1.247

# DEGs (Differentially Expressed Genes) Summary.
siggenes.table<-samr.compute.siggenes.table(samr.pac, delta.opt, samdata, delta.table, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
DEGs<-rbind(up,down) #bind high/low table rows
DEGs <- merge(DEGs,Annotation,by.x="Gene ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work
write.csv(DEGs, file="RMAcontrolbmp-kdbmpSAM.csv")

####### limma. #######
#kd
Infocont <- Info[Info$batch=="1",]
design <- model.matrix(~ -1 + factor(Infocont$condition,levels=unique(Infocont$condition)))#creates design matrix for linear model
colnames(design) <- unique(Infocont$condition)
design

#adding bmp
cont.matrix <- makeContrasts(control-controlbmp,levels=design)#creates contrast matrix 
#kd
cont.matrix <- makeContrasts(controlbmp-kdbmp,levels=design)
#kd
cont.matrix <- makeContrasts(control-knockdown,levels=design)

fit <- lmFit(KD.rma1, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

ModtPvals<-data.frame(PSID=rownames(fit2$p.value),modtpvals=fit2$p.value,row.names=NULL)
save(ModtPvals,file="ModtPvals.rda")

length(fit2$p.value[fit2$p.value<0.01])#find out how many genes pass p-value cutoffs
hist(fit2$p.value,breaks=50)

AllGenes<-topTable(fit2,number=nrow(KD.rma1),genelist=NULL,lfc=.59, sort.by="logFC",adjust.method="BH",p.value=0.05)
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)
# we want to bind gene annotation onto this data table
setwd("/Users/Kristen/Documents/GustafsonLab")
Annotation<-read.csv("hugeneinfo.csv",header=TRUE)#these files are available at Affymetrix.com for each array typ
# Merge results, sample means and annotation
AllGenes<-merge(AllGenes,Annotation,by.x="Probe.Set.ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work

write.csv(AllGenes, file="RMA.kdcontrol-controlbmpLimma.csv")
write.csv(AllGenes, file="RMAkdbmp-controlbmpLimma.csv")
write.csv(AllGenes, file="RMAcontrol-kdLimma.csv")

####### limma. 
#oe
Infocont <- Info[Info$batch=="2",]
design <- model.matrix(~ -1 + factor(Infocont$condition,levels=unique(Infocont$condition)))#creates design matrix for linear model
colnames(design) <- unique(Infocont$condition)
design

#adding bmp
cont.matrix <- makeContrasts(control-controlbmp,levels=design)#creates contrast matrix 
#oe
cont.matrix <- makeContrasts(controlbmp-oebmp,levels=design)
#oe
cont.matrix <- makeContrasts(control-overexpression,levels=design)

fit <- lmFit(OE.rma1, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

ModtPvals<-data.frame(PSID=rownames(fit2$p.value),modtpvals=fit2$p.value,row.names=NULL)
save(ModtPvals,file="ModtPvals.rda")

length(fit2$p.value[fit2$p.value<0.01])#find out how many genes pass p-value cutoffs
hist(fit2$p.value,breaks=50)

AllGenes<-topTable(fit2,number=nrow(OE.rma1),genelist=NULL,lfc=.59, sort.by="logFC",adjust.method="BH",p.value=0.05)
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)
# we want to bind gene annotation onto this data table
AllGenes<-merge(AllGenes,Annotation,by.x="Probe.Set.ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work

write.csv(AllGenes, file="RMA.oecontrol-controlbmpLimma.csv")
write.csv(AllGenes, file="RMAoebmp-controlbmpLimma.csv")
write.csv(AllGenes, file="RMAcontrol-oeLimma.csv")

######## TO DO
# double check FC are relative to control
### control vs oe, control vs kd
## write up quality control summary and pca summary
setwd("/Users/Kristen/Desktop/CEL Files")
ncol(OE.rma1)
plotMA(KD.rma1, main="MA Plot for KD")
pairs(KD.rma1, pch=".",main="Scatter plots",col=1:12) 
plotMA(OE.rma1, main="MA Plot for OE")
pairs(KD.rma1, pch=".",main="Scatter plots",col=1:12) 

## check on red/green genes in the paper for comparison. p-values
biomarkers <- c(8162373,8152512,8152617,8112971,8060850,8141140,8057506,8161919,7927631,8162394,8162388,7898693,7960947,7951271,7951309, 7963664)

setwd("/Users/Kristen/Documents/GustafsonLab")
Annotation<-read.csv("hugeneinfo.csv",header=TRUE)
Info <- fileinfo

####### limma. 
#kd
Infocont <- Info[Info$batch=="1",]
design <- model.matrix(~ -1 + factor(Infocont$condition,levels=unique(Infocont$condition)))#creates design matrix for linear model
colnames(design) <- unique(Infocont$condition)
design

#adding bmp
cont.matrix <- makeContrasts(controlbmp-control,levels=design)#creates contrast matrix 
#kd
cont.matrix <- makeContrasts(kdbmp-controlbmp,levels=design)
#kd
cont.matrix <- makeContrasts(knockdown-control,levels=design)

fit <- lmFit(KD.rma1, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

ModtPvals<-data.frame(PSID=rownames(fit2$p.value),modtpvals=fit2$p.value,row.names=NULL)

length(fit2$p.value[fit2$p.value<0.01])#find out how many genes pass p-value cutoffs
hist(fit2$p.value,breaks=50)

AllGenes<-topTable(fit2,number=nrow(KD.rma1),genelist=NULL, sort.by="logFC",adjust.method="BH")
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)

# Merge results, sample means and annotation
AllGenes<-merge(AllGenes,Annotation,by.x="Probe.Set.ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work

write.csv(AllGenes, file="RMA.kdcontrol-controlbmpLimmaALL.csv")
write.csv(AllGenes, file="RMAkdbmp-controlbmpLimmaALL.csv")
write.csv(AllGenes, file="RMAcontrol-kdLimmaALL.csv")

kd.biom <- AllGenes[AllGenes[,1] %in% biomarkers, ]
write.csv(kd.biom, file="biomarkers_controlbmpvcontrol.csv")

####### limma. ######
# oe
Infocont <- Info[Info$batch=="2",]
design <- model.matrix(~ -1 + factor(Infocont$condition,levels=unique(Infocont$condition)))#creates design matrix for linear model
colnames(design) <- unique(Infocont$condition)
design

# adding bmp
cont.matrix <- makeContrasts(controlbmp-control,levels=design)#creates contrast matrix 
# oe
cont.matrix <- makeContrasts(oebmp-controlbmp,levels=design)
# oe
cont.matrix <- makeContrasts(overexpression-control,levels=design)

fit <- lmFit(OE.rma1, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

ModtPvals<-data.frame(PSID=rownames(fit2$p.value),modtpvals=fit2$p.value,row.names=NULL)
length(fit2$p.value[fit2$p.value<0.01])#find out how many genes pass p-value cutoffs
hist(fit2$p.value,breaks=50)

AllGenes<-topTable(fit2,number=nrow(OE.rma1),genelist=NULL,lfc=0, sort.by="logFC",adjust.method="BH",p.value=1)
AllGenes<-data.frame(Probe.Set.ID=row.names(AllGenes),AllGenes)

# We want to bind gene annotation onto this data table
AllGenes<-merge(AllGenes,Annotation,by.x="Probe.Set.ID",by.y="Probe.Set")#both objects need a column named "Probe.Set.ID" for this to work

# Write files
write.csv(AllGenes, file="RMA.oecontrol-controlbmpLimmaALL.csv")
write.csv(AllGenes, file="RMAoebmp-controlbmpLimma0.5.csv")
write.csv(AllGenes, file="RMAcontrol-oeLimmaALL.csv")
oe.biom <- AllGenes[AllGenes[,1] %in% biomarkers, ]
write.csv(oe.biom, file="biomarkers_controlvoe.csv")