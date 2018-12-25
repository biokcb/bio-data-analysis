setwd("~/Desktop")
germ <- read.table(file="GSE696-GPL539_series_matrix.txt", skip = 84, stringsAsFactors = F, sep="\t", fill = T, row.names=1)
colnames(germ) <- germ[1,]
germ <- germ[2:20002,]
germ <- germ[,8:43]

platform <- read.table("GPL539-tbl-1.txt", stringsAsFactors = F, sep="\t",row.names = 1)
platform <- platform[,1,drop=F]
#germs <- merge(germ, platform) didn't work, takes too much time

germ <- germ[order(rownames(germ)),]
germ <- germ[2:20001,]
platform <- platform[order(rownames(platform)),,drop=F]
identical(rownames(germ),rownames(platform))
germs <- cbind(platform, germ)
germs2 <- cbind(platform, germ)


argos <- read.table("argolist2", stringsAsFactors = F)
geneID_list <- read.csv("c_elegans.WS230.geneIDs.txt", header=F, stringsAsFactors = FALSE)

for (i in 1:nrow(argos)){
  temp_gene<- argos[i,1]
  if (tolower(temp_gene)%in%geneID_list[,2]){
    argos[i,2] <- geneID_list[geneID_list[,2]%in%tolower(temp_gene),3]
  } else {argos[i,2]<-temp_gene}
}

#did not work for sago-1/wago-8(row#10), sago-2/WAGO-6(row#14), PPW-1/WAGO-7(#18), 
#PPW-2/WAGO-3(#23), alg-3, alg-4
missing<-geneID_list[geneID_list[,2]%in%c("sago-1","wago-8","sago-2","wago-6","ppw-1","wago-7","ppw-2","wago-3"),]
missing
argos[10,2]<-missing[4,3]
argos[14,2]<-missing[3,3]
argos[18,2]<-missing[1,3]
argos[23,2]<-missing[2,3]
argos[21,2] <- "T22B3.2"
argos[22,2] <- "ZK757.3"

germs <- germs[germs[,1]%in%argos[,2],]
henn1 <- germs2[germs2[,1]%in%"C02F5.6",]
pup3 <-germs2[germs2[,1]%in%"F59A3.9",]
plot(x,pup3[,2:37])

argomatch<-germs[,1]
setdiff(argos[,2],germs[,1])
rownames(germs)<-germs[,1]
germs<-germs[,2:37]

germs<-sapply(germs,as.numeric)
rownames(germs)<-argomatch

germs_avg<-matrix(NA, nrow = 24, ncol=12)
rownames(germs_avg)<-argomatch
for (i in 1:24){
  for (j in seq(1,36,3)){
    germs_avg[i,(j+2)/3] <- mean(c(germs[i,j],germs[i,j+1],germs[i,j+2]), na.rm=T)
  }
}

x <- rep(1:12,3)
x <- x[order(x)]
plot(1:13,c(seq(-2,4,0.5)), type="n")
#points(x,germs["T23D8.7",], pch=20, col=1)
#lines(1:12,germs_avg["T23D8.7",], lwd=3)

for (i in 1:24){
  lines(1:12,germs_avg[i,], col=i)
 # points(x,germs[i,], pch=20, col=i)
  text(12.3,germs_avg[i,12],labels=names(germs_avg[i,1]), cex=0.5, col=i)
}

plot(1:13,c(seq(-2,4,0.5)), type="n")
#points(x,germs["T23D8.7",], pch=20, col=1)
lines(1:12,germs_avg["T23D8.7",], lwd=3)
for (i in 9:16){
  lines(1:12,germs_avg[i,], col=i)
  text(12.3,germs_avg[i,12],labels=names(germs_avg[i,1]), cex=0.3)
}

plot(1:13,c(seq(-2,4,0.5)), type="n")
#points(x,germs["T23D8.7",], pch=20, col=1)
lines(1:12,germs_avg["T23D8.7",], lwd=3)
for (i in 17:24){
  lines(1:12,germs_avg[i,], col=i)
  text(12.3,germs_avg[i,12],labels=names(germs_avg[i,1]), cex=0.3)
}



plot(1:13,c(seq(-2,4,0.5)), type="n")
#points(x,germs["T23D8.7",], pch=20, col=1)
lines(1:12,germs_avg["T23D8.7",], lwd=3)
for (i in c(19,11,15,12,16,9,6,4,1)){
  lines(1:12,germs_avg[i,], col=i)
  text(12.3,germs_avg[i,12],labels=names(germs_avg[i,1]), cex=0.3,col=i)
}