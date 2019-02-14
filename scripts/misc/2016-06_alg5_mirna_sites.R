#6-9-16
setwd("~/Desktop/MontgomeryLab/")
mirnaSites <- read.table(file="ts_all_cel_mirnas_fixed.txt")

alg5_mirnas <- read.table("./new_smRNA/hpo24_mirnas.txt")

up <- read.table("TableS3_up.txt", stringsAsFactors = F)
down <- read.table("TableS3_down.txt", stringsAsFactors = F)
tables3 <- read.table("TableS3.txt")

alg5_mirna_sites <- mirnaSites[mirnaSites$miRNA_family_ID %in% alg5_mirnas$V1,]
sites<-0
siteNames<-''
for (i in 1:nrow(tables3)){
  for (j in 1:nrow(alg5_mirna_sites)){
  if (rownames(tables3[i,]) %in% alg5_mirna_sites[j,9]){
    sites <- sites+1
    siteNames <- c(siteNames,as.character(alg5_mirna_sites[j,2]))
  } else {sites <- sites+0
  siteNames <- c(siteNames)}
  }
  tables3[i,5] <- sites
  tables3[i,6]<- toString(siteNames[2:length(siteNames)])
  sites<-0
  siteNames<-''
}

tables3[,6] <- gsub("NA,","",tables3[,6])

up[,6] <- tables3[rownames(up),5]  
up[,7] <- tables3[rownames(up),6]  

down[,6] <- tables3[rownames(down),5]  
down[,7] <- tables3[rownames(down),6]  


write.table(up,"TableS3_up_2.txt", quote=F, sep="\t", col.names = F)
write.table(down,"TableS3_down_2.txt", quote=F, sep="\t", col.names = F)
