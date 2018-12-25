#tailing
setwd("~/Desktop/")

mut2 <- read.table("out_KB7.txt", sep="\t", header=FALSE, row.names = 1, stringsAsFactors = FALSE)
mut2_tails <- mut2[mut2[,2]>1,]
mut2_tails2 <-c()

for (i in 1:nrow(mut2_tails)){
  a <- unlist(strsplit(mut2_tails[i,3],""))
  length(a) <- nchar(mut2_tails[nrow(mut2_tails),3])
  mut2_tails2 <- rbind(mut2_tails2,a)
}

mut2_tails3 <- c()
mut2_tails4 <- matrix(NA, nrow(mut2_tails2), ncol(mut2_tails2))
mut2_tail_weights <- 1/(mut2_tails[,2]-1)
for (i in 1:ncol(mut2_tails2)){
  b <- c(nrow(mut2_tails2[mut2_tails2[,i]%in%"A",]),nrow(mut2_tails2[mut2_tails2[,i]%in%"C",]),nrow(mut2_tails2[mut2_tails2[,i]%in%"T",]),nrow(mut2_tails2[mut2_tails2[,i]%in%"G",]))
  length(b)<-4
  mut2_tails3 <- rbind(mut2_tails3,b)
  
  for (j in 1:nrow(mut2_tails2)){
    if (!is.na(mut2_tails2[j,i])){
      mut2_tails4[j,i] <- mut2_tail_weights[j]
    } else{mut2_tails4[j,i] <- 0}
  }
  }
colnames(mut2_tails3)<-c("A","C","T","G")
rownames(mut2_tails3)<-1:nrow(mut2_tails3)
barplot(t(mut2_tails3), beside=FALSE, legend=TRUE,col=c("palevioletred3","wheat","lightseagreen","grey40"), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,50), xlab="Position", ylab="Number")
hist((mut2_tails[,2]-1), breaks=0:45, xlim=c(1,45), main="mut-2 tailing", xlab="Length of trimmed sequence")
par(lab=c(10,10,7))

mut2_tails5 <- t(mut2_tails3)
for (i in 1:ncol(mut2_tails5)){
  c <- mut2_tails4[,i]
  mut2_tails5[1,i] <- sum(c[mut2_tails2[,i]%in%"A"])  
  mut2_tails5[2,i] <- sum(c[mut2_tails2[,i]%in%"C"]) 
  mut2_tails5[3,i] <- sum(c[mut2_tails2[,i]%in%"T"])
  mut2_tails5[4,i] <- sum(c[mut2_tails2[,i]%in%"G"]) 
}

barplot(wt_tails5, beside=TRUE, legend=TRUE,col=c("palevioletred3","wheat","lightseagreen","grey40"), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,10), xlab="Position", ylab="Number")
barplot(mut7_tails5["A",], legend=FALSE,col=sample(colors(),4), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,10), xlab="Position", ylab="Number")
barplot(mut7_tails5["T",], legend=FALSE,col=sample(colors(),4), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,10), xlab="Position", ylab="Number")
barplot(mut7_tails5["G",], legend=FALSE,col=sample(colors(),4), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,10), xlab="Position", ylab="Number")
barplot(mut7_tails5["C",], legend=FALSE,col=sample(colors(),4), main="Frequency of trimmed nucleotides", sub="mut-2(-)", ylim=c(0,10), xlab="Position", ylab="Number")

###############
mut7 <- read.table("out_KB8.txt", sep="\t", header=FALSE, row.names = 1, stringsAsFactors = FALSE)
mut7_tails <- mut7[mut7[,2]>1,]
mut7_tails2 <-c()

for (i in 1:nrow(mut7_tails)){
  a <- unlist(strsplit(mut7_tails[i,3],""))
  length(a) <- nchar(mut7_tails[nrow(mut7_tails),3])
  mut7_tails2 <- rbind(mut7_tails2,a)
}

mut7_tails3 <- c()
mut7_tails4 <- matrix(NA, nrow(mut7_tails2), ncol(mut7_tails2))
mut7_tail_weights <- 1/(mut7_tails[,2]-1)
for (i in 1:ncol(mut7_tails2)){
  b <- c(nrow(mut7_tails2[mut7_tails2[,i]%in%"A",]),nrow(mut7_tails2[mut7_tails2[,i]%in%"C",]),nrow(mut7_tails2[mut7_tails2[,i]%in%"T",]),nrow(mut7_tails2[mut7_tails2[,i]%in%"G",]))
  length(b)<-4
  mut7_tails3 <- rbind(mut7_tails3,b)
  
  for (j in 1:nrow(mut7_tails2)){
    if (!is.na(mut7_tails2[j,i])){
      mut7_tails4[j,i] <- mut7_tail_weights[j]
    } else{mut7_tails4[j,i] <- 0}
  }
}
colnames(mut7_tails3)<-c("A","C","T","G")
rownames(mut7_tails3)<-1:nrow(mut7_tails3)
barplot(t(mut7_tails3), beside=FALSE, legend=TRUE,col=c("palevioletred3","wheat","lightseagreen","grey40"), main="Frequency of trimmed nucleotides", sub="mut-7(-)", xlab="Position", ylab="Number")
hist(mut7_tails[,2]-1, breaks = 0:45,xlim=c(1,45), main="mut-7 tailing", xlab="Length of trimmed sequence")
hist(mut7_tails[,2]-1, breaks = 45,xlim=c(1,45),ylim=c(0,5), main="mut-7 tailing", xlab="Length of trimmed sequence")

mut7_tails5 <- t(mut7_tails3)
for (i in 1:ncol(mut7_tails5)){
  c <- mut7_tails4[,i]
  mut7_tails5[1,i] <- sum(c[mut7_tails2[,i]%in%"A"])  
  mut7_tails5[2,i] <- sum(c[mut7_tails2[,i]%in%"C"]) 
  mut7_tails5[3,i] <- sum(c[mut7_tails2[,i]%in%"T"])
  mut7_tails5[4,i] <- sum(c[mut7_tails2[,i]%in%"G"]) 
}

barplot(mut7_tails5,legend=TRUE,col=c("palevioletred3","wheat","lightseagreen","grey40"), main="Frequency of trimmed nucleotides", sub="mut-7(-)", ylim=c(0,20), xlab="Position", ylab="Number")



########
wt <- read.table("out_KB9.txt", sep="\t", header=FALSE, row.names = 1, stringsAsFactors = FALSE)
wt_tails <- wt[wt[,2]>1,]
wt_tails2 <-c()

for (i in 1:nrow(wt_tails)){
  a <- unlist(strsplit(wt_tails[i,3],""))
  length(a) <- nchar(wt_tails[nrow(wt_tails),3])
  wt_tails2 <- rbind(wt_tails2,a)
}

wt_tails3 <- c()
wt_tails4 <- matrix(NA, nrow(wt_tails2), ncol(wt_tails2))
wt_tail_weights <- 1/(wt_tails[,2]-1)
for (i in 1:ncol(wt_tails2)){
  b <- c(nrow(wt_tails2[wt_tails2[,i]%in%"A",]),nrow(wt_tails2[wt_tails2[,i]%in%"C",]),nrow(wt_tails2[wt_tails2[,i]%in%"T",]),nrow(wt_tails2[wt_tails2[,i]%in%"G",]))
  length(b)<-4
  wt_tails3 <- rbind(wt_tails3,b)
  
  for (j in 1:nrow(wt_tails2)){
    if (!is.na(wt_tails2[j,i])){
      wt_tails4[j,i] <- wt_tail_weights[j]
    } else{wt_tails4[j,i] <- 0}
  }
}
colnames(wt_tails3)<-c("A","C","T","G")
rownames(wt_tails3)<-1:nrow(wt_tails3)
barplot(t(wt_tails3),beside=FALSE, legend=TRUE,col=c("palevioletred3","wheat","lightseagreen","grey40"), main="Frequency of trimmed nucleotides", sub="wild type", ylim=c(0,40), xlab="Position", ylab="Number")
hist(wt_tails[,2]-1, breaks=0:45,xlim=c(1,45), main="wild type tailing", xlab="Length of trimmed sequence")


wt_tails5 <- t(wt_tails3)
for (i in 1:ncol(wt_tails5)){
  c <- wt_tails4[,i]
  wt_tails5[1,i] <- sum(c[wt_tails2[,i]%in%"A"])  
  wt_tails5[2,i] <- sum(c[wt_tails2[,i]%in%"C"]) 
  wt_tails5[3,i] <- sum(c[wt_tails2[,i]%in%"T"])
  wt_tails5[4,i] <- sum(c[wt_tails2[,i]%in%"G"]) 
}

barplot(wt_tails5, legend=TRUE,col=sample(colors(),4), main="Frequency of trimmed nucleotides", sub="wild type", ylim=c(0,10), xlab="Position", ylab="Number")
