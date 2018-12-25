###################### Kristen's commonly used functions ######################
#                                                                             #
# a bunch of custom functions I use in R to do data analysis and plotting     #
# there's probably many improvements to be made, but whatever I'm lazy        # 
#                                                                             #
#                                                                             #
#                                                                             #
# Last update: 2/3/17                                                         #
#                                                                             #
##############################################################################


# adding log2 tick marks to a generic plot. Uses base R plotting #
plot_log2_ticks<-function(min,max){

tick_sep <- seq(from=min, to=max, by=4)
tick_sep <- 2^tick_sep
for (i in 1:length(tick_sep)){
  #major ticks
  axis(1,at=log2(tick_sep[i]), labels=prettyNum(tick_sep[i],digits=2,format="f"),las=1,cex.axis=0.8)
  axis(2,at=log2(tick_sep[i]), labels=prettyNum(tick_sep[i],digits=2,format="f"),las=1,cex.axis=0.8)
}

for (i in 2:(length(tick_sep))){
  #add minor ticks
  minor_sep <- (tick_sep[i] - tick_sep[i-1])/8
  minor_ticks <- seq(from=tick_sep[i-1], to=tick_sep[i], by=minor_sep)
  axis(1,at=log2(minor_ticks[1:8]), labels=F, tck=-0.005)
  axis(2,at=log2(minor_ticks[1:8]), labels=F, tck=-0.005)
  axis(1, at=log2(minor_ticks[5]), tck=-0.009, labels=F)
  axis(2, at=log2(minor_ticks[5]), tck=-0.009, labels=F)
}
}

#creates scatter plots showing enrichment or depletion of small RNAs in mutants or IPs. Makes a PDF if desired
smrna_scatter <- function(sample1, sample2, classes, ylim, log2Trans = TRUE, sample1Name = "Wild type", sample2Name="mutant",pdfOpt=3, plotSeparate=FALSE, plotNumPerPage=c(1,1), pdfSize=c(10,10), highlight,highlightclass=3){
  #samples should be a vector of filenames containing smRNA read totals for each class to be plotted
  #if the reads should be log2 transformed, specify log2Trans as TRUE
  #classes are the descriptors for each class of small rna in the list, used for color coding purposes

  
  if (length(sample1)==length(sample2) & length(sample1)==length(classes)){
    sort(classes)
    #get variables (counts) from filenames
    for (i in 1:length(sample1)){
      temp <- read.table(file=sample1[i],sep="\t",stringsAsFactors = FALSE)
      temp <- temp[!duplicated(temp[,1]),]
      rownames(temp) <- temp[,1]
      temp <- temp[,2,drop=FALSE]
      
      temp2 <- read.table(file=sample2[i],sep="\t",stringsAsFactors = FALSE)
      temp2 <- temp2[!duplicated(temp2[,1]),]
      rownames(temp2) <- temp2[,1]
      temp2 <- temp2[,2,drop=FALSE]
      
      
      assign(paste0("sample1_",classes[i]), temp)
      assign(paste0("sample2_",classes[i]), temp2)
    }
  
    #place each into lists
    s1 <- mget(ls(pattern="sample1_"))
    s2 <- mget(ls(pattern="sample2_"))
    
    #log2 transform reads
    if (log2Trans){
      for (i in 1:length(s1)){
        s1[[i]]<-log2(s1[[i]])
      }
      for (i in 1:length(s2)){
        s2[[i]]<-log2(s2[[i]])
      }
    }
    
    #calculate plot boundaries
   
     if (missing(ylim)){
       am <- c(sapply(s1, function(x) min(x[!is.infinite(x)])),sapply(s2, function(x) min(x[!is.infinite(x)])),sapply(s1, function(x) max(x[!is.infinite(x)])),sapply(s2, function(x) max(x[!is.infinite(x)])))
      plMin <- round(min(am[!is.infinite(am)]))
      plMax <- round(max(am[!is.infinite(am)]))
    } else {
      plMin <- ylim[1]
      plMax <- ylim[2]
    }
    if (pdfOpt==1 | pdfOpt==3){
      if (plotSeparate){
        pdf(file=paste0(gsub("/","",sample2Name),"_","scatterPlots.pdf"),width=pdfSize[1],height=pdfSize[2])
        par(mfrow=plotNumPerPage,mar = c(4.5, 7, 2.5, 2.5), mgp = c(3, 1, 0),pty="s")
        #populate with points for each class
        for (i in 1:length(classes)){
          #create an empty plot window
          plot(plMin:plMax,plMin:plMax, type="n", yaxt='n', xaxt='n', ann=FALSE)
          points(s1[[i]][,1], s2[[i]][,1], pch=15,col="tomato1", cex=0.3)
          abline(h=log2(10), v=log2(10), lty=2, col="grey80")
          abline(0,1, col="grey80")
          abline(log2(0.5),1, col="grey80")
          abline(log2(2),1, col="grey80")
          if (i == highlightclass){
            if (length(highlight)>=1){
              points(s1[[i]][highlight,1],s2[[i]][highlight,1], pch=20,col="deepskyblue", cex=0.3)
            }
          }
          title(main=paste0("Enrichment of ",classes[i]," in ",sample2Name), xlab = paste0("Normalized reads in ",sample1Name),ylab = paste0("Normalized reads in ",sample2Name,"\n"),font.lab=2)
          plot_log2_ticks(plMin,plMax)
        } 
        dev.off()
      } else{
        pdf(file=paste0(sample2Name,"_","scatterPlots.pdf"),width=pdfSize[1],height=pdfSize[2])
        par(mfrow=plotNumPerPage,mar = c(4.5, 7, 2.5, 2.5), mgp = c(3, 1, 0),pty="s")
        #create an empty plot window
        plot(plMin:plMax,plMin:plMax, type="n", yaxt='n', xaxt='n', ann=FALSE)
    
        #populate with points for each class
          for (i in 1:length(classes)){
            colorlist = c("tomato1","palegreen","deeppink4","goldenrod1","skyblue4","violetred","brown","darkorange4")
            points(s1[[i]],s2[[i]], pch=15,col=colorlist[i], cex=0.3)
            title(main=paste0("Enrichment of ",classes[i]," in ",sample2Name), xlab = paste0("Normalized reads in ",sample1Name),ylab = paste0("Normalized reads in ",sample2Name,"\n"),font.lab=2)
          }
        abline(h=log2(10), v=log2(10), lty=2, col="grey80")
        abline(0,1, col="grey80")
        plot_log2_ticks(plMin,plMax)
        text(paste0("Depleted in ",sample2Name))
        if (i == highlightclass){
          if (length(highlight)>=1){
            points(s1[[i]][highlight,1],s2[[i]][highlight,1], pch=20,col="deepskyblue", cex=0.3)
          }
        }
        }
      } else {
        if (plotSeparate){
          par(mfrow=plotNumPerPage,mar = c(4.5, 7, 2.5, 2.5), mgp = c(3, 1, 0),pty="s")
          #populate with points for each class
          for (i in 1:length(classes)){
            #create an empty plot window
            plot(plMin:plMax,plMin:plMax, type="n", yaxt='n', xaxt='n', ann=FALSE)
            points(s1[[i]][,1], s2[[i]][,1], pch=15,col="tomato1", cex=0.3)
            abline(h=log2(10), v=log2(10), lty=2, col="grey80")
            abline(0,1, col="grey80")
            abline(log2(0.5),1, col="grey80")
            abline(log2(2),1, col="grey80")
            title(main=paste0("Enrichment of ",classes[i]," in ",sample2Name), xlab = paste0("Normalized reads in ",sample1Name),ylab = paste0("Normalized reads in ",sample2Name,"\n"),font.lab=2)
            plot_log2_ticks(plMin,plMax)
            if (i == highlightclass){
              if (length(highlight)>=1){
                points(s1[[i]][highlight,1],s2[[i]][highlight,1], pch=20,col="deepskyblue", cex=0.3)
              }
            }
            
          } 
         } else{
          par(mfrow=plotNumPerPage,mar = c(4.5, 7, 2.5, 2.5), mgp = c(3, 1, 0),pty="s")
          #create an empty plot window
          plot(plMin:plMax,plMin:plMax, type="n", yaxt='n', xaxt='n', ann=FALSE)
          
          #populate with points for each class
          for (i in 1:length(classes)){
            colorlist = c("tomato1","palegreen","deeppink4","goldenrod1","skyblue4","violetred","brown","darkorange4")
            points(s1[[i]][,1],s2[[i]][,1], pch=15,col=colorlist[i], cex=0.3)
          }
          title(main=paste0("Enrichment of smRNAs in ",sample2Name), xlab = paste0("Normalized reads in ",sample1Name),ylab = paste0("Normalized reads in ",sample2Name,"\n"),font.lab=2)
          legend(plMax,plMax,classes,colorlist[1:length(classes)])
          abline(h=log2(10), v=log2(10), lty=2, col="grey80")
          abline(0,1, col="grey80")
          plot_log2_ticks(plMin,plMax)
          text(paste0("Depleted in ",sample2Name))
         }
      }
  } else if (length(sample1)!=length(sample2)){
    print("The number of files in sample1 is not equal to number of files in sample2")
  } else if (length(sample1)!=length(classes) | length(sample2)!=length(classes)){
      print("The number of small RNA classes does not match the number of filenames given")
    }
}

#creates tables with RPM for each small RNA class
make_smrnatables <- function(sampleFiles, sampleNames, className, pdfOpt=3){
  smRNAtable <- c()
  for (i in 1:length(sampleFiles)){
    temp <- read.table(file=sampleFiles[i],sep="\t",stringsAsFactors = FALSE)
    temp <- temp[!duplicated(temp[,1]),]
    rownames(temp) <- temp[,1]
    temp <- temp[,2,drop=F]
    if (length(smRNAtable) > 0){
      smRNAtable <- cbind(smRNAtable, temp)
    } else {smRNAtable <- temp}
  }
  colnames(smRNAtable) <- sampleNames
  if (pdfOpt%in%c(1,3)){
    write.table(smRNAtable, file=paste0("smRNA_RPM_", className,".txt"), col.names = T, row.names = T, quote=F, sep="\t")
  }
  if (pdfOpt%in%c(2,3)){
    return(smRNAtable)
  }
}


#barplots/piecharts for class enrichment

class_plots <- function(sample1, sample2, classes,sampleNames, pdfSize=c(12,8)){
  if (length(sample1)==length(sample2) & length(sample1)==length(classes)){
    sort(classes)
    #get variables (counts) from filenames
    for (i in 1:length(sample1)){
      temp <- read.table(file=sample1[i],sep="\t",stringsAsFactors = FALSE)
      temp <- temp[!duplicated(temp[,1]),]
      rownames(temp) <- temp[,1]
      temp <- temp[,2]
      
      temp2 <- read.table(file=sample2[i],sep="\t",stringsAsFactors = FALSE)
      temp2 <- temp2[!duplicated(temp2[,1]),]
      rownames(temp2) <- temp2[,1]
      temp2 <- temp2[,2]
      
      assign(paste0("sample1_",classes[i]), temp)
      assign(paste0("sample2_",classes[i]), temp2)
    }
  } else if (length(sample1)!=length(classes) | length(sample2)!=length(classes)){
    print("The number of small RNA classes does not match the number of filenames given")
  }
    classvalues_1 <- mget(ls(pattern="sample1_"))
    classvalues_2 <- mget(ls(pattern="sample2_"))
    a_1 <- sapply(classvalues_1, sum)
    a_2 <- sapply(classvalues_2, sum)
    other_1 <- 1e6 - sum(a_1)
    other_2 <- 1e6 - sum(a_2)
    
    ratios <- c(a_1/a_2,other_1/other_2)
    for (i in 1:length(classes)+1){
      if (ratios[i]<1){
        ratios[i]<- -1/ratios[i]
      } 
    }
    
    pdf(file=paste0(gsub(":","",sampleNames[1]),"_","classPlots.pdf"),width=pdfSize[1],height=pdfSize[2])
    par(mfrow=c(1,3))
    barplot(ratios, horiz=T, names.arg = c(classes,"other"), axis.lty = 1, xlim=c(min(ratios)-1,max(ratios)+1),main=paste0(sampleNames[1]," vs ",sampleNames[2]))
    pie(c(a_1,other_1), main=sampleNames[1],labels=c(classes,"other"))
    pie(c(a_2,other_2), main=sampleNames[2],labels=c(classes,"other"))
    dev.off()
}

#a barplot generator of individually enriched or depleted smRNAs
smRNA_bars <- function(sample1, sample2, sampleNames, depletion=F, threshold=2, rpm=10, col=c("conservedcombined","conservedbyspecies","default"),choices=sample(colors(),3)){
  s1 <- read.table(sample1,sep="\t",stringsAsFactors = FALSE)
  s2 <- read.table(sample2,sep="\t",stringsAsFactors = FALSE)
  s1 <- s1[order(s1[,1]),]
  s2 <- s2[order(s2[,1]),]
  ratios <- data.frame(s1,s2)
  
  if (identical(ratios[,1],ratios[,3])){
    ratios <- ratios[!duplicated(ratios[,1]),]
    rownames(ratios) <- ratios[,1]
    ratios <- ratios[,c(2,4)]
    ratios[,3]<-ratios[,1]/ratios[,2]
  } else {print("smRNAs are missing or are mismatched from one sample")}
  
  ratios <- ratios[ratios[,1]>=10|ratios[,2]>=10,]
  ratios <- ratios[!is.nan(ratios[,3]) & !is.infinite(ratios[,3]) & !is.na(ratios[,3]),]
  
  if (depletion){
    for (i in 1:nrow(ratios)){
      if (ratios[i,3]<1){
        ratios[i,3]<- -1/ratios[i,3]
      }
    }
    ratiosT <- ratios[ratios[,3] < threshold,]
  } else {
    ratiosT <- ratios[ratios[,3] > threshold,]
  }
  if (col=="conservedbyspecies"){
    hsa <- read.table("~/Desktop/MontgomeryLab/mature_seed_hsa.txt", sep="\t", stringsAsFactors = F)
    dme <- read.table("~/Desktop/MontgomeryLab/mature_seed_dme.txt", sep="\t", stringsAsFactors = F)
    cel <- read.table("~/Desktop/MontgomeryLab/mirnas_v21_seeds.txt", sep="\t", stringsAsFactors = F)
    leg <- c("conserved in humans","conserved in drosophila", "not conserved")
    cel_dme <- cel[cel[,2]%in%dme[,2],]
    cel_hsa <- cel[cel[,2]%in%hsa[,2],]
    colors <- c()
    for (i in 1:nrow(ratiosT)){
      if (rownames(ratiosT[i,])%in%cel_hsa[,1]){
        colors <- c(colors, choices[1])
      } else if (rownames(ratiosT[i,])%in%cel_dme[,1]){
        colors <- c(colors, choices[2])
      } else {colors <- c(colors, choices[3])}
    } 
  } else if (col=="conservedcombined"){
    hsa <- read.table("~/Desktop/MontgomeryLab/mature_seed_hsa.txt", sep="\t", stringsAsFactors = F)
    dme <- read.table("~/Desktop/MontgomeryLab/mature_seed_dme.txt", sep="\t", stringsAsFactors = F)
    cel <- read.table("~/Desktop/MontgomeryLab/mirnas_v21_seeds.txt", sep="\t", stringsAsFactors = F)
    leg <- c("conserved", "not conserved")
    cel_dme <- cel[cel[,2]%in%dme[,2],]
    cel_hsa <- cel[cel[,2]%in%hsa[,2],]
    colors <- c()
    for (i in 1:nrow(ratiosT)){
      if (rownames(ratiosT[i,])%in%cel_hsa[,1] | rownames(ratiosT[i,])%in%cel_dme[,1]){
        colors <- c(colors, choices[1])
      } else {colors <- c(colors, choices[2])}
    } 
  }
  barplot(ratiosT[order(ratiosT[,3]),3], horiz=T, names.arg =gsub("cel-","",rownames(ratiosT[order(ratiosT[,3]),])), xlim=c(0,max(ratiosT[,3])),main=paste0(sampleNames[1]," vs ",sampleNames[2]), border=NA, las=1, cex.names = 0.5, col = colors[order(ratiosT[,3])])
  legend(150, 8, leg ,fill=choices[1:length(leg)], bty="n")
  print(choices)
}