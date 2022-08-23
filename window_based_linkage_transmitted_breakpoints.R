###################################################################################################################################################
#
# Title: Window based linkage 
# Description: Find pairs of sequences with potentially transmitted recombinants 
#
# Input files: (1) alignment to test (either positive controls or 164 A1/D recombinants )
#             (2) Genome_breakpoints from other script with subtype information on the above alignment 
#              
# Notes 
#### Positive controls 
#### Script modified from Donald Kwan -->NetworkByWndow.R 
###################################################################################################################################################
#find # Given a set of aligned [sequences]
# - subset the string by moving [window] and [shift]
# - for each window 
# - - calculate pairwise distance
# - - filter edges by [threshold] and [identity]
# - - return edgelist
# - combine the edgelist as one network by [method]
# 

refs<-read.dna("data/ALN/pol2011_fullG_all_HXB2coord_EDITED_TO_ALIGN_WITH_62_HISTORICAL_A_D_only.fasta", format="fasta")
rownames(refs)
refs<-read.dna("data/ALN/HARRIS_DOWLING_ALIGNED_8292.fasta", format='interleaved')
rownames(refs) # DOWLING_AF457051   DOWLING_AF457072      DOWLING_AF457060
refs<-refs[c("DOWLING_AF457060", "DOWLING_AF457072", "DOWLING_AF457051"),]
write.dna(refs, "data/ALN/DOWLINGCRF.fasta")
library(ape)
library(igraph)
library(rlang)
#outputfn.prefix <- "C:\\Users\\CEID\\Desktop\\NCC\\seq\\Routput\\byWindow\\byWindow_"

#outputfn.prefix <- "/Volumes/cmvm-1/RECOMBINATION/Network/"
outputfn.prefix <- getwd()
#Genome_breakpoints <- readRDS("~/Desktop/HIV_recombination_in_R/data/saved_objects/Genome_breakpoints_pdfs_duplicates.RDS") 

#write.dna(al164.nfl, "al164.nfl.fasta",format='fasta')
al.hiv <- read.dna("data/ALN/Best_alignment_22March.fasta", format='fasta')


### edit to rec-1 new names 
#https://rdrr.io/cran/phylotools/man/rename.fasta.html
#don't need this twice 
library(phylotools)
old_name <- get.fasta.name("./data/al164.nfl.fasta")
Genome_breakpoints<-readRDS("data/saved_objects/Genome_breakpoints_164_with_psudeo_labs.RDS")
ref2 <- Genome_breakpoints[,c(1,9)]
rename.fasta(infile = "./data/al164.nfl.fasta", ref_table = ref2, outfile = "data/ALN/164renamed.fasta")
#done 

Genome_breakpoints<-readRDS("data/saved_objects/Genome_breakpoints_164_with_psudeo_labs.RDS")
fasta.fn<-"data/ALN/164renamed.fasta"
#then rename in df
#Genome_breakpoints<-Genome_breakpoints_simulated
Genome_breakpoints$Name<-Genome_breakpoints$psudo_lab

windowsize <- 300
threshold <- 0.02
#threshold <- 0.02
input.fasta <- read.dna(fasta.fn, format = "fasta", as.character=T, as.matrix=T)
#input.fasta<-al.nfl
numSeq <- nrow(input.fasta)
all.edgelist <- data.frame(check.names=F)
window <- 1
numWindows <- seq(1,length(input.fasta[1,])%/%windowsize*windowsize ,by=windowsize)

join.edgelist<-data.frame(Source = 'blank', Target ='blank', Dist="blank", Window='blank',check.names = F)
#layout(matrix(1:30, 6,5,byrow=TRUE))
for (i in numWindows) {
  print(paste(Sys.time(), " Processing: ", as.character(window), "/" , as.character(length(numWindows)), sep=""))
  subseq <- as.DNAbin(input.fasta[,seq(i,i+windowsize-1)])
  pairwise.dist <- dist.dna(subseq, model = "TN93", variance = FALSE,
                            gamma = FALSE, pairwise.deletion = TRUE,
                            base.freq = NULL, as.matrix = TRUE)
  edgelist <- data.frame(check.names=F)
  lowertri <- lower.tri(pairwise.dist) #matrix of logicals 
  #hist(pairwise.dist, main=paste("window", window)) #shows the distribution of pairwise distances 
  #abline(v=mean(pairwise.dist), col="purple")
  
  for (j in which(pairwise.dist < threshold)) {
    if (lowertri[j]) {
      if (j %% numSeq == 0) { # Don't do this.. keep going up until you reached num sequences = j  
        edgelist<- rbind(edgelist, c(rownames(pairwise.dist)[j %/% numSeq], rownames(pairwise.dist)[numSeq], pairwise.dist[j]),stringsAsFactors=F)
      } else { #CREATE THE EDGE LIST IN THIS PART 
        edgelist<- rbind(edgelist, c(rownames(pairwise.dist)[j %/% numSeq+1], rownames(pairwise.dist)[j %% numSeq], pairwise.dist[j]),stringsAsFactors=F)
      }
    }
  }
  if(nrow(edgelist)!=0){
  edgelist$Window <- window
  colnames(edgelist) <- c("Source", "Target", "Dist", "Window")
  #write.csv(edgelist, paste(outputfn.prefix, "Window", as.character(window), "_edgelist.csv",sep=""))
  all.edgelist <- rbind(all.edgelist, edgelist)
  }
  #join.edgelist<- left_join(join.edgelist, edgelist, by=c("Source" = "Source", "Target" = "Target"))
  window <- window + 1
  
}

allg <- graph_from_data_frame(all.edgelist, directed=F)
plot(simplify(allg), edge.width=count_multiple(allg)) #all the windows combined? 
simg <- simplify(allg, edge.attr.comb=list(Window=function(x) length(x), Dist=function(x) mean(as.double(x))))
plot(simg, edge.width=E(simg)$Window, layout=layout_with_fr, vertex.size=4)
plot(simg, edge.width=E(simg)$Window, layout=layout_with_fr, vertex.size=4, labels=FALSE)
write.csv(all.edgelist, paste(outputfn.prefix, "edgelist5pc.csv", sep=""))



################################################################################
library(tidyr)
merged<-unite(all.edgelist, "pair", c("Source", "Target"), sep = ",", remove = TRUE)
merged<-merged[-2]
all.edgelist$Dist<-as.numeric(all.edgelist$Dist)
plot(tapply(all.edgelist$Dist, all.edgelist$Window, mean), main='Average pairwise distance by window, for SUBTYPE D', xlab='window', ylab='pdist')
plot(all.edgelist$Dist, all.edgelist$Window, main='Average pairwise distance by window, for CRF45', xlab='window', ylab='pdist')


tapply(all.edgelist$Dist, all.edgelist$Window, mean)
tail(merged)
#df4<-spread(merged, key="Window", value="Window", fill="")
df4<-spread(merged, key="Window", value="Window", fill=NA)

length(unique(df4$pair)) #101 unique pairs #205 with 2pc #9072 with 5% 
#melted<-melt(data = all.edgelist, id.vars = "id", measure.vars = c("blue", "red"))
tail(df4)
#remove rows with only one window 
dim(df4)
df4$count_na<-rowSums(is.na(df4[,-1]))
df4$windowcount<-27-df4$count_na
#df4$windowcount<-31-df4$count_na
df4$windowcount<-ncol(df4)-2-df4$count_na

layout(1)
hist(as.numeric(df4$windowcount), breaks=27, xlim=c(0,27),col='pink')
unique(df4$windowcount)
#dim(df4[df4$windowcount==1,] )#single windows account for 90 #190 at 2pc

head(df4)

par(mfrow=c(9,1))
#par(mai=c(1, 2.5, 1, 0.2)) ###CAN SEE LABELS NOW

test<-df4[df4$windowcount>=2,] 
write.csv( test, "test.csv")
test<-test[order(-test$windowcount),] #re-order them to put most boxes first 


dim(test)
#test<-test[1:9,]
crfs<-test[c(10,11,14),]

#test<-test[-c(10,11,14),]

order<-c(5, 1,2 ,6, 3,4,8,7,9,10,11,12,13,14, 15)
#order<-c(10, 1,2 ,3,11,8,12 , 15,14,17,18,19,5,6,9,13,4,7,16,20) #old order with 389 
test<-test[order,]
  
layout(matrix(1:15,5,3, byrow = TRUE))
seq1<-c()
seq2<-c()
c<-1
Genome_breakpoints$Name<-Genome_breakpoints$psudo_lab
for(t in 1:nrow(test)){
  name1<-strsplit(test[t,]$pair,",")[[1]][1]
  pair1<-Genome_breakpoints[Genome_breakpoints$Name==name1,]
  name2<-strsplit(test[t,]$pair,",")[[1]][2]
  pair2<-Genome_breakpoints[Genome_breakpoints$Name==name2,]
  recombs<-rbind(pair1,pair2)
  #draw
  seq1<-c(seq1, name1)
  seq2<-c(seq2, name2)
#  plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("num windows matching=", test[t,]$windowcount),
#       main="",ylab="",axes=F,cex.lab=1.5)
  #plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("seq1:", pair2$SubtypesDetailed, "\n" ,"seq2:", pair1$SubtypesDetailed),
  #     main="",ylab="",axes=F,cex.lab=1)
  #plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("seq1:", pair2$Name, "\n" ,"seq2:", pair1$Name),
  #     main="",ylab="",axes=F,cex.lab=1)
  
  plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste(pair2$Name, "\n" , pair1$Name),
       main="",ylab="",axes=F,cex.lab=1)
  Axis(side=1,labels=T,cex.axis=1)
  
  if(t=='1' | t=='4' | t=='7'){
    
  }else{
    
    text(200, -1.7, paste("Pair", c, sep=' '), cex=1.5)
    c<-c+1
  }
  
  for(i in 1:dim(recombs)[1]){
    BG1<- recombs[i,]
    subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
    brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
    if(length(brpnts)==0){ #this bit is for pure genomes without breakpoitns 
      brpnts<-c('0', '8287')
    } else{
    brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
    brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
    }
    for(y in 1:length(subs)){
      if(subs[y]=="A1"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=1, col='dodgerblue3')
      }
      if (subs[y]=="D"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=1, col='darkgoldenrod1')
      }
    }
  }
  #text(4000, 3, paste(name1, name2, sep="     and     "))
  #text(4000, 6, paste("windows mathcing=", test[t,]$windowcount))
  wins<-test[t,2:28]
  #wins<-prepend(wins,0)
  for(w in 1:length(wins)){
    if(!is.na(wins[w])){
      #rect(as.numeric(wins[w])*300, i-1.15, as.numeric(wins[w+1])*300, i+1, lwd=1, col='#FF003322')
      rect(w*300-300, i-1.15, w*300, i+1.15, lwd=1, border="dimgray")#, col='#FF003322')
      
    }
  }
}
#vesion wihtout the weird box thing 
for(t in 1:nrow(test)){
  name1<-strsplit(test[t,]$pair,",")[[1]][1]
  pair1<-Genome_breakpoints[Genome_breakpoints$Name==name1,]
  name2<-strsplit(test[t,]$pair,",")[[1]][2]
  pair2<-Genome_breakpoints[Genome_breakpoints$Name==name2,]
  recombs<-rbind(pair1,pair2)
  #draw
  seq1<-c(seq1, name1)
  seq2<-c(seq2, name2)
  #  plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("num windows matching=", test[t,]$windowcount),
  #       main="",ylab="",axes=F,cex.lab=1.5)
  #plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("seq1:", pair2$SubtypesDetailed, "\n" ,"seq2:", pair1$SubtypesDetailed),
  #     main="",ylab="",axes=F,cex.lab=1)
  #plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste("seq1:", pair2$Name, "\n" ,"seq2:", pair1$Name),
  #     main="",ylab="",axes=F,cex.lab=1)
  
  plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),xlab=paste(pair2$Name, "\n" , pair1$Name),
       main="",ylab="",axes=F,cex.lab=1)
  Axis(side=1,labels=T,cex.axis=1)
  
  #if(t=='1' | t=='4' | t=='7'){
  #  
  #}else{
  #  text(200, -1.7, paste("Pair", c, sep=' '), cex=1.5)
  #  
  #  c<-c+1
  #}
  
  for(i in 1:dim(recombs)[1]){
    BG1<- recombs[i,]
    subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
    brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
    if(length(brpnts)==0){ #this bit is for pure genomes without breakpoitns 
      brpnts<-c('0', '8287')
    } else{
      brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
      brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
    }
    for(y in 1:length(subs)){
      if(subs[y]=="A1"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=1, col='dodgerblue3')
      }
      if (subs[y]=="D"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=1, col='darkgoldenrod1')
      }
    }
  }
  #text(4000, 3, paste(name1, name2, sep="     and     "))
  #text(4000, 6, paste("windows mathcing=", test[t,]$windowcount))
  wins<-test[t,2:28]
  #wins<-prepend(wins,0)
  for(w in 1:length(wins)){
    if(!is.na(wins[w])){
      #rect(as.numeric(wins[w])*300, i-1.15, as.numeric(wins[w+1])*300, i+1, lwd=1, col='#FF003322')
      rect(w*300-300, i-1.15, w*300, i+1.15, lwd=1, border="dimgray")#, col='#FF003322')
      
    }
  }
}


Genome_breakpoints <- readRDS("data/saved_objects/Genome_breakpoints_pdfs_duplicates.RDS")

dim(test) #15 pairs 
removal<-c()
for(p in 1:length(test$pair)){
  remove<-strsplit(test$pair[p], ',')[[1]][[1]]
  removal<-c(removal, remove)
}
removal<-unique(removal) #13 genomes 
dim(Genome_breakpoints[Genome_breakpoints$Name %in% removal,]) #yes 13 
Genome_breakpoints <- Genome_breakpoints[!Genome_breakpoints$Name %in% removal, ]
dim(Genome_breakpoints) #454 + 12 = 466 (-1 same one twice in CRF)
saveRDS(Genome_breakpoints, "data/saved_objects/Genome_breakpoints_no_dup_no_transmitted_2_27_2pc.RDS")
