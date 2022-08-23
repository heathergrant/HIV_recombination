
###################################################################################################################################################
#
# Title: SCUEAL > Genome_breakpoints 
# Description: Raw SCUEAL files > useable R recombination object 
#
# Input files: (1) SCUEAL raw output
#              (2) Full genome alignment 
# Notes 
#previous versions (not git versioned) include recombinatino and recombination2.R (january 2019) and (subtype editing.R) 
# cleanUpScuealRawFile() written by Donald Kwan 
###################################################################################################################################################

library(seqinr)
library(tidyverse)
library(ape)
library(dplyr)
require(ggplot2)

cleanUpScuealRawFile <- function(scueal.output.fn) {
  SCUEAL <- read.table(scueal.output.fn, header=T, as.is=T, row.names=1, sep='\t')
  SCUEAL <- SCUEAL[, c(1,2,3,8)]
  SCUEAL$Name <- trimws(as.character(SCUEAL$Name), which="both")
  SCUEAL$Simplified.Subtype <- gsub(" intra-subtype recombinant", "", SCUEAL$Simplified.Subtype, fixed=T)
  SCUEAL$Simplified.Subtype <- gsub(" inter-subtype recombinant", "", SCUEAL$Simplified.Subtype, fixed=T)
  SCUEAL$Simplified.Subtype <- trimws(SCUEAL$Simplified.Subtype, which="both")
  SCUEAL$Simplified.Subtype <- gsub(" recombinant","", SCUEAL$Simplified.Subtype, fixed=T)
  
  SCUEAL$Subtype <- gsub(" intra-subtype recombinant", "", SCUEAL$Subtype, fixed=T)
  SCUEAL$Subtype <- gsub(" inter-subtype recombinant", "", SCUEAL$Subtype, fixed=T)
  SCUEAL$Subtype <- trimws(SCUEAL$Subtype, which="both") #trimws = trim white space - leading or trailing (in our case, 'both')
  SCUEAL$Subtype <- gsub(" recombinant","", SCUEAL$Subtype, fixed=T)
  
  SCUEAL$Breakpoints <- gsub("[(][0-9]+[-][0-9]+[)]", "", SCUEAL$Breakpoints)

  return(SCUEAL)
}

S1277<-cleanUpScuealRawFile("data/raw_SCUEAL/1277_new_alignment.fasta.out")
S1277$PDF<-rownames(S1277)
S603<-cleanUpScuealRawFile("data/raw_SCUEAL/603_for_subtyping.fasta.out")
S603$PDF<-as.numeric(rownames(S603))+as.numeric(1277)
S<-rbind(S1277, S603)
#S$PDF<-rownames(S) #so that numbers actually match to pdfs now. 
head(S)
###################################################################################################################################################
#
# Add Breakpoint number 
seqs <- read.fasta("data/ALN/Best_alignment_22March.fasta")
subty<-S
numBP <- c()
for (i in 1:nrow(subty)){
  nam <- subty$Name[i]
  seq1 <- seqs[[which(attr(seqs,"name")==nam)]]
  numBP <- c(numBP, length(which(seq1!="-")))
  numBP<-as.numeric(numBP)
}
S<-cbind(S, numBP) 
#Filter out only NFLG over 8000bp 
S<-S[which(S$numBP>8000),]
###################################################################################################################################################
  
###Loop which agressively edits the subtype
## And finds gaps over 100 bp
###And find the new breakpoints adjusted for gaps. 

subty<-S
newBrks <- c()
newSimpSub <- c()
gapLocs <- c()
subty2 <- subty
numBrkPt<-c()
#A.As<-c()
#D.As<-c()
#D.Ds<-c()
#A.Ds<-c()
#Inter.list<-c()
simpSub<-c()
#find_intra<-'no' #switches off some code 
for (i in 1:nrow(subty)){
  nam <- subty$Name[i]
  subt <- strsplit(subty$Subtype[i],",")[[1]]
  brk <- as.numeric(strsplit(subty$Breakpoints[i],";")[[1]])

  #replace any B (only B, not B/C/D...etc) with D
  subt[which(subt=="B")] <- "D"  #important to do this before erasing intra-SCUEAL.Subtype recomb!!
  subt[which(subt=="B/D" | subt=="D/B" | subt=="BF" | subt=="BG")] <- "D"
  
  #replace 01, 02, and A3 with A1 also A2 
  subt[which(subt=="01" | subt=="02" | subt=="A3"| subt=="A2")] <- "A1"
  #subt[which(subt=="01" | subt=="A3"| subt=="A2")] <- "A1"
  
  #replace the recombinants of these three with A1
  comp <- which(grepl("/",subt,fixed=T))
  for(kj in comp){
    tx <- strsplit(subt[kj], "/")[[1]]
    #tx <- tx[grep("A3|02|01", tx, invert=T)] #changed january 2019 to add in A1/A2 that had been missed 
    tx <- tx[grep("A3|A2|A1|02|01", tx, invert=T)]
    if(length(tx)==0){ #the complex only contained A3, 02, and 01
      tx <- "A1"
    }
    if(all(tx=="A1")){ #if after removing A3, 02, 01, all thats left is A1..
      tx <- "A1"
    }	
    subt[kj] <- paste(tx, collapse="/")
  }
  
  
  #subt <- gsub("[[:blank:]][(][0-9]+[[:blank:]][[:alpha:]]+[)]", "", subt)
  #intra<-which(grep("breakpoints",subt,fixed=T))
  ###added below march to turn A1 (2 breakpoints) into A1, A1, A1 
  # removed to get rid of 
  #if(find_intra=="TRUE"){
  if(length(subt)==1){
    if(grepl('breakpoints', subt)){
    subt<-gsub("A1", "A", subt)
    insidebracket<-str_extract(subt, "\\-*\\d+\\.*\\d*")
    subt<-gsub("A", "A1", subt)
    withoutbracket<-gsub("[[:blank:]][(][0-9]+[[:blank:]][[:alpha:]]+[)]", "", subt)
    subt<-rep(withoutbracket, as.numeric(insidebracket)+1)
    }
  }
  #}
  #brk <- numeric(0) #there are no breakpoints
  #  }
  
  #In sequences with inter-SCUEAL.Subtype recomb,
  #get rid of intra-SCUEAL.Subtype recomb (compress SCUEAL.Breakpoints so only beteween SCUEAL.Subtype change, not within)
  #those with only inter-SCUEAL.Subtype recomb will retain their SCUEAL.Breakpoints
  if(length(subt)>1){
    keep <- 1
    for(j in 2:length(subt)){
      if(subt[j]!=subt[j-1]){#if the subtype is not the same as the one before it 
        keep <- c(keep, j) #add the subtype to the keep list 
      }
    }#keep is the NUMBER of subtypes to keep.e.g. h D A1 C D makes 1 2, 3, 4, 5
    keepBr <- keep-1 # this makes 0,1,2,3,4 (sequentially)
    subt <- subt[keep]#index of the subtypes to keep
    #brk <- c(brk[keepBr],100000)
    brk <- brk[keepBr]#index of breakpoints to keep
  } #else {
  #it's a pure subtype - but want to get rid of "(1 breakpoints)" etc
  #subt
  #brk
  
  cleanSub <- unique(subt)
  cleanSub <- cleanSub[order(cleanSub)]
  simpSub<-c(simpSub, paste(cleanSub, collapse=','))
  totSubs<-c() #added?
  totSubs <- c(totSubs, paste(cleanSub,collapse=","))
  numBrkPt <- c(numBrkPt, length(brk))
  subty2[i,"Subtype"] <- paste(subt,collapse=",")
  #adjust SCUEAL.Breakpoints
  seq1 <- seqs[[which(attr(seqs,"name")==nam)]]
  
  newBrk <- c() #new breakpoints empty list 
  if(length(brk)!=0){
    #newBrk <- c() #new breakpoints empty list 
    for(br in brk){
      gaps <- sum(length(which(seq1[1:br]=="-")) + length(which(seq1[1:br]=="n")))#gaps before EACH breakpoint 
      
    #}
  #}
  #keep track of gaps
    findG <- gregexpr("[-]+", paste(seq1, collapse=""))[[1]]
    gps <- c()
    for(j in 1:length(findG)){
        strt <- findG[j]
        len <- attr(findG, "match.length")[j]
        if(len > 100){
        gps <- c(gps, paste(c(strt,"-",(strt+len-1)), collapse=""))
          }
         }
  newBrk <- c(newBrk, (br+gaps))
  #numBrkPt<-c() #added? 
  }
    }else{
    gps<-""
    gaps<-0
    #newBrk<-""
  }
  
  gapLocs <- c(gapLocs, paste(gps, collapse=";"))
  subty2[i,"Breakpoints"] <- paste(brk,collapse=";")
  newBrks <- c(newBrks, paste(newBrk, collapse=";"))
  
  #}else{
  #  subty2[i,"Breakpoints"] <- ""
  #}
  
  #our own simplified subtype rules
  #if (length(cleanSub) > 3){
  #  newSimpSub <- c(newSimpSub, "complex")
  #} else {
  #  if( length(grep("/", cleanSub)) > 0){
  #  } else {
  #    newSimpSub <- c(newSimpSub, "complex")
  #    newSimpSub <- c(newSimpSub, paste(cleanSub, collapse=","))
  #  }
  #}
  #A.A<-c()
  #A.D<-c()
  #D.A<-c()
  #D.D<-c()
  #Inter<-c()
  #if(length(newBrk)!=0){
  #if(length(newBrk)!=0){
  #  for(b in 1:length(newBrk)){
      #if(subt[b]=="A1" && subt[b+1]=="A1"){
      #  A.A<-c(A.A, newBrk[b])
      #}
      #if(subt[b]=="D" && subt[b+1]=="D"){
      #  D.D<-c(D.D, newBrk[b])
      #}
      #if(subt[b]=="A1" && subt[b+1]=="D"){
      #  A.D<-c(A.D, newBrk[b])
      #}
      #if(subt[b]=="D" && subt[b+1]=="A1"){
      #  D.A<-c(D.A, newBrk[b])
      #}
      #if(subt[b]!=subt[b+1]){
      #  Inter<-c(Inter, newBrk[b])
      #}
    #}
  #}
  #A.As<-c(A.As, paste(A.A, collapse=";"))
  #A.Ds<-c(A.Ds,paste(A.D, collapse=";"))
  #D.As<-c(D.As, paste(D.A, collapse=";"))
  #D.Ds<-c(D.Ds, paste(D.D, collapse=";"))
  #Inter.list<-c(Inter.list, paste(Inter, collapse = ";"))

}

#Genome_breakpoints <- as.data.frame(cbind(subty$Name, subty2$Subtype, simpSub, newBrks, A.As, A.Ds, D.As, D.Ds, Inter.list, numBrkPt, subty$numBP, gapLocs, subty$PDF))

Genome_breakpoints <- as.data.frame(cbind(subty$Name, subty2$Subtype, simpSub, newBrks, numBrkPt, subty$numBP, gapLocs, subty$PDF))

#colnames(Genome_breakpoints) <- c("Name", "Subtypes", "CleanSub","Breakpoints", "AA", "AD", "DA", "DD", "inter","num_breaks", "numBP", "gaps",
                                #"PDF")
colnames(Genome_breakpoints) <- c("Name", "Subtypes", "CleanSub","Breakpoints", "num_breaks", "numBP", "gaps",
                                  "PDF")
head(Genome_breakpoints)
tail(Genome_breakpoints)

Genome_breakpoints_all_lengths <- as.data.frame(cbind(subty$Name, subty2$Subtype, simpSub, newBrks, A.As, A.Ds, D.As, D.Ds, numBrkPt, subty$numBP, gapLocs, subty$PDF))

colnames(Genome_breakpoints_all_lengths) <- c("Name", "Subtypes", "CleanSub","Breakpoints", "AA", "AD", "DA", "DD", "num_breaks", "numBP", "gaps",
                                  "PDF")

Genome_breakpoints_all_lengths$numBP<-as.numeric(as.character(Genome_breakpoints_all_lengths$numBP))
hist(Genome_breakpoints_all_lengths$numBP)
max(Genome_breakpoints_all_lengths$numBP)
#############################
#Duplicate removal
duplicates_and_scueal <- read_csv("data/duplicates_and_scueal.csv")
duplicates<-duplicates_and_scueal$`Duplicate 2`
duplicates <-duplicates[1:23]#remove NAs at end 
str(duplicates)
dups<-duplicates#$Name
#Dup<-Genome_breakpoints[Genome_breakpoints$Name==(duplicates),]
Genome_breakpoints <- Genome_breakpoints[!Genome_breakpoints$Name %in% duplicates, ]
Genome_breakpoints_all_lengths <- Genome_breakpoints_all_lengths[!Genome_breakpoints_all_lengths$Name %in% duplicates, ]
#dups<- Genome_breakpoints[Genome_breakpoints$Name %in% duplicates, ] 5 duplicates were 1D, 3 A1, 1 A1/D. 
dim(Genome_breakpoints) #465 - removed the 5 duplicates that are over 8000base pairs. 
#saveRDS(Genome_breakpoints, "saved_objects/Genome_breakpoints.RDS")
saveRDS(Genome_breakpoints, "data/saved_objects/Genome_breakpoints_no_intra_no_dups_465.RDS")
saveRDS(Genome_breakpoints, "data/saved_objects/Genome_breakpoints_inc_inter_no_dupes.RDS")
saveRDS(Genome_breakpoints_all_lengths, "data/saved_objects/Genome_breakpoints_all_lengths.RDS")
#############################
