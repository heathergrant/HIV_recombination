###################################################################################################################################################
#
# Title: simulate intra recombination
# Description: 
#
# Input files: alignment / sequences with 0 breakpoints (as determined by SCUEAL)
#              
# Notes 
###################################################################################################################################################
al.dna<-read.dna("data/ALN/Best_alignment_22March.fasta", format='fasta')

rownames(al.dna)==""
#D1 PANGEA-H17240-UG2013-0 # D_UG_1998_505 match 
D1<-al.dna["PANGEA-H17240-UG2013-0",]
#D2  PANGEA-A21437-UG2007-0 # D_KE_2001 090  
D2<-al.dna["PANGEA-A21437-UG2007-0",]
#D3 PANGEA-F50477-UG2013-0 #D_YE 903
D3<-al.dna["PANGEA-F50477-UG2013-0",]
# then A1 KE 055 pdf 490 
A1<-al.dna["PANGEA-I31498-UG2013-0",]
#A1 KE 2000_ 089 pdf 854 
A2<-al.dna["PANGEA-O28713-UG2013-0",]
#then A1 ug 1998 pick pdf 33 
A3<-al.dna["PANGEA-A27274-UG2013-0",]

refs<-rbind.DNAbin(A1,A2, A3, D1, D2, D3)
short.ref.names<-c()
for(i in 1:3){
  short.ref.names<-c(short.ref.names, paste0("A", i))
}
for(i in 4:6){
  short.ref.names<-c(short.ref.names, paste0("D", i))
}

rownames(refs)<-short.ref.names
plot(nj(dist.dna(refs)))
make_intRA_recombs <- function(num.recombs) {
  #Genome_breakpoints_simulated<-data.frame(Name='', Subtypes='', Breakpoints='', SubtypesDetailed='')
  #if(whichsubtype=="A"){
  #  refs<-refs[1:3,]
  #}
  #if(whichsubtype=="D"){
  #  refs<-refs[4:6,]
  #}
  #refs<-refs[1:3,]# #subtype A 
  refs<-refs[4:6,]# subtype D
  simulated.DNA<-refs[1,]
  text.records<-data.frame()
  #recombs.DNA<-c()
  #recombs.records<-c()
  break.list<-seq(300,7800,300) #consider changing the minimum to 300 bp... somehow? 
  br.num<-seq(1,2,1) #consider ecdf() emprical distribution sampling - the actual breakpoint freqs. 
  #ref.list<-rbind.DNAbin(A1.1, A1.2, A1.3, D.1, D.2, D.3)
  #genomes<-sample(refs,1)
  for(i in 1:num.recombs){
    rando.breaks<-sample(br.num, 1)
    rando.locations.base<-sort(sample(break.list, rando.breaks))
    rando.locations<-c(1, rando.locations.base, (dim(refs)[2]+1))#length align +1 
    #random.seq<-refs[sample(1:nrow(refs), rando.breaks+1),] # 4 sequences with which to recombine I think? 
    #random.seqD<-refs[sample(1:3, 1),] #pick random D
    #random.seqA<-refs[sample(4:6, 1),] #pick random A 
    #choices<-c("A", "D")
    #firstchoice<-sample(choices, 1)
    #if(firstchoice=="A"){
    #  secondchoice<-"D"
    #}
    #if(firstchoice=="D"){
    #  secondchoice<-"A"
    #}
    #order<-rep(c(firstchoice, secondchoice), 10)

    random.seq<-refs[sample(1:nrow(refs), rando.breaks+1),]

    
    
    #create record for it 
    Name<-paste0("sim_recomb", i)
    Breakpoints<-paste0(rando.locations.base, collapse=";")
    SubtypesDetailed<-paste0(rownames(random.seq), collapse=',')
    Subtypes<-gsub('[0-9]+', '', SubtypesDetailed)
    text.record<-cbind(Name, Subtypes, Breakpoints, SubtypesDetailed)
    text.record<-data.frame(text.record)
    text.records<-rbind(text.records, text.record)
    #DNA
    seed<-refs[1,1]
    #Record 
    rownames(seed)<-Name
    record.dna<-seed
    record.dna<-record.dna[,-1] #blank dna record with name simulated_recomb1 etc 
    for(n in 1:(length(rando.locations)-1)){
      seqs<-rownames(random.seq)
      chunk<-random.seq[n,(rando.locations[n]):(rando.locations[n+1]-1)]
      rownames(chunk)<-Name
      record.dna<-cbind.DNAbin(record.dna,chunk)
    }
    simulated.DNA<-rbind.DNAbin(simulated.DNA, record.dna)
    #Breakpoints<-strsplit(Breakpoints,';')[[1]][-1]
    #Breakpoints<-Breakpoints[-length(Breakpoints)]
  }
  simulated.DNA<-simulated.DNA[-1,]#remove the starter I needed to create the dna bin object 
  #write.dna(simulated.DNA, "simulated.DNA.fasta")
  #return(simulated.DNA)
  #write.csv(text.records, "text.records.csv")
  #return(Genome_breakpoints_simulated)
  Genome_breakpoints_INTRA_simulated_D<<-text.records
  DNA.simulated_INTRA_D<<-simulated.DNA
}
make_intRA_recombs(10) 
Genome_breakpoints_INTRA_simulated_A
saveRDS(Genome_breakpoints_INTRA_simulated_A, "Genome_breakpoints_INTRA_simulated_A.R")
write.dna(file = "DNA.simulated_INTRA_A.fasta", DNA.simulated_INTRA_A, format='fasta')

Genome_breakpoints_INTRA_simulated_D
saveRDS(Genome_breakpoints_INTRA_simulated_D, "Genome_breakpoints_INTRA_simulated_D.R")
write.dna(file = "DNA.simulated_INTRA_D.fasta", DNA.simulated_INTRA_D, format='fasta')

#make 10 fasta files containing each 100 replicates of 10 simulated recombinants 
for(i in 1:nrow(DNA.simulated_INTRA_D)){
  seed<-DNA.simulated_INTRA_D[i,]
  rownames(seed)<-paste(rownames(seed), ".", "1", sep='')
  recombi<-seed
  for(r in 2:100){
    seed<-DNA.simulated_INTRA_D[i,]
    rownames(seed)<-paste(rownames(seed), ".", r, sep='')
    recombi<-rbind.DNAbin(recombi, seed)
  }
  write.dna(recombi, paste("simulated_intra_D/Inter_sim_recomb", i, ".fasta", sep=''), format='fasta')
  
}


#path='data/raw_SCUEAL/sim/stoppingcritera100/'
path="data/raw_SCUEAL/sim/IntraD_10by100/"
Genome_breakpoints_simulated<-Genome_breakpoints_INTRA_simulated_D
#path="data/raw_SCUEAL/sim/IntraA_10by100/"
#Genome_breakpoints_simulated<-Genome_breakpoints_INTRA_simulated_A
folders<-list.files(path)
folders<-folders[-11]#remove plots folder
#recipe for a 4 bar strip with a lefthand side column of space.   
l<-matrix(1:20 , 4,5, byrow = F)

print(l)
layout(l) 
#plot.
id2<-c()
records2<-c()
timescolum<-c()
for(f in 1:length(folders)){
  
  
  inscueal<-folders[f][which(grepl('1', folders))]
  pathnew<-(paste0(path, f, "/"))
  folderfiles<-list.files(paste0(pathnew))
  report<-folderfiles[which(grepl("out", folderfiles))] #contining out - the scueal report 
  seqs<-folderfiles[-which(grepl("out", folderfiles))] #not containing out - the fasta file 
  truth<-Genome_breakpoints_simulated[f,]
  recombsscueal<-cleanUpScuealRawFile_CI(paste0(pathnew, report),paste0(pathnew, seqs))
  
  #par(mfrow=c(2,1))
  #plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+1),main="",ylab="",axes=F,cex.lab=1)
  plot(0,type='n',xlim=c(0,8287),ylim=c(1,nrow(recombs)+2),xlab="",main="",ylab="",axes=F,cex.lab=1)
  
  mtext(paste("Simulated A1 intra-subtype", f), cex=0.9)
  #Axis(side=1,labels=T,cex.axis=1)
  recombs<-truth
  for(i in 1:dim(recombs)[1]){
    BG1<- recombs[i,]
    subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
    
    brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
    brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
    brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
    for(y in 1:length(subs)){
      if(subs[y]=="A"| subs[y]=="35"){
        #rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='orange')
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='orange') #make them thicker?? 
        #rect(BG1[y,1],(i+0.1),BG1[y,2],i+0.9,lwd=0.2,col='firebrick1')
        recttext(brpnts[y],(i+0.1),brpnts[y+1],i+0.9, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][y], textArgs = list(col='black', cex=1))
      }
      if (subs[y]=="D"| subs[y]=="10"){
        #rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='dodgerblue3')
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='dodgerblue3')
        recttext(brpnts[y],(i+0.1),brpnts[y+1],i+0.9, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][y], textArgs = list(col='black', cex=1))
      }
    }
  }
  #go over first text so it shows -  
  recttext(brpnts[1],(i+0.1),brpnts[2],i+0.9, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][1], textArgs = list(col='black', cex=1))
  br.list.other<-c()
  for (i in 1:nrow(recombsscueal)){
    brk <- strsplit(as.character(recombsscueal$DD[i]),";")
    br.list.other<-c(br.list.other, brk)
  }
  br.list.other<-as.numeric(unlist(br.list.other))
  hist(as.numeric(unlist(br.list.other)), breaks = seq(from = 0, to = refs.length, by = 100), col='red', main='', xlab = "", ylab='Number', cex.lab=1.6, ylim=c(0,100))
  
  br.list.inter<-c()
  for (i in 1:nrow(recombsscueal)){
    brk <- strsplit(as.character(recombsscueal$inter[i]),";")
    br.list.inter<-c(br.list.inter, brk)
  }
  br.list.inter<-as.numeric(unlist(br.list.inter))
  hist(as.numeric(unlist(br.list.inter)), breaks = seq(from = 0, to = refs.length, by = 100), col='black', add=TRUE)
  #print(br.list.inter)
  #answers<-unique(recombsscueal$Subtypes)
  #for(a in 1:length(answers)){
  #  text(paste(answers[a], "found", sum(recombsscueal$Subtypes==answers[a]), "times"), x=2000, y=50+20*a)

  #}
  times<-c()
  id<-c()
  records<-c()
  cleananswers<-unique(recombsscueal$CleanSub)
  for(a in 1:length(cleananswers)){
    #text(paste(cleananswers[a], "called", sum(recombsscueal$CleanSub==cleananswers[a]), "times"), x=2000, y=50+20*a)
    times<-c(times, sum(recombsscueal$CleanSub==cleananswers[a]))
  }
  id<-rep(paste("Simulated D intra-subtype", f), length(cleananswers))
  id2<-c(id2, id)
  records<-cleananswers
  timescolum<-c(timescolum, times)
  records2<-c(records2, as.character(records))


}

records <- as.data.frame(cbind(id2, records2, timescolum))
colnames(records)<-c("Recombinant ID", "Unique answers", "Frequency")
write.csv(records, "records.csv")

#try again with the other ones 
#this makes the plots so you can see 
folders<-list.files(path)
folders<-folders[-11]#remove plots folder
for(f in 1:length(folders)){
  
  #png(filename=paste0(path, "plots/plot", f, ".png"), width = 2500, height = 2000, units = "px")
  pdf(paste0(path, "plots/plot", f, ".pdf"), width = 50, height = 40)
  layout(l) 
  inscueal<-folders[f][which(grepl('1', folders))]
  pathnew<-(paste0(path, f, "/"))
  folderfiles<-list.files(paste0(pathnew))
  report<-folderfiles[which(grepl("out", folderfiles))] #contining out - the scueal report 
  seqs<-folderfiles[-which(grepl("out", folderfiles))] #not containing out - the fasta file 
  truth<-Genome_breakpoints_simulated[f,]
  recombsscueal<-cleanUpScuealRawFile_CI(paste0(pathnew, report),paste0(pathnew, seqs))
  
  #par(mfrow=c(2,1))
  plot(0,type='n',xlim=c(0,8287),ylim=c(-1,nrow(recombs)),xlab=paste("In silico recombinant", f),main="",ylab="",axes=F,cex.lab=5)
  Axis(side=1,labels=T,cex.axis=2)
  recombs<-truth
  for(i in 1:dim(recombs)[1]){
    BG1<- recombs[i,]
    subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
    
    brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
    brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
    brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
    for(y in 1:length(subs)){
      if(subs[y]=="A"| subs[y]=="35"){
        #rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='orange')
        rect(brpnts[y],(0.1),brpnts[y+1],10,lwd=0.3, col='orange') #make them thicker?? 
        #rect(BG1[y,1],(i+0.1),BG1[y,2],i+0.9,lwd=0.2,col='firebrick1')
        recttext(brpnts[y],(0.1),brpnts[y+1],10, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][y], textArgs = list(col='black', cex=5))
      }
      if (subs[y]=="D"| subs[y]=="10"){
        #rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='dodgerblue3')
        rect(brpnts[y],(0.1),brpnts[y+1],10,lwd=0.3, col='dodgerblue3')
        recttext(brpnts[y],(0.1),brpnts[y+1],10, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][y], textArgs = list(col='black', cex=5))
      }
    }
  }
  #sim1<-hist(as.numeric(unlist(br.list.other)), breaks = seq(from = 0, to = refs.length, by = 100), col='black', main='', xlab = "")
  
  #)R 
  recombs<-recombsscueal
  plot(0,type='n',xlim=c(0,8287),ylim=c(0,nrow(recombs)+2),xlab="100 SCUEAL replicates",main="",ylab="",axes=F,cex.lab=5)
  Axis(side=1,labels=T,cex.axis=2)
  
  for(i in 1:dim(recombs)[1]){
    BG1<- recombs[i,]
    subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
    #
    brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
    if(is_empty(brpnts)){
      if(subs=="A1"){
        rect(0,(i+0.15),(8287),i+0.85,lwd=0.5, col='orange')
      }
      if(subs=="D"){
        rect(0,(i+0.15),8287,i+0.85,lwd=0.5, col='dodgerblue3')
      }
      
    }
    else{
      
      brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
      brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
    }
    for(y in 1:length(subs)){
      
      if(subs[y]=="A1"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=0.5, col='orange')
        #rect(BG1[y,1],(i+0.1),BG1[y,2],i+0.9,lwd=0.2,col='firebrick1')
        #rect(brpnts[y+1], (i+0.1), (brpnts[y+1]+0.1), i+0.9, lwd=0, col='white')
      }
      if (subs[y]=="D"){
        rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=0.5, col='dodgerblue3')
      }
      if (subs[y]=="19"){
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='purple')
        #recttext(brpnts[y],(i+0.1),brpnts[y+1],i+0.9, strsplit(as.character(recombs$SubtypesDetailed), ',')[[1]][y], textArgs = list(col='black', cex=1.5))
      }
      if (subs[y]=="10"){
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='lightblue')
      }   
      if (subs[y]=="07"| subs[y]=="04"| subs[y]=="15"| subs[y]=="32"| subs[y]=="13/27"| subs[y]=="33/34"| subs[y]=="16"|subs[y]=="18"|subs[y]=="11"|subs[y]=="23"|subs[y]=="15/11"|subs[y]=="19"|subs[y]=="36"){
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='pink')
      }
      if (subs[y]=="BG/G/N/O" | subs[y]=="B/BF/D/F1/F2/K" | subs[y]=="N/O" | subs[y]=="BG/G/N/O"| subs[y]=="04/06/09/11/13/14/15/16/18/20/22/23/24/25/27/32/33/34/35/36/37/43/BG/G"| subs[y]=="15/16/22/33/34/35/36/37" | subs[y]=="14/20/23/24/BG/G") {
        rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='grey')
      }
      
    }
  }
  
  plot(0,type='n',xlim=c(0,37),ylim=c(0,5),axes=F, xlab="", ylab="")
  rect(1,3,2,4,col='orange')
  rect(5,3,6,4,col='yellow')
  rect(9,3,10,4,col='dodgerblue3')       
  rect(13,3,14,4,col='lightblue')   
  rect(17,3,18,4,col='pink') 
  rect(21,3,22,4,col='grey')       
  
  text(3,3.5,"A1",cex=5)
  text(7.5,3.5,"CRF35",cex=5)
  text(10.5,3.5,"D",cex=5)
  text(15,3.5,"CRF10",cex=5)
  text(19,3.5,"other CRF",cex=5)
  text(23, 3.5, "complex", cex=5)
  
  dev.off()
  
}
