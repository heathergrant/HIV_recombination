### 
###################################################################################################################################################
#
# Title: 
# Description: INTRA WITHIN SCUEAL VALIDATION  
#
# Input files: (1) 
#              
# Notes 
#  
#######################################################################################################################
####
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
mat<-as.matrix(dist.dna(refs))
heatmap.2(mat,cellnote=mat)
heatmap.2(as.matrix(dist.dna(refs)),cellnote=as.matrix(dist.dna(refs)))


make_inter_recombs <- function(num.recombs) {
  #Genome_breakpoints_simulated<-data.frame(Name='', Subtypes='', Breakpoints='', SubtypesDetailed='')
  simulated.DNA<-refs[1,]
  text.records<-data.frame()
  #recombs.DNA<-c()
  #recombs.records<-c()
  break.list<-seq(300,7800,300) #consider changing the minimum to 300 bp... somehow? 
  br.num<-seq(1,3,1) #consider ecdf() emprical distribution sampling - the actual breakpoint freqs. 
  #ref.list<-rbind.DNAbin(A1.1, A1.2, A1.3, D.1, D.2, D.3)
  #genomes<-sample(refs,1)
  for(i in 1:num.recombs){
    rando.breaks<-sample(br.num, 1)
    rando.locations.base<-sort(sample(break.list, rando.breaks))
    rando.locations<-c(1, rando.locations.base, (dim(refs)[2]+1))#length align +1 
    #random.seq<-refs[sample(1:nrow(refs), rando.breaks+1),] # 4 sequences with which to recombine I think? 
    #random.seqD<-refs[sample(1:3, 1),] #pick random D
    #random.seqA<-refs[sample(4:6, 1),] #pick random A 
    choices<-c("A", "D")
    firstchoice<-sample(choices, 1)
    if(firstchoice=="A"){
      secondchoice<-"D"
    }
    if(firstchoice=="D"){
      secondchoice<-"A"
    }
    order<-rep(c(firstchoice, secondchoice), 10)
    for(o in 1:((rando.breaks)+1)){
      if(order[o]=="A" && o==1){
        random.seq<-refs[sample(1:3, 1),]  #pick random A 
      }
      if(order[o]=="D" && o==1){
        random.seq<-refs[sample(4:6, 1),]  #pick random D 
      }
      if(order[o]=="A" && o!=1 ){
        random.seq<-rbind(random.seq, refs[sample(1:3, 1),] ) #pick random A 
      }
      if(order[o]=="D"&& o!=1 ){
        random.seq<-rbind(random.seq, refs[sample(4:6, 1),] ) #pick random D 
      }
    } 
    
    
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
  Genome_breakpoints_simulated<<-text.records
  DNA.simulated<<-simulated.DNA
}


make_inter_recombs(10)
DNA.simulated
Genome_breakpoints_simulated
write.FASTA(DNA.simulated, "Inter_sim_recombs_A1-D6.fasta")
saveRDS(Genome_breakpoints_simulated, "Genome_breakpoints_simulated_SCUEAL_INTER.R")
#Genome_breakpoints_simulated_desire<-readRDS("Genome_breakpoints_simulated_SCUEAL_desire.R")
image(DNA.simulated)

#make 10 fasta files containing each 100 replicates of 10 simulated recombinants 
for(i in 1:nrow(DNA.simulated)){
  seed<-DNA.simulated[i,]
  rownames(seed)<-paste(rownames(seed), ".", "1", sep='')
  recombi<-seed
  for(r in 2:100){
    seed<-DNA.simulated[i,]
    rownames(seed)<-paste(rownames(seed), ".", r, sep='')
    recombi<-rbind.DNAbin(recombi, seed)
  }
  write.dna(recombi, paste("simulated_inter/Inter_sim_recomb", i, ".fasta", sep=''), format='fasta')
  
}

#


