###################################################################################################################################################
#
# Title: Plotting functions 
# Description: Plotting of the stacked genomes to show their subtype along the genome 
#
# Input files: (1) R object - Genome_breakpoints - from other scripts 
#              (2) K means clustering results in another script (subtype_distance.R)
# Notes 
#
###################################################################################################################################################
Genome_breakpoints<-readRDS("....")


Group1_breakpoints <- Genome_breakpoints[Genome_breakpoints$Name %in% group1$Name, ]
Group2_breakpoints <- Genome_breakpoints[Genome_breakpoints$Name %in% group2$Name, ]
Group3_breakpoints <- Genome_breakpoints[Genome_breakpoints$Name %in% group3$Name, ]
Group4_breakpoints <- Genome_breakpoints[Genome_breakpoints$Name %in% group4$Name, ]
Group5_breakpoints <- Genome_breakpoints[Genome_breakpoints$Name %in% group5$Name, ]

head(recombs)

recombs<-rbind(Group1_breakpoints, Group2_breakpoints, Group3_breakpoints, Group4_breakpoints,Group5_breakpoints,Group6_breakpoints,Group7_breakpoints,Group8_breakpoints,Group9_breakpoints)
#recombs<-rbind(Group1_breakpoints, Group2_breakpoints, Group3_breakpoints)

par(mai=c(1, 1.1, 1, 0.2)) ###CAN SEE LABELS NOW
pdf("test.pdf")
#plot(1:10)
#dev.off()
recombs<-Group9_breakpoints
plot(0,type='n',xlim=c(0,8287),ylim=c(0,164+5),xlab="Genome position",main="",ylab="",axes=F,cex.lab=2)
plot(0,type='n',xlim=c(0,8287),ylim=c(0, dim(recombs)[2]+9),xlab="Genome position",main="",ylab="",axes=F,cex.lab=2)

Axis(side=1,labels=T,cex.axis=1)
#rect(0,(i+0.1),8289,i+0.9,lwd=0.3, col='white')
for(i in 1:dim(recombs)[1]){
  BG1<- recombs[i,]
  subs<-strsplit(as.character(BG1$Subtypes), ',')[[1]] #subtype list 
  
  brpnts<-strsplit(as.character(BG1$Breakpoints), ';')[[1]] #breakpoint list 
  brpnts<-prepend(brpnts, "0", before = 1) #add a 0 as first breakpoint 
  brpnts<-append(brpnts, "8287") #add 8287 as the final breakpont 
  for(y in 1:length(subs)){
    
    if(subs[y]=="A1"){
      rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=0.5, col='dodgerblue3')
      #rect(BG1[y,1],(i+0.1),BG1[y,2],i+0.9,lwd=0.2,col='firebrick1')
      #rect(brpnts[y+1], (i+0.1), (brpnts[y+1]+0.1), i+0.9, lwd=0, col='white')
    }
    if (subs[y]=="D"){
      rect(((brpnts[y])),(i+0.15),((brpnts[y+1])),i+0.85,lwd=0.5, col='darkgoldenrod1')
    }
    if (subs[y]=="C"){
      rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='red')
    }   
    if (subs[y]=="G"| subs[y]=="BG/G"){
      rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='pink')
    }
    if (subs[y]=="BG/G/N/O" | subs[y]=="B/BF/D/F1/F2/K" | subs[y]=="N/O" | subs[y]=="BG/G/N/O"| subs[y]=="A1,B/BF/D/F1/F2/K"| subs[y]=="A1,B/BF/D/F1/F2/K" | subs[y]=="BF/F1/F2/K") {
      rect(brpnts[y],(i+0.1),brpnts[y+1],i+0.9,lwd=0.3, col='grey')
    }
    text(9200, i+0.5, BG1$psudo_lab, cex=1.5)
    #text(-2000, end+nrow(Group4_breakpoints)-10, "Group1", cex=1.5)
  }
}

#dev.off()
#rect(#x start, Ystart, xend, Yend )
#rect(#left, bottom, right, top ), 

dim(recombs)

#side panel labels 
end<-1
rect(-300, end, -650, end+nrow(Group4_breakpoints), col='coral')
text(-2000, end+nrow(Group4_breakpoints)-10, "Group1", cex=1.5)
end<-end+nrow(Group4_breakpoints)

rect(-300, end, -650, end+nrow(Group2_breakpoints), col='lightpink')
text(-2000, end+nrow(Group2_breakpoints)-10, "Group2", cex=1.5)
end<-end+nrow(Group2_breakpoints)
rect(-300, end, -650, end+nrow(Group3_breakpoints), col='coral')
text(-2000, end+nrow(Group3_breakpoints)-10, "Group3", cex=1.5)
end<-end+nrow(Group3_breakpoints)
rect(-300, end, -650, end+nrow(Group1_breakpoints), col='lightpink')
text(-2000, end+nrow(Group1_breakpoints)-10, "Group4", cex=1.5)
end<-end+nrow(Group1_breakpoints)

rect(-300, end, -650, end+nrow(Group5_breakpoints), col='coral')
text(-2000, end+nrow(Group5_breakpoints)-10, "Group5", cex=1.5)




hmin=-3

# HIV gene map 
rect(0, hmin, 396, hmin+1, lwd=1, col="dodgerblue1") #gag 
rect(397, hmin, 1304, hmin+1, lwd=1, col="dodgerblue2") #gag 
rect(1305, hmin, 1500, hmin+1, lwd=1, col="dodgerblue3") #gag 

#rect(1302, hmin+3, 4317, hmin+0.1, lwd=1, col='darkorange') #pol 
rect(1302, hmin-2, 4317, hmin-1, lwd=1, col='darkorange') #pol 
rect(1472, hmin-2, 1769, hmin-1, lwd=1, col='darkorange1') #pol 
rect(1770, hmin-2, 3450, hmin-1, lwd=1, col='darkorange2') #pol 
rect(3451, hmin-2, 4317, hmin-1, lwd=1, col='darkorange3') #pol 

rect(4265, hmin, 4842, hmin+1, lwd=1, col='gray') #vif 
rect(4785, hmin-2, 5072, hmin-1, lwd=1, col = "dodgerblue1") #vpr 
rect(5057, hmin-1, 5273, hmin, lwd=1, col="magenta4")#tat1
rect(7264, hmin, 7305, hmin+1, lwd=1, col="magenta4")#tat2 
rect(5290, hmin-1, 5534, hmin, lwd=1, col= "brown")# vpu
rect(5198, hmin-2, 5273, hmin-1, lwd=1, col= "chartreuse4")# rev1
rect(7264, hmin-1, 7538, hmin, lwd=1, col= "chartreuse4")# rev2

rect(5455, hmin-2, 6644, hmin-1, lwd=1, col="coral1") #gp120 
rect(6645, hmin-2, 7680, hmin-1, lwd=1, col="coral2") #gp41
rect(7579, hmin, 8292, hmin+1, lwd=1, col="cyan4") #nef

##top too 
hmin=168
rect(0, hmin, 396, hmin+1, lwd=1, col="dodgerblue1") #gag 
rect(397, hmin, 1304, hmin+1, lwd=1, col="dodgerblue2") #gag 
rect(1305, hmin, 1500, hmin+1, lwd=1, col="dodgerblue3") #gag 

#rect(1302, hmin+3, 4317, hmin+0.1, lwd=1, col='darkorange') #pol 
rect(1302, hmin-2, 4317, hmin-1, lwd=1, col='darkorange') #pol 
rect(1472, hmin-2, 1769, hmin-1, lwd=1, col='darkorange1') #pol 
rect(1770, hmin-2, 3450, hmin-1, lwd=1, col='darkorange2') #pol 
rect(3451, hmin-2, 4317, hmin-1, lwd=1, col='darkorange3') #pol 

rect(4265, hmin, 4842, hmin+1, lwd=1, col='gray') #vif 
rect(4785, hmin-2, 5072, hmin-1, lwd=1, col = "dodgerblue1") #vpr 
rect(5057, hmin-1, 5273, hmin, lwd=1, col="magenta4")#tat1
rect(7264, hmin, 7305, hmin+1, lwd=1, col="magenta4")#tat2 
rect(5290, hmin-1, 5534, hmin, lwd=1, col= "brown")# vpu
rect(5198, hmin-2, 5273, hmin-1, lwd=1, col= "chartreuse4")# rev1
rect(7264, hmin-1, 7538, hmin, lwd=1, col= "chartreuse4")# rev2

rect(5455, hmin-2, 6644, hmin-1, lwd=1, col="coral1") #gp120 
rect(6645, hmin-2, 7680, hmin-1, lwd=1, col="coral2") #gp41
rect(7579, hmin, 8292, hmin+1, lwd=1, col="cyan4") #nef)
dev.off()

