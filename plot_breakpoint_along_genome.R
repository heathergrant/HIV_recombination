###################################################################################################################################################
#
# Title: Histograms of inter and intra subtype breakpoints 
# Description: above and below line histogram 
#
# Input files: (1) Genome_breakpoints R object created in other scripts 
#              
# Notes 
#  
###################################################################################################################################################
Genome_breakpoints <- readRDS("data/saved_objects/Genome_breakpoints_pdfs_duplicates.RDS")
#Genome_breakpoints <- readRDS("data/saved_objects/Genome_breakpoints_no_dup_no_transmitted_2_27_2pc.RDS")

#Genome_breakpoints <- readRDS("data/saved_objects/Genome_breakpoints_no_dup_no_transmitted.RDS")
#Firstly, number of breakpoints per genome 

hist(as.numeric(Genome_breakpoints$num_breaks), breaks=c(0,1,2,3,4,5,6,7), col='lightblue', xlab="Number of breakpoints per genome", main="Distribution of breakpoints per genome")
     
A<-Genome_breakpoints[Genome_breakpoints$CleanSub=="A1",]
D<-Genome_breakpoints[Genome_breakpoints$CleanSub=="D",]
A1.D<-Genome_breakpoints[Genome_breakpoints$CleanSub=="A1,D",]
pure<-rbind(A,D)
dim(pure)#225 
dim(A1.D)#164 # 151  is 13 ? YES 


###################################################################################################################################################
br.list.intraD<-c()
for (i in 1:nrow(D)){
  brk <- strsplit(as.character(D$DD[i]),";")
  br.list.intraD<-c(br.list.intraD, brk)
}
#for (i in 1:nrow(AD)){
#  brk <- strsplit(as.character(D$DD[i]),";")
#  br.list.intraD<-c(br.list.intraD, brk)
#}

br.list.intraD
br.list.intraD<-unlist(br.list.intraD)
br.list.intraD<-as.numeric(br.list.intraD)
br.list.intraD<-br.list.intraD[!is.na(br.list.intraD)]
###################################################################################################################################################
br.list.intraA<-c()
for (i in 1:nrow(A)){
  brk <- strsplit(as.character(A$AA[i]),";")
  br.list.intraA<-c(br.list.intraA, brk)
}
#for (i in 1:nrow(AD)){                             #remove the intra breakpoints within A1/D recombinants for normalization purposes 
#  brk <- strsplit(as.character(A$AA[i]),";")      
#  br.list.intraA<-c(br.list.intraA, brk)           
#}

br.list.intraA
br.list.intraA<-unlist(br.list.intraA)
br.list.intraA<-as.numeric(br.list.intraA)
br.list.intraA<-br.list.intraA[!is.na(br.list.intraA)]

###################################################################################################################################################
br.list.AD<-c()
for (i in 1:nrow(AD)){
  brk <- strsplit(as.character(A1.D$AD[i]),";")
  br.list.AD<-c(br.list.AD, brk)
}
for (i in 1:nrow(AD)){
  brk <- strsplit(as.character(A1.D$DA[i]),";")
  br.list.AD<-c(br.list.AD, brk)
}

br.list.AD
br.list.AD<-unlist(br.list.AD)
br.list.AD<-as.numeric(br.list.AD)
br.list.AD<-br.list.AD[!is.na(br.list.AD)]
###################################################################################################################################################
#normalizing
#nom.br.list.AD<-h1$counts/(146+165+100)
#nom.br.list.A<-h2$counts/(146+165)
#nom.br.list.D<-h3$counts/(100+165)

#plot(nom.br.list.A)
###################################################################################################################################################



complex<-Genome_breakpoints[Genome_breakpoints$CleanSub!="A1"   &Genome_breakpoints$CleanSub!="C"   & Genome_breakpoints$CleanSub!="D"  & Genome_breakpoints$CleanSub!="A1,D",]
#complex<-Genome_breakpoints[Genome_breakpoints$CleanSub=='C,D',]


br.list.other<-c()
for (i in 1:nrow(complex)){
  brk <- strsplit(as.character(complex$inter[i]),";")
  br.list.other<-c(br.list.other, brk)
}
br.list.other<-as.numeric(unlist(br.list.other))

#START FIGURE 1 
#recipe for a 4 bar strip with a lefthand side column of space.   
refs.length<-8293
win=300

#betw<-setdiff(br.list.8, br.list.8I)
#intra.list<-c(br.list.intraA, br.list.intraD)
###UPSIDE DOWN HIST PLOTS 
h1<-hist(as.numeric(unlist(br.list.AD)), breaks = seq(from = 0, to = refs.length, by = win), plot=FALSE)
#h2<-hist(as.numeric(unlist(intra.list)), breaks = seq(from = 0, to = refs.length, by = win), plot=FALSE)
h2<-hist(as.numeric(unlist(br.list.other)), breaks = seq(from = 0, to = refs.length, by = win), plot=FALSE)
#h3<-hist(as.numeric(unlist(br.list.intraD)), breaks = seq(from = 0, to = refs.length, by = win), plot=FALSE)
#
h4<-hist(as.numeric(rep(0, 27)), breaks = seq(from = 0, to = refs.length, by = win), plot=FALSE)
h4$counts[1]<-0
plot(h4) #blank histogram
h2$counts = - h2$counts
#h3$counts = - h3$counts
#h4$counts = -h4$counts
h1$density<-h1$counts/(dim(A1.D)[1])
h2$density<-h2$counts/nrow(complex) #143 #
h4$density<-h2$density#+h3$density
#h2$density<-h2$counts/(nrow(A)+nrow(D)) #As. 

hmax = max(h1$density)

hmin=min(h4$density)  #-0.16

X = c(h1$breaks, h2$breaks)#, h3$breaks)
xmax = max(X)
xmin = min(X)
par(mai=c(1, 2, 1, 0.2)) ###CAN SEE LABELS NOW
plot(h1, ylim=c(hmin, hmax), col="tomato", xlim=c(xmin, xmax), main="", xlab='Genome position', 
     freq = FALSE, cex.main=3, cex.lab=2, 
     ylab=expression("Normalized frequency\n(breakpoints/window/genome)"),
     yaxt='n',
     xaxt='n'
)
#mycol <- rgb(0, 0, 255, max = 255, alpha = 025, names = "blue50")
#lines(h3, freq=FALSE, col=mycol) #intraD #make transparent

#lines(h4,freq=FALSE, col="slategray4") #THIS IS ALLLLL combined 
lines(h2, freq=FALSE, col='slategray1') #the h2 A breaks normalised by all intra 


Axis(side=1,labels=T,cex.axis=1.5,at = seq(0, 8200, by = 300))  
Axis(side=2, labels=abs(seq(-0.4,0.2, by=0.2)), cex.axis=1.5, at=seq(-0.4,0.2, by=0.2))
#Axis(side=2, labels=abs(seq(0,0.2, by=0.05)), cex.axis=2, at=seq(0,0.2, by=0.05))
#legend("topleft", c("Inter-subtype (A/D)", "Intra-subtype (A/A)", "Intra-subtype (D/D)"), col=c("tomato", rgb(1,0,0,1/4), rgb(0,0,1,1/4)), lwd=10, cex=1.5)

legend("topright", c("Inter-subtype (A/D)", "Intra-subtype (A/A)", "Intra-subtype (D/D)"), col=c("tomato", "slategray1", "slategray4"), lwd=10, cex=1.3)

## rectangle annotations 

## rectangle annotations 
rect(0, hmin+0.1, 396, hmin+0.15, lwd=1) #, col="dodgerblue1") #gag 
rect(397, hmin+0.1, 1304, hmin+0.15, lwd=1) #, col="dodgerblue2") #gag 
rect(1305, hmin+0.1, 1500, hmin+0.15, lwd=1) #, col="dodgerblue3") #gag 

#rect(1302, hmin+0.05, 4317, hmin+0.1, lwd=1) #, col='darkorange') #pol 
rect(1302, hmin, 4317, hmin+0.05, lwd=1) #, col='darkorange') #pol 
rect(1472, hmin, 1769, hmin+0.05, lwd=1) #, col='darkorange1') #pol 
rect(1770, hmin, 3450, hmin+0.05, lwd=1) #, col='darkorange2') #pol 
rect(3451, hmin, 4317, hmin+0.05, lwd=1) #, col='darkorange3') #pol 

rect(4265, hmin+0.1, 4842, hmin+0.15, lwd=1) #, col='gray') #vif 
rect(4785, hmin, 5072, hmin+0.05, lwd=1) #, col = "dodgerblue1") #vpr 
rect(5057, hmin+0.05, 5273, hmin+0.1, lwd=1) #, col="magenta4")#tat1
rect(7264, hmin+0.1, 7305, hmin+0.15, lwd=1) #, col="magenta4")#tat2 
rect(5290, hmin+0.05, 5534, hmin+0.1, lwd=1) #, col= "brown")# vpu
rect(5198, hmin+0.0, 5273, hmin+0.05, lwd=1) #, col= "chartreuse4")# rev1
rect(7264, hmin+0.05, 7538, hmin+0.1, lwd=1) #, col= "chartreuse4")# rev2

rect(5455, hmin, 6644, hmin+0.05, lwd=1) #, col="coral1") #gp120 
rect(6645, hmin, 7680, hmin+0.05, lwd=1) #, col="coral2") #gp41
rect(7579, hmin+0.1, 8292, hmin+0.15, lwd=1) #, col="cyan4") #nef


text(200, (hmin+0.125),"p17", col="black", cex=1)
text(850, (hmin+0.125),"p24", col="black", cex=1)

text(1600, (hmin+0.025),"prot", col="black", cex=1)
text(2500, (hmin+0.025),"RT", col="black", cex=1)
text(3900, (hmin+0.025),"int", col="black", cex=1)
text(4500, (hmin+0.125),"vif", col="black", cex=1)
text(4900, (hmin+0.025),"vpr", col="black", cex=1)

#text(5150, (hmin+0.11), "tat", col='magenta4', cex=1)
#text(7150, (hmin+0.125), "tat", col='magenta4', cex=1)
text(5150, (hmin+0.075), "tat", col='black', cex=1)
text(7400, (hmin+0.125), "tat", col='black', cex=1)
#curveGrob(5150, hmin+0.11, 7200, hmin+0.125)

#text(5400, (hmin+0.108), "vpu", col='brown', cex=1)

#text(5250, (hmin-0.0125), "rev", col='chartreuse4', cex=1)
#text(7150, (hmin+0.09), "rev", col='chartreuse4', cex=1)
text(5400, (hmin+0.075), "vpu", col='black', cex=1)

text(5250, (hmin-0.0125), "rev", col='black', cex=1)
text(7400, (hmin+0.075), "rev", col='black', cex=1)
text(7900, (hmin+0.125),"nef", col="black", cex=1)

text(6100, (hmin+0.025),"g120", col="black", cex=1)
text(7100, (hmin+0.025),"gp41", col="black", cex=1)

####end figure 


#complex genomes show same patterns? 
complex<-Genome_breakpoints[Genome_breakpoints$CleanSub!="A1"   & Genome_breakpoints$CleanSub!="D"  & Genome_breakpoints$CleanSub!="A1,D" & Genome_breakpoints$CleanSub!="C" ,]
#complex<-Genome_breakpoints[Genome_breakpoints$CleanSub=="A1,C"  ,]
#complex<-Genome_breakpoints[Genome_breakpoints$CleanSub=="C"  ,]
#complex<-Genome_breakpoints[Genome_breakpoints$CleanSub=="A1,C,D"  ,]
#complex<-Genome_breakpoints[Genome_breakpoints$CleanSub=="A1,C,D"  ,]
#nrow(complex) #13 A1, C genomes 
#nrow(complex) #8 C genomes 
#nrow(complex) # 13 genomes NO breakpoints in envelope. 

br.list.other<-c()
for (i in 1:nrow(complex)){
  brk <- strsplit(as.character(complex$Breakpoints[i]),";")
  br.list.other<-c(br.list.other, brk)
}
br.list.other<-as.numeric(unlist(br.list.other))
h1<-hist(as.numeric(unlist(br.list.other)), breaks = seq(from = 0, to = refs.length, by = 300), col='tomato', main='non A, D or AD recombs, n=76', plot=F)
#yes. wow. 
#
nrow(complex) #76 genomes 
h1$density<-h1$counts/(76)
hmax = max(h1$density)
#hmin = min(h2$density+h3$density)
hmin=-0.1
#hmin = min(h2$counts)
#X = c(h1$breaks, h2$breaks, h3$breaks)
xmax = max(X)
xmin = min(X)
par(mai=c(1, 2, 1, 0.2)) ###CAN SEE LABELS NOW
#ylab=cat("Normalized frequency", "\n" ,"(breakpoints/window/genome)")
plot(h1, ylim=c(hmin, hmax), col="goldenrod4", xlim=c(xmin, xmax), main="", xlab='Genome position', 
     freq = FALSE, cex.main=3, cex.lab=1.5, 
     ylab=expression("Normalized frequency\n(breakpoints/window/genome)"),
     yaxt='n',
     xaxt='n'
)
Axis(side=1,labels=T,cex.axis=1.5,at = seq(0, 8200, by = 300))
Axis(side=2, labels=abs(seq(0,0.4, by=0.2)), cex.axis=1.5, at=seq(0,0.4, by=0.2))
#Axis(side=2, labels=abs(seq(0,0.2, by=0.05)), cex.axis=2, at=seq(0,0.2, by=0.05))
legend("topleft", c("All other genomes (non A, non D, non A1/D)"), col="goldenrod4", lwd=10, cex=1.5)



## rectangle annotations 
rect(0, hmin+0.04, 396, hmin+0.06, lwd=1) #, col="dodgerblue1") #gag 
rect(397, hmin+0.04, 1304, hmin+0.06, lwd=1) #, col="dodgerblue2") #gag 
rect(1305, hmin+0.04, 1500, hmin+0.06, lwd=1) #, col="dodgerblue3") #gag 

#rect(1302, hmin+0.05, 4317, hmin+0.1, lwd=1) #, col='darkorange') #pol 
rect(1302, hmin, 4317, hmin+0.02, lwd=1) #, col='darkorange') #pol 
rect(1472, hmin, 1769, hmin+0.02, lwd=1) #, col='darkorange1') #pol 
rect(1770, hmin, 3450, hmin+0.02, lwd=1) #, col='darkorange2') #pol 
rect(3451, hmin, 4317, hmin+0.02, lwd=1) #, col='darkorange3') #pol 

rect(4265, hmin+0.04, 4842, hmin+0.06, lwd=1) #, col='gray') #vif 
rect(4785, hmin, 5072, hmin+0.02, lwd=1) #, col = "dodgerblue1") #vpr 
rect(5057, hmin+0.02, 5273, hmin+0.04, lwd=1) #, col="magenta4")#tat1
rect(7264, hmin+0.04, 7305, hmin+0.06, lwd=1) #, col="magenta4")#tat2 
rect(5290, hmin+0.02, 5534, hmin+0.04, lwd=1) #, col= "brown")# vpu
rect(5198, hmin, 5273, hmin+0.02, lwd=1) #, col= "chartreuse4")# rev1
rect(7264, hmin+0.02, 7538, hmin+0.04, lwd=1) #, col= "chartreuse4")# rev2

rect(5455, hmin, 6644, hmin+0.02, lwd=1) #, col="coral1") #gp120 
rect(6645, hmin, 7680, hmin+0.02, lwd=1) #, col="coral2") #gp41
rect(7579, hmin+0.04, 8292, hmin+0.06, lwd=1) #, col="cyan4") #nef


text(200, (hmin+0.05),"p17", col="black", cex=1)
text(850, (hmin+0.05),"p24", col="black", cex=1)

text(1600, (hmin+0.01),"prot", col="black", cex=1)
text(2500, (hmin+0.01),"RT", col="black", cex=1)
text(3900, (hmin+0.01),"int", col="black", cex=1)
text(4500, (hmin+0.05),"vif", col="black", cex=1)
text(4900, (hmin+0.01),"vpr", col="black", cex=1)
text(7900, (hmin+0.05),"nef", col="black", cex=1)


text(6100, (hmin+0.01),"g120", col="black", cex=1)
text(7100, (hmin+0.01),"gp41", col="black", cex=1)

hmin=0.1
rect(0, hmin+0.6, 396, hmin+0.9, lwd=1 , col="dodgerblue1") #gag 
rect(397, hmin+0.6, 1304, hmin+0.9, lwd=1 , col="dodgerblue2") #gag 
rect(1305, hmin+0.6, 1500, hmin+0.9, lwd=1 , col="dodgerblue3") #gag 

rect(1302, hmin, 4317, hmin+0.3, lwd=1 , col='darkorange') #pol 
rect(1302, hmin, 4317, hmin+0.3, lwd=1 , col='darkorange') #pol 
rect(1472, hmin, 1769, hmin+0.3, lwd=1 , col='darkorange1') #pol 
rect(1770, hmin, 3450, hmin+0.3, lwd=1 , col='darkorange2') #pol 
rect(3451, hmin, 4317, hmin+0.3, lwd=1 , col='darkorange3') #pol 

rect(4265, hmin+0.6, 4842, hmin+0.9, lwd=1 , col='gray') #vif 
rect(4785, hmin, 5072, hmin+0.3, lwd=1 , col = "dodgerblue1") #vpr 
rect(5057, hmin+0.3, 5273, hmin+0.6, lwd=1 , col="magenta4")#tat1
rect(7264, hmin+0.6, 7305, hmin+0.9, lwd=1 , col="magenta4")#tat2 
rect(5290, hmin+0.3, 5534, hmin+0.6, lwd=1 , col= "brown")# vpu
rect(5198, hmin, 5273, hmin+0.3, lwd=1 , col= "chartreuse4")# rev1
rect(7264, hmin+0.3, 7538, hmin+0.6, lwd=1 , col= "chartreuse4")# rev2

rect(5455, hmin, 6644, hmin+0.3, lwd=1 , col="coral1") #gp120 
rect(6645, hmin, 7680, hmin+0.3, lwd=1 , col="coral2") #gp41
rect(7579, hmin+0.6, 8292, hmin+0.9, lwd=1 , col="cyan4") #nef


text(200, (hmin+0.75),"p17", col="black", cex=1.7)
text(850, (hmin+0.75),"p24", col="black", cex=1.7)

text(1600, (hmin+0.15),"prot", col="black", cex=1.5)
text(2500, (hmin+0.15),"RT", col="black", cex=1.7)
text(3900, (hmin+0.15),"int", col="black", cex=1.7)
text(4500, (hmin+0.75),"vif", col="black", cex=1.7)
text(4900, (hmin+0.15),"vpr", col="black", cex=1.7)
text(7900, (hmin+0.75),"nef", col="black", cex=1.7)


text(6100, (hmin+0.15),"g120", col="black", cex=1.7)
text(7100, (hmin+0.15),"gp41", col="black", cex=1.7)
text(5400, hmin+0.45, "vpu", cex=1.5)

#text(5400, hmin+0.75, "rev1 and tat1", cex=1.5)
text(5150, hmin+0.45, "tat", cex=1.5)
text(7400, hmin+0.75, "tat", cex=1.5)
text(7400, hmin+0.45, "rev", cex=1.4)
# #plot gp41 out 


###coreelatino
cor.test(h1$counts, h2$counts)
cor.test(h1$density, h2$density)
cor.test(br.list.other, br.list.AD)
####


p17start<-1 ; p17end<-396 #gag rf1 
p24start<-397 ; p24end<-1304 #gag rf1 
p67start<-1305 ; p67end<-1500 #gag rf1 
#
polstart<-1302 ; polend<-1471 #pol rf3 
protstart<-1472 ; protend<-1769 #pol rf3 
RTstart<-1770 ; RTend<-3450 #pol rf3 
INstart<-3451 ; INend<-4317 #pol rf3 
#
vifstart<-4265 ; vifend<- 4842 #rf1
vprstart<-4785 ; vprend<-5072 #rf3 
tat1start<-5057 ; tat1end<-5273  #rf2
tat2start<-7264 ; tat2end<-7305  #rf1
vpustart<- 5290 ; vpuend<- 5534 #5,290 -> 5,534
rev1start<-5198 ; rev1end<- 5273 #rf3
rev2start<-7264 ; rev2end<- 7538 #rf2
gp120start<-5455 ; gp120end<-6645 #rf3
gp41start<-6646 ; gp41end<-7572 #rf3
nefstart<-7579 ; nefend<-8292 #rf1
#hmin=-0.2



