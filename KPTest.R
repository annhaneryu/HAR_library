###Goal: Fit model to human/chimp ratios of RNA/DNA per HAR, test for cis and trans effects

##setup
#set directory
setwd("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/data/")
load("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/KPQuant.RData")

#these are non-permutation HARs
#h1 is RNA/DNA ratios for human, library 1 only (not logged)
#c1 is RNA/DNA ratios for chimp, library 1 only (not logged)
#lib1 is log2((human RNA/DNA ratio)/(chimp RNA/DNA ratio)), library 1 only
#useData is indicator of library 1 samples

##expressed HARs
#retain if human or chimp active vs 75th percentile of negative controls for that sample
ncut<-apply(negratios[,useData],2,function(x) quantile(x,0.75))
hexpress<-t(t(h1)>ncut)
cexpress<-t(t(c1)>ncut)
#plot
pdf("expressed.pdf")
barplot(hexpress,main="Human RNA/DNA > 0.75'tile of Negatives",horiz=TRUE,oma=c(1,3,1,1),las=2,cex.names=0.5)
barplot(cexpress,main="Chimp RNA/DNA > 0.75'tile of Negatives",horiz=TRUE,oma=c(1,3,1,1),las=2,cex.names=0.5)
heatmap(cbind(apply(hexpress,1,sum),apply(cexpress,1,sum)),labCol=c("Human","Chimp"),cexCol=0.8,labRow=NA,Colv=NA,Rowv=NA,col=cvBWO,scale="none",main="Number of samples with HAR expressed")
heatmap(matrix(as.numeric(hexpress),nr=nrow(hexpress)),scale="none",Colv=NA,Rowv=NA,main="Human RNA/DNA (expressed = yellow, not = red)",oma=c(2,1,3,1),labRow=NA,labCol=colnames(hexpress))
heatmap(matrix(as.numeric(cexpress),nr=nrow(cexpress)),scale="none",Colv=NA,Rowv=NA,main="Chimp RNA/DNA (expressed = yellow, not = red",labRow=NA,labCol=colnames(cexpress))
dev.off()
#choose which ones to filter out
#check how many not expressed in any condition
rownames(hexpress)[apply(hexpress,1,sum)==0]
rownames(cexpress)[apply(cexpress,1,sum)==0]
rownames(hexpress)[apply(hexpress+cexpress,1,sum)==0]
#distribution of replicates with expression
table(apply(hexpress,1,sum))
table(apply(cexpress,1,sum))
table(apply(cexpress+hexpress,1,sum))
#filter out the ones that are not expressed in any sample (human or chimp)
# i.e., <75th percentile of negative controls for human and chimp RNA/DNA
expressed<-lib1[rownames(hexpress)[apply(hexpress+cexpress,1,sum)>0],]

##fit linear model
library(limma)
#first lump together cell species and stages, adjust for batch effects
# and test for HARs with hg/pt ratios significantly different from zero
harDesign0<-model.matrix(~passage+weeks+unique, data=IDs)
fit0<-lmFit(expressed,harDesign0,na.action="na.exclude")
fit0<-eBayes(fit0)
int0<-topTable(fit0,coef=1,number="all")
#note: not sure we want to put in batch effects since human/chimp cancels these out.
#now look at significance of batch effects
pass<-topTable(fit0,coef=2,number="all")
wks<-topTable(fit0,coef=3,number="all")
uni<-topTable(fit0,coef=4,number="all")
#adjusted pvalue distributions for human vs. chimp seq and batch effects
pdf("adjp.pdf")
boxplot(cbind(int0$adj.P.Val,pass$adj.P.Val,wks$adj.P.Val,uni$adj.P.Val),names=c("seq","passage","week","UMI"))
heatmap(cbind(int0[rownames(lib1),"adj.P.Val"],pass[rownames(lib1),"adj.P.Val"],wks[rownames(lib1),"adj.P.Val"],uni[rownames(lib1),"adj.P.Val"]),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Adjusted Pvals",labRow=NA,labCol=c("seq","passage","week","UMI"))
dev.off()
pdf("ranks.pdf")
qqplot(int1[rownames(lib1),"adj.P.Val"],int0[rownames(lib1),"adj.P.Val"],xlab="Without batch effects",ylab="With batch effects",main="QQ-Plot of Adjusted Pvalues")
plot(log10(int1[rownames(lib1),"adj.P.Val"]),log10(int0[rownames(lib1),"adj.P.Val"]),xlab="Without batch effects",ylab="With batch effects",main="Log10 Adjusted Pvalues")
dev.off()
#conclusion: batch effects are not necessary with ratio set up

#redo linear model no batch effects included
#this model tests H0: log2(human:chimp)=0 across all conditions
harDesign1 <- model.matrix(~1, data=IDs)
fit1<-lmFit(expressed,harDesign1,na.action="na.exclude")
fit1<-eBayes(fit1)
out1<-topTable(fit1,coef=1,number="all")
#look at differentially active HARs
#default sort is on p-value, look at top 50
out1[1:50,]
#sort on fold change, look at top 50
out1[rev(order(abs(int1$logFC))),][1:50,]
#Cyagen HARs are still hits without batch effects?
cy<-c("2xHAR.133","2xHAR.548","2xHAR.356","2xHAR.308")
out1[cy,]
#all are significant, 2xHAR.548 is biggest fold change
#grab average RNA/DNA ratios for human and chimp
h1means<-apply(h1[rownames(out1),],1,mean)
c1means<-apply(c1[rownames(out1),],1,mean)
#write to a file
write.table(cbind(out1,h1means,c1means),"DiffExprHARs.tab",sep="\t",quote=F)

#model expression in each cell type and species of origin
#separtely test for differential expression of human vs. chimp alleles in NPCs
# and the same in GPCs
#also HARs with log2(human:chimp) different in human vs. chimp cells
# and HARs with log2(human:chimp) different in NPC vs. GPC stage
IDs$combo<-interaction(IDs$stage,IDs$species)
harDesign2<-model.matrix(~0+combo,IDs)
colnames(harDesign2)<-levels(IDs$combo)
fit2<-lmFit(expressed,harDesign2,na.action="na.exclude")
cm2<-makeContrasts(
    npc = npc.hg+npc.pt,
    gpc = gpc.hg+gpc.pt,
    spp = npc.hg+gpc.hg-(npc.pt+gpc.pt),
    stg = npc.hg+npc.pt-(gpc.hg+gpc.pt),
    levels = levels(IDs$combo))
fit2<-contrasts.fit(fit2,contrasts=cm2)
fit2<-eBayes(fit2)
sum(topTableF(fit2,number="all")$adj.P.Val<0.05) #401 have at least one sig contrast

#npc only
npchits<-topTable(fit2,coef="npc",number="all")
sum(npchits$adj.P.Val<0.05) #356
#compare to overall results
out1[rownames(npchits[npchits$adj.P.Val<0.05,]),]
npchits[rownames(out1[out1$adj.P.Val<0.05,]),]
sum(out1[rownames(npchits[npchits$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #338
sum(npchits[rownames(out1[out1$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #338

#gpc only
gpchits<-topTable(fit2,coef="gpc",number="all")
sum(gpchits$adj.P.Val<0.05) #342
#compare to overall results
out1[rownames(gpchits[gpchits$adj.P.Val<0.05,]),]
gpchits[rownames(out1[out1$adj.P.Val<0.05,]),]
sum(out1[rownames(gpchits[gpchits$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #325
sum(gpchits[rownames(out1[out1$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #325
#compare to NPCs
sum(gpchits[rownames(npchits[npchits$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #269
sum(npchits[rownames(gpchits[gpchits$adj.P.Val<0.05,]),]$adj.P.Val<0.05) #269

#species of origin of cell line effects
topTable(fit2,coef="spp",number=20)
#only 2xHAR.518 and HAR51 have significant species or origin effects

#NPC vs. GPC effects
topTable(fit2,coef="stg",number=20)
#only 2xHAR.319, HAR5, 2xHAR.28, 2xHAR.238, and 2xHAR.2 have significant NPC vs. GPC

write.table(cbind(topTableF(fit2,number="all")[rownames(out1),],topTable(fit2,coef="npc",number="all")[rownames(out1),"adj.P.Val"],topTable(fit2,coef="gpc",number="all")[rownames(out1),"adj.P.Val"],topTable(fit2,coef="spp",number="all")[rownames(out1),"adj.P.Val"],topTable(fit2,coef="stg",number="all")[rownames(out1),"adj.P.Val"]),"Contrasts.tab",sep="\t",quote=F)

#save image
save.image("KPTest.RData")
