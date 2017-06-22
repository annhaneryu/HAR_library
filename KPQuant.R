###Goal: Quantify RNA/DNA ratio for each HAR and each species, then take human/chimp ratio 

##set directory
setwd("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/data/")
load("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/KPQuant.RData")


#this function reads in counts for each oligo from an RNA or DNA file
#modified from Sean's to drop terms that later cancel out and not compute ratio yet
getHARdata2 <- function (path,d1f,d2f,d3f,r1f,r2f,r3f) {
	# read in DNA and RNA files for all 3 replicates
	d1 <- read.delim(paste(path,d1f,sep="/"),header=F)
	d2 <- read.delim(paste(path,d2f,sep="/"),header=F)
	d3 <- read.delim(paste(path,d3f,sep="/"),header=F)
	r1 <- read.delim(paste(path,r1f,sep="/"),header=F)
	r2 <- read.delim(paste(path,r2f,sep="/"),header=F)
	r3 <- read.delim(paste(path,r3f,sep="/"),header=F)
	
	# merge data by replicate
	data.1 <- merge(r1,d1,by="V1")[,-3]
	data.2 <- merge(r2,d2,by="V1")[,-3]
	data.3 <- merge(r3,d3,by="V1")[,-3]
	colnames(data.1) <- c("barcodes","X","Y","name")
	colnames(data.2) <- c("barcodes","X","Y","name")
	colnames(data.3) <- c("barcodes","X","Y","name")

	# get category of enhancer fragment
	data.1$category <- sapply(strsplit(as.character(data.1$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
	data.2$category <- sapply(strsplit(as.character(data.2$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
	data.3$category <- sapply(strsplit(as.character(data.3$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })

	# get counts for each fragment across all RNA and DNA barcodes
	res.1 <- as.data.frame(t(sapply(unique(data.1$category),FUN=function(x) { sel <- which(data.1$category == x); c(sum(data.1$X[sel]),sum(data.1$Y[sel]),length(sel)) } )))
	res.2 <- as.data.frame(t(sapply(unique(data.2$category),FUN=function(x) { sel <- which(data.2$category == x); c(sum(data.2$X[sel]),sum(data.2$Y[sel]),length(sel)) } )))
	res.3 <- as.data.frame(t(sapply(unique(data.3$category),FUN=function(x) { sel <- which(data.3$category == x); c(sum(data.3$X[sel]),sum(data.3$Y[sel]),length(sel)) } )))
	res.1$name <- rownames(res.1)
	res.2$name <- rownames(res.2)
	res.3$name <- rownames(res.3)
	
	#return object
	v <- c()
	v$res1 <- res.1
	v$res2 <- res.2
	v$res3 <- res.3
	return(v)
}

##read in data and compute RNA/DNA ratio per oligo
wtc.npc.l1 <- getHARdata2(getwd(),"WTC-NPC-1-DNA.tsv","WTC-NPC-2-DNA.tsv","WTC-NPC-3-DNA.tsv","WTC-NPC-1-RNA.tsv","WTC-NPC-2-RNA.tsv","WTC-NPC-3-RNA.tsv")
pt2a.npc.l1 <- getHARdata2(getwd(),"Pt2a-NPC-1-DNA.tsv","Pt2a-NPC-2-DNA.tsv","Pt2a-NPC-3-DNA.tsv","Pt2a-NPC-1-RNA.tsv","Pt2a-NPC-2-RNA.tsv","Pt2a-NPC-3-RNA.tsv")
wtc.gpc.l1 <- getHARdata2(getwd(),"WTc-1-DNA.tsv","WTc-2-DNA.tsv","WTc-3-DNA.tsv","WTc-1-RNA.tsv","WTc-2-RNA.tsv","WTc-3-RNA.tsv")
wtc.gpc.l2 <- getHARdata2(getwd(),"WTC-GPC-1-DNA.tsv","WTC-GPC-2-DNA.tsv","WTC-GPC-3-DNA.tsv","WTC-GPC-1-RNA.tsv","WTC-GPC-2-RNA.tsv","WTC-GPC-3-RNA.tsv")
pt2a.gpc.l1 <- getHARdata2(getwd(),"Pt2A-1-DNA.tsv","Pt2A-2-DNA.tsv","Pt2A-3-DNA.tsv","Pt2A-1-RNA.tsv","Pt2A-2-RNA.tsv","Pt2A-3-RNA.tsv")
pt2a.gpc.l2 <- getHARdata2(getwd(),"Pt2A-GPC-1-DNA.tsv","Pt2A-GPC-2-DNA.tsv","Pt2A-GPC-3-DNA.tsv","Pt2A-GPC-1-RNA.tsv","Pt2A-GPC-2-RNA.tsv","Pt2A-GPC-3-RNA.tsv")
pt5c.gpc.l1 <- getHARdata2(getwd(),"Pt5C-GPC-1-DNA.tsv","Pt5C-GPC-2-DNA.tsv","Pt5C-GPC-3-DNA.tsv","Pt5C-GPC-1-RNA.tsv","Pt5C-GPC-2-RNA.tsv","Pt5C-GPC-3-RNA.tsv")
pt5c.npc.l1 <- getHARdata2(getwd(),"Pt5C-NPC-1-DNA.tsv","Pt5C-NPC-2-DNA.tsv","Pt5C-NPC-3-DNA.tsv","Pt5C-NPC-1-RNA.tsv","Pt5C-NPC-2-RNA.tsv","Pt5C-NPC-3-RNA.tsv")
hs1.gpc.l1 <- getHARdata2(getwd(),"HS1-11-GPC-1-DNA.tsv","HS1-11-GPC-2-DNA.tsv","HS1-11-GPC-3-DNA.tsv","HS1-11-GPC-1-RNA.tsv","HS1-11-GPC-2-RNA.tsv","HS1-11-GPC-3-RNA.tsv")
hs1.npc.l1 <- getHARdata2(getwd(),"HS1-11-NPC-1-DNA.tsv","HS1-11-NPC-2-DNA.tsv","HS1-11-NPC-3-DNA.tsv","HS1-11-NPC-1-RNA.tsv","HS1-11-NPC-2-RNA.tsv","HS1-11-NPC-3-RNA.tsv")

##organize annotations of oligos
# get all possible fragment names from all samples and replicates
uid <- data.frame(sort(unique(c(
	wtc.gpc.l1$res1$name,
	wtc.gpc.l1$res2$name,
	wtc.gpc.l1$res3$name,
	wtc.gpc.l2$res1$name,
	wtc.gpc.l2$res2$name,
	wtc.gpc.l2$res3$name,
	pt2a.gpc.l1$res1$name,
	pt2a.gpc.l1$res2$name,
	pt2a.gpc.l1$res3$name,
	pt2a.gpc.l2$res1$name,
	pt2a.gpc.l2$res2$name,
	pt2a.gpc.l2$res3$name,
	wtc.npc.l1$res1$name,
	wtc.npc.l1$res1$name,
	wtc.npc.l1$res1$name,
	pt2a.npc.l1$res1$name,
	pt2a.npc.l1$res1$name,
	pt2a.npc.l1$res1$name,
	hs1.gpc.l1$res1$name,
	hs1.gpc.l1$res2$name,
	hs1.gpc.l1$res3$name,
	pt5c.npc.l1$res1$name,
	pt5c.npc.l1$res2$name,
	pt5c.npc.l1$res3$name
	))))
colnames(uid)[1] <- "name"
un<-nrow(uid)

## oligo metadata from names
#HAR/other names
aID <- apply(uid,1,function(a) unlist(strsplit(as.character(a[1]),split="__"))[1] )
#species (1=human, 2=chimp, 0=control)
species <- rep(0,un)
species[grep("species=human",as.character(uid$name))] <- 1
species[grep("species=chimp",as.character(uid$name))] <- 2
#cut (i.e., oligo number for that HAR)
cut <- apply(uid,1,function(a) unlist(strsplit(unlist(strsplit(as.character(a[1]),split="cut="))[2],split="_"))[1])
cut[is.na(cut)] <- 0
#permutations
perm <- rep(0,un)
perm[grep(as.character(uid$name),pattern="permutation")] <- 1
#library 2
lib2 <- rep(0,un)
lib2[grep(as.character(uid$name),pattern="HARE5")] <- 1
lib2[perm==1] <- 1
#controls
pos <- neg <- control <- rep(0,un)
control[grep("CONTROL",as.character(uid$name))] <-1
pos[grep("POSITIVE_CONTROL",as.character(uid$name))] <- 1
neg[grep("NEGATIVE_CONTROL",as.character(uid$name))] <- 1
vpos<-c("brain_vista_pos_ctrl10__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__", "ENCODE_pos_ctrl5__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__")
vneg<-c("N1_K27me3_neg_ctrl42__cut=1_of_2__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__", "ENCODE_neg_ctrl42__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__")
#compile into a dataframe
info <- data.frame(name=uid$name,aID,species,cut,perm,lib2,control,pos,neg)

##combine oligo-level RNA counts across samples (type=1)
# same for DNA counts (type=2) or number of barcodes (type=3)
mergeSample<-function(info,counts,sampname,type=1){
      un<-nrow(info)
      # match oligo name to all the rows of data with that name in each replicate
      #res1
      x<-match(info$name,counts$res1$name)
      w<-which(!is.na(x))
      # grab counts
      r1<-rep(NA,un)
      r1[w]<-counts$res1[x[w],type]
      #res2
      x<-match(info$name,counts$res2$name)
      w<-which(!is.na(x))
      # grab counts
      r2<-rep(NA,un)
      r2[w]<-counts$res2[x[w],type]
      #res3
      x<-match(info$name,counts$res3$name)
      w<-which(!is.na(x))
      # grab counts
      r3<-rep(NA,un)
      r3[w]<-counts$res3[x[w],type]
      #output
      r<-cbind(r1,r2,r3)
      colnames(r)<-paste(sampname,colnames(r),sep=".")
      return(r)
}
#run on each sample
rna<-dna<-codes<-NULL
rownames(rna)<-rownames(dna)<-rownames(codes)<-info$name
#RNA first
rna<-cbind(rna,mergeSample(info,wtc.npc.l1,"wtc.npc.l1",1))
rna<-cbind(rna,mergeSample(info,pt2a.npc.l1,"pt2a.npc.l1",1))
rna<-cbind(rna,mergeSample(info,wtc.gpc.l1,"wtc.gpc.l1",1))
rna<-cbind(rna,mergeSample(info,wtc.gpc.l2,"wtc.gpc.l2",1))
rna<-cbind(rna,mergeSample(info,pt2a.gpc.l1,"pt2a.gpc.l1",1))
rna<-cbind(rna,mergeSample(info,pt2a.gpc.l2,"pt2a.gpc.l2",1))
rna<-cbind(rna,mergeSample(info,pt5c.npc.l1,"pt5c.npc.l1",1))
rna<-cbind(rna,mergeSample(info,hs1.gpc.l1,"hs1.gpc.l1",1))
rna<-cbind(rna,mergeSample(info,hs1.npc.l1,"hs1.npc.l1",1))
rna<-cbind(rna,mergeSample(info,pt5c.gpc.l1,"pt5c.gpc.l1",1))
#DNA
dna<-cbind(dna,mergeSample(info,wtc.npc.l1,"wtc.npc.l1",2))
dna<-cbind(dna,mergeSample(info,pt2a.npc.l1,"pt2a.npc.l1",2))
dna<-cbind(dna,mergeSample(info,wtc.gpc.l1,"wtc.gpc.l1",2))
dna<-cbind(dna,mergeSample(info,wtc.gpc.l2,"wtc.gpc.l2",2))
dna<-cbind(dna,mergeSample(info,pt2a.gpc.l1,"pt2a.gpc.l1",2))
dna<-cbind(dna,mergeSample(info,pt2a.gpc.l2,"pt2a.gpc.l2",2))
dna<-cbind(dna,mergeSample(info,pt5c.npc.l1,"pt5c.npc.l1",2))
dna<-cbind(dna,mergeSample(info,hs1.gpc.l1,"hs1.gpc.l1",2))
dna<-cbind(dna,mergeSample(info,hs1.npc.l1,"hs1.npc.l1",2))
dna<-cbind(dna,mergeSample(info,pt5c.gpc.l1,"pt5c.gpc.l1",2))
#Barcodes
codes<-cbind(codes,mergeSample(info,wtc.npc.l1,"wtc.npc.l1",3))
codes<-cbind(codes,mergeSample(info,pt2a.npc.l1,"pt2a.npc.l1",3))
codes<-cbind(codes,mergeSample(info,wtc.gpc.l1,"wtc.gpc.l1",3))
codes<-cbind(codes,mergeSample(info,wtc.gpc.l2,"wtc.gpc.l2",3))
codes<-cbind(codes,mergeSample(info,pt2a.gpc.l1,"pt2a.gpc.l1",3))
codes<-cbind(codes,mergeSample(info,pt2a.gpc.l2,"pt2a.gpc.l2",3))
codes<-cbind(codes,mergeSample(info,pt5c.npc.l1,"pt5c.npc.l1",3))
codes<-cbind(codes,mergeSample(info,hs1.gpc.l1,"hs1.gpc.l1",3))
codes<-cbind(codes,mergeSample(info,hs1.npc.l1,"hs1.npc.l1",3))
codes<-cbind(codes,mergeSample(info,pt5c.gpc.l1,"pt5c.gpc.l1",3))


## compile HAR level RNA, DNA and barcode counts
# using human and chimp oligos per HAR, no library 2
uaID <- unique(sort(info$aID))
keep <- 1:length(uaID)*0
for (i in 1:length(uaID)) {
	ws <- which(as.character(info$aID)==uaID[i] & info$lib2==0)
	if (any(info$species[ws]>0)) {keep[i] <- 1}
}
uaID.f <- uaID[which(keep==1)]
uun<-length(uaID.f)
nn<-ncol(rna)
hrna<-hdna<-hcodes<-crna<-cdna<-ccodes<-matrix(0,nr=uun,nc=nn,dimnames=list(as.character(uaID.f),colnames(rna)))
noligos<-matrix(0,nr=uun,nc=2,dimnames=list(as.character(uaID.f),c("human","chimp")))
for (i in 1:uun) {
        #identify all HAR fragments that are human sequence and not a permutation or control
        w <- which(as.character(info$aID)==uaID.f[i] & info$species==1 & info$perm==0)
        noligos[as.character(uaID.f[i]),"human"] <- length(w)
        # get sum of counts (rna, dna, barcodes) over oligos per HAR 	
	hrna[as.character(uaID.f[i]),] <- apply(matrix(rna[w,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
        hdna[as.character(uaID.f[i]),] <- apply(matrix(dna[w,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
        hcodes[as.character(uaID.f[i]),] <- apply(matrix(codes[w,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
        #identify all HAR fragments that are chimp sequence and not a permutation or control
        ww <- which(as.character(info$aID)==uaID.f[i] & info$species==2 & info$perm==0)
        noligos[as.character(uaID.f[i]),"chimp"] <- length(ww)
        # get sum of counts (rna, dna, barcodes) over oligos per HAR 	
	crna[as.character(uaID.f[i]),] <- apply(matrix(rna[ww,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
        cdna[as.character(uaID.f[i]),] <- apply(matrix(dna[ww,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
        ccodes[as.character(uaID.f[i]),] <- apply(matrix(codes[ww,],ncol=nn),2,function(a) if (any(is.finite(a))) {sum(a,na.rm=T)} else {NA})
}
#compile RNA/DNA ratios for validated controls
vposrna<-rna[which(as.character(info$name)%in%vpos),]
vnegrna<-rna[which(as.character(info$name)%in%vneg),]
vposdna<-dna[which(as.character(info$name)%in%vpos),]
vnegdna<-dna[which(as.character(info$name)%in%vneg),]
vposratios<-vposrna/vposdna
vnegratios<-vnegrna/vnegdna
#RNA/DNA ratios for ALL controls
posrna<-rna[which(info$pos==1),]
negrna<-rna[which(info$neg==1),]
posdna<-dna[which(info$pos==1),]
negdna<-dna[which(info$neg==1),]
posratios<-posrna/posdna
negratios<-negrna/negdna

##compute RNA/DNA ratios
#human
hratios<-hrna/hdna
#chimp
cratios<-crna/cdna
#human/chimp RNA/DNA ratio
ratios<-hratios/cratios
#ratios[is.na(ratios)]<-0

##select samples to include for downstream testing
#keep library 1 (l1) and drop l2 samples
useData <- grep("l1",colnames(ratios))
#check that all HARs have DNA reads
summary(hdna[,useData])
summary(cdna[,useData])
#looks good - minimum is hundreds of DNA reads
#human and chimp RNA/DNA ratios for library 1 only
h1<-hratios[,useData]
c1<-cratios[,useData]
#also take log base 2 so ratios are centered at zero
lib1 <- log2(ratios[,useData])

##annotate library 1 samples
#start with info from colnames: cell line, cell stage, and replicate
library(plyr)
IDs <- data.frame(matrix(unlist(strsplit(as.vector(colnames(lib1)),"\\.")),ncol=4,byrow=TRUE))[,-3]
IDs <- rename(IDs, c("X1" = "line", "X2" = "stage", "X4" = "replicate"))
#add cell species
species <- c(rep("hg",3),rep("pt",3),rep("hg",3),rep("pt",6),rep("hg",6),rep("pt",3))
IDs <- data.frame(cbind(species),IDs)
#add age (in weeks and passage) and unique counts from batch effects file (manually edited from Hane's xls)
batch <- read.csv("HAR_batchEffect_Info_4.14.17.csv",header=TRUE)
#double check that order of rows is the same as in IDs
#data.frame(batch[batch$rna,],IDs)[,c("line","line.1","stage","stage.1","replicate","replicate.1")]
#merge in weeks and passage
IDs <- data.frame(IDs,batch[batch$rna,c("passage","weeks","unique")])
#set the baselines (i.e., intercept of the linear model will be for this set of values)
IDs$species <- relevel(IDs$species, ref= "hg")
IDs$stage <- relevel(IDs$stage, ref= "npc")

##plots
setwd("/Users/kpollard/Desktop/hane/")
#colors defined
getColGrad <- function(a,b,x) {
	rgb(
		round((b[1]-a[1])*c(1:x)/x+a[1],2),
		round((b[2]-a[2])*c(1:x)/x+a[2],2),
		round((b[3]-a[3])*c(1:x)/x+a[3],2)
	)
}
whiteColorVector <- c(255,255,255)/255
redColorVector <- c(255,0,0)/255
cvWR <- getColGrad(whiteColorVector,redColorVector,100)
orangeColorVector <- c(255,155,0)/255
blueColorVector <- c(0,25.5,255)/255
cvBW <- getColGrad(blueColorVector,whiteColorVector,100)
cvWO <- getColGrad(whiteColorVector,orangeColorVector,100)
cvBWO <- c(cvBW,cvWO)
#qqplots of replicates with controls
pdf("replicates.pdf")
par(mfrow=c(2,3))
for(x in unique(IDs$combo)){
    hdat<-hratios[,useData][,which(IDs$combo==x)]
    cdat<-cratios[,useData][,which(IDs$combo==x)]
    ndat<-vnegratios[,useData][,which(IDs$combo==x)]
    nndat<-negratios[,useData][,which(IDs$combo==x)]
    pdat<-vposratios[,useData][,which(IDs$combo==x)]
    ppdat<-posratios[,useData][,which(IDs$combo==x)]
    subIDs<-IDs[which(IDs$combo==x),]
    for(l in unique(subIDs$line)){
        #cat(x,l,"\n")
        hsubdat<-hdat[,subIDs$line==l]
        csubdat<-cdat[,subIDs$line==l]
        nsubdat<-ndat[,subIDs$line==l]
        nnsubdat<-nndat[,subIDs$line==l]
        psubdat<-pdat[,subIDs$line==l]
        ppsubdat<-ppdat[,subIDs$line==l]
        plot(0:8,0:8,type="n",main=paste(x,"-",l),xlab="r1",ylab="r2")
        points(hsubdat[,1],hsubdat[,2],col="blue",pch=".")
        points(csubdat[,1],csubdat[,2],col="orange",pch=".")
        points(ppsubdat[,1],ppsubdat[,2],col="light green",pch="+")
        points(nnsubdat[,1],nnsubdat[,2],col="pink",pch="-")
        points(psubdat[,1],psubdat[,2],col="dark green",pch="+")
        points(nsubdat[,1],nsubdat[,2],col="red",pch="-")
        points(quantile(nnsubdat[,1],0.75),quantile(nnsubdat[,2],0.75),col="purple",pch=2)
        plot(0:8,0:8,type="n",main=paste(x,"-",l),xlab="r1",ylab="r3")
        points(hsubdat[,1],hsubdat[,3],col="blue",pch=".")
        points(csubdat[,1],csubdat[,3],col="orange",pch=".")
        points(ppsubdat[,1],ppsubdat[,3],col="light green",pch="+")
        points(nnsubdat[,1],nnsubdat[,3],col="pink",pch="-")
        points(psubdat[,1],psubdat[,3],col="dark green",pch="+")
        points(nsubdat[,1],nsubdat[,3],col="red",pch="-")
        points(quantile(nnsubdat[,1],0.75),quantile(nnsubdat[,3],0.75),col="purple",pch=2)
        plot(0:8,0:8,type="n",main=paste(x,"-",l),xlab="r2",ylab="r3")
        points(hsubdat[,2],hsubdat[,3],col="blue",pch=".")
        points(csubdat[,2],csubdat[,3],col="orange",pch=".")
        points(ppsubdat[,2],ppsubdat[,3],col="light green",pch="+")
        points(nnsubdat[,2],nnsubdat[,3],col="pink",pch="-")
        points(psubdat[,2],psubdat[,3],col="dark green",pch="+")
        points(nsubdat[,2],nsubdat[,3],col="red",pch="-")
        points(quantile(nnsubdat[,2],0.75),quantile(nnsubdat[,3],0.75),col="purple",pch=2)
    }
}
dev.off()
#human and chimp RNA/DNA for dynamic range
pdf("rna_dna.pdf")
heatmap(log2(hratios),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Human log2(RNA/DNA) - Lib2 Included",labRow=NA)
heatmap(log2(hratios[,useData]),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Human log2(RNA/DNA) - No Lib2",labRow=NA)
heatmap(log2(cratios),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Chimp log2(RNA/DNA) - Lib2 Included",labRow=NA)
heatmap(log2(cratios[,useData]),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Chimp log2(RNA/DNA) - No Lib2",labRow=NA)
boxplot(hratios,main="Human RNA/DNA",pars=list(oma=c(1,3,1,1),las=2,cex.axis=0.5),horizontal=TRUE)
boxplot(cratios,main="Chimp RNA/DNA",pars=list(oma=c(1,3,1,1),las=2,cex.axis=0.5),horizontal=TRUE)
dev.off()
#human/chimp ratio of RNA/DNA values
pdf("ratios.pdf")
heatmap(log2(ratios),col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Library 2 Included",labRow=NA)
heatmap(lib1,col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Library 1 Only",labRow=NA)
boxplot(log2(ratios),main="log2(Human RNA/DNA:Chimp RNA/DNA",pars=list(oma=c(1,3,1,1),las=2,cex.axis=0.5),horizontal=TRUE)
par(mfrow=c(2,3))
for(i in 1:30){hist(ratios[,i],main=colnames(ratios)[i])}
dev.off()
#number of barcodes
pdf("barcodes.pdf")
heatmap(hcodes,col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Human - with lib2",labRow=NA)
heatmap(ccodes,col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Chimp - with lib2",labRow=NA)
heatmap(hcodes-ccodes,col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Human-Chimp (with lib2)",labRow=NA)
dev.off()
#number of oligos
pdf("noligos.pdf")
heatmap(noligos,col=cvBWO,scale="none",Colv=NA,Rowv=NA,main="Number of oligos",labRow=NA,labCol=colnames(noligos))
plot(noligos[,"human"],noligos[,"chimp"],main="Number of oligos",xlab="human",ylab="chimp")
dev.off()

#save image
save.image("KPQuant.RData")


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
DiffExprHARs <-cbind(out1,h1means,c1means)

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

##Plot all expressed HARs in a sorted, clustered heatmap by human v chimp seq (z-scores) for all replicates and cell lines
library(gplots)
library(preprocessCore)

#make list of expressed HARs
expressed_HARs <-rownames(DiffExprHARs)

#only take lib1 ratios for both human and chimp RNA/DNA matrices
lib1_hratios <- hratios[,which(!grepl("l2", colnames(hratios)))]
lib1_cratios <- cratios[,which(!grepl("l2", colnames(cratios)))]

#remove/reset outliers with RNA/DNA ratios > 8
#lib1_hratios[which(lib1_hratios>8)] <- 8
#lib1_cratios[which(lib1_cratios>8)] <- 8

#take log2() of the ratios to center around zero
log2_hratios <-log2(lib1_hratios)
log2_cratios <-log2(lib1_cratios)

#quantile normalization
#normalized_human_data <- normalize.quantiles(data.matrix(human_data), copy=TRUE)

#merge normalized human and chimp matrices
HumanandChimp <- merge(human_data, chimp_data, by='row.names', suffixes=c(".h", ".c"))
rownames(HumanandChimp) <- HumanandChimp[,c(1)]
HumanandChimp <- HumanandChimp[,-c(1)]

#filter merged human and chimp matrix by expressedHARs
dfiltered <- HumanandChimp[which(rownames(HumanandChimp) %in% expressed_HARs),]

#remove batch effects in merged human and chimp matrix for all expressed HARs
merged_batch <- c(rep("A",6),rep("B",9),rep("A",3),rep("B",6),rep("A",6),rep("B",9),rep("A",3),rep("B",6))
merged_data <- removeBatchEffect(dfiltered, merged_batch)


##calculate z-scores for all expressed HARs matrix (mean subtraction)
centered_data = scale(merged_data, scale=T) # center rows, mean substracted


#perform spearman correlation for clustering and heatmap plotting for all expressed HARs
cr = cor(centered_data, method='spearman')
gene_dist = dist(centered_data, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
#myheatcol = redgreen(255)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), s=.5, v=.75, start=0.1, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
quantile.range <- quantile(centered_data, probs = seq(0, 1, 0.001))
palette.breaks <- seq(quantile.range["2.0%"], quantile.range["98.0%"], 0.05)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(length(palette.breaks) - 1)
pdf("diffExprHARs_byHandCseq_heatmap.pdf")
heatmap.2(centered_data, symkey=T, symm=T, symbreaks=T, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=FALSE, col = color.palette, breaks= palette.breaks,,RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,4,.5), margins=c(12,5), labRow=FALSE)
dev.off()


##Plot all sig diff expressed HARs in clustered heatmap by human v chimp seq (z-scores) for all replicates and cell lines
library(gplots)
library(preprocessCore)

#make list of sigHARs
sigHARs <- rownames(DiffExprHARs[which(DiffExprHARs$adj.P.Val<0.01),])

#only take lib1 ratios for both human and chimp RNA/DNA matrices
lib1_hratios <- hratios[,which(!grepl("l2", colnames(hratios)))]
lib1_cratios <- cratios[,which(!grepl("l2", colnames(cratios)))]

#take log2() of the ratios to center around zero
log2_hratios <-log2(lib1_hratios)
log2_cratios <-log2(lib1_cratios)

#quantile normalization
#normalized_human_data <- normalize.quantiles(data.matrix(human_data), copy=TRUE)

#merge normalized human and chimp matrices
HumanandChimp <- merge(log2_hratios, log2_cratios, by='row.names', suffixes=c(".h", ".c"))
rownames(HumanandChimp) <- HumanandChimp[,c(1)]
HumanandChimp <- HumanandChimp[,-c(1)]

#filter merged human and chimp matrix by sigHARs
sig_dfiltered <- HumanandChimp[which(rownames(HumanandChimp) %in% sigHARs),]

#remove batch effects in merged human and chimp matrix for all sig diff HARs
merged_batch <- c(rep("A",6),rep("B",9),rep("A",3),rep("B",6),rep("A",6),rep("B",9),rep("A",3),rep("B",6))
sig_merged_data <- removeBatchEffect(sig_dfiltered, merged_batch)

##calculate z-scores for all sig diff HARs matrix (mean subtraction)
sig_centered_data = scale(sig_merged_data, scale=T) # center rows, mean substracted

#perform spearman correlation for clustering and heatmap plotting for all sig diff HARs
sig_cr = cor(sig_centered_data, method='spearman')
sig_gene_dist = dist(sig_centered_data, method='euclidean')
sig_hc_genes = hclust(sig_gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
#myheatcol = redgreen(255)
sig_gene_partition_assignments <- cutree(as.hclust(sig_hc_genes), k=6);
sig_partition_colors = rainbow(length(unique(sig_gene_partition_assignments)), s=.5, v=.75, start=0.1, end=0.95)
sig_gene_colors = partition_colors[sig_gene_partition_assignments]
sig_quantile.range <- quantile(sig_centered_data, probs = seq(0, 1, 0.001))
sig_palette.breaks <- seq(sig_quantile.range["2.0%"], sig_quantile.range["98.0%"], 0.05)
sig_color.palette  <- colorRampPalette(c("blue", "white", "red"))(length(sig_palette.breaks) - 1)
pdf("sig_diffExprHARs_byHandCseq_heatmap.pdf")
heatmap.2(sig_centered_data, symkey=T, symm=T, symbreaks=T, dendrogram='row', Rowv=as.dendrogram(sig_hc_genes), Colv=FALSE, col = sig_color.palette, breaks= sig_palette.breaks,,RowSideColors=sig_gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,4,.5), margins=c(12,5), labRow=FALSE)
dev.off()




