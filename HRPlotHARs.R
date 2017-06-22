##Goal: plot human and chimp seq enhancer scores for all expressed HARs in all cell lines

##setup
#set directory
setwd("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/data/")
load("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/data/KPQuant.RData")
load("/Users/haneryu/Desktop/Dropbox/HAR_MPRA_Hane_March2017/data/KPTest.RData")


##Plot all sig diff expressed HARs in clustered heatmap by human v chimp seq (z-scores) for all replicates and cell lines
library(gplots)
library(preprocessCore)

#make list of sigHARs
sigHARs <- rownames(DiffExprHARs[which(DiffExprHARs$adj.P.Val<0.01),])

#only take lib1 ratios for both (human RNA/DNA)/(chimp RNA/DNA) matrices
ratios <- ratios[,which(!grepl("l2", colnames(ratios)))]

#take log2() of the ratios to center around zero
log2_ratios <-log2(ratios)

#filter log2(human/chimp) matrix by sigHARs
sig_dfiltered <- log2_ratios[which(rownames(log2_ratios) %in% sigHARs),]

#remove batch effects in merged human and chimp matrix for all sig diff HARs
#merged_batch <- c(rep("A",6),rep("B",9),rep("A",3),rep("B",6))
#sig_merged_data <- removeBatchEffect(sig_dfiltered, merged_batch)

##calculate z-scores for all sig diff HARs matrix (mean subtraction)
sig_centered_data = scale(sig_dfiltered, scale=T) # center rows, mean substracted

#perform spearman correlation for clustering and heatmap plotting for all sig diff HARs
sig_cr = cor(sig_centered_data, method='spearman')
sig_gene_dist = dist(sig_centered_data, method='euclidean')
sig_hc_genes = hclust(sig_gene_dist, method='complete')
hc_samples = hclust(as.dist(1-sig_cr), method="complete") # cluster conditions
#color bar on side
sig_gene_partition_assignments <- cutree(as.hclust(sig_hc_genes), k=6);
sig_partition_colors = rainbow(length(unique(sig_gene_partition_assignments)), s=.5, v=.75, start=0.1, end=0.95)
sig_gene_colors = partition_colors[sig_gene_partition_assignments]

sig_quantile.range <- quantile(sig_centered_data, probs = seq(0, 1, 0.001))
sig_palette.breaks <- seq(sig_quantile.range["2.0%"], sig_quantile.range["98.0%"], 0.05)
sig_color.palette  <- colorRampPalette(c("blue", "white", "red"))(length(sig_palette.breaks) - 1)
pdf("sig_diffExprHARs_byHandCseq_heatmap.pdf")
heatmap.2(sig_centered_data, symkey=T, symm=T, symbreaks=T, dendrogram='row', Rowv=as.dendrogram(sig_hc_genes), Colv=FALSE, col = sig_color.palette, breaks= sig_palette.breaks,,RowSideColors=sig_gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,4,.5), margins=c(12,5), labRow=FALSE)
dev.off()


##Plot all expressed HARs in a sorted, clustered heatmap by human v chimp seq (z-scores) for all replicates and cell lines
library(gplots)
library(preprocessCore)

#make list of expressed HARs
expressed_HARs <-rownames(DiffExprHARs)

#only take lib1 ratios for both human and chimp RNA/DNA matrices
lib1_hratios <- hratios[,which(!grepl("l2", colnames(hratios)))]
lib1_cratios <- cratios[,which(!grepl("l2", colnames(cratios)))]
#try including negative controls in the data plotted in heatmap?


#take log2() of the ratios to center around zero
log2_hratios <-log2(lib1_hratios)
log2_cratios <-log2(lib1_cratios)

#quantile normalization
#normalized_human_data <- normalize.quantiles(data.matrix(human_data), copy=TRUE)

#merge normalized human and chimp matrices
HumanandChimp <- merge(log2_hratios, log2_cratios, by='row.names', suffixes=c(".h", ".c"))
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
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), s=.5, v=.75, start=0.1, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
quantile.range <- quantile(centered_data, probs = seq(0, 1, 0.001))
palette.breaks <- seq(quantile.range["2.0%"], quantile.range["98.0%"], 0.05)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(length(palette.breaks) - 1)
pdf("diffExprHARs_byHandCseq_heatmap.pdf")
heatmap.2(centered_data, symkey=T, symm=T, symbreaks=T, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=FALSE, col = color.palette, breaks= palette.breaks,,RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,4,.5), margins=c(12,5), labRow=FALSE)
dev.off()

