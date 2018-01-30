library(tidyr) #for spread
library(gplots) #for heatmap.2
library(RColorBrewer)
library(bigrquery)
#library(cluster) #for daisy, gower
#read GSEA summary and plot heatmap of tumor-type vs. pathways 
# (cell value = Normalized enrichment score)

#read excel files to compile summary

gseaDF = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/MYCN_GO_MolecularFunction/MYCN_GSEAResultsSummary_GO_MolecularFunction_PositiveCorrelations.tsv',header = T,sep = "\t")
#gseaDF$NES = as.numeric(levels(gseaDF$NES))[gseaDF$NES]
#scatter plot NES vs. FDR
par(mfrow=c(2,2))
plot(density(gseaDF$NES,na.rm = T),main="Normalized Enrichment Score")
plot(density(gseaDF$Size),main="Gene set sizes (# of genes)")
plot(gseaDF$NES,gseaDF$Size,col="blue",main="Normalized Enrichment Score\nvs.\nGene set size")
plot(gseaDF$NES,gseaDF$qvalue,col="blue",main="Normalized Enrichment Score\nvs.\nFDR q-value",xlab="Normalized Enrichment Score",ylab="FDR q-value")


#boxplot NES per tumor-type 
par(mfrow=c(1,1))
bxp = boxplot(gseaDF$NES~gseaDF$Study,main="Distribution of NES (GO Biological Processes)\n(per tumor type)",las=2)

gseaNESWide = gseaDF[,c('Study','Name','NES')]
gseaNESWide = spread(data=gseaNESWide,Name,value = NES,fill=NA,drop=F,sep=NULL)
rownames(gseaNESWide) = gseaNESWide[[1]] #study names become row names
gseaNESWide = gseaNESWide[,c(-1)] #remove redundant study name column

gseaQvalWide = gseaDF[,c('Study','Name','qvalue')]
gseaQvalWide = spread(data=gseaQvalWide,Name,value = qvalue,fill=NA,drop=F,sep=NULL)
rownames(gseaQvalWide) = gseaQvalWide[[1]] #study names become row names
gseaQvalWide = gseaQvalWide[,c(-1)] #remove redundant study name column

#intersect with gene sets retained from reduce_overlap
reducedOverlapGeneSets = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/MYCN_GO_MolecularFunction/ReducedOverlap/MYCN_ReducedOverlap_GeneSets.tsv',sep = "\t",header = T)
gseaNESWide = gseaNESWide[,intersect(colnames(gseaNESWide),reducedOverlapGeneSets$term)]
gseaReducedOverlapQValue = gseaQvalWide[,intersect(colnames(gseaQvalWide),reducedOverlapGeneSets$term)]
write.table(t(gseaReducedOverlapQValue),file='/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYC/MYC_GO_MolecularFunction/ReducedOverlap/MYC_ReducedOverlap_QValues_PerGeneSetPerStudy.tsv',row.names = T,col.names = T,sep = "\t")
#All tumor-types have atleast one NA gene set..so cannot filter along that dimension

#show gene sets with top 20 mean enrichment score; show NAs as gray cells
#first, order gene sets by mean and median separately
#!!!NOTE: sortign order changes for positively and negatively correlated gene sets

#for positive correlations
#gseaNESWide_OrderedByMeanNES = gseaNESWide[,order(colMeans(gseaNESWide,na.rm = T),decreasing = T)]
gseaNESWide_OrderedByMedianNES = gseaNESWide[,order(apply(gseaNESWide,2,median,na.rm=T),decreasing = T)]


#for negative correlations
#gseaNESWide_OrderedByMeanNES = gseaNESWide[,order(colMeans(gseaNESWide,na.rm = T),decreasing = F)]
#gseaNESWide_OrderedByMedianNES = gseaNESWide[,order(apply(gseaNESWide,2,median,na.rm=T),decreasing = F)]

#normalize mean ES to range 0-1
minMedianNES = min(gseaNESWide_OrderedByMedianNES[!is.na(gseaNESWide_OrderedByMedianNES)])
maxMedianNES = max(gseaNESWide_OrderedByMedianNES[!is.na(gseaNESWide_OrderedByMedianNES)])
gseaNESWide_OrderedByMedianNES = (gseaNESWide_OrderedByMedianNES - minMedianNES )/ (maxMedianNES - minMedianNES)


#heatmap with top 25 gene sets ordered by mean NES
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
#strip 'GO_' prefix from gene set names 
colnames(gseaNESWide_OrderedByMedianNES) = substring(colnames(gseaNESWide_OrderedByMedianNES),first = 4)
###replace NAS with -1
gseaNESWide_OrderedByMedianNES[is.na(gseaNESWide_OrderedByMedianNES)]=-1
#cluster both GS as well as tumor types
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4)
lwid = c(0.5,4,2)
#with row labels
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMedianNES[,1:100])),cexRow=0.7,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "MYCN\nPositively Correlated Gene Sets\n(Median Normalized Enrichment Score)",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.2,key.title = 'Score')
#without row labels
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMedianNES[,1:100])),cexRow=0.7,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "MYC\nPositively Correlated Gene Sets\n(Median Normalized Enrichment Score)",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.2,key.title = 'Score',labRow = "")

#row and column clustering for MYCN
mycnColors = c("gray",colorRampPalette(colors = c("white","red"))(1000))
mycnBreaks = c(-1,seq(0,1,length.out = 1001))
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMedianNES[,1:100])),cexRow=0.7,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "MYCN\nPositively Correlated Gene Sets\n(Median Normalized Enrichment Score)",na.color = "gray",symbreaks = F,density.info="none",keysize = 0.2,key.title = 'Score',breaks = mycnBreaks,col=mycnColors,scale = "none")
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMedianNES[,1:100])),cexRow=0.7,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "MYCN\nPositively Correlated Gene Sets\n(Median Normalized Enrichment Score)",na.color = "gray",symbreaks = F,density.info="none",keysize = 0.2,key.title = 'Score',breaks = mycnBreaks,col=mycnColors,scale = "none",labRow = "")

#cluster just tumor types, gene sets ordered by mean enrichment score top to bottom
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4)
lwid = c(0.2,4,2)
tumorClust = hclust(dist(gseaNESWide_OrderedByMedianNES))
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMedianNES[,1:100])),cexRow=0.8,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,Rowv = FALSE,Colv = as.dendrogram(tumorClust),dendrogram = "column",trace = "none",main = "MYCN\nPositively Correlated Gene Sets\n(Median Normalized Enrichment Score)",na.color = "gray",symbreaks = F,col=scaleblue,density.info="none",keysize = 0.4,key.title = "Score")

#################################################################
###MYC ~ miRNA EXPRESSION CORRELATION############################
#################################################################
billingProject = 'isb-cgc-04-0007' 
mirnaEXP = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/miRNA/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',sep = ",",header = T)
mirnaEXP = mirnaEXP[mirnaEXP$Correction=='Corrected',]
rownames(mirnaEXP) = mirnaEXP$Genes
mirnaEXP = subset(mirnaEXP,select=-c(Genes,Correction))
colnames(mirnaEXP) = gsub('.','-',substr(colnames(mirnaEXP),1,15),fixed = T)
mirnaEXP = t(mirnaEXP)
mirnaLog2EXP = log2(mirnaEXP+1)
#read myc log2 exp
mycEXPQuery = "SELECT SampleBarcode , Study, log2_count  FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol = 'MYC'"
mycEXP = query_exec(mycEXPQuery,project =billingProject,max_pages = Inf,use_legacy_sql = F)

#retain only overlapping samples between the 2 data sets
overlapSamples = intersect(rownames(mirnaLog2EXP),as.vector(mycEXP$SampleBarcode))

mirnaLog2EXP = mirnaLog2EXP[overlapSamples,]
mycLog2EXP = mycEXP[match(overlapSamples,mycEXP$SampleBarcode),]

######pan-can correlations
pairwiseCor = cor(mirnaLog2EXP,mycLog2EXP$log2_count,method = "spearman",use = "pairwise.complete.obs")
pairwiseCorPval = apply(mirnaLog2EXP,2,function(x) (cor.test(x,mycLog2EXP$log2_count,method = "spearman",use="pairwise.complete.obs",exact = F))$p.value)

plot(density(-log10(pairwiseCorPval)))

#######Study specific correlations
pairwiseCorPerStudy = data.frame()
pairwiseCorPvalPerStudy = data.frame()
for(thisStudy in unique(mycLog2EXP$Study))
{
  thisStudyMYC = mycLog2EXP$log2_count[mycLog2EXP$Study == thisStudy] 
  thisStudyMIR = mirnaLog2EXP[mycLog2EXP$Study == thisStudy,]
  
  thisCor = cor(thisStudyMIR,thisStudyMYC,method = "spearman",use = "pairwise.complete.obs")
  pairwiseCorPerStudy = rbind(pairwiseCorPerStudy,t(thisCor[,1])) 
  rownames(pairwiseCorPerStudy)[nrow(pairwiseCorPerStudy)] = thisStudy
  
  thisCorPval = apply(thisStudyMIR,2,function(x) (cor.test(x,thisStudyMYC,method = "spearman",use="pairwise.complete.obs",exact = F))$p.value)
  pairwiseCorPvalPerStudy = rbind(pairwiseCorPvalPerStudy,as.vector(thisCorPval))
  rownames(pairwiseCorPvalPerStudy)[nrow(pairwiseCorPvalPerStudy)] = thisStudy
}
colnames(pairwiseCorPvalPerStudy) = colnames(mirnaLog2EXP)
#median pvalue per miRNA
medianPval = apply(pairwiseCorPvalPerStudy,2,median,na.rm=T)
#heatmap with 100 miRNAs with least median p-values
least100MedianPval = pairwiseCorPvalPerStudy[,order(medianPval,decreasing = F)[1:100]]
least100MedianPval = least100MedianPval[order(rownames(least100MedianPval)),]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.9,4.7)
lwid = c(0.5,4,1.5)
colorScaleWhiteRed = colorRampPalette(colors = c("white","red"))(1000)
heatmap.2(-log10(as.matrix(t(least100MedianPval))),density="none",col=colorScaleWhiteRed,Colv=F,Rowv=F,dendrogram = "none",trace = "none",key.title="-log10(p-value)",main="Correlation between MYC expression\nand\nmiRNA expression",lmat = lmat,lhei = lhei,lwid = lwid,cexRow = 0.8)

######for tumor-type specific enrichment analysis, write sorted miRNA lists
###to separate tumor-specific files
for (thisStudy in rownames(pairwiseCorPvalPerStudy))
{
  sortedMIR = rownames(t(pairwiseCorPvalPerStudy[thisStudy,order(pairwiseCorPvalPerStudy[thisStudy,])]))
  write.table(sortedMIR,paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/miRNA/",thisStudy,"_miRNASortedByPValue.tsv",sep = ""),sep = "\n",row.names = F,col.names = F)
}


##########################GSEA Reduce Overlap########################################
library(GOplot)
library(tidyr)
gsPVal = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/MYCN_GO_MolecularFunction/MYCN_GSEAResultsSummary_GO_MolecularFunction_PositiveCorrelations.tsv',sep = "\t",header = T)
gsPVal = gsPVal[,c('Study','Name','qvalue')]
gsPVal = spread(gsPVal,key=Name,value = qvalue)
rownames(gsPVal) = gsPVal$Study
gsPVal = subset(gsPVal,select=-c(Study))
gsPVal['MedianPval',] = apply(gsPVal,2,median,na.rm=T)
gsPVal = as.data.frame(t(gsPVal))
gsPVal[,'Term'] = rownames(gsPVal)
gseaResults = gsPVal[,c('Term','MedianPval')]
gseaResults[,'Category'] = 'MF'

#get gene list per gene set
geneList = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/c5.mf.v6.0.GeneList.tsv',sep = "\t",header = F)  
gseaResults[,'Genes'] = geneList$V2[match(gseaResults$Term,geneList$V1)]  
rownames(gseaResults) = NULL 
##Note: some gene lists are NA 
gseaResults = gseaResults[,c('Category','Term','Genes','MedianPval')]
colnames(gseaResults) = c('category','term','genes','adj_pval')  


#read ranked gene list per tumor type
#compute median correlation score per gene across all tumor types
rnkFilenames = list.files('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/',pattern = "*rnk",recursive = F,include.dirs = F,full.names = F)
corPerGenePerStudy = data.frame()
for(filename in rnkFilenames)
{
  thisTumorGeneList = read.table(paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/",filename,sep = ""),sep = "\t",header = F)
  if(ncol(corPerGenePerStudy)==0)
  {
    corPerGenePerStudy = thisTumorGeneList
    rownames(corPerGenePerStudy) = corPerGenePerStudy[,'V1']
    #corPerGenePerStudy = corPerGenePerStudy[,-c(1)]
  }
  else
  {
    corPerGenePerStudy[,ncol(corPerGenePerStudy)+1] = thisTumorGeneList[match(rownames(corPerGenePerStudy),thisTumorGeneList[,'V1']),'V2']
  }
}
corPerGenePerStudy = corPerGenePerStudy[,-c(1)]
corPerGenePerStudy[,'MedianCor'] = apply(corPerGenePerStudy,1,median,na.rm=T)
medianCorPerGene = corPerGenePerStudy
medianCorPerGene[,'ID'] = rownames(medianCorPerGene)
medianCorPerGene = medianCorPerGene[,c('ID','MedianCor')]
colnames(medianCorPerGene) = c('ID','logFC')
rownames(medianCorPerGene) = NULL

circ = circle_dat(gseaResults,medianCorPerGene)
reducedCirc = reduce_overlap(circ,overlap = 0.75)
write.table(reducedCirc,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/GSEA1/MYCN/MYCN_GO_MolecularFunction/MYCN_ReducedOverlap_GeneSets.tsv',sep = "\t",col.names = T,row.names = F)
