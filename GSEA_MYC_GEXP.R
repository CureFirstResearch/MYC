library(tidyr) #for spread
library(gplots) #for heatmap.2
library(RColorBrewer)
#library(cluster) #for daisy, gower
#read GSEA summary and plot heatmap of tumor-type vs. pathways 
# (cell value = Normalized enrichment score)

gseaDF = read.table('/Users/varsha/Documents/SEngine/MYC/GSEA1/GSEAResultsSummary_AllTumorTypes.tsv',header = T,sep = "\t")

#scatter plot NES vs. FDR
par(mfrow=c(2,2))
plot(density(gseaDF$NES),main="Normalized Enrichment Score")
plot(density(gseaDF$Size),main="Gene set sizes (# of genes)")
plot(gseaDF$NES,gseaDF$Size,col="blue",main="Normalized Enrichment Score\nvs.\nGene set size")
plot(gseaDF$NES,gseaDF$qValue,col="blue",main="Normalized Enrichment Score\nvs.\nFDR q-value",xlab="Normalized Enrichment Score",ylab="FDR q-value")


#boxplot NES per tumor-type 
par(mfrow=c(1,1))
boxplot(gseaDF$NES~gseaDF$Study,main="Distribution of NES\n(per tumor type)",las=2)

gseaNES = gseaDF[,c('Study','GSName','NES')]
gseaNESWide = spread(data=gseaNES,GSName,value = NES,fill=NA,drop=F,sep=NULL)
rownames(gseaNESWide) = gseaNESWide[[1]] #study names become row names
gseaNESWide = gseaNESWide[,c(-1)] #remove redundant study name column

#eliminate gene sets that have NA values for any tumor-type
#only 8 gene sets with no NAs!!!
gseaFilteredGeneSets = gseaNESWide[,colSums(is.na(gseaNESWide))==0]
#heatmap(as.matrix(gseaFilteredGeneSets))

#All tumor-types have atleast one NA gene set..so cannot filter along that dimension

#show gene sets with top 20 mean enrichment score; show NAs as gray cells
#first, order gene sets by mean and median separately
gseaNESWide_OrderedByMeanNES = gseaNESWide[,order(colMeans(gseaNESWide,na.rm = T),decreasing = T)]
gseaNESWide_OrderedByMedianNES = gseaNESWide[,order(apply(gseaNESWide,2,median,na.rm=T),decreasing = T)]


#heatmap with top 50 gene sets ordered by mean NES
scaleblue<-colorRampPalette(colors=c("white","blue"))(1000)
#strip 'GO_' prefix from gene set names 
colnames(gseaNESWide_OrderedByMeanNES) = substring(colnames(gseaNESWide_OrderedByMeanNES),first = 4)
#cluster both GS as well as tumor types
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1.5,4)
lwid = c(0.5,4,2)
#heatmap.2(as.matrix(t(gseaNESWide_OrderedByMeanNES[,1:25])),lmat = lmat,trace = "none",main = "Positively Correlated Gene Sets\n(NormalizedEnrichmentScore)",na.color = "gray",cexRow = 0.5,symbreaks = F,col=scaleblue,density.info="none",keysize = 0.8,key.title = "Score")

#cluster just tumor types
heatmap.2(as.matrix(t(gseaNESWide_OrderedByMeanNES[,1:25])),lmat = lmat,lhei = lhei,lwid = lwid,Rowv = FALSE,Colv = TRUE,dendrogram = "column",trace = "none",main = "Positively Correlated Gene Sets\n(Normalized Enrichment Score)",na.color = "gray",cexRow = 0.8,symbreaks = F,col=scaleblue,density.info="none",keysize = 0.6,key.title = "Score")
