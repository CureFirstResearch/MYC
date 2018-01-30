library(bigrquery)
library(lattice)
library(tidyr)
library(reshape2)
library(ggplot2)
#get gexp values for MYC, MYCL, MYCN for all tumor type samples 
#NOTE: looking at all sample types for now (solid, normal etc.)

billingProject = 'isb-cgc-04-0007' 
querySql <- paste("SELECT trmSampleBarcode as SampleBarcode, Study , Symbol , normalized_count , log(normalized_count+1,2) as log2_count",
                  "FROM (",
                  "(SELECT SUBSTR(SampleBarcode,1,LENGTH(SampleBarcode )-1) as trmSampleBarcode, Study, Symbol , normalized_count",  
                  "FROM `isb-cgc-01-0008.Annotated.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_24Jan2017_annot`",  
                  "WHERE Symbol IN ('MYC','MYCL','MYCN')) as gexp",
                  "JOIN (",
                  "SELECT SAMPLE_BARCODE", 
                  "FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist` ) as wl",
                  "ON SAMPLE_BARCODE = trmSampleBarcode)",sep=" ")

mycGexpDF <- query_exec(query = querySql,project = billingProject,use_legacy_sql = F)
#write.table(mycGexpDF,file="/Users/varsha/Documents/SEngine/MYC/MYC_GEXP.tsv",sep = "\t")
#density plots
#mycGexpDF = read.table("/Users/varsha/Documents/SEngine/MYC/MYC_GEXP.tsv",header = TRUE,sep = "\t")
par(mfrow=c(3,1))
plot(density(mycGexpDF$log2_count),main='MYC/MYCL/MYCN\nGEXP Distribution\n(Pan-Cancer)',xlab="") 
plot(density(mycGexpDF[mycGexpDF$Symbol=='MYC',]$log2_count),main="MYC",xlab="")
#plot(density(mycGexpDF[mycGexpDF$Symbol=='MYCL',]$log2_count),main="MYCL",xlab="")
plot(density(mycGexpDF[mycGexpDF$Symbol=='MYCN',]$log2_count),main="MYCN",xlab = "log2(Normalized count+1)")

#per tumor-type boxplots
par(mfrow=c(3,1))
boxplot(mycGexpDF$log2_count~mycGexpDF$Study,main=paste('MYC/MYCL/MYCN\nGEXP Distribution\n(Per Tumor Type)',sep=" "),ylab="",las=2)

for (gene in c('MYC','MYCN'))
{
  
  gexpThisGene = mycGexpDF[mycGexpDF$Symbol==gene,c(1,2,5)]
  boxplot(gexpThisGene$log2_count~gexpThisGene$Study,main=gene,ylab="",las=2)
}

#mean vs. variance (per tumor type)
# querySql = paste("SELECT Symbol ,Study,avg(log2_count) as avgGEXP,variance(log2_count) as varGEXP", 
#                    "FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_WhiteList`",  
#                    "GROUP BY Symbol ,Study",sep = " ")
# resultDF = query_exec(querySql,project = billingProject,useLegacySql = FALSE)
# write.table(resultDF,file="/Users/varsha/Documents/SEngine/MYC/MYC_GEXPMeanVsVariance_TumorSpecific.tsv",sep = "\t")

resultDF = read.table(file = "/Users/varsha/Documents/SEngine/MYC/MYC_GEXPMeanVsVariance_TumorSpecific.tsv",header = TRUE,sep="\t")
xyplot(resultDF$avgGEXP/resultDF$varGEXP~as.factor(resultDF$Study),groups = resultDF$Symbol,auto.key = TRUE,cex=1.5,xlab='Study',ylab='mean(CNV)/var(CNV)',scales = list(x=list(rot=45)))
xyplot(resultDF$avgGEXP/resultDF$varGEXP~as.factor(resultDF$Study),groups = resultDF$Symbol,auto.key = list(space="top",just=0.95),cex=1.5,xlab='Study',ylab='mean(GEXP)/var(GEXP)',pch=c('o','+','*')[as.factor(resultDF$Symbol)],scales = list(x=list(rot=45)))

#correlation between MYC expression and expression of all ~30k genes
#the following query was run to compute spearman and pearson correlation between MYC expression and GEXP of all genes
# 
# SELECT Symbol ,mycGene,CORR(mycExp,geneExp) as pearson_corr, CORR(myc_rank,expr_rank) as spearman_corr, count(*) as N
# FROM (
#   SELECT SampleBarcode ,Symbol, geneExp,
#   RANK() OVER (PARTITION BY Symbol,mycGene ORDER BY geneExp ASC) AS expr_rank,
#   mycGene,mycExp,
#   RANK() OVER (PARTITION BY Symbol,mycGene ORDER BY mycExp ASC) AS myc_rank
#   FROM (
#     SELECT gexp.SampleBarcode ,mycGene,mycExp,Symbol,log2_count as geneExp 
#     FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_WhiteList` as gexp
#     JOIN
#     (SELECT SampleBarcode, Symbol as mycGene,log2_count  as mycExp   
#     FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_WhiteList` 
#     WHERE Symbol = 'MYC') as mycGEXP
#     ON gexp.SampleBarcode = mycGEXP.SampleBarcode
#     WHERE gexp.log2_count IS NOT NULL AND mycExp IS NOT NULL))
# GROUP BY Symbol, mycGene
# HAVING Symbol <> 'MYC' AND (pearson_corr IS NOT NULL OR spearman_corr IS NOT NULL)

####Computed correlation scores are saved to a local file to avoid rerunning BigQueries
querySql = "SELECT * FROM `isb-cgc-04-0007.MYC.MYC_AllGene_ExpCorr`"
resultDF = query_exec(querySql,project = billingProject,useLegacySql = FALSE)
resultDF = resultDF[!is.na(resultDF$Symbol),]  #remove rows where gene is NA/null
write.table(resultDF,file = "/Users/varsha/Documents/SEngine/MYC/MYC_AllGene_ExpCorr.tsv",sep = "\t")

corrDF = read.table("/Users/varsha/Documents/SEngine/MYC/MYC_AllGene_ExpCorr.tsv",header = TRUE)
par(mfrow=c(1,3))
plot(density(corrDF$spearman_corr),main="Spearman Correlation")
plot(density(corrDF$pearson_corr),main="Pearson Correlation")
plot(corrDF$pearson_corr,corrDF$spearman_corr)

#get top 10 genes based on correlation of their expression with MYC GEXP
topGenes_PosCorr = corrDF$Symbol[order(corrDF$spearman_corr,decreasing = TRUE)[1:10]]
topPosCorr = corrDF$spearman_corr[order(corrDF$spearman_corr,decreasing = TRUE)[1:10]]
topGenes_NegCorr = corrDF$Symbol[order(corrDF$spearman_corr,decreasing = FALSE)[1:10]]
topNegCorr = corrDF$spearman_corr[order(corrDF$spearman_corr,decreasing = FALSE)[1:10]]

#get expression values for these genes and for MYC - for plotting
# inList = paste(paste(shQuote(topGenes_PosCorr),collapse = ","),paste(shQuote(topGenes_NegCorr),collapse = ","),"'MYC'",sep = ',') #get exp values for MYC as well
# querySql = paste("SELECT SampleBarcode ,Symbol ,geneExp , mycExp",
#                  "FROM `isb-cgc-04-0007.MYC.MYC_AllGene_ExpRanks` ",
#                  "WHERE Symbol  IN (",inList,")",sep = " ")
# resultDF = query_exec(querySql,project = billingProject,max_pages = Inf,useLegacySql = FALSE)
# write.table(resultDF,file="/Users/varsha/Documents/SEngine/MYC/MYC_GEXP_GEXPForTopCorrelatedGenes.tsv",sep="\t")
#plot with and without log2 transformation
resultDF = read.table(file="/Users/varsha/Documents/SEngine/MYC/MYC_GEXP_GEXPForTopCorrelatedGenes.tsv",sep="\t",header = TRUE)

par(mfrow=c(2,5))
#top 10 positive correlations
for(i in c(1:10))
{
  thisDF = resultDF[resultDF$Symbol==as.vector(topGenes_PosCorr)[i],]
  plot(thisDF$mycExp,thisDF$geneExp,main=as.vector(topGenes_PosCorr)[i],xlab=paste(expression(rho),'=',round(topPosCorr[i],digits = 4)),ylab="Normalized Count")
  lines(lowess(thisDF$mycExp,thisDF$geneExp), col="blue",lwd=3)
}
#top 10 negative correlations
par(mfrow=c(2,5))
for(i in c(1:10))
{
  thisDF = resultDF[resultDF$Symbol==as.vector(topGenes_NegCorr)[i],]
  plot(thisDF$mycExp,thisDF$geneExp,main=as.vector(topGenes_NegCorr)[i],xlab = paste(expression(rho),'=',round(topNegCorr[i],digits = 4)),ylab = "Normalized Count")
  lines(lowess(thisDF$mycExp,thisDF$geneExp), col="red",lwd=3)
}


#TODO:tumor-specific correlation
# QUERY
# SELECT Study,Symbol ,mycGene,CORR(mycExp,geneExp) as pearson_corr, CORR(myc_rank,expr_rank) as spearman_corr, count(*) as N
# FROM (
#   SELECT SampleBarcode ,Study,Symbol, geneExp,
#   RANK() OVER (PARTITION BY Symbol,mycGene,Study  ORDER BY geneExp ASC) AS expr_rank,
#   mycGene,mycExp,
#   RANK() OVER (PARTITION BY Symbol,mycGene,Study ORDER BY mycExp ASC) AS myc_rank
#   FROM (
#     SELECT gexp.SampleBarcode,Study,mycGene,mycExp,Symbol,log2_count as geneExp 
#     FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_WhiteList` as gexp
#     JOIN
#     (SELECT SampleBarcode, Symbol as mycGene,log2_count  as mycExp   
#     FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_WhiteList` 
#     WHERE Symbol = 'MYC') as mycGEXP
#     ON gexp.SampleBarcode = mycGEXP.SampleBarcode
#     WHERE gexp.log2_count IS NOT NULL AND mycExp IS NOT NULL AND Symbol <> '?'))
# GROUP BY Symbol, mycGene, Study 
# HAVING Symbol <> 'MYC' AND (NOT IS_NAN(pearson_corr) OR NOT IS_NAN(spearman_corr))

# querySql = paste("SELECT Study,Symbol,spearman_corr,N",
#                  "FROM `isb-cgc-04-0007.MYC.MYC_AllGene_ExpCorr_PerStudy`",sep = " ")
#  
# resultDF = query_exec(querySql,project = billingProject,max_pages = Inf,useLegacySql = FALSE)
# write.table(resultDF,file = "/Users/varsha/Documents/SEngine/MYC/MYC_AllGene_ExpCorr_PerStudy.tsv",sep="\t")

resultDF = read.table(file="/Users/varsha/Documents/SEngine/MYC/MYC_AllGene_ExpCorr_PerStudy.tsv",header = TRUE,sep="\t")
par(mfrow = c(1,1))
boxplot(resultDF$spearman_corr~resultDF$Study,main=paste('Myc GEXP ~ All GEXP Correlation\n(Per Tumor Type)',sep=" "),ylab="Spearman Correlation",las=2)


#get top 10 positive and negative tumor-specific correlations
topPosCorr = resultDF[order(resultDF$spearman_corr,decreasing = TRUE)[1:10],]
topNegCorr =  resultDF[order(resultDF$spearman_corr,decreasing = FALSE)[1:10],]

#get expression values for these genes and for MYC - for plotting
# 
inList = paste(paste(shQuote(topPosCorr$Symbol),collapse = ","),paste(shQuote(topNegCorr$Symbol),collapse = ","),sep = ',') #get exp values for MYC as well
querySql = paste("SELECT SampleBarcode,Study, Symbol,geneExp,mycExp",
                     "FROM `isb-cgc-04-0007.MYC.MYC_AllGene_ExpRanks_PerStudy`", 
                     "WHERE Symbol IN (",inList,")",sep = " ")
resultDF = query_exec(querySql,project = billingProject,max_pages = Inf,useLegacySql = FALSE)

# CAUTION: The above result contains gexp values for all studies for each gene in the inList. 
# We want only the values for the study in which a gene had high correlation with MYC exp.

#plot top positive correlations
par(mfrow=c(2,5))
for (rn in 1:nrow(topPosCorr))
{
  thisExp = resultDF[resultDF$Study==topPosCorr$Study[rn] & resultDF$Symbol==topPosCorr$Symbol[rn],]
  plot(thisExp$mycExp,thisExp$geneExp,main=paste(topPosCorr$Symbol[rn],'|',topPosCorr$Study[rn]),xlab=paste(expression(rho),'=',round(as.vector(topPosCorr$spearman_corr[rn]),digits = 4)),ylab="log2(normalized_count+1)")
  lines(lowess(thisExp$mycExp,thisExp$geneExp), col="blue",lwd=3)
}

#plot top negative correlations
par(mfrow=c(2,5))
for (rn in 1:nrow(topNegCorr))
{
  thisExp = resultDF[resultDF$Study==topNegCorr$Study[rn] & resultDF$Symbol==topNegCorr$Symbol[rn],]
  plot(thisExp$mycExp,thisExp$geneExp,main=paste(topNegCorr$Symbol[rn],'|',topNegCorr$Study[rn]),xlab=paste(expression(rho),'=',round(as.vector(topNegCorr$spearman_corr[rn]),digits = 4)),ylab="log2(normalized_count+1)")
  lines(lowess(thisExp$mycExp,thisExp$geneExp), col="red",lwd=3)
}

for(study in levels(temp$Study))
{
  thisStudy = temp[temp$Study==study,]
  rankList = thisStudy[order(thisStudy$spearman_corr,decreasing = T),c('Symbol','spearman_corr')]
  write.table(rankList,file=paste("/Users/varsha/Documents/SEngine/MYC/GSEA1/",study,"_MYCExp_Corr.rnk",sep=""),sep="\t",row.names = F,col.names = F)
}


#temp code

studies = unique(exp$Study)
for(thisStudy in studies)
{
  rankedList = exp[exp$Study==thisStudy,c('Symbol','spearman_corr')]
  rankedList = rankedList[order(rankedList$spearman_corr,decreasing = T),]
  write.table(rankedList,file = paste(thisStudy,"_MYCNCorr_GOBP_RankedList.rnk",sep=""),sep = "\t",row.names = F,quote = F,col.names = F)
}

####PMN correlations
pmn = c('MYC','MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')
querySql = paste("SELECT SampleBarcode, Study, Symbol, log2_count", 
                 "FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist`",   
                 "WHERE Symbol IN ('MYC','MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')")
pmnGexpDF <- query_exec(querySql,project = billingProject, max_pages = Inf,use_legacy_sql = F)
pmnGexpDF_Wide = pmnGexpDF[,c('SampleBarcode','Symbol','log2_count')]
pmnGexpDF_Wide = aggregate(log2_count ~ SampleBarcode + Symbol,data = pmnGexpDF_Wide,mean)
pmnGexpDF_Wide = spread(pmnGexpDF_Wide, key = Symbol, value = log2_count)
rownames(pmnGexpDF_Wide) = pmnGexpDF_Wide$SampleBarcode
pmnGexpDF_Wide = subset(pmnGexpDF_Wide,select = -c(SampleBarcode))
pmnGexpDF_Wide[,'Study'] = pmnGexpDF$Study[match(rownames(pmnGexpDF_Wide),pmnGexpDF$SampleBarcode)]

######PMN GEXP clustering - October 17,2017############################
sampleVariance = apply(pmnGexpDF_Wide,1,var)
pmnGexpForHeatmap = pmnGexpDF_Wide
pmnGexpForHeatmap_Top50Percent = pmnGexpForHeatmap[sampleVariance>=quantile(sampleVariance,c(0.25)),]

pmnGexp_LongForDensityPlot = pmnGexpForHeatmap_Top50Percent
pmnGexp_LongForDensityPlot[,'SampleBarcode'] = rownames(pmnGexp_LongForDensityPlot)
pmnGexp_LongForDensityPlot = melt(pmnGexp_LongForDensityPlot,id.vars='SampleBarcode')
colnames(pmnGexp_LongForDensityPlot) = c('SampleBarcode','Gene','log2RPM')
#check gexp distribution per gene - any bimodality?
ggplot(pmnGexp_LongForDensityPlot,aes(log2RPM)) +
  geom_density() +
  facet_wrap(~Gene)

#for every PMN gene, top 25 percent samples get labeled as highly expressed, bottom 25% as low expressed
binGexp <- function(gexpVector)
{
  gexpBinarized = gexpVector
  gexpBinarized[gexpVector>=quantile(gexpVector,c(0.75))] = 1
  gexpBinarized[gexpVector<quantile(gexpVector,c(0.75))] = 0
  gexpBinarized[gexpVector<quantile(gexpVector,c(0.25))] = -1
  return(gexpBinarized)
}

pmnGexpForHeatmap_Top50Percent_Binarized = apply(pmnGexpForHeatmap_Top50Percent,2,function(x) binGexp(x))

myColors = c("blue","white","red")
#heirarchical clustering
pheatmap(as.matrix(pmnGexpForHeatmap_Top50Percent_Binarized),cluster_rows = T,cluster_cols = T,color = myColors,scale = "none",show_rownames = F,legend = T,legend_breaks = c(-1,0,1),legend_labels = c("Low","","High"),fontsize = 8,fontsize_col = 12,main = "Gene Expression in Proximal MYC Network\n(Binarized values, Heirarchical Clustering)")
#ordering wrt binarized MYC
pheatmap(as.matrix(pmnGexpForHeatmap_Top50Percent_Binarized[order(pmnGexpForHeatmap_Top50Percent_Binarized[,'MYC'],decreasing = T),]),cluster_rows = F,cluster_cols = F,color = myColors,scale = "none",show_rownames = F,legend = T,legend_breaks = c(-1,0,1),legend_labels = c("Low","","High"),fontsize = 8,fontsize_col = 12,main = "Gene Expression in Proximal MYC Network\n(Binarized values, MYC ordered)")
#ordering wrt continuous MYC
pheatmap(as.matrix(pmnGexpForHeatmap_Top50Percent[order(pmnGexpForHeatmap_Top50Percent[,'MYC'],decreasing = T),]),cluster_rows = F,cluster_cols = F,scale = "none",show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = "Gene Expression in Proximal MYC Network\n(Continuous values, MYC ordered)")
#heirarchical,  continuous 
pheatmap(as.matrix(pmnGexpForHeatmap_Top50Percent),cluster_rows = T,cluster_cols = T,scale = "none",show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = "Gene Expression in Proximal MYC Network\n(Continuous values)")

######################################################
# corPerStudy = data.frame()
# for(thisStudy in unique(pmnGexpDF_Wide$Study))
# {
#   print(thisStudy)
#   thisGexp = pmnGexpDF_Wide[pmnGexpDF_Wide$Study==thisStudy,1:13]
#   thisCor = cor(thisGexp,use="pairwise.complete.obs",method = "spearman")
#   thisCorLong = melt(thisCor)
#   colnames(thisCorLong) = c('Gene1','Gene2','Spearman')
#   thisCorLong[,'Study'] = thisStudy
#   corPerStudy = rbind(corPerStudy,thisCorLong)
# }
# write.table(corPerStudy,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/GEXP_Correlations/PMN_GEXPCorrelation_PerStudy.tsv',sep = "\t",col.names = T,row.names = F)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

corPerStudy = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/GEXP_Correlations/PMN_GEXPCorrelation_PerStudy.tsv',sep = "\t",header = T)
pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/GEXP_Correlations/HeatmapsPerStudy.pdf',onefile = T)
for(thisStudy in unique(pmnGexpDF$Study))
{
  thisStudyCorr = corPerStudy[as.vector(corPerStudy$Study)==thisStudy,c('Gene1','Gene2','Spearman')]
  thisStudyCorr_Wide = spread(thisStudyCorr,key = Gene2,value = Spearman)
  rownames(thisStudyCorr_Wide) = thisStudyCorr_Wide$Gene1
  thisStudyCorr_Wide = subset(thisStudyCorr_Wide,select = -c(Gene1))
  thisStudyCorr_Wide = reorder_cormat(thisStudyCorr_Wide)
  thisStudyCorr_Wide[upper.tri(thisStudyCorr_Wide,diag = T)] = NA
  thisStudyCorr_Long = melt(as.matrix(thisStudyCorr_Wide),na.rm = T)
  colnames(thisStudyCorr_Long) = c('Gene1','Gene2','Spearman')
  thisPlot = ggplot(data = thisStudyCorr_Long, aes(Gene1, Gene2, fill = Spearman)) +
             geom_tile(color = "white") +
             scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
             labs(title = paste(thisStudy,"\nProximal MYC Network")) +         
             theme_minimal() + 
             coord_fixed() +
             geom_text(aes(Gene1, Gene2,label = round(Spearman,2)), color = "black", size = 4) +
             theme(
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              legend.justification = c(1, 0),
              legend.position = c(0.4, 0.7),
              legend.direction = "horizontal",
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
              axis.text.y = element_text(size=12))+
              guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
             
  print(thisPlot)
}
dev.off()


#########Preparing miRNA data for Andrea Ventura#########################
billingProject = 'isb-cgc-04-0007' 
querySql <- paste("SELECT trmSampleBarcode as SampleBarcode, Study , Symbol , normalized_count , log(normalized_count+1,2) as log2_count",
                  "FROM (",
                  "(SELECT SUBSTR(SampleBarcode,1,LENGTH(SampleBarcode )-1) as trmSampleBarcode, Study, Symbol , normalized_count",  
                  "FROM `isb-cgc-01-0008.Annotated.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_24Jan2017_annot`",  
                  "WHERE Symbol IN ('MYC','MYCL','MYCN')) as gexp",
                  "JOIN (",
                  "SELECT SAMPLE_BARCODE", 
                  "FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist` ) as wl",
                  "ON SAMPLE_BARCODE = trmSampleBarcode)",sep=" ")

mycGexpDF <- query_exec(query = querySql,project = billingProject,use_legacy_sql = F)

mirnaEXP = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/miRNA/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',sep = ",",header = T)
mirnaEXP = mirnaEXP[as.vector(mirnaEXP$Correction)=="Corrected",]
rownames(mirnaEXP) = mirnaEXP$Genes
mirnaEXP = subset(mirnaEXP,select = -c(Genes,Correction))
colnames(mirnaEXP) = gsub(pattern = '.',replacement = '-',colnames(mirnaEXP),fixed = T)
colnames(mirnaEXP) = substr(colnames(mirnaEXP),1,15)
mirnaEXP = t(mirnaEXP)

#find overlap between mycGexpDF and mirnaEXP
sampleOverlap = intersect(rownames(mirnaEXP),mycGexpDF$SampleBarcode)

mycGexpDF = mycGexpDF[mycGexpDF$SampleBarcode %in% sampleOverlap & mycGexpDF$Symbol=="MYC",]
mycGexpDF = aggregate(log2_count~SampleBarcode+Study,data = mycGexpDF,mean)
rownames(mycGexpDF) = mycGexpDF$SampleBarcode
mycGexpDF = mycGexpDF[sampleOverlap,]
mirnaEXP = mirnaEXP[sampleOverlap,]
mirnaEXP = as.data.frame(mirnaEXP)
##append myc columns to mirna DF
mirnaEXP[,'MYCLog2Count'] = mycGexpDF$log2_count
mirnaEXP[,'Study'] = mycGexpDF$Study
write.table(mirnaEXP,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/miRNA/miRNAEXP_MYCGEXP.tsv',sep = "\t",col.names = T,row.names = T)
