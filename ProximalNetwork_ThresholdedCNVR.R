#re-analyze with relative copy numbers and thresholded AMP/DEL status - provided by Andy Cherniack

library(reshape) #for melt
library(ggplot2)
library(bigrquery)
library(tidyr) #for spread
library(gplots) #for heatmap.2
library(effsize) #for cohen.d
billingProject = 'isb-cgc-04-0007'

#######read thresholded AMP/DEL status
pancanThresholdedCNVR = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/ThresholdedCopyNumber/syn5049520_all_thresholded.by_genes_whitelisted.tsv',sep = "\t",header = T)
#retain only myc netwrok genes
mycNW_GeneList = c('MYC','MYCN','MYCL','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')
pancanThresholdedCNVR = pancanThresholdedCNVR[pancanThresholdedCNVR$Gene.Symbol %in% mycNW_GeneList,]
#remove locus ID and cytoband
pancanThresholdedCNVR = pancanThresholdedCNVR[,-c(2,3)]
#assign row names
rownames(pancanThresholdedCNVR) = pancanThresholdedCNVR[,'Gene.Symbol']
pancanThresholdedCNVR = pancanThresholdedCNVR[,-c(1)]
#trim colnames to be sample barcodes instead of full aliquot barcodes
colnames(pancanThresholdedCNVR) = substring(colnames(pancanThresholdedCNVR),1,15)
colnames(pancanThresholdedCNVR) = gsub('.','-',colnames(pancanThresholdedCNVR),fixed = T)
#transpose so that samples are along rows
pancanThresholdedCNVR = as.data.frame(t(pancanThresholdedCNVR))


########read pancan pathway whitelist 
sampleWhitelist = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/PanCanPathwayWhitelist.tsv',sep="\t",header = T)
sampleCountPerStudy = as.data.frame(table(sampleWhitelist$DISEASE))
names(sampleCountPerStudy) = c('Study','Count')

########filter to keep only whitelisted samples
pancanThresholdedCNVR = pancanThresholdedCNVR[as.vector(sampleWhitelist$SAMPLE_BARCODE),]

########plot distribution per gene 
#first, compute counts per category within a gene
countsPerGene = apply(pancanThresholdedCNVR, 2, FUN = function(x) table(x))
#convert into percentage of samples within each gene
percsPerGene = as.data.frame(countsPerGene*100/colSums(countsPerGene))
percsPerGene[,'Category'] = rownames(percsPerGene)
percsPerGene$Category = factor(percsPerGene$Category,levels = c("2","1","0","-1","-2"))
#next, melt DF into long form
percsPerGene_Long = melt(percsPerGene,id = c('Category'))
names(percsPerGene_Long) = c('ThresholdedCopyNumber','Gene','PercentSamples')

#######stacked bar plot
#TODO - blue to red color scale
ggplot(data=percsPerGene_Long, aes(x=Gene, y=PercentSamples, fill=ThresholdedCopyNumber)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))


#######stacked bar plots, facet_wrapped by tumor-type
#first, add Study as a column
pancanThresholdedCNVR[,'Study'] = sampleWhitelist[match(rownames(pancanThresholdedCNVR),sampleWhitelist$SAMPLE_BARCODE),'DISEASE']
#melt into long format
pancanThresholdedCNVR_Long = melt(pancanThresholdedCNVR,id=c('Study'))
names(pancanThresholdedCNVR_Long) = c('Study','Gene','ThresholdedCopyNumber')
#compute counts per Study, per gene, per copy number
countsPerGenePerStudy = aggregate(rep(1,nrow(pancanThresholdedCNVR_Long)),by=pancanThresholdedCNVR_Long[,c('Study','Gene','ThresholdedCopyNumber')],sum)
names(countsPerGenePerStudy) = c('Study','Gene','ThresholdedCopyNumber','Count')
countsPerGenePerStudy[,'TotalNThisStudy'] = sampleCountPerStudy[match(countsPerGenePerStudy$Study,sampleCountPerStudy$Study),'Count']
countsPerGenePerStudy[,'PercentSamples'] = countsPerGenePerStudy$Count*100/countsPerGenePerStudy$TotalNThisStudy
#convert copy number to string for clearer colors in ggplot
countsPerGenePerStudy$ThresholdedCopyNumber = factor(countsPerGenePerStudy$ThresholdedCopyNumber,levels = c(2,1,0,-1,-2))
ggplot(data=countsPerGenePerStudy, aes(x=Gene, y=PercentSamples, fill=ThresholdedCopyNumber)) +
  geom_bar(stat="identity") +
  facet_wrap(~Study) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5.5)) 
  
####heatmap of percent of sample amplified per network gene per tumor-type
#for first heatmap, look at both 1 and 2 copy gains
pancanBinaryAMP = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$ThresholdedCopyNumber)) > 0,c('Study','Gene','Count')]
pancanBinaryAMP = aggregate(pancanBinaryAMP$Count,by = pancanBinaryAMP[,c('Study','Gene')],sum)
names(pancanBinaryAMP) = c('Study','Gene','Count')
#spread for heatmap
pancanBinaryAMP = spread(pancanBinaryAMP,Gene,Count)
rownames(pancanBinaryAMP) = pancanBinaryAMP$Study
pancanBinaryAMP = pancanBinaryAMP[,-c(1)]
#Add sample count per Study as column
pancanBinaryAMP[,'TotalN']= sampleCountPerStudy[match(rownames(pancanBinaryAMP),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanBinaryAMP['PanCan',] = colSums(pancanBinaryAMP,na.rm = T)
#compute percentages
pancanBinaryAMP_Perc = pancanBinaryAMP*100/pancanBinaryAMP$TotalN
pancanBinaryAMP_Perc = pancanBinaryAMP_Perc[,-c(ncol(pancanBinaryAMP_Perc))]
pancanBinaryAMP_Perc[is.na(pancanBinaryAMP_Perc)] = 0
pancanBinaryAMP_Perc = as.matrix(pancanBinaryAMP_Perc)
pancanBinaryAMP_Perc = pancanBinaryAMP_Perc[,order(colnames(pancanBinaryAMP_Perc))]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
heatmap.2(pancanBinaryAMP_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryAMP_Perc==0,NA,round(pancanBinaryAMP_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples Amplified\n(1 or more copy gains)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.8)
#for second heatmap, look at only 2 copy gains
pancanBinaryAMP = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$ThresholdedCopyNumber)) > 1,c('Study','Gene','Count')]
pancanBinaryAMP = aggregate(pancanBinaryAMP$Count,by = pancanBinaryAMP[,c('Study','Gene')],sum)
names(pancanBinaryAMP) = c('Study','Gene','Count')
#spread for heatmap
pancanBinaryAMP = spread(pancanBinaryAMP,Gene,Count)
rownames(pancanBinaryAMP) = pancanBinaryAMP$Study
pancanBinaryAMP = pancanBinaryAMP[,-c(1)]
#Add sample count per Study as column
pancanBinaryAMP[,'TotalN']= sampleCountPerStudy[match(rownames(pancanBinaryAMP),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanBinaryAMP['PanCan',] = colSums(pancanBinaryAMP,na.rm = T)
#compute percentages
pancanBinaryAMP_Perc = pancanBinaryAMP*100/pancanBinaryAMP$TotalN
pancanBinaryAMP_Perc = pancanBinaryAMP_Perc[,-c(ncol(pancanBinaryAMP_Perc))]
pancanBinaryAMP_Perc[is.na(pancanBinaryAMP_Perc)] = 0
pancanBinaryAMP_Perc = as.matrix(pancanBinaryAMP_Perc)
pancanBinaryAMP_Perc = pancanBinaryAMP_Perc[,order(colnames(pancanBinaryAMP_Perc))]

lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
heatmap.2(pancanBinaryAMP_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryAMP_Perc==0,NA,round(pancanBinaryAMP_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples Amplified\n(2 or more copy gains)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.8)


####heatmap of percent of samples with deletions per network gene per tumor-type
#for first heatmap, look at both 1 and 2 copy losses
pancanBinaryDEL = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$ThresholdedCopyNumber)) < 0,c('Study','Gene','Count')]
pancanBinaryDEL = aggregate(pancanBinaryDEL$Count,by = pancanBinaryDEL[,c('Study','Gene')],sum)
names(pancanBinaryDEL) = c('Study','Gene','Count')
#spread for heatmap
pancanBinaryDEL = spread(pancanBinaryDEL,Gene,Count)
rownames(pancanBinaryDEL) = pancanBinaryDEL$Study
pancanBinaryDEL = pancanBinaryDEL[,-c(1)]
#Add sample count per Study as column
pancanBinaryDEL[,'TotalN']= sampleCountPerStudy[match(rownames(pancanBinaryDEL),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanBinaryDEL['PanCan',] = colSums(pancanBinaryDEL,na.rm = T)
#compute percentages
pancanBinaryDEL_Perc = pancanBinaryDEL*100/pancanBinaryDEL$TotalN
pancanBinaryDEL_Perc = pancanBinaryDEL_Perc[,-c(ncol(pancanBinaryDEL_Perc))]
pancanBinaryDEL_Perc[is.na(pancanBinaryDEL_Perc)] = 0
pancanBinaryDEL_Perc = as.matrix(pancanBinaryDEL_Perc)
pancanBinaryDEL_Perc = pancanBinaryDEL_Perc[,order(colnames(pancanBinaryDEL_Perc))]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
scaleBlue<-colorRampPalette(colors=c("white","blue"))(1000)
heatmap.2(pancanBinaryDEL_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryDEL_Perc==0,NA,round(pancanBinaryDEL_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples with Deletions\n(1 or more copy losses)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.8)
#for second heatmap, look at only 2 copy gains
pancanBinaryDEL = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$ThresholdedCopyNumber)) < -1,c('Study','Gene','Count')]
pancanBinaryDEL = aggregate(pancanBinaryDEL$Count,by = pancanBinaryDEL[,c('Study','Gene')],sum)
names(pancanBinaryDEL) = c('Study','Gene','Count')
#spread for heatmap
pancanBinaryDEL = spread(pancanBinaryDEL,Gene,Count)
rownames(pancanBinaryDEL) = pancanBinaryDEL$Study
pancanBinaryDEL = pancanBinaryDEL[,-c(1)]
#Add sample count per Study as column
pancanBinaryDEL[,'TotalN']= sampleCountPerStudy[match(rownames(pancanBinaryDEL),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanBinaryDEL['PanCan',] = colSums(pancanBinaryDEL,na.rm = T)
#compute percentages
pancanBinaryDEL_Perc = pancanBinaryDEL*100/pancanBinaryDEL$TotalN
pancanBinaryDEL_Perc = pancanBinaryDEL_Perc[,-c(ncol(pancanBinaryDEL_Perc))]
pancanBinaryDEL_Perc[is.na(pancanBinaryDEL_Perc)] = 0
pancanBinaryDEL_Perc = as.matrix(pancanBinaryDEL_Perc)
pancanBinaryDEL_Perc = pancanBinaryDEL_Perc[,order(colnames(pancanBinaryDEL_Perc))]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
heatmap.2(pancanBinaryDEL_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryDEL_Perc==0,NA,round(pancanBinaryDEL_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples with Deletions\n(2 copy losses)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.8)


###########MUTATIONS
#Proximal netwrok MAF Analysis
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.ProximalNetwork_MAF`"
mafDF <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#count of samples per tumor type that have mutations in proximal network genes
mutCountPerStudy = table(mafDF$Study,mafDF$Hugo_Symbol)
mutCountPerStudy = as.data.frame.matrix(mutCountPerStudy)
mutCountPerStudy['UVM',] = rep(0,ncol(mutCountPerStudy))
#add total sample count per study as column
mutCountPerStudy[,'TotalN'] = sampleCountPerStudy$Count[match(rownames(mutCountPerStudy),sampleCountPerStudy$Study)]
#add PanCan row
mutCountPerStudy['PanCan',] = colSums(mutCountPerStudy,na.rm = T)

#compute percentages
mutPercPerStudy = data.matrix(mutCountPerStudy)*100/as.numeric(mutCountPerStudy$TotalN)
#drop TotalN column
mutPercPerStudy = mutPercPerStudy[,1:ncol(mutPercPerStudy)-1]
heatmap.2(mutPercPerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mutPercPerStudy==0,NA,round(mutPercPerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nCoding mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

#########GET GEXP values for various boxplots to be plotted below
##1. MYC network GEXP
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol IN ('MYC','MYCN','MYCL','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')"
mycNetworkGEXP <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
mycNetworkGEXP_Wide = mycNetworkGEXP[,c('SampleBarcode','Symbol','normalized_count')]
#take average of duplicate sample-gene records
mycNetworkGEXP_Wide = aggregate(mycNetworkGEXP_Wide,by = list(SampleBarcode = mycNetworkGEXP_Wide$SampleBarcode,Symbol=mycNetworkGEXP_Wide$Symbol),median,na.rm=T)
mycNetworkGEXP_Wide = mycNetworkGEXP_Wide[,c(1,2,5)]
#spread mycNetworkGEXP into 2D sample X gene matrix
mycNetworkGEXP_Wide = spread(mycNetworkGEXP_Wide,key = Symbol,value = normalized_count)
rownames(mycNetworkGEXP_Wide) = mycNetworkGEXP_Wide$SampleBarcode
mycNetworkGEXP_Wide = mycNetworkGEXP_Wide[,-c(1)]
#compute average MYC network expression per sample
mycNetworkGEXP_Wide[,'MedianNetworkGEXP'] = apply(mycNetworkGEXP_Wide,1,function(x) median(x,na.rm = T))
mycNetworkGEXP_Wide[,'Log2MedianNetworkGEXP'] = log2(mycNetworkGEXP_Wide[,'MedianNetworkGEXP']+1)

##2. OTHER PATHWAYS
#analyse for GO_HELICASE_ACTIVITY, RNA Polymerase activity, Translation elongation, Translation Initiation, Chemokine_Activity
pathwayGenes = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/GO_TranslationElongation_GeneList.csv',sep = ',',header = T)

#get expression data for these genes
inList = paste(shQuote(pathwayGenes$GO_TRANSLATION_ELONGATION_FACTOR_ACTIVITY),collapse = ",")
querySql <- paste("SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol IN (",inList,")",sep="")
pathwayGEXP <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

pathwayGEXP_Wide = pathwayGEXP[,c('SampleBarcode','Symbol','normalized_count')]
pathwayGEXP_Wide = aggregate(pathwayGEXP_Wide,by = list(SampleBarcode = pathwayGEXP_Wide$SampleBarcode,Symbol=pathwayGEXP_Wide$Symbol),median,na.rm=T)
pathwayGEXP_Wide = pathwayGEXP_Wide[,c(1,2,5)]
#spread pathwayGEXP into 2D sample X gene matrix
pathwayGEXP_Wide = spread(pathwayGEXP_Wide,key = Symbol,value = normalized_count)
rownames(pathwayGEXP_Wide) = pathwayGEXP_Wide$SampleBarcode
pathwayGEXP_Wide = pathwayGEXP_Wide[,-c(1)]
#compute average MYC network expression per sample
pathwayGEXP_Wide[,'MedianNetworkGEXP'] = apply(pathwayGEXP_Wide,1,function(x) median(x,na.rm = T))
pathwayGEXP_Wide[,'Log2MedianNetworkGEXP'] = log2(pathwayGEXP_Wide[,'MedianNetworkGEXP']+1)


########START LOW THRESHOLD ANALYSES- Combine CNVR and MUTATIONS, #include all copy gains and losses and binarize
pancanBinaryAlt = pancanThresholdedCNVR[,-c(ncol(pancanThresholdedCNVR))] #remove Study column
pancanBinaryAlt[pancanBinaryAlt!=0] = 1 #include all copy gains and losses and binarize
#OR mutation data into pancanBinaryAlt
for(index in 1:nrow(mafDF))
{
  pancanBinaryAlt[rownames(pancanBinaryAlt)==mafDF$SAMPLE_BARCODE[index],colnames(pancanBinaryAlt)==mafDF$Hugo_Symbol[index]] = 1
}
#heatmap of alterations per network gene per tumor-type
#add Study column
pancanBinaryAlt[,'Study'] = pancanThresholdedCNVR$Study
#melt and aggregate sample counts per gene per study
pancanBinaryAlt_Long = melt(pancanBinaryAlt,id='Study')
names(pancanBinaryAlt_Long) = c('Study','Gene','Altered')
countAlteredPerGenePerStudy = aggregate(pancanBinaryAlt_Long$Altered,by=pancanBinaryAlt_Long[,c('Study','Gene')],sum)
countAlteredPerGenePerStudy_Wide = spread(countAlteredPerGenePerStudy,Gene,x)
rownames(countAlteredPerGenePerStudy_Wide) = countAlteredPerGenePerStudy_Wide$Study
countAlteredPerGenePerStudy_Wide = countAlteredPerGenePerStudy_Wide[,-c(1)]
#add TotalN as column
countAlteredPerGenePerStudy_Wide[,'TotalN'] = sampleCountPerStudy$Count[match(rownames(countAlteredPerGenePerStudy_Wide),sampleCountPerStudy$Study)]
#add PanCan row
countAlteredPerGenePerStudy_Wide['PanCan',] = colSums(countAlteredPerGenePerStudy_Wide)
#compute percentage
percAlteredPerGenePerStudy = countAlteredPerGenePerStudy_Wide*100/countAlteredPerGenePerStudy_Wide$TotalN
percAlteredPerGenePerStudy = percAlteredPerGenePerStudy[,-c(ncol(percAlteredPerGenePerStudy))]
percAlteredPerGenePerStudy = as.matrix(percAlteredPerGenePerStudy)
heatmap.2(percAlteredPerGenePerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(percAlteredPerGenePerStudy==0,NA,round(percAlteredPerGenePerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nAMP/DEL/Coding Mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

#############LOW THRESHOLD - Bar graph showing percentage alterations in MYC network per tumor-type (only +/-2 CNVR counted)
pancanBinaryAlt[,'NWAltered'] = ifelse(rowSums(pancanBinaryAlt[,c(1:ncol(pancanBinaryAlt)-1)],na.rm = T) != 0,1,0)
countAltered_AllNW_PerStudy = aggregate(pancanBinaryAlt$NWAltered,by=list(pancanBinaryAlt$Study),sum)
names(countAltered_AllNW_PerStudy) = c('Study','Count')
#add TotalN columnm
countAltered_AllNW_PerStudy[,'TotalN'] = sampleCountPerStudy$Count[match(countAltered_AllNW_PerStudy$Study,sampleCountPerStudy$Study)]
#compute perc
countAltered_AllNW_PerStudy[,'Perc'] = countAltered_AllNW_PerStudy$Count*100/countAltered_AllNW_PerStudy$TotalN
#order Study levels by Perc
countAltered_AllNW_PerStudy$Study = factor(countAltered_AllNW_PerStudy$Study,levels = as.vector(countAltered_AllNW_PerStudy$Study[order(countAltered_AllNW_PerStudy$Perc)]))
#barplot
ggplot(data=countAltered_AllNW_PerStudy, aes(x=Study, y=Perc)) +
  geom_bar(stat="identity") +
  labs(y="Percentage of Sample Altered",title='MYC Network Alterations') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7),plot.title = element_text(hjust = 0.5)) 

###Boxplots SummaryAlt vs. GEXP
#summarize binary Alt columns into WildType, MYCAltered, AnyButMYCAltered
summaryAlt = function(x) {
  #print(x['MYC'])
  #return()
  if(all(x==F,na.rm = T))
  {
    return('NoAlter')
  }else if(x['MYC']==F)
  {
    return('AnyButMYC')
  }else
  {
    return('MYCAlt')
  }
  
}
summaryAltPerSample = apply(pancanBinaryAlt[,c(1:13)],MARGIN = 1,FUN = summaryAlt)
pancanBinaryAlt[,'SummaryAlt'] = summaryAltPerSample 
pancanBinaryAlt[,'Log2MYCExp'] = log2(mycNetworkGEXP_Wide$MYC[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)
pancanBinaryAlt[,'Log2MYCMedianNWExp'] = mycNetworkGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]
pancanBinaryAlt[,'Log2MedianPathwayExp'] = pathwayGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(pathwayGEXP_Wide))]
pancanBinaryAlt = pancanBinaryAlt[!is.na(pancanBinaryAlt$Study),]
give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
#1.low threshold alterations, MYC gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title=""))
#2. low threshold alterations, MYC gexp, per Study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title=""))
#3. low threshold alterations, MYC median network gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCMedianNWExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median MYC Network Expression",title='MYC Network Expression') +
  guides(fill=guide_legend(title=""))
#4. low threshold alterations, MYC median network gexp, per Study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCMedianNWExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median MYC Network Expression",title='MYC Network Expression') +
  guides(fill=guide_legend(title=""))
#5. low threshold alterations, median pathway gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MedianPathwayExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO CHEMOKINE ACTIVITY') +
  guides(fill=guide_legend(title=""))
#6. low threshold alterations, median pathway gexp, per study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MedianPathwayExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO CHEMOKINE ACTIVITY') +
  guides(fill=guide_legend(title=""))



########START HIGH THRESHOLD analyses- Combine CNVR and MUTATIONS, #include only +/-2 gains and losses and binarize
pancanBinaryAlt = pancanThresholdedCNVR[,-c(ncol(pancanThresholdedCNVR))] #remove Study column
pancanBinaryAlt[pancanBinaryAlt==1 | pancanBinaryAlt==-1] = 0 
pancanBinaryAlt[pancanBinaryAlt==2 | pancanBinaryAlt==-2] = 1
#OR mutation data into pancanBinaryAlt
for(index in 1:nrow(mafDF))
{
  pancanBinaryAlt[rownames(pancanBinaryAlt)==mafDF$SAMPLE_BARCODE[index],colnames(pancanBinaryAlt)==mafDF$Hugo_Symbol[index]] = 1
}
#heatmap of alterations per network gene per tumor-type
#add Study column
pancanBinaryAlt[,'Study'] = pancanThresholdedCNVR$Study
#melt and aggregate sample counts per gene per study
pancanBinaryAlt_Long = melt(pancanBinaryAlt,id='Study')
names(pancanBinaryAlt_Long) = c('Study','Gene','Altered')
countAlteredPerGenePerStudy = aggregate(pancanBinaryAlt_Long$Altered,by=pancanBinaryAlt_Long[,c('Study','Gene')],sum)
countAlteredPerGenePerStudy_Wide = spread(countAlteredPerGenePerStudy,Gene,x)
rownames(countAlteredPerGenePerStudy_Wide) = countAlteredPerGenePerStudy_Wide$Study
countAlteredPerGenePerStudy_Wide = countAlteredPerGenePerStudy_Wide[,-c(1)]
#add TotalN as column
countAlteredPerGenePerStudy_Wide[,'TotalN'] = sampleCountPerStudy$Count[match(rownames(countAlteredPerGenePerStudy_Wide),sampleCountPerStudy$Study)]
#add PanCan row
countAlteredPerGenePerStudy_Wide['PanCan',] = colSums(countAlteredPerGenePerStudy_Wide)
#compute percentage
percAlteredPerGenePerStudy = countAlteredPerGenePerStudy_Wide*100/countAlteredPerGenePerStudy_Wide$TotalN
percAlteredPerGenePerStudy = percAlteredPerGenePerStudy[,-c(ncol(percAlteredPerGenePerStudy))]
percAlteredPerGenePerStudy = as.matrix(percAlteredPerGenePerStudy)
heatmap.2(percAlteredPerGenePerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(percAlteredPerGenePerStudy==0,NA,round(percAlteredPerGenePerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nAMP/DEL/Coding Mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

#############HIGH THRESHOLD - Bar graph showing percentage alterations in MYC network per tumor-type (only +/-2 CNVR counted)
pancanBinaryAlt[,'NWAltered'] = ifelse(rowSums(pancanBinaryAlt[,c(1:ncol(pancanBinaryAlt)-1)],na.rm = T) != 0,1,0)
countAltered_AllNW_PerStudy = aggregate(pancanBinaryAlt$NWAltered,by=list(pancanBinaryAlt$Study),sum)
names(countAltered_AllNW_PerStudy) = c('Study','Count')
#add TotalN columnm
countAltered_AllNW_PerStudy[,'TotalN'] = sampleCountPerStudy$Count[match(countAltered_AllNW_PerStudy$Study,sampleCountPerStudy$Study)]
#compute perc
countAltered_AllNW_PerStudy[,'Perc'] = countAltered_AllNW_PerStudy$Count*100/countAltered_AllNW_PerStudy$TotalN
#order Study levels by Perc
countAltered_AllNW_PerStudy$Study = factor(countAltered_AllNW_PerStudy$Study,levels = as.vector(countAltered_AllNW_PerStudy$Study[order(countAltered_AllNW_PerStudy$Perc)]))
#barplot
ggplot(data=countAltered_AllNW_PerStudy, aes(x=Study, y=Perc)) +
  geom_bar(stat="identity") +
  labs(y="Percentage of Sample Altered",title='MYC Network Alterations') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7),plot.title = element_text(hjust = 0.5)) 

###############HIGH THRESHOLD - MYC NW alterations vs. GEXP
#summarize binary Alt columns into WildType, MYCAltered, AnyButMYCAltered
summaryAlt = function(x) {
  #print(x['MYC'])
  #return()
  if(all(x==F,na.rm = T))
  {
    return('NoAlter')
  }else if(x['MYC']==F)
  {
    return('AnyButMYC')
  }else
  {
    return('MYCAlt')
  }
  
}
summaryAltPerSample = apply(pancanBinaryAlt[,c(1:13)],MARGIN = 1,FUN = summaryAlt)
pancanBinaryAlt[,'SummaryAlt'] = summaryAltPerSample 
pancanBinaryAlt[,'Log2MYCExp'] = log2(mycNetworkGEXP_Wide$MYC[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)
pancanBinaryAlt[,'Log2MYCMedianNWExp'] = mycNetworkGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]
pancanBinaryAlt[,'Log2MedianPathwayExp'] = pathwayGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(pathwayGEXP_Wide))]
pancanBinaryAlt = pancanBinaryAlt[!is.na(pancanBinaryAlt$Study),]
#compute anova and effect sizes
#compute one-way ANOVA
panCanModel = conover.test::conover.test(pancanBinaryAlt$Log2MedianPathwayExp,pancanBinaryAlt$SummaryAlt,method = "bh")
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_NoAlt)
print(Any_MYC)
print(MYC_NoAlt)

give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
#1.high threshold alterations, MYC gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title=""))
#2. high threshold alterations, MYC gexp, per Study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title=""))
#3. high threshold alterations, MYC median network gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCMedianNWExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median MYC Network Expression",title='MYC Network Expression') +
  guides(fill=guide_legend(title=""))
#4. high threshold alterations, MYC median network gexp, per Study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCMedianNWExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median MYC Network Expression",title='MYC Network Expression') +
  guides(fill=guide_legend(title=""))
#5. high threshold alterations, median pathway gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MedianPathwayExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO CHEMOKINE ACTIVITY') +
  guides(fill=guide_legend(title=""))
#6. high threshold alterations, median pathway gexp, per study
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MedianPathwayExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO CHEMOKINE ACTIVITY') +
  guides(fill=guide_legend(title=""))

