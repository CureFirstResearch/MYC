#re-analyze with focal copy number variation - provided by Andy Cherniack

library(ggplot2)
library(bigrquery)
library(tidyr) #for spread
library(reshape2)
library(gplots) #for heatmap.2
library(effsize) #for cohen.d
library(ggrepel)
library(ggsignif)
library(pheatmap)
billingProject = 'isb-cgc-04-0007'

#######read thresholded AMP/DEL status
pancanFocalCNVR = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/ISAR-GISTIC_focal_myc_pathway_values.txt',sep = "\t",header = T)
#pancanFocalCNVR = pancanFocalCNVR[,c('sample','tumor_type','ploidy','MAX_focal_alt','MGA_focal_alt','MLXIPL_focal_alt','MLXIP_focal_alt','MLX_focal_alt','MNT_focal_alt','MXD1_focal_alt','MXD3_focal_alt','MXD4_focal_alt','MXI1_focal_alt','MYCL1_focal_alt','MYCN_focal_alt','MYC_focal_alt')]

#remove study and ploidy
#pancanFocalCNVR = pancanFocalCNVR[,-c(2,3)]
#assign row names
rownames(pancanFocalCNVR) = pancanFocalCNVR[,'sample']
pancanFocalCNVR = pancanFocalCNVR[,-c(1)]

#trim rownames to be sample barcodes instead of full aliquot barcodes
rownames(pancanFocalCNVR) = substring(rownames(pancanFocalCNVR),1,15)

########read pancan pathway whitelist 
sampleWhitelist = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/PanCanPathwayWhitelist.tsv',sep="\t",header = T)


########filter to keep only whitelisted samples
pancanFocalCNVR = pancanFocalCNVR[as.vector(sampleWhitelist$SAMPLE_BARCODE),]
#remove samples with all NAs
pancanFocalCNVR = pancanFocalCNVR[rowSums(is.na(pancanFocalCNVR))!=ncol(pancanFocalCNVR),]

#reduce sampleWhitelist to only those samples that have CN data
sampleWhitelist = sampleWhitelist[sampleWhitelist$SAMPLE_BARCODE %in% rownames(pancanFocalCNVR),]
sampleCountPerStudy = as.data.frame(table(sampleWhitelist$DISEASE))
names(sampleCountPerStudy) = c('Study','Count')

########clustering of pan-cancer samples based on PMN CNVR#########################
##filter out any samples with all 0's in the matrix
focalCNVRForHeatmap = pancanFocalCNVR[rowSums(pancanFocalCNVR==0)!=ncol(pancanFocalCNVR),]
annRow = as.data.frame(as.vector(sampleWhitelist$DISEASE[match(rownames(focalCNVRForHeatmap),sampleWhitelist$SAMPLE_BARCODE)]))
rownames(annRow) = rownames(focalCNVRForHeatmap)
colnames(annRow) = c("Study")

pheatmap(as.matrix(focalCNVRForHeatmap),cluster_rows = T,cluster_cols = T,color = colorRampPalette(c("blue","white","red"))(100),breaks = c(seq(min(focalCNVRForHeatmap),0,length.out = 51),seq(max(focalCNVRForHeatmap)/100,max(focalCNVRForHeatmap),length.out = 50)),scale = "none",annotation_row = annRow,show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = "Focal Copy Numbers in Proximal MYC Network\n(Continuous values, Heirarchical clustering)")

pheatmap(as.matrix(focalCNVRForHeatmap[order(focalCNVRForHeatmap$MYC,decreasing = T),]),cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("blue","white","red"))(100),breaks = c(seq(min(focalCNVRForHeatmap),0,length.out = 51),seq(max(focalCNVRForHeatmap)/100,max(focalCNVRForHeatmap),length.out = 50)),scale = "none",annotation_row = annRow,show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = "Focal Copy Numbers in Proximal MYC Network\n(Continuous values, MYC ordered)")

#categorize : >0 = AMP, <0 = DEL
pancanFocalCNVR[pancanFocalCNVR>0] = 1
pancanFocalCNVR[pancanFocalCNVR<0] = -1
focalCNVRForHeatmap = pancanFocalCNVR[rowSums(pancanFocalCNVR==0)!=ncol(pancanFocalCNVR),]
annRow = as.data.frame(as.vector(sampleWhitelist$DISEASE[match(rownames(focalCNVRForHeatmap),sampleWhitelist$SAMPLE_BARCODE)]))
rownames(annRow) = rownames(focalCNVRForHeatmap)
colnames(annRow) = c("Study")
myColors = c("blue","white","red")
pheatmap(as.matrix(focalCNVRForHeatmap),cluster_rows = T,cluster_cols = T,color = myColors,scale = "none",annotation_row = annRow,show_rownames = F,legend = T,legend_breaks = c(-1,0,1),legend_labels = c("DEL","","AMP"),fontsize = 8,fontsize_col = 12,main = "Focal Copy Numbers in Proximal MYC Network\n(Binarized values, Heirarchical Clustering)")

pheatmap(as.matrix(focalCNVRForHeatmap[order(focalCNVRForHeatmap$MYC,decreasing = T),]),cluster_rows = F,cluster_cols = F,color = myColors,scale = "none",annotation_row = annRow,show_rownames = F,legend = T,legend_breaks = c(-1,0,1),legend_labels = c("DEL","","AMP"),fontsize = 8,fontsize_col = 12,main = "Focal Copy Numbers in Proximal MYC Network\n(Binarized values, MYC ordered)")

###########analysis of samples with just one PMN gene altered############################################
###get gexp data ready for subsetting
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

pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/PanCancerClustering/ExclusiveGeneAlterations.pdf',onefile = T)
for(thisGene in colnames(pancanFocalCNVR))
{
  thisGeneAMP = pancanFocalCNVR[pancanFocalCNVR[,thisGene]==1,]
  thisGeneOnlyAMP = thisGeneAMP[rowSums(thisGeneAMP!=0)==1,]
  thisGeneAMPSamples = as.data.frame(rownames(thisGeneOnlyAMP))
  colnames(thisGeneAMPSamples) = c('SampleBarcode')
  thisGeneAMPSamples[,'Study'] = sampleWhitelist$DISEASE[match(thisGeneAMPSamples$SampleBarcode,sampleWhitelist$SAMPLE_BARCODE)]
  thisGeneAMPSamples[,'Subtype'] = sampleWhitelist$SUBTYPE[match(thisGeneAMPSamples$SampleBarcode,sampleWhitelist$SAMPLE_BARCODE)]
  
  p=ggplot(thisGeneAMPSamples,aes(x=Study)) +
    geom_bar(aes(y = (..count..)/sum(..count..))) +
    scale_y_continuous(labels=scales::percent) +
    labs(y="Percent of samples",title=paste("Samples with only",thisGene,'amplified\nDistribution across tumor types\nTotal N =',nrow(thisGeneAMPSamples))) +
    theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))
  print(p)
  #PMN gexp profiles for these samples
  thisGexp = pmnGexpDF_Wide[rownames(pmnGexpDF_Wide) %in% thisGeneAMPSamples$SampleBarcode,]
  annRow = as.data.frame(thisGexp$Study)
  rownames(annRow) = rownames(thisGexp)
  colnames(annRow) = c("Study")
  p=pheatmap(as.matrix(thisGexp[,1:13]),cluster_rows = T,cluster_cols = T,scale = "none",annotation_row = annRow,show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = paste("Gene Expression in Proximal MYC Network\n(With only",thisGene,"amplified)"))
  print(p)
  #deletion
  thisGeneDEL = pancanFocalCNVR[pancanFocalCNVR[,thisGene]==-1,]
  thisGeneOnlyDEL = thisGeneDEL[rowSums(thisGeneDEL!=0)==1,]
  thisGeneDELSamples = as.data.frame(rownames(thisGeneOnlyDEL))
  colnames(thisGeneDELSamples) = c('SampleBarcode')
  thisGeneDELSamples[,'Study'] = sampleWhitelist$DISEASE[match(thisGeneDELSamples$SampleBarcode,sampleWhitelist$SAMPLE_BARCODE)]
  thisGeneDELSamples[,'Subtype'] = sampleWhitelist$SUBTYPE[match(thisGeneDELSamples$SampleBarcode,sampleWhitelist$SAMPLE_BARCODE)]
  
  p=ggplot(thisGeneDELSamples,aes(x=Study)) +
    geom_bar(aes(y = (..count..)/sum(..count..))) +
    scale_y_continuous(labels=scales::percent) +
    labs(y="Percent of samples",title=paste("Samples with only",thisGene,'deleted\nDistribution across tumor types\nTotal N =',nrow(thisGeneDELSamples))) +
    theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))
  print(p)
  #PMN gexp profiles for these samples
  thisGexp = pmnGexpDF_Wide[rownames(pmnGexpDF_Wide) %in% thisGeneDELSamples$SampleBarcode,]
  annRow = as.data.frame(thisGexp$Study)
  rownames(annRow) = rownames(thisGexp)
  colnames(annRow) = c("Study")
  p=pheatmap(as.matrix(thisGexp[,1:13]),cluster_rows = T,cluster_cols = T,scale = "none",annotation_row = annRow,show_rownames = F,legend = T,fontsize = 8,fontsize_col = 12,main = paste("Gene Expression in Proximal MYC Network\n(With only",thisGene,"deleted)"))
  print(p)
}
dev.off()
########plot distribution per gene ###########################
#first, compute counts per category within a gene
countsPerGene = apply(pancanFocalCNVR, 2, FUN = function(x) table(x))
#convert into percentage of samples within each gene
percsPerGene = as.data.frame(countsPerGene*100/colSums(countsPerGene))
percsPerGene[,'Category'] = rownames(percsPerGene)
percsPerGene$Category = factor(percsPerGene$Category,levels = c("1","0","-1"))
#next, melt DF into long form
percsPerGene_Long = melt(percsPerGene,id = c('Category'))
names(percsPerGene_Long) = c('FocalCopyNumber','Gene','PercentSamples')
percsPerGene_Long$FocalCopyNumber = as.vector(percsPerGene_Long$FocalCopyNumber)
#######stacked bar plot
percsPerGene_Long$FocalCopyNumber[percsPerGene_Long$FocalCopyNumber=="1"] = "AMP"
percsPerGene_Long$FocalCopyNumber[percsPerGene_Long$FocalCopyNumber=="-1"] = "DEL"
percsPerGene_Long$FocalCopyNumber[percsPerGene_Long$FocalCopyNumber=="0"] = "DUP"
percsPerGene_Long$FocalCopyNumber = factor(percsPerGene_Long$FocalCopyNumber)
group.colors <- c("AMP" = "red","DUP" ="gray", "DEL" = "blue")
pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/Figures/FigureS1/FocalCNVRDistribution_PerGene_PanCan.pdf',width = 9)
ggplot(data=percsPerGene_Long[percsPerGene_Long$FocalCopyNumber!="DUP",], aes(x=Gene, y=PercentSamples, fill=FocalCopyNumber)) +
  geom_bar(stat="identity",position = position_dodge()) +
  theme(legend.title = element_text(size = 12),legend.text = element_text(size = 10),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 12),axis.title = element_text(size = 14)) + 
  scale_fill_manual(values = group.colors)
dev.off()

#######stacked bar plots, facet_wrapped by tumor-type
#first, add Study as a column
pancanFocalCNVR[,'Study'] = sampleWhitelist[match(rownames(pancanFocalCNVR),sampleWhitelist$SAMPLE_BARCODE),'DISEASE']
#melt into long format
pancanFocalCNVR_Long = melt(pancanFocalCNVR,id=c('Study'))
names(pancanFocalCNVR_Long) = c('Study','Gene','FocalCopyNumber')
#compute counts per Study, per gene, per copy number
countsPerGenePerStudy = aggregate(rep(1,nrow(pancanFocalCNVR_Long)),by=pancanFocalCNVR_Long[,c('Study','Gene','FocalCopyNumber')],sum)
names(countsPerGenePerStudy) = c('Study','Gene','FocalCopyNumber','Count')
countsPerGenePerStudy[,'TotalNThisStudy'] = sampleCountPerStudy[match(countsPerGenePerStudy$Study,sampleCountPerStudy$Study),'Count']
countsPerGenePerStudy[,'PercentSamples'] = countsPerGenePerStudy$Count*100/countsPerGenePerStudy$TotalNThisStudy
#convert copy number to string for clearer colors in ggplot
countsPerGenePerStudy$FocalCopyNumber = factor(countsPerGenePerStudy$FocalCopyNumber,levels = c("1","0","-1"))
group.colors <- c("1" = "red","0" ="gray", "-1" = "blue")
ggplot(data=countsPerGenePerStudy, aes(x=Gene, y=PercentSamples, fill=FocalCopyNumber)) +
  geom_bar(stat="identity") +
  facet_wrap(~Study) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5.5)) +
  scale_fill_manual(values = group.colors)
  
####heatmap of percent of sample amplified per network gene per tumor-type
pancanBinaryAMP = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$FocalCopyNumber)) > 0,c('Study','Gene','Count')]
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
#heatmap.2
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.8,4.5)
lwid = c(0.5,4,1.5)
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
heatmap.2(pancanBinaryAMP_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryAMP_Perc==0,NA,round(pancanBinaryAMP_Perc,digits = 1)),notecex = 1.5,notecol="black",sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples Amplified\n(Focal Copy Gain)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.5)
#with tumor-types clustered
heatmap.2(pancanBinaryAMP_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryAMP_Perc==0,NA,round(pancanBinaryAMP_Perc,digits = 1)),notecex = 1.5,notecol="black",sepcolor = "cyan",Rowv=TRUE,Colv = FALSE,dendrogram="row",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples Amplified\n(Focal Copy Gain)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.5)

#ToDo : ggplot heatmap

####heatmap of percent of samples with deletions per network gene per tumor-type
pancanBinaryDEL = countsPerGenePerStudy[as.numeric(as.vector(countsPerGenePerStudy$FocalCopyNumber)) < 0,c('Study','Gene','Count')]
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
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.8,4.5)
lwid = c(0.5,4,1.5)
scaleBlue<-colorRampPalette(colors=c("white","blue"))(1000)
heatmap.2(pancanBinaryDEL_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryDEL_Perc==0,NA,round(pancanBinaryDEL_Perc,digits = 1)),notecex = 1.5,notecol="black",sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples with Deletions\n(Focal Copy Loss)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.5)
#with tumor-types clustered
heatmap.2(pancanBinaryDEL_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanBinaryDEL_Perc==0,NA,round(pancanBinaryDEL_Perc,digits = 1)),notecex = 1.5,notecol="black",sepcolor = "cyan",Rowv=TRUE,Colv = FALSE,dendrogram="row",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples with Deletions\n(Focal Copy Loss)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.5)

###########MUTATIONS
#Proximal netwrok MAF Analysis
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.ProximalNetwork_MAF`"
mafDF <- query_exec(querySql, project=billingProject,use_legacy_sql = F,max_pages = Inf)

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
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.8,4.5)
lwid = c(0.5,4,1.5)
scaleGreen<-colorRampPalette(colors=c("white","darkgreen"))(1000)
heatmap.2(mutPercPerStudy,breaks = seq(0,100,length.out = 1001),Rowv=F,Colv=F,dendrogram='none',notecex = 1.5,cellnote = ifelse(mutPercPerStudy==0,NA,round(mutPercPerStudy,digits = 1)),notecol="black",cexRow=1.5,cexCol = 1.5,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nCoding mutations",na.color = "gray",symbreaks = F,col=scaleGreen,density.info="none",keysize = 0.5)
#with tumor-types clustered
heatmap.2(mutPercPerStudy,breaks = seq(0,100,length.out = 1001),Rowv=T,Colv=F,dendrogram='row',notecex = 1.5,cellnote = ifelse(mutPercPerStudy==0,NA,round(mutPercPerStudy,digits = 1)),notecol="black",cexRow=1.5,cexCol = 1.5,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nCoding mutations",na.color = "gray",symbreaks = F,col=scaleGreen,density.info="none",keysize = 0.5)

#########GET GEXP values for various boxplots to be plotted below
##1. MYC network GEXP
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol IN ('MYC','MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')"
mycNetworkGEXP <- query_exec(querySql, project = billingProject,use_legacy_sql = F,max_pages = Inf)
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
#MYC GEXP violin plots per study for manuscript
mycGEXP = mycNetworkGEXP[mycNetworkGEXP$Symbol %in% c('MYC','MYCN'),]
pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/Figures/Figure5D/MYCGexp.pdf',width = 11,height = 5)
ggplot(mycGEXP,aes(x=Study,y=log2_count)) +
  geom_boxplot() +
  facet_wrap(~Symbol,nrow = 2) +
  theme(axis.text.x = element_text(angle = 90,size=10),plot.title = element_text(size=10),axis.title = element_text(size = 14)) +
  labs(y="log2 normalized count")
dev.off()
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


########Combine CNVR and MUTATIONS, #include all copy gains and losses and binarize
pancanBinaryAlt = pancanFocalCNVR[,-c(ncol(pancanFocalCNVR))] #remove Study column
##save MYCL to look at just MYCL amplifications
pancanMYC = pancanBinaryAlt[,'MYCL1']
pancanBinaryAlt[pancanBinaryAlt!=0] = 1 #include all copy gains and losses and binarize
#OR mutation data into pancanBinaryAlt
for(index in 1:nrow(mafDF))
{
  pancanBinaryAlt[rownames(pancanBinaryAlt)==mafDF$SAMPLE_BARCODE[index],colnames(pancanBinaryAlt)==mafDF$Hugo_Symbol[index]] = 1
}
###restore MYC column to have just AMPlifications
pancanBinaryAlt$MYCL1 = pancanMYC
pancanBinaryAlt$MYCL1[pancanBinaryAlt$MYCL1 <0] = 0
#heatmap of alterations per network gene per tumor-type
#add Study column
pancanBinaryAlt[,'Study'] = pancanFocalCNVR$Study
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
heatmap.2(percAlteredPerGenePerStudy,breaks = seq(0,100,length.out = 1001),Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(percAlteredPerGenePerStudy==0,NA,round(percAlteredPerGenePerStudy,digits = 1)),notecex=1.5,cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nAMP/DEL/Coding Mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.5)

#############Bar graph showing percentage alterations in MYC network per tumor-type 
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
  scale_fill_grey() +
  labs(y="Percentage of Sample Altered",title='MYC Network Alterations') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size = 14),plot.title = element_text(hjust = 0.5,size=18)) 

###Boxplots SummaryAlt vs. GEXP
#summarize binary Alt columns into WildType, MYCAltered, AnyButMYCAltered
summaryAlt = function(x) {
  #print(x['MYC'])
  #return()
  if(all(x==F,na.rm = T))
  {
    return('NoAlter')
  }else if(x['MYCL1']==F)
  {
    return('AnyButMYCL')
  }else
  {
    return('MYCL_AMP')
  }
  
}
summaryAltPerSample = apply(pancanBinaryAlt[,c(1:13)],MARGIN = 1,FUN = summaryAlt)
pancanBinaryAlt[,'SummaryAlt'] = summaryAltPerSample 
pancanBinaryAlt[,'Log2MYCExp'] = log2(mycNetworkGEXP_Wide$MYC[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)
pancanBinaryAlt[,'Log2MYCLExp'] = log2(mycNetworkGEXP_Wide$MYCL1[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)
pancanBinaryAlt[,'Log2MYCNExp'] = log2(mycNetworkGEXP_Wide$MYCN[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)

pancanBinaryAlt[,'Log2MYCMedianNWExp'] = mycNetworkGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]
pancanBinaryAlt[,'Log2MedianPathwayExp'] = pathwayGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(pathwayGEXP_Wide))]
pancanBinaryAlt = pancanBinaryAlt[!is.na(pancanBinaryAlt$Study),]
#compute one-way ANOVA
panCanModel = conover.test::conover.test(pancanBinaryAlt$Log2MYCLExp,pancanBinaryAlt$SummaryAlt,method = "bh")
Any_MYC = cohen.d(pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="AnyButMYCL"],pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="MYCL_AMP"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="AnyButMYCL"],pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pancanBinaryAlt$Log2MYCLExp[pancanBinaryAlt$SummaryAlt=="MYCL_AMP"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_NoAlt)
print(Any_MYC)
print(MYC_NoAlt)
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
#1.alterations, MYC gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCLExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  geom_signif(comparisons=list(c("AnyButMYCL", "MYCL_AMP")), annotations="**", y_position = 19, vjust=0.4,textsize = 15) +
  geom_signif(comparisons=list(c("MYCL_AMP", "NoAlter")), annotations="***", y_position = 19, vjust=0.4,textsize = 15) +
  geom_signif(comparisons=list(c("AnyButMYCL", "NoAlter")), annotations="*", y_position = 20, vjust=0.4,textsize = 15) +
  stat_summary(fun.data = give.n, geom = "text",size=7) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position="bottom",
        legend.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5,size = 20)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYCL Expression') +
  guides(fill=guide_legend(title="",keyheight = 2))
#2. alterations, MYC gexp, per Study
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
#3.alterations, MYC median network gexp, pancan
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
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO TRANSLATION ELONGATION FACTOR ACTIVITY') +
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
  labs(x="MYC Pathway Alterations",y="Log2 Median Pathway Expression",title='GO TRANSLATION ELONGATION FACTOR ACTIVITY') +
  guides(fill=guide_legend(title=""))

#compute one-way ANOVA
panCanModel = conover.test::conover.test(pancanBinaryAlt$Log2MedianPathwayExp,pancanBinaryAlt$SummaryAlt,method = "bh")
panCanModel$P.adjusted
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pancanBinaryAlt$Log2MedianPathwayExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_NoAlt)
print(Any_MYC)
print(MYC_NoAlt)


###################################
#MUTUAL EXCLUSIVITY - Method 1 - genome wide background distribution
wgFocal = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/GISTIC.focal_data_by_genes.conf_95.txt',sep="\t",header = T)

wgFocal = subset(wgFocal,select = -c(Locus.ID,Cytoband))
rownames(wgFocal) = wgFocal$Gene.Symbol
wgFocal = wgFocal[,-c(1)]
colnames(wgFocal) = gsub('.','-',substr(colnames(wgFocal),1,15),fixed = T)

#retain only whitelisted samples
keepList = colnames(wgFocal)[colnames(wgFocal) %in% as.vector(sampleWhitelist$SAMPLE_BARCODE)]
wgFocal = wgFocal[,keepList]

#binarize focal alterations
wgFocal[wgFocal!=0] = 1

#transpose so that columns are genes
wgFocal = t(wgFocal)


#OR with MAF data
querySql = "SELECT SAMPLE_BARCODE,Hugo_Symbol,Mutated FROM `isb-cgc-04-0007.MYC.PanCan_MAF_JOIN_Whitelist`"
mafDF = query_exec(querySql,project = billingProject,max_pages = Inf,use_legacy_sql = F)
mafDF_Wide = spread(mafDF,key=Hugo_Symbol,value=Mutated)
rownames(mafDF_Wide) = mafDF_Wide$SAMPLE_BARCODE
mafDF_Wide = mafDF_Wide[rownames(wgFocal),]

mafDF_Wide = mafDF_Wide[,-c(1)]
mafDF_Wide[!is.na(mafDF_Wide)] = 1
mafDF_Wide[is.na(mafDF_Wide)] = 0

#Alter = Focal OR MAF
wgAlter = wgFocal
for(i in 1:ncol(mafDF_Wide))
{
  print(i)
  thisGene = colnames(mafDF_Wide)[i]
  if(thisGene %in% colnames(wgAlter))
  {
    wgAlter[,thisGene] = wgAlter[,thisGene] | mafDF_Wide[,thisGene]
  }
}

geneOverlap = intersect(colnames(wgFocal),colnames(mafDF_Wide))
wgAlter[,geneOverlap] = wgFocal[,geneOverlap] | mafDF_Wide[,geneOverlap]
#
write.table(wgAlter,gzfile('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCanWhitelist_AllGenes_BinaryAltState_FocalCNVR_OR_MAF.tsv.gz'),sep = "\t",row.names = T)

#for every gene-MYC pair, count N for GeneNMycN, GeneYMycN, GeneNMycY, GeneYMycY
wgAlter = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCanWhitelist_AllGenes_BinaryAltState_FocalCNVR_OR_MAF.tsv',sep = "\t",header = T)
countsTables = lapply(apply(wgAlter,2,function(x) list(table(x,wgAlter[,'MYC']))),'[[',1)
countsTables = as.data.frame(countsTables)
countsTables = countsTables[,grep('Freq',colnames(countsTables),fixed = T)]
rownames(countsTables) = c('MycN_GeneN','MycN_GeneY','MycY_GeneN','MycY_GeneY')
colnames(countsTables) = gsub('.Freq','',colnames(countsTables),fixed = T)
countsTables['mutualXPerc',] = countsTables['MycY_GeneN',]*100/(countsTables['MycY_GeneN',]+countsTables['MycY_GeneY',])
countsTables = t(countsTables)
countsTables = data.frame(countsTables)
countsTables = countsTables[rownames(countsTables) != 'MYC',]
#1. Pan-Can, All genes vs. MYC mutualX density plot

ggplot(data.frame(countsTables),aes(x=mutualXPerc)) + 
  geom_density() +
  geom_point(data=countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text_repel(data=countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],],aes(label=rownames(countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],]),angle=45,hjust=0,vjust=0),size=2.5,y=0.012,show.legend = F) +
  geom_point(data=countsTables[order(countsTables$mutualXPerc,decreasing = F)[1],],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text_repel(data=countsTables[order(countsTables$mutualXPerc,decreasing = F)[1],],aes(label=rownames(countsTables[order(countsTables$mutualXPerc,decreasing = F)[1],]),angle=45,hjust=0,vjust=0),size=2.5,y=0.012,show.legend = F) +
  geom_point(data=countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text_repel(data=countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),],aes(label=rownames(countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),]),angle=45,hjust=0,vjust=0),size=2.5,y=0.012,show.legend = F) +
  labs(x='P(Gene is not altered|MYC is altered)',title='Mutually Exclusive Coding Mutation\n relative to MYC Alterations') + 
  theme(plot.title = element_text(hjust = 0.5))

######BOX PLOTS : MYC GEXP vs. top 10 MYC-mutually exclusive genes
top10MXGenes = rownames(countsTables)[order(countsTables$mutualXPerc,decreasing = T)[1:10]]
mycGEXP = mycNetworkGEXP_Wide[rownames(wgAlter),'MYC']
effectDF = matrix(data = NA,nrow = length(top10MXGenes),ncol = 6)
for(ind in 1:length(top10MXGenes))
{
  thisGene = top10MXGenes[ind]
  MYC_GeneAltState = wgAlter[,c(thisGene,'MYC')]
  MYC_GeneAltState[,'Category'] = ifelse((MYC_GeneAltState[,thisGene]==0 & MYC_GeneAltState[,'MYC']==0),'NoneAltered',ifelse((MYC_GeneAltState[,thisGene]==1 & MYC_GeneAltState[,'MYC']==0),paste(thisGene,'Altered',sep = ""),ifelse((MYC_GeneAltState[,thisGene]==0 & MYC_GeneAltState[,'MYC']==1),'MYCAltered','BothAltered')))
  MYC_GeneAltState[,'Category'] = factor(MYC_GeneAltState[,'Category'],levels = c('MYCAltered',paste(thisGene,'Altered',sep = ""),'BothAltered','NoneAltered'))
  MYC_GeneAltState[,'log2MYC'] = log2(mycGEXP+1)

  ##compute effect sizes
  None_Both = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="NoneAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="BothAltered"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  None_MYC = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="NoneAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="MYCAltered"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  None_Gene = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="NoneAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category==paste(thisGene,"Altered",sep = "")],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  Both_MYC = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="BothAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="MYCAltered"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  Both_Gene = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="BothAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category==paste(thisGene,"Altered",sep = "")],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  MYC_Gene = cohen.d(MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category=="MYCAltered"],MYC_GeneAltState$log2MYC[MYC_GeneAltState$Category==paste(thisGene,"Altered",sep = "")],pooled = T,paired = F,na.rm = T,hedges.correction = T)
  
  effectDF[ind,] = c(paste(None_Both$estimate,'(',None_Both$magnitude,')',sep = ""),paste(None_MYC$estimate,'(',None_MYC$magnitude,')',sep = ""),paste(None_Gene$estimate,'(',None_Gene$magnitude,')',sep = ""),paste(Both_MYC$estimate,'(',Both_MYC$magnitude,')',sep = ""),paste(Both_Gene$estimate,'(',Both_Gene$magnitude,')',sep = ""),paste(MYC_Gene$estimate,'(',MYC_Gene$magnitude,')',sep = ""))

}
rownames(effectDF) = top10MXGenes
write.table(effectDF,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/Top10MutualXGenes_MYCGEXP_Boxplots_EffectSizes.tsv',sep = "\t",row.names = T)
#boxplot - run once per gene
thisGene = "PIK3CA"
pdf(paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/',thisGene,'_MYCBoxplots.pdf',sep = ""))

MYC_GeneAltState = wgAlter[,c(thisGene,'MYC')]
MYC_GeneAltState[,'Category'] = ifelse((MYC_GeneAltState[,thisGene]==0 & MYC_GeneAltState[,'MYC']==0),'NoneAltered',ifelse((MYC_GeneAltState[,thisGene]==1 & MYC_GeneAltState[,'MYC']==0),paste(thisGene,'Altered',sep = ""),ifelse((MYC_GeneAltState[,thisGene]==0 & MYC_GeneAltState[,'MYC']==1),'MYCAltered','BothAltered')))
MYC_GeneAltState[,'Category'] = factor(MYC_GeneAltState[,'Category'],levels = c('MYCAltered',paste(thisGene,'Altered',sep = ""),'BothAltered','NoneAltered'))
MYC_GeneAltState[,'log2MYC'] = log2(mycGEXP+1)
thisPlot = ggplot(MYC_GeneAltState, aes(Category,log2MYC,fill=Category)) + 
  geom_boxplot() + 
  geom_signif(comparisons=list(c("MYCAltered", paste(thisGene,"Altered",sep = ""))), annotations = "", y_position = 16, vjust=0.4) +
  geom_signif(comparisons=list(c(paste(thisGene,"Altered",sep = ""), "BothAltered")), annotations = "", y_position = 16, vjust=0.4) +
  geom_signif(comparisons=list(c("BothAltered", "NoneAltered")), annotations = "", y_position = 16, vjust=0.4) +
  geom_signif(comparisons=list(c(paste(thisGene,"Altered",sep = ""), "NoneAltered")), annotations = "", y_position = 17, vjust=0.4) +
  geom_signif(comparisons=list(c("MYCAltered", "BothAltered")), annotations = "", y_position = 18, vjust=0.4) +
  geom_signif(comparisons=list(c("MYCAltered", "NoneAltered")), annotations = "", y_position = 19, vjust=0.4) +
  stat_summary(fun.data = give.n, geom = "text",size=7) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position="bottom",
        legend.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5,size = 20)) +
  labs(x=paste("MYC-",thisGene,"Alterations",sep = ""),y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title="",keyheight = 2))
print(thisPlot)
dev.off()
#### DATA FOR ONCOPRINT #1. MYC Network
#1. MAF
querySql = "SELECT trmSampleBarcode, Hugo_Symbol, HGVSp_Short, Variant_Classification FROM `isb-cgc-04-0007.MYC.PanCanWhitelist_MYCNetwork_MAF`"
mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)
#2. Focal
pancanFocalCNVR_ForOncoprint = pancanFocalCNVR
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA
pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'

#3. retain only those samples from MAF data that also have CN data
mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% rownames(pancanFocalCNVR),]

#4. rbind CN and MAF
alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)

#5. Add samples that have no events
uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
noEventSamples = rownames(pancanFocalCNVR)[!(rownames(pancanFocalCNVR) %in% uniqueAlterSamples)]

write.table(alterForOncoprint,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCan_MYCNetwork_AlterationsForOncoprint.tsv',sep = "\t", row.names = F)
write.table(noEventSamples,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/NoEventSamplesForOncoprint.tsv',sep = "\n", row.names = F)

#### DATA FOR ONCOPRINT #2. Top 10 mutually exclusive genes
#1. MAF
inList = paste(paste(shQuote(top10MXGenes),shQuote('MYC'),sep=","),collapse = ",")
querySql = paste("SELECT trmSampleBarcode,Hugo_Symbol,HGVSp_Short,Variant_Classification",
                 "FROM (",
                 "SELECT SUBSTR(Tumor_SampleBarcode,1,length(Tumor_SampleBarcode)-1) as trmSampleBarcode, Study, Hugo_Symbol, HGVSp_Short , Variant_Classification", 
                 "FROM `isb-cgc-01-0008.Filtered.MC3_MAF_V5_coding`", 
                 "WHERE Variant_Classification in ('Splice_Site','Missense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site','Nonsense_Mutation','In_Frame_Del','Frame_Shift_Del')",
                 "AND Hugo_Symbol IN (",inList,"))",
                 "WHERE trmSampleBarcode In (SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`)",sep = " ")      
 
mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)
#2. wgFocal
top10MXGenes = append(top10MXGenes,'MYC')
pancanFocalCNVR_ForOncoprint = wgFocal
pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[,top10MXGenes]
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA
pancanFocalCNVR_ForOncoprint = as.data.frame(pancanFocalCNVR_ForOncoprint)
pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'

#3. retain only those samples from MAF data that also have CN data
mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% rownames(pancanFocalCNVR),]

#4. rbind CN and MAF
alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)

#5. Add samples that have no events
uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
noEventSamples = rownames(wgFocal)[!(rownames(wgFocal) %in% uniqueAlterSamples)]

write.table(alterForOncoprint,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCan_Top10MXGenes_AlterationsForOncoprint.tsv',sep = "\t", row.names = F)
write.table(noEventSamples,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCan_Top10MXGenes_NoEventSamplesForOncoprint.tsv',sep = "\n", row.names = F)

##MUTUAL EXCLUSIVITY - Method 1 per tumor type
wgAlter = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCanWhitelist_AllGenes_BinaryAltState_FocalCNVR_OR_MAF.tsv',sep = "\t",header=T)
wgAlter[,'Study'] = sampleWhitelist$DISEASE[match(rownames(wgAlter),sampleWhitelist$SAMPLE_BARCODE)]
ggPlots = list()
nMYCAltered = list()
showGenes = list()
pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/TumorTypeSpecific_WholeGenome_MYCMutualXDistribution_RepelText_6In1.pdf',onefile = T)
#mutualXPerStudy = data.frame()
#for(thisStudy in unique(as.vector(wgAlter$Study)))
for(thisStudy in c('OV','ESCA','LUSC','STAD','LUAD','BRCA'))
{
  print(thisStudy)
  thisDF = wgAlter[as.vector(wgAlter$Study)==thisStudy,]
  print(head(thisDF[,1:6]))
  #remove Study column from thisDF
  thisDF = thisDF[,-c(ncol(thisDF))]
  nMYCAltered = c(nMYCAltered,sum(thisDF$MYC))
  #for every gene-MYC pair, count N for GeneNMycN, GeneYMycN, GeneNMycY, GeneYMycY
  print('lapply')
  countsTables = lapply(apply(thisDF,2,function(x) list(table(x,thisDF[,'MYC']))),'[[',1)
  countsTables = as.data.frame(countsTables)
  countsTables = countsTables[,grep('Freq',colnames(countsTables),fixed = T)]
  rownames(countsTables) = c('MycN_GeneN','MycY_GeneN','MycN_GeneY','MycY_GeneY')
  colnames(countsTables) = gsub('.Freq','',colnames(countsTables),fixed = T)
  countsTables['mutualXPerc',] = countsTables['MycY_GeneN',]*100/(countsTables['MycY_GeneN',]+countsTables['MycY_GeneY',])
  countsTables = t(countsTables)
  countsTables = data.frame(countsTables)
  countsTables = countsTables[rownames(countsTables) != 'MYC',]
  #print('rbinding')
  #mutualXPerStudy = rbind(mutualXPerStudy,countsTables$mutualXPerc)
  #rownames(mutualXPerStudy)[nrow(mutualXPerStudy)] = thisStudy
  #write.table(countsTables,file=paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/',thisStudy,'_MYCvsAllGenes_MutualXCounts.tsv',sep = ""),sep = "\t",row.names = T)
  #1. All genes vs. MYC mutualX density plot
  thisShowGenes = countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],]
  thisShowGenes= rbind(thisShowGenes,countsTables[order(countsTables$mutualXPerc,decreasing = F)[1],])
  thisShowGenes = rbind(thisShowGenes,countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),])
  showGenes[[length(showGenes)+1]] = thisShowGenes
  ggPlots = c(ggPlots,list(ggplot(data.frame(countsTables),aes(x=mutualXPerc)) + 
    geom_density(aes(color='blue'),show.legend = F) + 
    geom_point(data=thisShowGenes,aes(x=mutualXPerc),y=0.01,col='blue'))) #+
    #geom_text_repel(data=showGenes,aes(label=rownames(showGenes),angle=45),size=3,y=0.012,show.legend = F) +
    #geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[match(thisStudy,c('OV','ESCA'))]]))) +
    #labs(title=thisStudy,x="",y="") + 
    #theme(plot.title = element_text(hjust = 0.5,size = 0.5))))

  #print(thisPlot)
  
}
grid.arrange(ggPlots[[1]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[1]],aes(label=rownames(showGenes[[1]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="OV"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[1]]))),
             ggPlots[[2]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[2]],aes(label=rownames(showGenes[[2]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="ESCA"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[2]]))),
             ggPlots[[3]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[3]],aes(label=rownames(showGenes[[3]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="LUSC"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[3]]))),
             ggPlots[[4]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[4]],aes(label=rownames(showGenes[[4]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="STAD"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[4]]))),
             ggPlots[[5]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[5]],aes(label=rownames(showGenes[[5]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="LUAD"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[5]]))),
             ggPlots[[6]]+xlim(c(0,100))+geom_text_repel(data=showGenes[[6]],aes(label=rownames(showGenes[[6]]),angle=45),size=3,y=0.012,show.legend = F)+geom_text(aes(x=90,y=Inf,vjust=1.5,label="BRCA"))+geom_text(aes(x=-Inf,y=Inf,hjust=-0.5,vjust=1.5,label=paste('N =',nMYCAltered[[6]]))),
             ncol=1,
             top='Mutual Exclusivity with MYC\n(Relative Focal CNVR OR Coding Mutations)',
             left='Density',
             bottom='P(Gene not altered | MYC altered)')
dev.off()
colnames(mutualXPerStudy) = rownames(countsTables)
#for ggplot violin plot, melt mutualXPerStudy into long form
mutualXPerStudy[,'Study'] = rownames(mutualXPerStudy)
mutualXPerStudy_Long = melt(mutualXPerStudy,id.vars = 'Study')
colnames(mutualXPerStudy_Long) = c('Study','Gene','MutualXPercent')
mutualXPerStudy_Long_6Studies = mutualXPerStudy_Long[mutualXPerStudy_Long$Study %in% c('OV','ESCA','LUSC','STAD','LUAD','BRCA'),]
mutualXPerStudy_Long_6Studies$Study = as.factor(mutualXPerStudy_Long_6Studies$Study)
mutualXPerStudy_Long_6Studies$MutualXPercent = as.numeric(mutualXPerStudy_Long_6Studies$MutualXPercent)
mutualXPerStudy_Long$Study = as.factor(mutualXPerStudy_Long$Study)
mutualXPerStudy_Long$MutualXPercent = as.numeric(mutualXPerStudy_Long$MutualXPercent)
#max gene labels
showGenes = countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],]
#showGenes= rbind(showGenes,countsTables[order(countsTables$mutualXPerc,decreasing = F)[1],])
#showGenes = rbind(showGenes,countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),])
ggplot(mutualXPerStudy_Long_6Studies, aes(x=1,y=MutualXPercent,fill=Study)) + 
  geom_violin(trim = F) +
  geom_boxplot(width=0.1) +
  facet_wrap(~Study,ncol = 3) +
  labs(title="Percent Exclusivity With MYC",y = "P(Gene is not altered|MYC is altered)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Dark2") +
#  scale_x_discrete(limits=c("OV", "BRCA", "STAD",'ESCA','LUSC','LUAD'))
########MUTUAL EXCLUSIVITY - Method 2 - DISCOVER
#write.table(pancanBinaryAlt,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCan_MYCNetwork_BinaryAlt.tsv',sep = "\t",row.names = T,col.names = T)
library(discover)
wgAlter = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCanWhitelist_AllGenes_BinaryAltState_FocalCNVR_OR_MAF.tsv',sep = "\t",header=T)
wgAlterStudy = as.vector(sampleWhitelist$DISEASE[match(rownames(wgAlter),sampleWhitelist$SAMPLE_BARCODE)])
wgAlter = as.matrix(t(wgAlter))
events = discover.matrix(wgAlter)
save(events,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanDiscoverMatrix.rda')
load('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanDiscoverMatrix.rda')


# gene.groups is a grouping vector that is 1 in the entry corresponding to TP53, and 2 for all other genes
gene.groups <- ifelse(rownames(events) == "MYC", 1, 2)
result <- pairwise.discover.test(events, gene.groups,alternative = "greater")
save(result, file = "/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanDiscoverCoOccurrenceResult.rda")

####DISCOVER RESULTS AND ONCOPRINTS##################
load("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanDiscoverCoOccurrenceResult.rda")
resultDF = data.frame()
resultDF = result$p.values['MYC',]
resultDF = rbind(resultDF,result$q.values['MYC',])
rownames(resultDF) = c('pValue','qValue')
resultDF = as.data.frame(t(resultDF))
resultDF = resultDF[order(resultDF$qValue,decreasing = F),]
resultDF_FDR1 = resultDF[resultDF$qValue<=0.01,]
resultDF_FDR1 = resultDF_FDR1[order(resultDF_FDR1$qValue,decreasing = F),]
resultDF_FDR1 = resultDF_FDR1[!is.na(resultDF_FDR1$qValue),]
write.table(resultDF,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanResults_CoOccurrenceWithMYC_SortedByQValue.tsv',sep = "\t",row.names = T,col.names = T)
resultDF_FDR1 = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCanResults_MutualXWithMYC_FDR1.tsv',sep = "\t",header = T)

#### DATA FOR ONCOPRINT. Top 10 mutually exclusive genes
top10MXGenes = rownames(resultDF_FDR1)[1:10]
#1. MAF

inList = paste(paste(shQuote(top10MXGenes),shQuote('MYC'),sep=","),collapse = ",")
querySql = paste("SELECT trmSampleBarcode,Hugo_Symbol,HGVSp_Short,Variant_Classification",
                 "FROM (",
                 "SELECT SUBSTR(Tumor_SampleBarcode,1,length(Tumor_SampleBarcode)-1) as trmSampleBarcode, Study, Hugo_Symbol, HGVSp_Short , Variant_Classification", 
                 "FROM `isb-cgc-01-0008.Filtered.MC3_MAF_V5_coding`", 
                 "WHERE Variant_Classification in ('Splice_Site','Missense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site','Nonsense_Mutation','In_Frame_Del','Frame_Shift_Del')",
                 "AND Hugo_Symbol IN (",inList,"))",
                 "WHERE trmSampleBarcode In (SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`)",sep = " ")      

mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,use_legacy_sql = FALSE,max_pages = Inf)
colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)
#2. wgFocal
top10MXGenes = append(top10MXGenes,'MYC')
pancanFocalCNVR_ForOncoprint = as.data.frame(t(wgFocal)) #wgFocal not binarized, and has genes along rows

pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[,top10MXGenes]
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA

pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'

#3. retain only those samples from MAF data that also have CN data
mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% colnames(wgFocal),]

#4. rbind CN and MAF
alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)

#5. Add samples that have no events
uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
noEventSamples = colnames(wgFocal)[!(colnames(wgFocal) %in% uniqueAlterSamples)]

write.table(alterForOncoprint,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCan_Top10MXGenes_AlterationsForOncoprint.tsv',sep = "\t", row.names = F)
write.table(noEventSamples,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCan_Top10MXGenes_NoEventSamplesForOncoprint.tsv',sep = "\n", row.names = F)

###DISCOVER per tumor-type
#write tumor-type specific data in separate files to run as separate DISCOVER jobs
#library(openxlsx)
#wb = createWorkbook("TumorSpecificDF")
for(thisStudy in unique(wgAlterStudy))
{
  print(thisStudy)
  # thisStudyDF = read.table(paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'_WGAlter.tsv',sep = ""),sep = "\t",header = T)
  # thisStudyDF = as.matrix(t(thisStudyDF))
  # events = discover.matrix(thisStudyDF)
  # save(events,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'DiscoverMatrix.rda',sep = ""))
  # 
  #load(paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'DiscoverMatrix.rda',sep = ""))
  
  # gene.groups is a grouping vector that is 1 in the entry corresponding to MYC, and 2 for all other genes
  #gene.groups <- ifelse(rownames(events) == "MYC", 1, 2)
  #result <- pairwise.discover.test(events, gene.groups,alternative = "greater")
  #save(result, file = paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/",thisStudy,"DiscoverCoOccurrenceResult.rda",sep = ""))
  
  load(paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/",thisStudy,"DiscoverCoOccurrenceResult.rda",sep = ""))
  
  resultDF = data.frame()
  resultDF = result$p.values['MYC',]
  resultDF = rbind(resultDF,result$q.values['MYC',])
  rownames(resultDF) = c('pValue','qValue')
  resultDF = as.data.frame(t(resultDF))
  resultDF = resultDF[order(resultDF$qValue,decreasing = F),]
  #addWorksheet(wb,sheetName = thisStudy)
  #writeData(wb,sheet = thisStudy,x = resultDF,rowNames = T,colNames = T)
  
  #DATA FOR tumor-type specific oncoprints
  resultDF_FDR1 = resultDF[resultDF$qValue<=0.01,]
  resultDF_FDR1 = resultDF_FDR1[!is.na(resultDF_FDR1$qValue),]
  geneList = rownames(resultDF_FDR1)

  if(length(geneList)>0)
  {
    geneList = append(geneList,'MYC')
    #1. MAF

    inList = paste(shQuote(geneList),collapse = ",")
    querySql = paste("SELECT trmSampleBarcode,Hugo_Symbol,HGVSp_Short,Variant_Classification",
                     " FROM (",
                     " SELECT SUBSTR(Tumor_SampleBarcode,1,length(Tumor_SampleBarcode)-1) as trmSampleBarcode, Study, Hugo_Symbol, HGVSp_Short , Variant_Classification",
                     " FROM `isb-cgc-01-0008.Filtered.MC3_MAF_V5_coding`",
                     " WHERE Variant_Classification in ('Splice_Site','Missense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site','Nonsense_Mutation','In_Frame_Del','Frame_Shift_Del')",
                     " AND Hugo_Symbol IN (",inList,")"," AND Study = '",thisStudy,"')",
                     " WHERE trmSampleBarcode In (SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`)",sep = "")

    mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,use_legacy_sql = FALSE,max_pages = Inf)
    colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
    mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)

    #2. wgFocal

    pancanFocalCNVR_ForOncoprint = as.data.frame(t(wgFocal)) #wgFocal not binarized, and has genes along rows

    pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[,geneList]
    pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[wgAlterStudy==thisStudy,]
    pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
    pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
    pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA

    pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
    pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
    colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
    pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'

    #3. retain only those samples from MAF data that also have CN data
    mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% colnames(wgFocal),]

    #4. rbind CN and MAF
    alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)

    #5. Add samples that have no events
    uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
    thisStudyAllSamples = colnames(wgFocal)[wgAlterStudy==thisStudy]
    noEventSamples = thisStudyAllSamples[!(thisStudyAllSamples %in% uniqueAlterSamples)]

    write.table(alterForOncoprint,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'_SignificantGenes_AlterationsForOncoprint.tsv',sep = ""),sep = "\t", row.names = F)
    write.table(noEventSamples,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'_SignificantGenes_NoEventSamplesForOncoprint.tsv',sep=""),sep = "\n", row.names = F)

  }
}
saveWorkbook(wb,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/AllTumorTypeCoOccurrenceResults.xlsx',overwrite = T)







###################################DISCOVER : MYC amplification vs. mutations on all genes######################
wgFocal = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/GISTIC.focal_data_by_genes.conf_95.txt',sep="\t",header = T)

wgFocal = subset(wgFocal,select = -c(Locus.ID,Cytoband))
rownames(wgFocal) = wgFocal$Gene.Symbol
wgFocal = wgFocal[,-c(1)]
colnames(wgFocal) = gsub('.','-',substr(colnames(wgFocal),1,15),fixed = T)

#retain only whitelisted samples
keepList = colnames(wgFocal)[colnames(wgFocal) %in% as.vector(sampleWhitelist$SAMPLE_BARCODE)]
wgFocal = wgFocal[,keepList]

#binarize focal alterations
wgFocal[wgFocal>0] = 1 #retain only amplifications
wgFocal[wgFocal<0] = 0

#transpose so that columns are genes
wgFocal = t(wgFocal)

mycAMP = wgFocal[,'MYC']

#OR with MAF data
querySql = "SELECT SAMPLE_BARCODE,Hugo_Symbol,Mutated FROM `isb-cgc-04-0007.MYC.PanCan_MAF_JOIN_Whitelist`"
mafDF = query_exec(querySql,project = billingProject,max_pages = Inf,use_legacy_sql = F)
mafDF_Wide = spread(mafDF,key=Hugo_Symbol,value=Mutated)
rownames(mafDF_Wide) = mafDF_Wide$SAMPLE_BARCODE
mafDF_Wide = mafDF_Wide[rownames(wgFocal),]

mafDF_Wide = mafDF_Wide[,-c(1)]
mafDF_Wide[!is.na(mafDF_Wide)] = 1
mafDF_Wide[is.na(mafDF_Wide)] = 0

#let MYC be just amplifications
#run DISCOVER for MYC amplifications vs. all MAF - PanCan
library(discover)
wgAlter = mafDF_Wide
wgAlter[,'MYC'] = mycAMP
#get tabular counts for every MYC-gene pair
countsPerPair = apply(wgAlter,2,function(x) table(wgAlter[,'MYC'],x))
countsPerPairDF = as.data.frame(countsPerPair)
countsPerPairDF = countsPerPairDF[,grep('Freq',colnames(countsPerPairDF))]
rownames(countsPerPairDF) = c('MycN_GeneN','MycY_GeneN','MycN_GeneY','MycY_GeneY')
colnames(countsPerPairDF) = gsub('.Freq','',colnames(countsPerPairDF),fixed = T)
countsPerPairDF = t(countsPerPairDF)
#construct DISCOVER event matrix
wgAlter = as.matrix(t(wgAlter))
events = discover.matrix(wgAlter)
save(events,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanDiscoverMatrix_MYCAMP_AllGeneMAF.rda')
load('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanDiscoverMatrix_MYCAMP_AllGeneMAF.rda')

# gene.groups is a grouping vector that is 1 in the entry corresponding to TP53, and 2 for all other genes
gene.groups <- ifelse(rownames(events) == "MYC", 1, 2)
result <- pairwise.discover.test(events, gene.groups,alternative = "greater")
save(result, file = "/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanDiscoverCoOccurrenceResult_MYCAMP_AllGeneMAF.rda")

load("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanDiscoverMutualXResult_MYCAMP_AllGeneMAF.rda")
resultDF = data.frame()
resultDF = result$p.values['MYC',]
resultDF = rbind(resultDF,result$q.values['MYC',])
rownames(resultDF) = c('pValue','qValue')
resultDF = as.data.frame(t(resultDF))
resultDF = resultDF[order(resultDF$qValue,decreasing = F),]
resultDF = cbind(resultDF,countsPerPairDF[match(rownames(resultDF),rownames(countsPerPairDF)),])
resultDF_FDR1 = resultDF[resultDF$qValue<=0.01,]
resultDF_FDR1 = resultDF_FDR1[order(resultDF_FDR1$qValue,decreasing = F),]
resultDF_FDR1 = resultDF_FDR1[!is.na(resultDF_FDR1$qValue),]
write.table(resultDF,'/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanResults_CoOccurrence_MYCAMP_AllGeneMAF_SortedByQValue.tsv',sep = "\t",row.names = T,col.names = T)

#####Pan-Can DISCOVER Oncoprints - MYC AMP vs. all-gene MAF
#### DATA FOR ONCOPRINT. Top 10 mutually exclusive genes
# resultDF = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/PanCanResults_MutualX_MYCAMP_AllGeneMAF_SortedByQValue.tsv', sep = "\t",header = T)
# top10MXGenes = rownames(resultDF_FDR1)[1:10]
# #1. MAF
# 
# inList = paste(paste(shQuote(top10MXGenes),shQuote('MYC'),sep=","),collapse = ",")
# querySql = paste("SELECT trmSampleBarcode,Hugo_Symbol,HGVSp_Short,Variant_Classification",
#                  "FROM (",
#                  "SELECT SUBSTR(Tumor_SampleBarcode,1,length(Tumor_SampleBarcode)-1) as trmSampleBarcode, Study, Hugo_Symbol, HGVSp_Short , Variant_Classification", 
#                  "FROM `isb-cgc-01-0008.Filtered.MC3_MAF_V5_coding`", 
#                  "WHERE Variant_Classification in ('Splice_Site','Missense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site','Nonsense_Mutation','In_Frame_Del','Frame_Shift_Del')",
#                  "AND Hugo_Symbol IN (",inList,"))",
#                  "WHERE trmSampleBarcode In (SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`)",sep = " ")      
# 
# mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,use_legacy_sql = FALSE,max_pages = Inf)
# colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
# mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)
# #2. wgFocal
# top10MXGenes = append(top10MXGenes,'MYC')
# pancanFocalCNVR_ForOncoprint = as.data.frame(t(wgFocal)) #wgFocal not binarized, and has genes along rows
# 
# pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[,top10MXGenes]
# pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
# pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
# pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA
# 
# pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
# pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
# colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
# pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'
# 
# #3. retain only those samples from MAF data that also have CN data
# mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% colnames(wgFocal),]
# 
# #4. rbind CN and MAF
# alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)
# 
# #5. Add samples that have no events
# uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
# noEventSamples = colnames(wgFocal)[!(colnames(wgFocal) %in% uniqueAlterSamples)]
# 
# write.table(alterForOncoprint,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCan_Top10MXGenes_AlterationsForOncoprint.tsv',sep = "\t", row.names = F)
# write.table(noEventSamples,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/PanCan_Top10MXGenes_NoEventSamplesForOncoprint.tsv',sep = "\n", row.names = F)


###DISCOVER MYC AMP vs. all-gene MAF : per tumor-type
#write tumor-type specific data in separate files to run as separate DISCOVER jobs
#library(openxlsx)
#wb = createWorkbook("TumorSpecificDF")
wgAlter = t(wgAlter)
wgAlterStudy = as.vector(sampleWhitelist$DISEASE[match(rownames(wgAlter),as.vector(sampleWhitelist$SAMPLE_BARCODE))])
for(thisStudy in unique(wgAlterStudy))
{
  print(thisStudy)
  #thisStudyDF = wgAlter[wgAlterStudy==thisStudy,]
  #write.table(thisStudyDF,paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/',thisStudy,'_WGAlter.tsv',sep = ""),sep = "\t",col.names = T,row.names = T)
  # thisStudyDF = read.table(paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/',thisStudy,'_WGAlter.tsv',sep = ""),sep = "\t",header = T)
  #thisStudyDF = as.matrix(t(thisStudyDF))
  #events = discover.matrix(thisStudyDF)
  #save(events,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/',thisStudy,'DiscoverMatrix.rda',sep = ""))
  # 
  #load(paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/',thisStudy,'DiscoverMatrix.rda',sep = ""))
  
  #gene.groups is a grouping vector that is 1 in the entry corresponding to MYC, and 2 for all other genes
  #gene.groups <- ifelse(rownames(events) == "MYC", 1, 2)
  #result <- pairwise.discover.test(events, gene.groups,alternative = "less")
  #save(result, file = paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/",thisStudy,"DiscoverMutualXResult.rda",sep = ""))
  
  load(paste("/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/",thisStudy,"DiscoverMutualXResult.rda",sep = ""))
  resultDF = data.frame()
  resultDF = result$p.values['MYC',]
  resultDF = rbind(resultDF,result$q.values['MYC',])
  rownames(resultDF) = c('pValue','qValue')
  resultDF = as.data.frame(t(resultDF))
  resultDF = resultDF[order(resultDF$qValue,decreasing = F),]
  # addWorksheet(wb,sheetName = thisStudy)
  # writeData(wb,sheet = thisStudy,x = resultDF,rowNames = T,colNames = T)

  #DATA FOR tumor-type specific oncoprints
  # resultDF_FDR1 = resultDF[resultDF$qValue<=0.01,]
  # resultDF_FDR1 = resultDF_FDR1[!is.na(resultDF_FDR1$qValue),]
  # geneList = rownames(resultDF_FDR1)
  # 
  # if(length(geneList)>0)
  # {
  #   #geneList = append(geneList,'MYC')
  #   #1. MAF
  # 
  #   inList = paste(shQuote(geneList),collapse = ",")
  #   querySql = paste("SELECT trmSampleBarcode,Hugo_Symbol,HGVSp_Short,Variant_Classification",
  #                    " FROM (",
  #                    " SELECT SUBSTR(Tumor_SampleBarcode,1,length(Tumor_SampleBarcode)-1) as trmSampleBarcode, Study, Hugo_Symbol, HGVSp_Short , Variant_Classification",
  #                    " FROM `isb-cgc-01-0008.Filtered.MC3_MAF_V5_coding`",
  #                    " WHERE Variant_Classification in ('Splice_Site','Missense_Mutation','Frame_Shift_Ins','In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site','Nonsense_Mutation','In_Frame_Del','Frame_Shift_Del')",
  #                    " AND Hugo_Symbol IN (",inList,")"," AND Study = '",thisStudy,"')",
  #                    " WHERE trmSampleBarcode In (SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`)",sep = "")
  # 
  #   mafDF_ForOncoPrint = query_exec(querySql, project=billingProject,use_legacy_sql = FALSE,max_pages = Inf)
  #   colnames(mafDF_ForOncoPrint) = c('Sample','Gene','Alteration','Type')
  #   mafDF_ForOncoPrint$Type = gsub('_Mutation','',mafDF_ForOncoPrint$Type,fixed = T)
  # 
  #   #2. MYC Focal amplification
  # 
  #   pancanFocalCNVR_ForOncoprint = as.data.frame(t(wgFocal)) #wgFocal not binarized, and has genes along rows
  # 
  #   pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[,geneList]
  #   pancanFocalCNVR_ForOncoprint = pancanFocalCNVR_ForOncoprint[wgAlterStudy==thisStudy,]
  #   pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint>0] = 'AMP'
  #   pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint<0] = 'HOMDEL'
  #   pancanFocalCNVR_ForOncoprint[pancanFocalCNVR_ForOncoprint==0] = NA
  # 
  #   pancanFocalCNVR_ForOncoprint[,'Sample'] = rownames(pancanFocalCNVR_ForOncoprint)
  #   pancanFocalCNVR_ForOncoprint_Long = melt(pancanFocalCNVR_ForOncoprint,id.vars = 'Sample',na.rm = T)
  #   colnames(pancanFocalCNVR_ForOncoprint_Long) = c('Sample','Gene','Alteration')
  #   pancanFocalCNVR_ForOncoprint_Long[,'Type'] = 'CNA'
  # 
  #   #3. retain only those samples from MAF data that also have CN data
  #   mafDF_ForOncoPrint = mafDF_ForOncoPrint[mafDF_ForOncoPrint$Sample %in% colnames(wgFocal),]
  # 
  #   #4. rbind CN and MAF
  #   alterForOncoprint = rbind(pancanFocalCNVR_ForOncoprint_Long,mafDF_ForOncoPrint)
  # 
  #   #5. Add samples that have no events
  #   uniqueAlterSamples = unique(as.vector(alterForOncoprint$Sample))
  #   thisStudyAllSamples = colnames(wgFocal)[wgAlterStudy==thisStudy]
  #   noEventSamples = thisStudyAllSamples[!(thisStudyAllSamples %in% uniqueAlterSamples)]
  # 
  #   write.table(alterForOncoprint,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'_SignificantGenes_AlterationsForOncoprint.tsv',sep = ""),sep = "\t", row.names = F)
  #   write.table(noEventSamples,file = paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/TumorTypeSpecificData/',thisStudy,'_SignificantGenes_NoEventSamplesForOncoprint.tsv',sep=""),sep = "\n", row.names = F)
  # 
  # }
}
saveWorkbook(wb,file = '/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/DISCOVER/MYCAMP_AllGeneMAF/TumorTypeSpecificData/AllTumorTypeMutualXResults.xlsx',overwrite = T)



#######################BAR PLOTS : whole genome median alt count given MYC
library(reshape2)
medianAltCountGivenMYC = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/MutualX/FocalCNVR/MYCvsWholeGenomeAlter_DataForBarplotsPerStudy.tsv',sep = "\t",header = T,stringsAsFactors = F)
medianAltCountGivenMYC_Long = melt(medianAltCountGivenMYC,id.vars = c('Study','MYCAltered'))
colnames(medianAltCountGivenMYC_Long) = c('Study','MYCAltered','GenomicAlteration','MedianCount')
medianAltCountGivenMYC_Long$GenomicAlteration = as.vector(medianAltCountGivenMYC_Long$GenomicAlteration)
medianAltCountGivenMYC_Long$MYCAltered[medianAltCountGivenMYC_Long$MYCAltered==1] = 'TRUE'
medianAltCountGivenMYC_Long$MYCAltered[medianAltCountGivenMYC_Long$MYCAltered==0] = 'FALSE'
medianAltCountGivenMYC_Long$GenomicAlteration[medianAltCountGivenMYC_Long$GenomicAlteration=='MedianGenesAltered'] = 'CNVR'
medianAltCountGivenMYC_Long$GenomicAlteration[medianAltCountGivenMYC_Long$GenomicAlteration=='MedianGenesMutated'] = 'MUT'
medianAltCountGivenMYC_Long[,'MedianPercent'] = medianAltCountGivenMYC_Long$MedianCount*100/24204
ggplot(medianAltCountGivenMYC_Long,aes(x=GenomicAlteration,y=MedianPercent,fill=MYCAltered)) +
  geom_bar(stat="identity",position="dodge") + 
  facet_wrap(~Study)
ggplot(medianAltCountGivenMYC_Long[medianAltCountGivenMYC_Long$GenomicAlteration=='CNVR',],aes(x=GenomicAlteration,y=MedianPercent,fill=MYCAltered)) +
  geom_bar(stat="identity",position="dodge") + 
  facet_wrap(~Study)
ggplot(medianAltCountGivenMYC_Long[medianAltCountGivenMYC_Long$GenomicAlteration=='MUT',],aes(x=GenomicAlteration,y=MedianPercent,fill=MYCAltered)) +
  geom_bar(stat="identity",position="dodge") + 
  facet_wrap(~Study)


#######################BAR GRAPH PER PMN GENE####################################

library(openxlsx)
library(reshape2)
library(ggplot2)
geneList = c('MYC','MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL')
pdf('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PMN_AmpDelMutPercentages_Barplots.pdf',onefile = T,width = 9)
for(thisGene in geneList)
{
  thisGeneData = read.xlsx('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PMN_FocalAmpDelMutPercentages_ForBarplots.xlsx',sheet = thisGene)
  thisGeneData_Long = thisGeneData[,c('Study','AMP','DEL','MUT')]
  thisGeneData_Long = melt(thisGeneData_Long,id.vars=1)
  colnames(thisGeneData_Long) = c('Study','Alteration','Percentage')
  #reorder levels of Study factor so that PanCan is last
  thisGeneData_Long$Study = factor(thisGeneData_Long$Study,levels = c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM','PanCan'))
  thisPlot = ggplot(thisGeneData_Long,aes(x=Study,y=Percentage)) +
    geom_bar(stat = "identity",position = "dodge",aes(fill=Alteration)) +
    theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust=0.5)) +
    labs(title=thisGene,y="Percentage Samples")
  print(thisPlot)
}
dev.off()  

