#re-analyze with relative copy numbers and thresholded AMP/DEL status - provided by Andy Cherniack

library(reshape) #for melt
library(ggplot2)
library(bigrquery)
library(tidyr) #for spread
library(gplots) #for heatmap.2
library(effsize) #for cohen.d
billingProject = 'isb-cgc-04-0007'

#######read thresholded AMP/DEL status
pancanRelativeCNVR = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/RelativeCopyNumberCorrectedForPloidy/tcga_pancan_myc_pathway_CN_relativescores.txt',sep = "\t",header = T)
pancanRelativeCNVR = pancanRelativeCNVR[,c('sample','tumor_type','ploidy','MAX_r','MGA_r','MLXIPL_r','MLXIP_r','MLX_r','MNT_r','MXD1_r','MXD3_r','MXD4_r','MXI1_r','MYCL1_r','MYCN_r','MYC_r')]

#assign row names
rownames(pancanRelativeCNVR) = pancanRelativeCNVR[,'sample']
pancanRelativeCNVR = pancanRelativeCNVR[,-c(1)]

#trim rownames to be sample barcodes instead of full aliquot barcodes
rownames(pancanRelativeCNVR) = substring(rownames(pancanRelativeCNVR),1,15)
colnames(pancanRelativeCNVR) = gsub('_r','',colnames(pancanRelativeCNVR))
########read pancan pathway whitelist 
sampleWhitelist = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/PanCanPathwayWhitelist.tsv',sep="\t",header = T)

########filter to keep only whitelisted samples
pancanRelativeCNVR = pancanRelativeCNVR[as.vector(sampleWhitelist$SAMPLE_BARCODE),]
#remove samples with all NAs (missing ploidy call, hence missing relative copy number data)
pancanRelativeCNVR = pancanRelativeCNVR[rowSums(is.na(pancanRelativeCNVR[3:15]))!=13,]

#reduce sampleWhitelist to only those samples that have CN data
sampleWhitelist = sampleWhitelist[sampleWhitelist$SAMPLE_BARCODE %in% rownames(pancanRelativeCNVR),]
sampleCountPerStudy = as.data.frame(table(sampleWhitelist$DISEASE))
names(sampleCountPerStudy) = c('Study','Count')
##########NOTE: the DF still contains some NAs

########plot distribution per gene 

par(mfrow=c(1,1))

#per study distribution
pancanRelativeCNVR = pancanRelativeCNVR[,-c(2)] #remove ploidy
pancanRelativeCNVR_Long = melt(pancanRelativeCNVR,id.vars = 'tumor_type')
colnames(pancanRelativeCNVR_Long) = c('Study','Gene','RelativeCopyNumber')
pancanRelativeCNVR_Long = pancanRelativeCNVR_Long[!is.na(pancanRelativeCNVR_Long$Study),]
par(mar(4,2,4,1))
ggplot(pancanRelativeCNVR_Long, aes(Gene,RelativeCopyNumber)) + 
  geom_violin(scale = "count",draw_quantiles = c(0.5),trim = F) +
  theme(axis.text.x = element_text(size = 12,angle = 90,hjust = 1,vjust = 0.5),axis.title = element_text(size=14),plot.title = element_text(size=18,hjust = 0.5)) +
  geom_hline(yintercept = c(0.5,1,1.5),linetype="dashed",colour=c("blue","orange","red")) +
  labs(title="Pan-Cancer Broad Copy Number Variation")


ggplot(pancanRelativeCNVR_Long, aes(Gene,RelativeCopyNumber)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  guides(fill=guide_legend(title=""))

####START LOW THRESHOLD. RelCN >=1.5 is AMP , RelCN<=0.5 is DEL
####LOW THRESHOLD : heatmap of percent of sample amplified per network gene per tumor-type
pancanAMP = pancanRelativeCNVR_Long
pancanAMP$RelativeCopyNumber[pancanAMP$RelativeCopyNumber<1.5] = 0
pancanAMP$RelativeCopyNumber[pancanAMP$RelativeCopyNumber>=1.5] = 1

pancanAMPCount = aggregate(pancanAMP$RelativeCopyNumber,by = pancanAMP[,c('Study','Gene')],sum,na.rm=T)
names(pancanAMPCount) = c('Study','Gene','Count')
#spread for heatmap
pancanAMPCount_Wide = spread(pancanAMPCount,Gene,Count)
rownames(pancanAMPCount_Wide) = pancanAMPCount_Wide$Study
pancanAMPCount_Wide = pancanAMPCount_Wide[,-c(1)]
#Add sample count per Study as column
pancanAMPCount_Wide[,'TotalN']= sampleCountPerStudy[match(rownames(pancanAMPCount_Wide),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanAMPCount_Wide['PanCan',] = colSums(pancanAMPCount_Wide,na.rm = T)
#compute percentages
pancanAMP_Perc = pancanAMPCount_Wide*100/pancanAMPCount_Wide$TotalN
pancanAMP_Perc = pancanAMP_Perc[,-c(ncol(pancanAMP_Perc))]
pancanAMP_Perc[is.na(pancanAMP_Perc)] = 0
pancanAMP_Perc = as.matrix(pancanAMP_Perc)
pancanAMP_Perc = pancanAMP_Perc[,order(colnames(pancanAMP_Perc))]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.8,4.5)
lwid = c(0.5,4,1.5)
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
heatmap.2(pancanAMP_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,notecol="black",cellnote = ifelse(pancanAMP_Perc==0,NA,round(pancanAMP_Perc,digits = 1)),notecex = 1.5,sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples Amplified\n(Relative Copy Number >=1.5)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.5)

####LOW THRESHOLD : heatmap of percent of samples with deletions per network gene per tumor-type
pancanDEL = pancanRelativeCNVR_Long
pancanDEL$RelativeCopyNumber[pancanDEL$RelativeCopyNumber>0.5] = NA
pancanDEL$RelativeCopyNumber[pancanDEL$RelativeCopyNumber<=0.5] = 1
pancanDEL$RelativeCopyNumber[is.na(pancanDEL$RelativeCopyNumber)] = 0

pancanDELCount = aggregate(pancanDEL$RelativeCopyNumber,by = pancanDEL[,c('Study','Gene')],sum,na.rm=T)
names(pancanDELCount) = c('Study','Gene','Count')
#spread for heatmap
pancanDELCount_Wide = spread(pancanDELCount,Gene,Count)
rownames(pancanDELCount_Wide) = pancanDELCount_Wide$Study
pancanDELCount_Wide = pancanDELCount_Wide[,-c(1)]
#Add sample count per Study as column
pancanDELCount_Wide[,'TotalN']= sampleCountPerStudy[match(rownames(pancanDELCount_Wide),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanDELCount_Wide['PanCan',] = colSums(pancanDELCount_Wide,na.rm = T)
#compute percentages
pancanDEL_Perc = pancanDELCount_Wide*100/pancanDELCount_Wide$TotalN
pancanDEL_Perc = pancanDEL_Perc[,-c(ncol(pancanDEL_Perc))]
pancanDEL_Perc[is.na(pancanDEL_Perc)] = 0
pancanDEL_Perc = as.matrix(pancanDEL_Perc)
pancanDEL_Perc = pancanDEL_Perc[,order(colnames(pancanDEL_Perc))]
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(0.8,4.5)
lwid = c(0.5,4,1.5)
scaleBlue<-colorRampPalette(colors=c("white","blue"))(1000)
heatmap.2(pancanDEL_Perc,breaks = seq(0,100,length.out = 1001),lmat=lmat,lhei=lhei,lwid=lwid,notecol="black",cellnote = ifelse(pancanDEL_Perc==0,NA,round(pancanDEL_Perc,digits = 1)),notecex = 1.5,sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5),trace = "none",main = "Percentage of samples Deleted\n(Relative Copy Number <= 0.5)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.5)
###################################

####START HIGH THRESHOLD. RelCN >=2 is AMP , RelCN<0.5 (= 0) is DEL
####HIGH THRESHOLD : heatmap of percent of sample amplified per network gene per tumor-type
pancanAMP = pancanRelativeCNVR_Long
pancanAMP$RelativeCopyNumber[pancanAMP$RelativeCopyNumber<2] = 0
pancanAMP$RelativeCopyNumber[pancanAMP$RelativeCopyNumber>=2] = 1

pancanAMPCount = aggregate(pancanAMP$RelativeCopyNumber,by = pancanAMP[,c('Study','Gene')],sum,na.rm=T)
names(pancanAMPCount) = c('Study','Gene','Count')
#spread for heatmap
pancanAMPCount_Wide = spread(pancanAMPCount,Gene,Count)
rownames(pancanAMPCount_Wide) = pancanAMPCount_Wide$Study
pancanAMPCount_Wide = pancanAMPCount_Wide[,-c(1)]
#Add sample count per Study as column
pancanAMPCount_Wide[,'TotalN']= sampleCountPerStudy[match(rownames(pancanAMPCount_Wide),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanAMPCount_Wide['PanCan',] = colSums(pancanAMPCount_Wide,na.rm = T)
#compute percentages
pancanAMP_Perc = pancanAMPCount_Wide*100/pancanAMPCount_Wide$TotalN
pancanAMP_Perc = pancanAMP_Perc[,-c(ncol(pancanAMP_Perc))]
pancanAMP_Perc[is.na(pancanAMP_Perc)] = 0
pancanAMP_Perc = as.matrix(pancanAMP_Perc)
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
heatmap.2(pancanAMP_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanAMP_Perc==0,NA,round(pancanAMP_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples Amplified\n(Relative Copy Number >=2)",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.8)

####HIGH THRESHOLD : heatmap of percent of samples with deletions per network gene per tumor-type
pancanDEL = pancanRelativeCNVR_Long
pancanDEL$RelativeCopyNumber[pancanDEL$RelativeCopyNumber>0] = NA
pancanDEL$RelativeCopyNumber[pancanDEL$RelativeCopyNumber<0.5] = 1
pancanDEL$RelativeCopyNumber[is.na(pancanDEL$RelativeCopyNumber)] = 0

pancanDELCount = aggregate(pancanDEL$RelativeCopyNumber,by = pancanDEL[,c('Study','Gene')],sum,na.rm=T)
names(pancanDELCount) = c('Study','Gene','Count')
#spread for heatmap
pancanDELCount_Wide = spread(pancanDELCount,Gene,Count)
rownames(pancanDELCount_Wide) = pancanDELCount_Wide$Study
pancanDELCount_Wide = pancanDELCount_Wide[,-c(1)]
#Add sample count per Study as column
pancanDELCount_Wide[,'TotalN']= sampleCountPerStudy[match(rownames(pancanDELCount_Wide),sampleCountPerStudy$Study),'Count']
#Add PanCan row
pancanDELCount_Wide['PanCan',] = colSums(pancanDELCount_Wide,na.rm = T)
#compute percentages
pancanDEL_Perc = pancanDELCount_Wide*100/pancanDELCount_Wide$TotalN
pancanDEL_Perc = pancanDEL_Perc[,-c(ncol(pancanDEL_Perc))]
#pancanDEL_Perc[is.na(pancanDEL_Perc)] = 0
pancanDEL_Perc = as.matrix(pancanDEL_Perc)
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1,4.5)
lwid = c(0.5,4,1.5)
scaleBlue<-colorRampPalette(colors=c("gray","blue"))(1000)
heatmap.2(pancanDEL_Perc,lmat=lmat,lhei=lhei,lwid=lwid,cellnote = ifelse(pancanDEL_Perc==0,NA,round(pancanDEL_Perc,digits = 1)),sepcolor = "cyan",Rowv=FALSE,Colv = FALSE,dendrogram="none",cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5),trace = "none",main = "Percentage of samples Deleted\n(Relative Copy Number = 0)",na.color = "gray",col = scaleBlue,symbreaks = F,density.info="none",keysize = 0.8)
#############################
###########MUTATIONS
#Proximal network MAF Analysis
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
heatmap.2(mutPercPerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mutPercPerStudy==0,NA,round(mutPercPerStudy,digits = 1)),notecex = 1.5,cexRow=1.5,cexCol = 1.5,adjCol = c(0.8,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nCoding mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.5)

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
pathwayGenes = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/GO_RNAPolymeraseActivity_GeneList.csv',sep = ',',header = T)

#get expression data for these genes
inList = paste(shQuote(pathwayGenes$GO_RNA_POLYMERASE_ACTIVITY),collapse = ",")
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
pancanBinaryAlt = subset(pancanRelativeCNVR,select = -c(tumor_type)) #remove Study column
#include all copy gains and losses and binarize, CN>=1.5 and CN<=0.5
pancanBinaryAlt[pancanBinaryAlt<=0.5] = 2
pancanBinaryAlt[pancanBinaryAlt<1.5] = 0
pancanBinaryAlt[pancanBinaryAlt>=1.5] = 1
#OR mutation data into pancanBinaryAlt
for(index in 1:nrow(mafDF))
{
  pancanBinaryAlt[rownames(pancanBinaryAlt)==mafDF$SAMPLE_BARCODE[index],colnames(pancanBinaryAlt)==mafDF$Hugo_Symbol[index]] = 1
}
#heatmap of alterations per network gene per tumor-type
#add Study column
pancanBinaryAlt[,'Study'] = pancanRelativeCNVR$tumor_type
#melt and aggregate sample counts per gene per study
pancanBinaryAlt_Long = melt(pancanBinaryAlt,id='Study')
names(pancanBinaryAlt_Long) = c('Study','Gene','Altered')
countAlteredPerGenePerStudy = aggregate(pancanBinaryAlt_Long$Altered,by=pancanBinaryAlt_Long[,c('Study','Gene')],sum,na.rm=T)
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
percAlteredPerGenePerStudy = percAlteredPerGenePerStudy[,order(colnames(percAlteredPerGenePerStudy))]
heatmap.2(percAlteredPerGenePerStudy,breaks = seq(0,100,length.out = 1001),Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(percAlteredPerGenePerStudy==0,NA,round(percAlteredPerGenePerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nAMP/DEL/Coding Mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

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
  labs(y="Percentage of Sample Altered",title='MYC Network Alterations') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7),plot.title = element_text(hjust = 0.5)) 

###Boxplots SummaryAlt vs. GEXP
#summarize binary Alt columns into WildType, MYCAltered, AnyButMYCAltered
summaryAlt = function(x) {
  #print(x['MYC'])
  #return()
  if(!any(is.na(x)) & all(x==F,na.rm = T)) # no NAs, all F including MYC
  {
    return('NoAlter')
  }
  if(!is.na(x['MYC']) & x['MYC']==F & sum(x,na.rm = T)!=0) #MYC is F, something else is T
  {
    return('AnyButMYC')
  }else if(!is.na(x['MYC']) & x['MYC']==T)
  {
    return('MYCAlt')
  }else
  {
    return(NA)
  }
  
}
summaryAltPerSample = apply(pancanBinaryAlt[,c(1:13)],MARGIN = 1,FUN = summaryAlt)
pancanBinaryAlt[,'SummaryAlt'] = summaryAltPerSample 
pancanBinaryAlt = pancanBinaryAlt[!is.na(pancanBinaryAlt$SummaryAlt),]
pancanBinaryAlt[,'Log2MYCExp'] = log2(mycNetworkGEXP_Wide$MYC[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]+1)
pancanBinaryAlt[,'Log2MYCMedianNWExp'] = mycNetworkGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(mycNetworkGEXP_Wide))]
pancanBinaryAlt[,'Log2MedianPathwayExp'] = pathwayGEXP_Wide$Log2MedianNetworkGEXP[match(rownames(pancanBinaryAlt),rownames(pathwayGEXP_Wide))]
pancanBinaryAlt = pancanBinaryAlt[!is.na(pancanBinaryAlt$Study),]

#compute one-way ANOVA
panCanModel = conover.test::conover.test(pancanBinaryAlt$Log2MYCMedianNWExp,pancanBinaryAlt$SummaryAlt,method = "bh")
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="AnyButMYC"],pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="NoAlter"],pancanBinaryAlt$Log2MYCMedianNWExp[pancanBinaryAlt$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_NoAlt)
print(Any_MYC)
print(MYC_NoAlt)
###
give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
#1.alterations, MYC gexp, pancan
ggplot(pancanBinaryAlt, aes(SummaryAlt,Log2MYCExp,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="MYC Pathway Alterations",y="Log2 normalized read count",title='MYC Expression') +
  guides(fill=guide_legend(title=""))
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

######################################
##MUTUAL EXCLUSIVITY - METHOD 1 : distribution across whole genome
###################################
#MUTUAL EXCLUSIVITY - Method 1 - genome wide background distribution
wgRelative = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/RelativeCopyNumberCorrectedForPloidy/ABSOLUTE.relative_gene_scores.txt',sep="\t",header = T)

wgRelative = subset(wgRelative,select = -c(tumor_type))
rownames(wgRelative) = wgRelative$sample
wgRelative = wgRelative[,-c(1)]
colnames(wgRelative) = gsub('.','-',substr(colnames(wgRelative),1,15),fixed = T)

#transpose so that samples are along columns
wgRelative = t(wgRelative)
colnames(wgRelative) = substr(colnames(wgRelative),1,15)
#retain only whitelisted samples
keepList = colnames(wgRelative)[colnames(wgRelative) %in% as.vector(sampleWhitelist$SAMPLE_BARCODE)]
wgRelative = wgRelative[,keepList]

#remove rows with all NAs
wgRelative1 = wgRelative[rowSums(is.na(wgRelative))!=ncol(wgRelative),]

#binarize relative alterations
wgRelative1[wgRelative1<=0.5] = 2
wgRelative1[wgRelative1<1.5] = 0
wgRelative1[wgRelative1>=1.5] = 1

#transpose so that columns are genes
wgRelative1 = t(wgRelative1)



#OR with MAF data
querySql = "SELECT SAMPLE_BARCODE,Hugo_Symbol,Mutated FROM `isb-cgc-04-0007.MYC.PanCan_MAF_JOIN_Whitelist`"
mafDF = query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
mafDF_Wide = spread(mafDF,key=Hugo_Symbol,value=Mutated)
rownames(mafDF_Wide) = mafDF_Wide$SAMPLE_BARCODE

#subset rows to match wgRelative1
mafDF_Wide = mafDF_Wide[rownames(wgRelative1),]

mafDF_Wide = mafDF_Wide[,-c(1)]
mafDF_Wide[!is.na(mafDF_Wide)] = 1
mafDF_Wide[is.na(mafDF_Wide)] = 0

#Alter = Relative CN OR MAF
wgAlter = wgRelative1

geneOverlap = intersect(colnames(wgRelative1),colnames(mafDF_Wide))
wgAlter[,geneOverlap] = wgRelative1[,geneOverlap] | mafDF_Wide[,geneOverlap]
#
write.table(wgAlter,gzfile('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/RelativeCopyNumberCorrectedForPloidy/PanCanWhitelist_AllGenes_BinaryAltState_RelativeCNVR_OR_MAF.tsv.gz'),sep = "\t",row.names = T)

#for every gene-MYC pair, count N for GeneNMycN, GeneYMycN, GeneNMycY, GeneYMycY
countsTables = lapply(apply(wgAlter,2,function(x) list(table(x,wgAlter[,'MYC']))),'[[',1)
countsTables = as.data.frame(countsTables)
countsTables = countsTables[,grep('Freq',colnames(countsTables),fixed = T)]
rownames(countsTables) = c('MycN_GeneN','MycY_GeneN','MycN_GeneY','MycY_GeneY')
colnames(countsTables) = gsub('.Freq','',colnames(countsTables),fixed = T)
countsTables['mutualXPerc',] = countsTables['MycY_GeneN',]*100/(countsTables['MycY_GeneN',]+countsTables['MycY_GeneY',])
countsTables = t(countsTables)
countsTables = data.frame(countsTables)
countsTables = countsTables[rownames(countsTables) != 'MYC',]
#1. Pan-Can, All genes vs. MYC mutualX density plot

ggplot(data.frame(countsTables),aes(x=mutualXPerc)) + 
  geom_density() + 
  geom_point(data=countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text(data=countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],],aes(label=rownames(countsTables[order(countsTables$mutualXPerc,decreasing = T)[1:3],]),angle=45,hjust=0,vjust=0),size=2.5,y=0.012,show.legend = F,check_overlap = T) +
  geom_point(data=countsTables[order(countsTables$mutualXPerc,decreasing = F)[1:5],],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text(data=countsTables[order(countsTables$mutualXPerc,decreasing = F)[1:5],],aes(label=rownames(countsTables[order(countsTables$mutualXPerc,decreasing = F)[1:5],]),angle=45,hjust=0,vjust=0),size=2.5,y=0.012,show.legend = F,check_overlap = T) +
  geom_point(data=countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),],aes(x=mutualXPerc),y=0.01,col='blue') +
  geom_text(data=countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),],aes(label=rownames(countsTables[c('MYCN','MYCL1','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL'),]),angle=45,hjust=0,vjust=0),check_overlap = T,size=2.5,y=0.012,show.legend = F) +
  labs(x='P(Gene is not altered|MYC is altered)',title='Mutual Exclusivity with MYC\n(Relative Broad CNVR OR Coding Mutations)') + 
  theme(plot.title = element_text(hjust = 0.5))

