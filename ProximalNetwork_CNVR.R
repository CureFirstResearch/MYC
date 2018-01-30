#MYC proximal network analysis
library(bigrquery)
library(lattice)
library(tidyr)
library(gplots) #for heatmap.2
library(ggplot2)
library(corrplot)
library(effsize)
#get GISTIC estimates

billingProject = 'isb-cgc-04-0007' 

querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.ProximalNetwork_GISTIC`"
mycCNVRDF <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
mycCNVRDF_Wide = spread(mycCNVRDF,key=Gene_Symbol,value = GISTIC_Calls)

#plot distribution of GISTIC values per gene
boxplot(mycCNVRDF_Wide[,3:15],las=2,main='Distribution of GISTIC values per gene')
abline(h=c(1,1.5,2),col=c('yellow','orange','red'),lty=2,lwd=2)
ggplot(mycCNVRDF,aes(x=GISTIC_Calls,fill=Gene_Symbol)) + geom_density(alpha=0.5) + ylim(0,5)

#heatmap of GISTIC calls 
#cluster both genes as well as tumor types
rownames(mycCNVRDF_Wide) = mycCNVRDF_Wide[,1]
mycCNVRDF_Wide = data.matrix(mycCNVRDF_Wide[,3:ncol(mycCNVRDF_Wide)],rownames.force = T)

par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1.5,4)
lwid = c(0.5,4,2)
heatmap.2(mycCNVRDF_Wide,cexRow=0.5,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "MYC Proximal Network CNVR\n(GISTIC)",na.color = "gray",symbreaks = F,density.info="none",keysize = 0.8)


#calculate percentage of samples with GISTIC value >1.5 or <-0.5 per gene per tumor type
querySql = paste("SELECT Study, Gene_Symbol, SUM(GISTIC_Calls >0 AND GISTIC_Calls <= 1) as AMP1,",
                 "SUM(GISTIC_Calls >1 AND GISTIC_Calls <=1.5) as AMP1.5,",
                 "SUM(GISTIC_Calls> 1.5 AND GISTIC_Calls <= 2) as AMP2,",
                 "SUM(GISTIC_Calls > 2) as AMPgt2,",
                 "SUM(GISTIC_Calls <0 AND GISTIC_Calls >=-0.5) as DEL0.5,",
                 "SUM(GISTIC_Calls <-0.5 AND GISTIC_Calls >=-1) as DEL1,",
                 "SUM(GISTIC_Calls<-1 AND GISTIC_Calls >= -1.5) as DEL1.5,",
                 "SUM(GISTIC_Calls < -1.5 AND GISTIC_Calls >= -2) as DEL2,",
                 "SUM(GISTIC_Calls < -2) as DELlt2, COUNT(*) AS TOTAL", 
                 "FROM [isb-cgc-04-0007:MYC.ProximalNetwork_GISTIC]", 
                 "GROUP BY Study , Gene_Symbol",sep=" ")
mycCNVRCounts <- query_exec(querySql, project=billingProject,useLegacySql = TRUE,max_pages = Inf)

#sum counts over all tumor types for each gene
mycCNVRPanCanCounts = aggregate(mycCNVRCounts[,3:ncol(mycCNVRCounts)],by = list(Gene = mycCNVRCounts$Gene_Symbol),sum)
rownames(mycCNVRPanCanCounts) = mycCNVRPanCanCounts[,1]
mycCNVRPanCanCounts = mycCNVRPanCanCounts[,2:ncol(mycCNVRPanCanCounts)]
mycCNVRPanCanPerc = data.matrix(mycCNVRPanCanCounts[,1:8]*100/mycCNVRPanCanCounts$TOTAL,rownames.force = T)
mycCNVRPanCanCounts = data.matrix(mycCNVRPanCanCounts[,1:8],rownames.force = T)

#heatmap of PanCan percentages and counts
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1.5,4)
lwid = c(0.5,4,2)
scaleRed<-colorRampPalette(colors=c("red","white"))(1000)
mycCNVRPanCanPerc = mycCNVRPanCanPerc[,c(8,7,6,5,1,2,3,4)]
mycCNVRPanCanCounts = mycCNVRPanCanCounts[,c(8,7,6,5,1,2,3,4)]
heatmap.2(mycCNVRPanCanPerc,cellnote = round(mycCNVRPanCanPerc,digits = 1),Colv=F,dendrogram="row",sepcolor = "cyan",labCol = c('<-2','(-2,-1.5)','(-1.5,-1)','(-1,0)','(0,1)','(1,1.5)','(1.5,2)','>2'),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nper GISTIC range\nper Gene",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.8)
#heatmap.2(mycCNVRPanCanCounts,cellnote = mycCNVRPanCanCounts,Colv=F,dendrogram="row",sepcolor = "cyan",labCol = c('<-2','(-2,-1.5)','(-1.5,-1)','(-1,0)','(0,1)','(1,1.5)','(1.5,2)','>2'),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Counts of samples with GISTIC > 1",na.color = "gray",col = scaleRed,symbreaks = F,density.info="none",keysize = 0.8)


#heatmap showing all tumor types, counts of samples with >1.5 GISTIC estimated value
mycCNVR_AMP = mycCNVRCounts[,c('Study','Gene_Symbol','AMP2','AMPgt2','TOTAL')]
mycCNVR_AMP[,'AMPgt1.5'] = mycCNVR_AMP$AMP2 + mycCNVR_AMP$AMPgt2
mycCNVR_AMP = mycCNVR_AMP[,c('Study','Gene_Symbol','AMPgt1.5','TOTAL')]
#aggregate across all tumor types
ampPanCan = aggregate(mycCNVR_AMP[,c('AMPgt1.5','TOTAL')],by = list(Gene_Symbol=mycCNVR_AMP$Gene_Symbol),sum)
#add 'PanCan' label and append to mycCNVR_AMP for heatmap
ampPanCan[,'Study'] = 'PanCan'
#reorder columns before you can append
ampPanCan = ampPanCan[,c('Study','Gene_Symbol','AMPgt1.5','TOTAL')]
mycCNVR_AMP = rbind(mycCNVR_AMP,ampPanCan)
#compute percentages
mycCNVR_AMP[,'Perc'] = mycCNVR_AMP$AMPgt1.5*100/mycCNVR_AMP$TOTAL
mycCNVR_AMP = mycCNVR_AMP[,c('Study','Gene_Symbol','Perc')]
#spread into 2D matrix
mycCNVR_AMP_Wide = spread(mycCNVR_AMP,key = Gene_Symbol,value = Perc)
rownames(mycCNVR_AMP_Wide) = mycCNVR_AMP_Wide[,1]
mycCNVR_AMP_Wide = data.matrix(mycCNVR_AMP_Wide[,2:ncol(mycCNVR_AMP_Wide)],rownames.force = T)
#reorder rows so that PanCan is the last row
target = rownames(mycCNVR_AMP_Wide)
target = target[target != 'PanCan']
target = append(target,'PanCan') 
mycCNVR_AMP_Wide = mycCNVR_AMP_Wide[match(target,rownames(mycCNVR_AMP_Wide)),]  

#heatmap for amplifications
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)

#cluster both genes as well as tumor types
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1.5,4)
lwid = c(0.5,4,2)

heatmap.2(mycCNVR_AMP_Wide,cellnote = ifelse(mycCNVR_AMP_Wide==0,NA,round(mycCNVR_AMP_Wide,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nGISTIC > 1.5",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.8)
#one without row clustering so that PanCan can be shown in the last row of the heatmap
heatmap.2(mycCNVR_AMP_Wide,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mycCNVR_AMP_Wide==0,NA,round(mycCNVR_AMP_Wide,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nGISTIC > 1.5",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.4)


#spread counts for DEL cases GISTIC < -0.5
mycCNVR_DEL = mycCNVRCounts[,c('Study','Gene_Symbol','DEL1','DEL1_5','DEL2','TOTAL')]
mycCNVR_DEL[,'DELgt0.5'] = mycCNVR_DEL$DEL1 + mycCNVR_DEL$DEL1_5 + mycCNVR_DEL$DEL2
mycCNVR_DEL = mycCNVR_DEL[,c('Study','Gene_Symbol','DELgt0.5','TOTAL')]
#aggregate across all tumor types
delPanCan = aggregate(mycCNVR_DEL[,c('DELgt0.5','TOTAL')],by = list(Gene_Symbol=mycCNVR_DEL$Gene_Symbol),sum)
#add 'PanCan' label and append to mycCNVR_DEL for heatmap
delPanCan[,'Study'] = 'PanCan'
#reorder columns before you can append
delPanCan = delPanCan[,c('Study','Gene_Symbol','DELgt0.5','TOTAL')]
mycCNVR_DEL = rbind(mycCNVR_DEL,delPanCan)
#compute percentages
mycCNVR_DEL[,'Perc'] = mycCNVR_DEL$DELgt0.5*100/mycCNVR_DEL$TOTAL
mycCNVR_DEL = mycCNVR_DEL[,c('Study','Gene_Symbol','Perc')]
#spread into 2D matrix
mycCNVR_DEL_Wide = spread(mycCNVR_DEL,key = Gene_Symbol,value = Perc)
rownames(mycCNVR_DEL_Wide) = mycCNVR_DEL_Wide[,1]
mycCNVR_DEL_Wide = data.matrix(mycCNVR_DEL_Wide[,2:ncol(mycCNVR_DEL_Wide)],rownames.force = T)
#reorder rows so that PanCan is the last row
target = rownames(mycCNVR_DEL_Wide)
target = target[target != 'PanCan']
target = append(target,'PanCan') 
mycCNVR_DEL_Wide = mycCNVR_DEL_Wide[match(target,rownames(mycCNVR_DEL_Wide)),]  

#heatmap for deletions
scaleBlue<-colorRampPalette(colors=c("white","blue"))(1000)

#cluster both genes as well as tumor types
par(mfrow=c(1,1))
lmat = rbind(c(0,3,4),c(2,1,0))
lhei=c(1.5,4)
lwid = c(0.5,4,2)

heatmap.2(mycCNVR_DEL_Wide,cellnote = ifelse(mycCNVR_DEL_Wide==0,NA,round(mycCNVR_DEL_Wide,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nGISTIC < -1",na.color = "gray",symbreaks = F,col=scaleBlue,density.info="none",keysize = 0.6)
#one without row clustering so that PanCan can be shown in the last row of the heatmap
heatmap.2(mycCNVR_DEL_Wide,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mycCNVR_DEL_Wide==0,NA,round(mycCNVR_DEL_Wide,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nGISTIC < -0.5",na.color = "gray",symbreaks = F,col=scaleBlue,density.info="none",keysize = 0.6)


#Proximal netwrok MAF Analysis
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.ProximalNetwork_MAF`"
mafDF <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#count of samples per tumor type that have mutations in proximal network genes
mutCountPerStudy = table(mafDF$Study,mafDF$Hugo_Symbol)
PanCan = colSums(mutCountPerStudy)
mutCountPerStudy = rbind(mutCountPerStudy,PanCan)

#get total sample counts per study
querySql <- "SELECT DISEASE,count(*) as N FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist` GROUP BY DISEASE"
totalCountPerStudy <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
totalCountPanCan = sum(totalCountPerStudy$N)
totalCountPerStudy = rbind(totalCountPerStudy,c('PanCan',totalCountPanCan))

#add as column to mutCount
mutCountPerStudy = data.frame(mutCountPerStudy)
mutCountPerStudy[,'TotalN'] = totalCountPerStudy$N[match(rownames(mutCountPerStudy),totalCountPerStudy$DISEASE)]
#compute percentages
mutPercPerStudy = data.matrix(mutCountPerStudy)*100/as.numeric(mutCountPerStudy$TotalN)
#drop TotalN column
mutPercPerStudy = mutPercPerStudy[,1:ncol(mutPercPerStudy)-1]
heatmap.2(mutPercPerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mutPercPerStudy==0,NA,round(mutPercPerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nCoding mutations",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

#combining AMP, DEL and MUT events
#mycCNVRDF_Wide must be a matrix from line 24
mycSigAlter = (mycCNVRDF_Wide > 1.5) | (mycCNVRDF_Wide < -0.5) #AMP and DEL
#add mutations
for(index in 1:nrow(mafDF))
{
  mycSigAlter[rownames(mycSigAlter)==mafDF$SAMPLE_BARCODE[index],colnames(mycSigAlter)==mafDF$Hugo_Symbol[index]] = T
}
#add Study column
sampleStudy = unique(mycCNVRDF[,c(1,2)])

mycSigAlter = data.frame(mycSigAlter)
############################################################
###MUTUAL EXCLUSIVITY  - METHOD 1
###correlation between network genes
#since these are binary variables (alteration or no alteration), cor() will compute phi constant
#corMat = cor(mycSigAlter)

#diag(corMat) = NA
#corrplot(corMat,method = "circle",col = 'cyan',tl.col = "black",cl.pos = "n",type="upper",diag=F,addCoef.col = T,number.cex = 0.8)
#corrplot(corMat,type = "upper",na.label = "square",na.label.col = "gray",cl.lim = c(min(corMat),max(corMat)))
#for each pair of genes, compute ratio of mutually exclusive events to all true events
chiSquared = matrix(data = NA,nrow=ncol(mycSigAlter),ncol=ncol(mycSigAlter))
chiPValue = matrix(data = NA,nrow=ncol(mycSigAlter),ncol=ncol(mycSigAlter))
pairWiseCorr = matrix(data = NA,nrow=ncol(mycSigAlter),ncol=ncol(mycSigAlter))
rownames(chiSquared) = colnames(mycSigAlter)
colnames(chiSquared) = colnames(mycSigAlter)
rownames(chiPValue) = colnames(mycSigAlter)
colnames(chiPValue) = colnames(mycSigAlter)
rownames(pairWiseCorr) = colnames(mycSigAlter)
colnames(pairWiseCorr) = colnames(mycSigAlter)

for(colIndex in 1:ncol(chiSquared))
{
  for(rowIndex in colIndex:nrow(chiSquared))
  {
    
    pairTable = table(mycSigAlter[,rowIndex],mycSigAlter[,colIndex])
    #ignore F-F cell of table
    #1st subscript is row value, 2nd subscript is column value
    totalN = pairTable['TRUE','TRUE'] + pairTable['TRUE','FALSE'] + pairTable['FALSE','TRUE']
    totalRowT = pairTable['TRUE','FALSE']+pairTable['TRUE','TRUE']
    totalColT = pairTable['TRUE','TRUE']+pairTable['FALSE','TRUE']
    totalRowF = pairTable['FALSE','TRUE']
    totalColF = pairTable['TRUE','FALSE']
    
    expTT = (totalRowT*totalColT)/totalN
    expTF = (totalRowT*totalColF)/totalN
    expFT = (totalRowF*totalColT)/totalN
    
    term1 = ifelse(expTT==pairTable['TRUE','TRUE'],0,((expTT - pairTable['TRUE','TRUE'])^2)/expTT)
    term2 = ifelse(expTF==pairTable['TRUE','FALSE'],0,((expTF - pairTable['TRUE','FALSE'])^2)/expTF)
    term3 = ifelse(expFT==pairTable['FALSE','TRUE'],0,((expFT - pairTable['FALSE','TRUE'])^2)/expFT)
    
    chiSquared[rowIndex,colIndex] = term1+term2+term3
    chiSquared[colIndex,rowIndex] = term1+term2+term3
    chiPValue[rowIndex,colIndex] = pchisq(chiSquared[rowIndex,colIndex],df = 1,lower.tail = F)
    chiPValue[colIndex,rowIndex] = pchisq(chiSquared[rowIndex,colIndex],df = 1,lower.tail = F)
    
    x = cor.test(as.numeric(mycSigAlter[mycSigAlter[,rowIndex]==TRUE | mycSigAlter[,colIndex]==TRUE,rowIndex]),as.numeric(mycSigAlter[mycSigAlter[,rowIndex]==TRUE | mycSigAlter[,colIndex]==TRUE,colIndex]),method="spearman",exact=T)
    pairWiseCorr[rowIndex,colIndex] = pairWiseCorr[colIndex,rowIndex] = x$estimate
    
    
  }
}
#diag(xPercPairwiseMat) = NA
corrplot(-log10(chiPValue),is.corr = F,type="upper",diag = F,p.mat = chiPValue,sig.level = 0.0001,insig = "pch",pch = "*",pch.col = "red",title="MYC Network\nPairwise Chi-squared Test",method = "circle",col = "cyan",tl.col = "black",cl.pos = "n",addCoef.col = T,number.cex = 0.8,tl.cex = 0.7)
par(mar=c(4,4,4,4))
diag(pairWiseCorr) = 1
corrplot(pairWiseCorr,type="upper",diag=F,cl.lim = c(min(pairWiseCorr),max(pairWiseCorr)),tl.col = "black")
######
mycSigAlter[,'Study'] = sampleStudy$Study[match(rownames(mycSigAlter),sampleStudy$SampleBarcode)]
mycSigAlterCountsPerStudy = aggregate(mycSigAlter[,1:13],by = list(Study=mycSigAlter$Study),sum)
rownames(mycSigAlterCountsPerStudy) = mycSigAlterCountsPerStudy[,1]
mycSigAlterCountsPerStudy = mycSigAlterCountsPerStudy[,-c(1)]


#add PanCan row
PanCan = colSums(mycSigAlterCountsPerStudy)
mycSigAlterCountsPerStudy = rbind(mycSigAlterCountsPerStudy,PanCan)
rownames(mycSigAlterCountsPerStudy)[nrow(mycSigAlterCountsPerStudy)] = 'PanCan'
#add TotalN column (assuming it already has PanCan count from lines 194-198 above)
mycSigAlterCountsPerStudy[,'TotalN'] = totalCountPerStudy$N[match(rownames(mycSigAlterCountsPerStudy),totalCountPerStudy$DISEASE)]

#compute percentage
mycSigAlterCountsPerStudy = data.matrix(mycSigAlterCountsPerStudy)*100/as.numeric(mycSigAlterCountsPerStudy$TotalN)
#drop TotalN column
mycSigAlterCountsPerStudy = mycSigAlterCountsPerStudy[,1:ncol(mycSigAlterCountsPerStudy)-1]
heatmap.2(mycSigAlterCountsPerStudy,Rowv=F,Colv=F,dendrogram='none',cellnote = ifelse(mycSigAlterCountsPerStudy==0,NA,round(mycSigAlterCountsPerStudy,digits = 1)),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "Percentage of samples\nAMP/DEL/Mutation",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.6)

#######################################
##compute pairwise stats just for UCEC
mycSigAlterUCEC = mycSigAlter[mycSigAlter$Study=='UCEC',]
mycSigAlterUCEC = mycSigAlterUCEC[,1:ncol(mycSigAlterUCEC)-1]
chiSquared = matrix(data = NA,nrow=ncol(mycSigAlterUCEC),ncol=ncol(mycSigAlterUCEC))
chiPValue = matrix(data = NA,nrow=ncol(mycSigAlterUCEC),ncol=ncol(mycSigAlterUCEC))
pairWiseCorr = matrix(data = NA,nrow=ncol(mycSigAlterUCEC),ncol=ncol(mycSigAlterUCEC))
rownames(chiSquared) = colnames(mycSigAlterUCEC)
colnames(chiSquared) = colnames(mycSigAlterUCEC)
rownames(chiPValue) = colnames(mycSigAlterUCEC)
colnames(chiPValue) = colnames(mycSigAlterUCEC)
for(colIndex in 1:ncol(chiSquared))
{
  for(rowIndex in colIndex:nrow(chiSquared))
  {
    
    pairTable = table(mycSigAlterUCEC[,rowIndex],mycSigAlterUCEC[,colIndex])
    #ignore F-F cell of table
    #1st subscript is row value, 2nd subscript is column value
    totalN = pairTable['TRUE','TRUE'] + pairTable['TRUE','FALSE'] + pairTable['FALSE','TRUE']
    totalRowT = pairTable['TRUE','FALSE']+pairTable['TRUE','TRUE']
    totalColT = pairTable['TRUE','TRUE']+pairTable['FALSE','TRUE']
    totalRowF = pairTable['FALSE','TRUE']
    totalColF = pairTable['TRUE','FALSE']
    
    expTT = (totalRowT*totalColT)/totalN
    expTF = (totalRowT*totalColF)/totalN
    expFT = (totalRowF*totalColT)/totalN
    
    term1 = ifelse(expTT==pairTable['TRUE','TRUE'],0,((expTT - pairTable['TRUE','TRUE'])^2)/expTT)
    term2 = ifelse(expTF==pairTable['TRUE','FALSE'],0,((expTF - pairTable['TRUE','FALSE'])^2)/expTF)
    term3 = ifelse(expFT==pairTable['FALSE','TRUE'],0,((expFT - pairTable['FALSE','TRUE'])^2)/expFT)
    
    chiSquared[rowIndex,colIndex] = term1+term2+term3
    chiSquared[colIndex,rowIndex] = term1+term2+term3
    chiPValue[rowIndex,colIndex] = pchisq(chiSquared[rowIndex,colIndex],df = 1,lower.tail = F)
    chiPValue[colIndex,rowIndex] = pchisq(chiSquared[rowIndex,colIndex],df = 1,lower.tail = F)
    
    x = cor.test(as.numeric(mycSigAlterUCEC[mycSigAlterUCEC[,rowIndex]==TRUE | mycSigAlterUCEC[,colIndex]==TRUE,rowIndex]),as.numeric(mycSigAlterUCEC[mycSigAlterUCEC[,rowIndex]==TRUE | mycSigAlterUCEC[,colIndex]==TRUE,colIndex]),method="spearman",exact=T)
    pairWiseCorr[rowIndex,colIndex] = pairWiseCorr[colIndex,rowIndex] = x$estimate
    
    
  }
}
#diag(xPercPairwiseMat) = NA
corrplot(chiSquared,is.corr = F,type="upper",diag = F,p.mat = chiPValue,sig.level = 0.05,insig = "pch",pch = "*",pch.col = "red",title="MYC Network\nPairwise Chi-squared Test\nUCEC",method = "circle",col = "cyan",tl.col = "black",cl.pos = "n",addCoef.col = T,number.cex = 0.8,tl.cex = 0.7)
#######################################

#association between MYC expression and alterations in MYC pathway
querySql <- "SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol = 'MYC'"
mycGEXP <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#add as column to mycSigAlter
mycSigAlter[mycGEXP$SampleBarcode,'MYC_GEXP'] = mycGEXP$log2_count

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
summaryAltPerSample = apply(mycSigAlter[,c(1:13)],MARGIN = 1,FUN = summaryAlt)
mycSigAlter[,'SummaryAlt'] = summaryAltPerSample 
mycSigAlter = mycSigAlter[!is.na(mycSigAlter$Study),]

#compute one-way ANOVA
panCanModel = conover.test::conover.test(mycSigAlter$MYC_GEXP,mycSigAlter$SummaryAlt,method = "bh")
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="AnyButMYC"],mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="AnyButMYC"],mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="NoAlter"],mycSigAlter$MYC_GEXP[mycSigAlter$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_MYC)
print(Any_NoAlt)
print(MYC_NoAlt)

#do one PanCan boxplot
par(mfrow=c(1,1))

give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
ggplot(mycSigAlter, aes(SummaryAlt,MYC_GEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations") +
  guides(fill=guide_legend(title=""))


#per study
ggplot(mycSigAlter, aes(SummaryAlt,MYC_GEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations") +
  guides(fill=guide_legend(title=""))
  
#######################################

#association between average MYC network expression and alterations in MYC pathway
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
mycNetworkGEXP_Wide[,'AvgNetworkGEXP'] = rowSums(mycNetworkGEXP_Wide,na.rm = T)/ncol(mycNetworkGEXP_Wide)
mycNetworkGEXP_Wide[,'AvgNetworkGEXP'] = log2(mycNetworkGEXP_Wide[,'AvgNetworkGEXP']+1)
#add SummaryAlt and Study columns from mycSigAlter to mycNetworkGEXP 
mycNetworkGEXP_Wide[,'SummaryAlt'] = mycSigAlter[rownames(mycNetworkGEXP_Wide),'SummaryAlt']
mycNetworkGEXP_Wide[,'Study'] = mycSigAlter[rownames(mycNetworkGEXP_Wide),'Study']
mycNetworkGEXP_Wide= mycNetworkGEXP_Wide[!is.na(mycNetworkGEXP_Wide$SummaryAlt),]

#compute one-way ANOVA
panCanModel = conover.test::conover.test(mycNetworkGEXP_Wide$AvgNetworkGEXP,mycNetworkGEXP_Wide$SummaryAlt,method = "bh")
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="AnyButMYC"],mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="AnyButMYC"],mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="NoAlter"],mycNetworkGEXP_Wide$AvgNetworkGEXP[mycNetworkGEXP_Wide$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_MYC)
print(Any_NoAlt)
print(MYC_NoAlt)

######do one PanCan boxplot
par(mfrow=c(1,1))

give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
ggplot(mycNetworkGEXP_Wide, aes(SummaryAlt,AvgNetworkGEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations","Average MYC Network Expression") +
  guides(fill=guide_legend(title=""))


#per study
ggplot(mycNetworkGEXP_Wide, aes(SummaryAlt,AvgNetworkGEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations","Average MYC Network Expression") +
  guides(fill=guide_legend(title=""))


######################
#How do MYC network alterations affect expression of downstream pathways
#analyse for GO_HELICASE_ACTIVITY, RNA Polymerase activity, Translation elongation, Translation Initiation, Chemokine_Activity
helicaseGenes = read.table('/Users/varsha/Documents/SEngine/MYC/ProximalNetwork/GO_CHEMOKINE_ACTIVITY_GeneList.csv',sep = ',',header = T)

#get expression data for these genes
inList = paste(shQuote(helicaseGenes$GO_CHEMOKINE_ACTIVITY),collapse = ",")
querySql <- paste("SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_GEXP_JOIN_Whitelist` WHERE Symbol IN (",inList,")",sep="")
helicaseGEXP <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

helicaseGEXP_Wide = helicaseGEXP[,c('SampleBarcode','Symbol','normalized_count')]
helicaseGEXP_Wide = aggregate(helicaseGEXP_Wide,by = list(SampleBarcode = helicaseGEXP_Wide$SampleBarcode,Symbol=helicaseGEXP_Wide$Symbol),median,na.rm=T)
helicaseGEXP_Wide = helicaseGEXP_Wide[,c(1,2,5)]
#spread mycNetworkGEXP into 2D sample X gene matrix
helicaseGEXP_Wide = spread(helicaseGEXP_Wide,key = Symbol,value = normalized_count)
rownames(helicaseGEXP_Wide) = helicaseGEXP_Wide$SampleBarcode
helicaseGEXP_Wide = helicaseGEXP_Wide[,-c(1)]
#compute average MYC network expression per sample
helicaseGEXP_Wide[,'AvgNetworkGEXP'] = rowSums(helicaseGEXP_Wide,na.rm = T)/ncol(helicaseGEXP_Wide)
helicaseGEXP_Wide[,'AvgNetworkGEXP'] = log2(helicaseGEXP_Wide[,'AvgNetworkGEXP']+1)
#add SummaryAlt and Study columns from mycSigAlter to mycNetworkGEXP 
helicaseGEXP_Wide[,'SummaryAlt'] =mycSigAlter[rownames(helicaseGEXP_Wide),'SummaryAlt'] 
helicaseGEXP_Wide[,'Study'] =mycSigAlter[rownames(helicaseGEXP_Wide),'Study']
helicaseGEXP_Wide= helicaseGEXP_Wide[!is.na(helicaseGEXP_Wide$SummaryAlt),]

#reorder levels of SummaryAlt
#helicaseGEXP_Wide$SummaryAlt = factor(helicaseGEXP_Wide$SummaryAlt,levels=c('NoAlter','MYCAlt','AnyButMYC'))
#compute one-way ANOVA
panCanModel = conover.test::conover.test(helicaseGEXP_Wide$AvgNetworkGEXP,helicaseGEXP_Wide$SummaryAlt,method = "bh")
#compute pairwise effect sizes (Hedge's g to account for different sample sizes)
Any_MYC = cohen.d(helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="AnyButMYC"],helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
Any_NoAlt = cohen.d(helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="AnyButMYC"],helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="NoAlter"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
MYC_NoAlt = cohen.d(helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="NoAlter"],helicaseGEXP_Wide$AvgNetworkGEXP[helicaseGEXP_Wide$SummaryAlt=="MYCAlt"],pooled = T,paired = F,na.rm = T,hedges.correction = T)
print(Any_MYC)
print(Any_NoAlt)
print(MYC_NoAlt)
######do one PanCan boxplot
par(mfrow=c(1,1))

give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
ggplot(helicaseGEXP_Wide, aes(SummaryAlt,AvgNetworkGEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations",y="Mean GO_TRANSLATION_ELONGATION_FACTOR_ACTIVITY Expression") +
  guides(fill=guide_legend(title=""))


#per study
ggplot(helicaseGEXP_Wide, aes(SummaryAlt,AvgNetworkGEXP,fill=SummaryAlt)) + 
  geom_boxplot() + 
  facet_wrap(~Study,ncol=6) +
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  labs(x="MYC Pathway Alterations",y="Mean GO_TRANSLATION_ELONGATION_FACTOR_ACTIVITY Expression") +
  guides(fill=guide_legend(title=""))


###############################
##barplot showing percentage of samples altered per MYC network gene per study
x = mycSigAlter[mycSigAlter$SummaryAlt!='NoAlter',]
x = as.data.frame(table(x$Study))
x[,'TotalN'] = totalCountPerStudy$N[match(x$Var1,totalCountPerStudy$DISEASE)] 
mode(x$TotalN) = "numeric"
x[,'Perc'] = x$Freq*100/x$TotalN

x$Var1 = factor(x$Var1,levels = x$Var1[order(x$Perc,decreasing = T)])

ggplot(x, aes(x = Var1 , y = Perc)) + 
  geom_bar(stat = "identity") + 
  ylim(0,100) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),plot.title = element_text(hjust = 0.5)) + 
  labs(x="Study",y="% samples",title="Alterations in MYC Proximal Network")


##################
#MNT GISTIC distribution
ggplot(as.data.frame(mycCNVRDF_Wide), aes(MNT)) + 
  geom_density() + 
  labs(title='MNT GISTIC estimates',x="GISTIC copy number estimate") + 
  theme(plot.title = element_text(hjust = 0.5))

###################################
###MUTUAL EXCLUSIVITY - Method 2
###get background distribution of co-alteration of MYC with all ~20K genes
###then see how likely co-mutation with MYC network genes is.
querySql = "SELECT * FROM `isb-cgc-04-0007.MYC.PanCan_MYC_AllGene_Count_MutualXAndCoccur_Whitelist`"
myc_AllGene_MutualCount <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
#convert counts into percentages per study
myc_AllGene_MutualPerc = myc_AllGene_MutualCount
myc_AllGene_MutualPerc[,c(3:6)] = myc_AllGene_MutualPerc[,c(3:6)]*100/rowSums(myc_AllGene_MutualPerc[,c(3:6)])
#percentages given MYC is TRUE
myc_AllGene_Mutual_PercGivenMYC = myc_AllGene_MutualCount[,c('Study','Gene_Symbol','MYCTrueGENEFalse','NCoMutated')]
myc_AllGene_Mutual_PercGivenMYC[,c(3,4)] = myc_AllGene_Mutual_PercGivenMYC[,c(3,4)]*100/rowSums(myc_AllGene_Mutual_PercGivenMYC[,c(3,4)])
myc_AllGene_Mutual_PercGivenMYC[,'CoMut_Zscore'] = ave(myc_AllGene_Mutual_PercGivenMYC$NCoMutated,as.factor(myc_AllGene_Mutual_PercGivenMYC$Study),FUN = function(x) (x-median(x,na.rm = T))/mad(x,na.rm = T) )
myc_AllGene_Mutual_PercGivenMYC[,'DifMut_Zscore'] = ave(myc_AllGene_Mutual_PercGivenMYC$MYCTrueGENEFalse,as.factor(myc_AllGene_Mutual_PercGivenMYC$Study),FUN = function(x) (x-median(x,na.rm = T))/mad(x,na.rm = T) )

# ggplot(data=myc_AllGene_Mutual_PercGivenMYC,aes(MYCTrueGENEFalse))+
#   geom_density()+ 
#   geom_point(data=myc_AllGene_Mutual_PercGivenMYC[myc_AllGene_Mutual_PercGivenMYC$Gene_Symbol=="TP53",],aes(x=MYCTrueGENEFalse,y=0.01,colour=Study)) + 
#   geom_point(data=myc_AllGene_Mutual_PercGivenMYC[myc_AllGene_Mutual_PercGivenMYC$Gene_Symbol=="TP53",],aes(x=NCoMutated,y=0.02,colour=Study))

#look at just one tumor type, eg. OV
study="UCEC"
mycNWGenes = c('PIK3CA','PTEN','TP53','MYC','MYCN','MYCL','MAX','MGA','MXD1','MXI1','MXD3','MXD4','MNT','MLX','MLXIP','MLXIPL','MXD2','MondoA','ChREBP')
studyDF = myc_AllGene_Mutual_PercGivenMYC[myc_AllGene_Mutual_PercGivenMYC$Study==study,]
#add z-scores as a column
study_mycNW_DF = studyDF[studyDF$Gene_Symbol %in% mycNWGenes,] #for plotting network genes on density plot

#distribution of CoMutationZscores (for a single study)
ggplot(data=studyDF,aes(CoMut_Zscore))+
  geom_density() +
  geom_point(data=study_mycNW_DF,aes(x=CoMut_Zscore,y=0.01,colour=Gene_Symbol)) +
  labs(title=paste(study,"- MYC Network"),x="Co-Alteration z-score") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = study_mycNW_DF$CoMut_Zscore, y = 0.03, label = study_mycNW_DF$Gene_Symbol,angle=90,size=2)

#distribution of MYCTrueGENEFalse_Zscore (for a singel study)
ggplot(data=studyDF,aes(DifMut_Zscore))+
  geom_density() +
  geom_point(data=study_mycNW_DF,aes(x=DifMut_Zscore,y=0.01,colour=Gene_Symbol)) +
  labs(title=paste(study,"- MYC Network"),x="Mutual Exclusivity z-score") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = study_mycNW_DF$DifMut_Zscore, y = 0.03, label = study_mycNW_DF$Gene_Symbol,angle=90,size=2)

#heatmap of CoMutation  z-scores
mycNW_CoMutZscores = myc_AllGene_Mutual_PercGivenMYC[myc_AllGene_Mutual_PercGivenMYC$Gene_Symbol %in% mycNWGenes,c('Study','Gene_Symbol','CoMut_Zscore')]
mycNW_CoMutZscores = spread(mycNW_CoMutZscores,key = Gene_Symbol,value = CoMut_Zscore)
rownames(mycNW_CoMutZscores) = mycNW_CoMutZscores$Study
mycNW_CoMutZscores = mycNW_CoMutZscores[,-c(1)]
mycNW_CoMutZscores = as.matrix(mycNW_CoMutZscores)
#remove studies where all z-scores are Nan
mycNW_CoMutZscores = mycNW_CoMutZscores[rowSums(is.finite(mycNW_CoMutZscores))>1,]

heatmap.2(mycNW_CoMutZscores,Rowv=F,Colv=F,dendrogram='none',cellnote = round(mycNW_CoMutZscores,digits = 1),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), trace = "none",main = "MYC co-alteration z-scores",na.color = "gray",symbreaks = T,col=bluered(1000),density.info="none",key.par = )

#heatmap of MYCTrueGENEFalse z-scores
mycNW_DifMutZscores = myc_AllGene_Mutual_PercGivenMYC[myc_AllGene_Mutual_PercGivenMYC$Gene_Symbol %in% mycNWGenes,c('Study','Gene_Symbol','DifMut_Zscore')]
mycNW_DifMutZscores = spread(mycNW_DifMutZscores,key = Gene_Symbol,value = DifMut_Zscore)
rownames(mycNW_DifMutZscores) = mycNW_DifMutZscores$Study
mycNW_DifMutZscores = mycNW_DifMutZscores[,-c(1)]
mycNW_DifMutZscores = as.matrix(mycNW_DifMutZscores)
#remove studies where all z-scores are Nan
mycNW_DifMutZscores = mycNW_DifMutZscores[rowSums(is.finite(mycNW_DifMutZscores))>1,]

heatmap.2(mycNW_DifMutZscores,Rowv=F,Colv=F,dendrogram='none',cellnote = round(mycNW_DifMutZscores,digits = 1),cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), trace = "none",main = "MYC Mutual Exclusivity z-scores",na.color = "gray",symbreaks = T,col=bluered(1000),density.info="none",key.par = )
####################################
#MUTUAL EXCLUSIVITY Method 3 : Permutations
#see mycNWAlterations_Permutations.R