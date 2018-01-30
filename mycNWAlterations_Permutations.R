#library(picante) #for randomizeMatrix
#library(foreach) #for parallel for loop
#library(doParallel) #for parallel for loop
library(tidyr)
library(picante) #for randomizeMatrix
library(foreach) #for parallel for loop
library(doParallel) #for parallel for loop

#INPUT FILE = columns sample, myc network genes, Study
myc_NWGene_BinaryAlt = read.table('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/PanCan_MYCNetwork_BinaryAlt.tsv',sep="\t",header=T) 
# args = commandArgs(trailing=T)
# print(args)
# thisStudy = args[1]
# nPerm = as.numeric(args[2])
# print(nPerm)
thisStudy='ACC'
nPerm=100

myc_NWGene_BinaryAlt_ThisStudy = myc_NWGene_BinaryAlt[myc_NWGene_BinaryAlt$Study==thisStudy,1:ncol(myc_NWGene_BinaryAlt)-1]
rownames(myc_NWGene_BinaryAlt_ThisStudy) = myc_NWGene_BinaryAlt_ThisStudy$SAMPLE
myc_NWGene_BinaryAlt_ThisStudy = myc_NWGene_BinaryAlt_ThisStudy[,-c(1)]
#myc_NWGene_BinaryAlt = subset(myc_NWGene_BinaryAlt_Mat,select=-c(MYC,MYCL,MYCN))

print(head(myc_NWGene_BinaryAlt_ThisStudy))
print(colSums(myc_NWGene_BinaryAlt_ThisStudy))
mutualX_CountsOriginal = as.vector(apply(myc_NWGene_BinaryAlt_ThisStudy,MARGIN = 2,FUN = function(x) table(x,myc_NWGene_BinaryAlt_ThisStudy[,'MYC'])['0','1']))
coAlter_CountsOriginal = as.vector(apply(myc_NWGene_BinaryAlt_ThisStudy,MARGIN = 2,FUN = function(x) ifelse(sum(x)!=0,table(x,myc_NWGene_BinaryAlt_ThisStudy[,'MYC'])['1','1'],0)))

print(mutualX_CountsOriginal)
print(coAlter_CountsOriginal)
#permute matrix 10k times, get MYC vs. Gene x contingency counts each time and save them
mutualX_CountPerPermutation = matrix(nrow=nPerm,ncol=ncol(myc_NWGene_BinaryAlt_ThisStudy))
#coAlter_CountPerPermutation = matrix(nrow=nPerm,ncol=ncol(myc_NWGene_BinaryAlt_Mat)-1)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer
registerDoParallel(cl)

finalCounts <- foreach(i=1:nPerm,.combine = rbind,.packages = "picante",.inorder = F) %dopar%
{
    permMat = randomizeMatrix(myc_NWGene_BinaryAlt_ThisStudy,null.model = "independentswap",iterations = 10000)
    thisPermMutualXCounts = as.vector(apply(permMat,MARGIN = 2,FUN = function(x) table(x,permMat[,'MYC'])['0','1']))
 #   thisPermCoAltCounts = as.vector(apply(permMat[,c(1:ncol(permMat)-1)],MARGIN = 2,FUN = function(x) table(x,permMat[,'MYC'])['1','1']))
    thisPermMutualXCounts
    #list(thisPermMutualXCounts,thisPermCoAltCounts)

}

stopCluster(cl)
#print(finalCounts)
mutualX_CountPerPermutation = finalCounts
#coAlter_CountPerPermutation = do.call(rbind,finalCounts[2,])
colnames(mutualX_CountPerPermutation) = colnames(myc_NWGene_BinaryAlt_ThisStudy)
#colnames(coAlter_CountPerPermutation) = colnames(myc_NWGene_BinaryAlt_Mat)[1:ncol(myc_NWGene_BinaryAlt_Mat)-1]

print(head(mutualX_CountPerPermutation))
#print(head(coAlter_CountPerPermutation))

mutualX_pVal = matrix(nrow=1,ncol=ncol(mutualX_CountPerPermutation))
#coAlter_pVal = rep(-1,ncol(coAlter_CountPerPermutation)) 
for(i in c(1:ncol(mutualX_CountPerPermutation)))
{
    mutualX_pVal[1,i] = (length(mutualX_CountPerPermutation[mutualX_CountPerPermutation[,i] >= mutualX_CountsOriginal[i],i])+1)/nPerm
    #coAlter_pVal[i] = length(coAlter_CountPerPermutation[coAlter_CountPerPermutation[,i] >= coAlter_CountsOriginal[i],i])/nPerm
}
colnames(mutualX_pVal) = colnames(mutualX_CountPerPermutation)
print(mutualX_pVal)
#print(coAlter_pVal)
write.table(mutualX_pVal,file=paste('/Users/varsha/Dropbox (SEngine)/VarshaAnalysis/MYC/ProximalNetwork/FocalCNVR/',thisStudy,'.tsv',sep=""),sep="\t",row.names=F)
