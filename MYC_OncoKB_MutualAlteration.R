library(httr)
library(bigrquery)

billingProject = 'isb-cgc-04-0007'

response = GET("http://oncokb.org/api/v1/genes")
respContent = content(response,"parse")
oncoKBGeneList = list()
for(index in 1:length(respContent$data))
{
  oncoKBGeneList[index] = respContent$data[[index]]$hugoSymbol
  
}
oncoKBGeneList = unlist(oncoKBGeneList)


#get CNVR events for all these genes for all PanCan samples.
#GISTIC >1.5 or GISTIC <-0.5
geneList = paste(shQuote(oncoKBGeneList),collapse = ",")
sqlQuery = "SELECT * FROM `isb-cgc-04-0007.MYC.OncoKB_GISTIC_Significant`"
mycCNVRDF <- query_exec(sqlQuery, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#get coding mutations
sqlQuery = "SELECT * FROM `isb-cgc-04-0007.MYC.OncoKB_MAF`"
mycMAF = query_exec(sqlQuery, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#populate binary matrix with AMP/DEL/MUT
#get data freeze sample list
sqlQuery = 'SELECT SAMPLE_BARCODE FROM `isb-cgc-04-0007.MYC.PanCanPathwayWhitelist`'
sampleListDF <- query_exec(sqlQuery, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
sampleList = as.vector(sampleListDF$SAMPLE_BARCODE)
#create empty matrix
altBinary = matrix(data=NA,nrow = length(sampleList),ncol=length(oncoKBGeneList))
rownames(altBinary) = sampleList
colnames(altBinary) = oncoKBGeneList

#loop through mycCNVRDF and populate altBinary
for(index in 1:nrow(altBinary))
{
  print(index)
  thisSampleCNVR = mycCNVRDF[grep(pattern = rownames(altBinary)[index],x = mycCNVRDF$SampleBarcode),c('Gene_Symbol','GISTIC_Calls')]
  altBinary[index,thisSampleCNVR$Gene_Symbol] = T
}

#loop through mycMAF to populate altBinary
for(index in 1:nrow(altBinary))
{
  print(index)
  thisSampleMAF = mycMAF[grep(pattern = rownames(altBinary)[index],x = mycMAF$SAMPLE_BARCODE),c('SAMPLE_BARCODE','Hugo_Symbol')]
  altBinary[index,thisSampleMAF$Hugo_Symbol] = T
}

#write to file because the for loops above are slow
write.table(altBinary,file = "/Users/varsha/Documents/SEngine/MYC/OncoKB_AMP_DEL_MUT_Binary.tsv",sep = "\t",row.names = T,col.names = T)

altBinary = data.frame(altBinary)
altBinary[is.na(altBinary)] = F
###correlation between network genes
#since these are binary variables (alteration or no alteration), cor() will compute phi constant

chiSquared = matrix(data = NA,nrow=ncol(altBinary),ncol=1)
chiPValue = matrix(data = NA,nrow=ncol(altBinary),ncol=1)
pairWiseCorr = matrix(data = NA,nrow=ncol(altBinary),ncol=1)
rownames(chiSquared) = colnames(altBinary)
colnames(chiSquared) = 'Chi-sq'
rownames(chiPValue) = colnames(altBinary)
colnames(chiPValue) = 'pvalue'
rownames(pairWiseCorr) = colnames(altBinary)
colnames(pairWiseCorr) = 'spearman'

for(colIndex in 1:nrow(chiSquared))
{
  #for(rowIndex in colIndex:nrow(chiSquared))
  #{
    
  pairTable = table(altBinary[,'MYC'],altBinary[,rownames(chiSquared)[colIndex]])
  #ignore F-F cell of table
  #1st subscript is row value, 2nd subscript is column value
  FT = ifelse(ncol(pairTable)==2,pairTable[1,2],0)
  TT = ifelse(ncol(pairTable)==2,pairTable[2,2],0)
  TF = pairTable[2,1]
  totalN = FT+TT+TF
  totalRowT = TF+TT
  totalColT = TT+FT
  totalRowF = FT
  totalColF = TF
    
  expTT = (totalRowT*totalColT)/totalN
  expTF = (totalRowT*totalColF)/totalN
  expFT = (totalRowF*totalColT)/totalN
    
  term1 = ifelse(expTT==TT,0,((expTT - TT)^2)/expTT)
  term2 = ifelse(expTF==TF,0,((expTF - TF)^2)/expTF)
  term3 = ifelse(expFT==FT,0,((expFT - FT)^2)/expFT)
    
  chiSquared[colIndex,1] = term1+term2+term3
   
  chiPValue[colIndex,1] = pchisq(chiSquared[colIndex,1],df = 1,lower.tail = F)
  
  x = cor.test(as.numeric(altBinary[altBinary[,'MYC']==TRUE | altBinary[,rownames(chiSquared)[colIndex]]==TRUE,'MYC']),as.numeric(altBinary[altBinary[,'MYC']==TRUE | altBinary[,rownames(chiSquared)[colIndex]]==TRUE,rownames(chiSquared)[colIndex]]),method="spearman",exact=T)
  pairWiseCorr[colIndex,1] = x$estimate
  
}
chiSquared = data.frame(chiSquared)
rownames(chiSquared)[order(chiSquared$Chi.sq,decreasing = T)[1:15]]
chiPValue = data.frame(chiPValue)
rownames(chiPValue)[order(chiPValue$pvalue,decreasing = F)[1:15]]
pairWiseCorr = data.frame(pairWiseCorr)
rownames(pairWiseCorr)[order(pairWiseCorr$spearman,decreasing = F)[1:15]]
