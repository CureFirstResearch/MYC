library(bigrquery)
library(lattice)
library(tidyr)

#get GISTIC calls for MYC, MYCL, MYCN for all tumor type samples 
#NOTE: looking at all sample types for now (solid, normal etc.)
#TODO: get this data from PanCan mirror 
#TODO: do JOINS on AliquotBarcodes or SampleBarcodes instead of ParticipantBarcodes

billingProject = 'isb-cgc-04-0007' 
# querySql <- paste("SELECT ParticipantBarcode,Study,Gene_Symbol,GISTIC_Calls",      
#                   "FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]",
#                   "WHERE Gene_Symbol IN ('MYC','MYCL','MYCN')",sep=" ")
# 
# gisticDF <- query_exec(querySql, project=billingProject)
# write.table(gisticDF,file="/Users/varsha/Documents/SEngine/MYC/MYC_GISTIC_Calls.tsv",sep = "\t")
#density plots
gisticDF = read.table("/Users/varsha/Documents/SEngine/MYC/MYC_GISTIC_Calls.tsv",header = TRUE,sep = "\t")
par(mfrow=c(4,1))
plot(density(gisticDF$GISTIC_Calls),main='MYC/MYCL/MYCN\nCopy Number Distribution\n(Pan-Cancer)',xlab="") 
plot(density((gisticDF[gisticDF$Gene_Symbol=='MYC',])$GISTIC_Calls),main=paste('MYC Copy Number Distribution\n(Pan-Cancer)',sep=" "),xlab="")
plot(density((gisticDF[gisticDF$Gene_Symbol=='MYCL',])$GISTIC_Calls),main=paste('MYCL Copy Number Distribution\n(Pan-Cancer)',sep=" "),xlab="")
plot(density((gisticDF[gisticDF$Gene_Symbol=='MYCN',])$GISTIC_Calls),main=paste('MYCN Copy Number Distribution\n(Pan-Cancer)',sep=" "),xlab = "GISTIC Copy Number estimation")

#per tumor-type boxplots
par(mfrow=c(4,1))
boxplot(gisticDF$GISTIC_Calls~gisticDF$Study,main=paste('MYC/MYCL/MYCN\nCopy Number Distribution\n(Per Tumor Type)',sep=" "),ylab="",las=2)

for (gene in c('MYC','MYCL','MYCN'))
{

  gisticThisGene = gisticDF[gisticDF$Gene_Symbol==gene,c(1,2,4)]
  boxplot(gisticThisGene$GISTIC_Calls~gisticThisGene$Study,main=paste(gene,'Copy Number Distribution\n(Per Tumor Type)',sep=" "),ylab="",las=2)
}

#mean vs. variance (per tumor type)
# querySql = paste("SELECT Gene_Symbol,Study,avg(GISTIC_Calls) as avgGISTIC,variance(GISTIC_Calls) as varGISTIC,", 
#                   "FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]",
#                   "WHERE Gene_Symbol IN ('MYC','MYCL','MYCN')",
#                   "GROUP BY Gene_Symbol,Study",sep = " ")
# resultDF = query_exec(querySql,project = billingProject)
# write.table(resultDF,file="/Users/varsha/Documents/SEngine/MYC/MYC_MeanVsVariance_TumorSpecific.tsv",sep = "\t")

resultDF = read.table(file = "/Users/varsha/Documents/SEngine/MYC/MYC_MeanVsVariance_TumorSpecific.tsv",header = TRUE,sep="\t")
xyplot(resultDF$avgGISTIC/resultDF$varGISTIC~as.factor(resultDF$Study),groups = resultDF$Gene_Symbol,auto.key = TRUE,cex=1.5,xlab='Study',ylab='mean(CNV)/var(CNV)',scales = list(x=list(rot=45)))
xyplot(resultDF$avgGISTIC/resultDF$varGISTIC~as.factor(resultDF$Study),groups = resultDF$Gene_Symbol,auto.key = list(space="top",just=0.95),cex=1.5,xlab='Study',ylab='mean(CNV)/var(CNV)',pch=c('o','+','*')[as.factor(resultDF$Gene_Symbol)],scales = list(x=list(rot=45)))

#correlation between MYC copy number variation and expression of all ~30k genes
#the following query was run to compute spearman and pearson correlation between MYC CNVR and GEXP of all genes

# SELECT HGNC_gene_symbol,mycGene,CORR(GISTIC_Calls,log2_count) as pearson_corr, CORR(myc_rank,expr_rank) as spearman_corr, count(*) as N
# FROM (
#   SELECT barcode,HGNC_gene_symbol,log2_count,
#   RANK() OVER (PARTITION BY HGNC_gene_symbol,mycGene ORDER BY log2_count ASC) AS expr_rank,
#   mycGene,GISTIC_Calls,
#   RANK() OVER (PARTITION BY HGNC_gene_symbol,mycGene ORDER BY GISTIC_Calls ASC) AS myc_rank,
#   FROM (
#     SELECT gexp.SampleBarcode as barcode,mycGene,GISTIC_Calls,HGNC_gene_symbol,log2(normalized_count +1) as log2_count 
#     FROM [isb-cgc:tcga_201607_beta.mRNA_UNC_HiSeq_RSEM] as gexp
#     JOIN EACH
#     (SELECT SampleBarcode,Gene_Symbol as mycGene,GISTIC_Calls   
#     FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]
#     WHERE Gene_Symbol == 'MYC') as mycCNV
#     ON gexp.SampleBarcode = mycCNV.SampleBarcode))
# GROUP BY HGNC_gene_symbol,mycGene
# HAVING pearson_corr IS NOT NULL OR spearman_corr IS NOT NULL
# ORDER BY spearman_corr

####Computed correlation scores are saved to a local file to avoid rerunning BigQueries
# querySql = "SELECT * FROM [isb-cgc-04-0007:MYC.MYC_CNV_GEXP_CORR]"
# resultDF = query_exec(querySql,project = billingProject)
# resultDF = resultDF[!is.na(resultDF$HGNC_gene_symbol),]  #remove rows where gene is NA/null
# write.table(resultDF,file = "/Users/varsha/Documents/SEngine/MYC/MYC_CNV_AllGEXP_Corr.tsv",sep = "\t")

corrDF = read.table("/Users/varsha/Documents/SEngine/MYC/MYC_CNV_AllGEXP_Corr.tsv",header = TRUE)
par(mfrow=c(1,3))
plot(density(corrDF$spearman_corr),main="Spearman Correlation")
plot(density(corrDF$pearson_corr),main="Pearson Correlation")
plot(corrDF$pearson_corr,corrDF$spearman_corr)

#get top 10 genes based on correlation of their expression with MYC CNVR
topGenes_PosCorr = corrDF$HGNC_gene_symbol[order(corrDF$spearman_corr,decreasing = TRUE)[1:10]]
topPosCorr = corrDF$spearman_corr[order(corrDF$spearman_corr,decreasing = TRUE)[1:10]]
topGenes_NegCorr = corrDF$HGNC_gene_symbol[order(corrDF$spearman_corr,decreasing = FALSE)[1:10]]
topNegCorr = corrDF$spearman_corr[order(corrDF$spearman_corr,decreasing = FALSE)[1:10]]

#get expression values for these genes and GISTIC estimates for MYC - for plotting
# inList = paste(paste(shQuote(topGenes_PosCorr),collapse = ","),paste(shQuote(topGenes_NegCorr),collapse = ","),"'MYC'",sep = ',') #get exp values for MYC as well
# querySql = paste("SELECT gexp.SampleBarcode, HGNC_gene_symbol,normalized_count,mycGene,GISTIC_Calls",
#                  "FROM",
#                  "(SELECT SampleBarcode,HGNC_gene_symbol,normalized_count", 
#                  "FROM [isb-cgc:tcga_201607_beta.mRNA_UNC_HiSeq_RSEM]", 
#                  "WHERE HGNC_gene_symbol IN (",inList,")) as gexp", 
#                  "JOIN EACH",
#                  "(SELECT SampleBarcode,Gene_Symbol as mycGene,GISTIC_Calls",   
#                  "FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]",
#                  "WHERE Gene_Symbol == 'MYC') as mycCNV",
#                  "ON gexp.SampleBarcode = mycCNV.SampleBarcode",sep = " ")
# resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)
# write.table(resultDF,file="/Users/varsha/Documents/SEngine/MYC/GEXPForTopCorrelatedGenes.tsv",sep="\t")
#plot with and without log2 transformation
resultDF = read.table(file="/Users/varsha/Documents/SEngine/MYC/GEXPForTopCorrelatedGenes.tsv",sep="\t",header = TRUE)
#1. transform DF from long to wide format
gexpLongDF = resultDF[,1:3]
#take mean when multiple aliquots of the same sample occur
gexpLongDF = aggregate(normalized_count~gexp_SampleBarcode+HGNC_gene_symbol,data = gexpLongDF,FUN = mean)
#transform to wide matrix format
gexpWideDF = spread(data = gexpLongDF,HGNC_gene_symbol,value = normalized_count,sep = NULL)
rownames(gexpWideDF) = gexpWideDF$gexp_SampleBarcode

mycLongDF = resultDF[,c(1,4,5)]
mycLongDF = aggregate(GISTIC_Calls~gexp_SampleBarcode+mycGene,data = mycLongDF,FUN = mean)
mycWideDF = spread(data = mycLongDF,mycGene,value = GISTIC_Calls,sep = NULL)
rownames(mycWideDF) = mycWideDF$gexp_SampleBarcode
mycWideDF = mycWideDF[rownames(gexpWideDF),] #reorder rows to match gexpWideDF 

par(mfrow=c(2,5))
#top 10 positive correlations
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr)[i]],main=as.vector(topGenes_PosCorr)[i],xlab=paste(expression(rho),'=',round(topPosCorr[i],digits = 4)),ylab="Normalized Count")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr)[i]]), col="blue",lwd=3)
}
#top 10 negative correlations
par(mfrow=c(2,5))
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr)[i]],main=as.vector(topGenes_NegCorr)[i],xlab = paste(expression(rho),'=',round(topNegCorr[i],digits = 4)),ylab = "Normalized Count")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr)[i]]), col="red",lwd=3)
}

#try log2 plots
gexpWideDF = log2(gexpWideDF[2:22]+1)
par(mfrow=c(2,5))
#top 10 positive correlations
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr)[i]],main=as.vector(topGenes_PosCorr)[i],xlab=paste(expression(rho),'=',round(topPosCorr[i],digits = 4)),ylab="log2(normalized_count+1)")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr)[i]]), col="blue",lwd=3)
}
#top 10 negative correlations
par(mfrow=c(2,5))
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr)[i]],main=as.vector(topGenes_NegCorr)[i],xlab = paste(expression(rho),'=',round(topNegCorr[i],digits = 4)),ylab = "log2(normalized_count+1)")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr)[i]]), col="red",lwd=3)
}

#TODO:tumor-specific correlation
# QUERY
# SELECT gexp.Study,HGNC_gene_symbol,mycGene,CORR(GISTIC_Calls,log2_count) as pearson_corr, CORR(myc_rank,expr_rank) as spearman_corr, count(*) as N
# FROM (
#   SELECT barcode,gexp.Study,HGNC_gene_symbol,log2_count,
#   RANK() OVER (PARTITION BY HGNC_gene_symbol,mycGene ORDER BY log2_count ASC) AS expr_rank,
#   mycGene,GISTIC_Calls,
#   RANK() OVER (PARTITION BY HGNC_gene_symbol,mycGene ORDER BY GISTIC_Calls ASC) AS myc_rank,
#   FROM (
#     SELECT gexp.SampleBarcode as barcode,gexp.Study,mycGene,GISTIC_Calls,HGNC_gene_symbol,log2(normalized_count +1) as log2_count 
#     FROM [isb-cgc:tcga_201607_beta.mRNA_UNC_HiSeq_RSEM] as gexp
#     JOIN EACH
#     (SELECT SampleBarcode,Study,Gene_Symbol as mycGene,GISTIC_Calls   
#     FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]
#     WHERE Gene_Symbol == 'MYC') as mycCNV
#     ON gexp.SampleBarcode = mycCNV.SampleBarcode))
# GROUP BY HGNC_gene_symbol,mycGene,gexp.Study
# HAVING pearson_corr IS NOT NULL OR spearman_corr IS NOT NULL
# ORDER BY spearman_corr

# querySql = paste("SELECT gexp_Study,HGNC_gene_symbol,spearman_corr,N",
#                 "FROM [isb-cgc-04-0007:MYC.MYC_CNVR_GEXP_CORR_PerTumor]",sep = " ")
# 
# resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)
# write.table(resultDF,file = "/Users/varsha/Documents/SEngine/MYC/TumorSpecificCorrelations.tsv",sep="\t")

resultDF = read.table(file="/Users/varsha/Documents/SEngine/MYC/TumorSpecificCorrelations.tsv",header = TRUE,sep="\t")
boxplot(resultDF$spearman_corr~resultDF$gexp_Study,main=paste('Myc CNVR ~ All GEXP Correlation\n(Per Tumor Type)',sep=" "),ylab="Spearman Correlation",las=2)

# query to get top 10 correlatons per tumor-type
# SELECT gexp_Study,HGNC_gene_symbol ,spearman_corr 
# FROM (
#   SELECT gexp_Study, HGNC_gene_symbol, spearman_corr , row_number() over (partition by gexp_Study order by spearman_corr DESC) as rn
#   from [isb-cgc-04-0007:MYC.MYC_CNVR_GEXP_CORR_PerTumor] )
# WHERE rn <=10

#get top correlations for UVM
resultDF_UVM = resultDF[resultDF$gexp_Study=="UVM",]
topGenes_PosCorr_UVM = resultDF_UVM$HGNC_gene_symbol[order(resultDF_UVM$spearman_corr,decreasing = TRUE)[1:10]]
topPosCorr_UVM = resultDF_UVM$spearman_corr[order(resultDF_UVM$spearman_corr,decreasing = TRUE)[1:10]]
topGenes_NegCorr_UVM = resultDF_UVM$HGNC_gene_symbol[order(resultDF_UVM$spearman_corr,decreasing = FALSE)[1:10]]
topNegCorr_UVM = resultDF_UVM$spearman_corr[order(resultDF_UVM$spearman_corr,decreasing = FALSE)[1:10]]

#get expression values for these genes and GISTIC estimates for MYC - for plotting
#NOTICE the WHERE clause that selects only UVM samples
#  inList = paste(paste(shQuote(topGenes_PosCorr_UVM),collapse = ","),paste(shQuote(topGenes_NegCorr_UVM),collapse = ","),sep = ',') #get exp values for MYC as well
#  querySql = paste("SELECT gexp.SampleBarcode, HGNC_gene_symbol,normalized_count,mycGene,GISTIC_Calls",
#                    "FROM",
#                    "(SELECT SampleBarcode,HGNC_gene_symbol,normalized_count", 
#                    "FROM [isb-cgc:tcga_201607_beta.mRNA_UNC_HiSeq_RSEM]", 
#                    "WHERE HGNC_gene_symbol IN (",inList,") AND Study=='UVM') as gexp", 
#                    "JOIN EACH",
#                    "(SELECT SampleBarcode,Gene_Symbol as mycGene,GISTIC_Calls",   
#                    "FROM [isb-cgc-01-0008:Filtered.all_CNVR_data_by_gene_filtered]",
#                    "WHERE Gene_Symbol == 'MYC' AND Study=='UVM') as mycCNV",
#                    "ON gexp.SampleBarcode = mycCNV.SampleBarcode",sep = " ")
# resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)
# write.table(resultDF,file="/Users/varsha/Documents/SEngine/MYC/GEXPForTopCorrelatedGenes_UVM.tsv",sep="\t")

resultDF = read.table(file="/Users/varsha/Documents/SEngine/MYC/GEXPForTopCorrelatedGenes_UVM.tsv",sep="\t",header = TRUE)
#1. transform DF from long to wide format
gexpLongDF = resultDF[,1:3]
#take mean when multiple aliquots of the same sample occur
gexpLongDF = aggregate(normalized_count~gexp_SampleBarcode+HGNC_gene_symbol,data = gexpLongDF,FUN = mean)
#transform to wide matrix format
gexpWideDF = spread(data = gexpLongDF,HGNC_gene_symbol,value = normalized_count,sep = NULL)
rownames(gexpWideDF) = gexpWideDF$gexp_SampleBarcode

mycLongDF = resultDF[,c(1,4,5)]
mycLongDF = aggregate(GISTIC_Calls~gexp_SampleBarcode+mycGene,data = mycLongDF,FUN = mean)
mycWideDF = spread(data = mycLongDF,mycGene,value = GISTIC_Calls,sep = NULL)
rownames(mycWideDF) = mycWideDF$gexp_SampleBarcode
mycWideDF = mycWideDF[rownames(gexpWideDF),] #reorder rows to match gexpWideDF 

#try log2 plots
gexpWideDF = log2(gexpWideDF[2:21]+1)
par(mfrow=c(2,5))
#top 10 positive correlations
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr_UVM)[i]],main=as.vector(topGenes_PosCorr_UVM)[i],xlab=paste(expression(rho),'=',round(topPosCorr_UVM[i],digits = 4)),ylab="log2(normalized_count+1)")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_PosCorr_UVM)[i]]), col="blue",lwd=3)
}
#top 10 negative correlations
par(mfrow=c(2,5))
for(i in c(1:10))
{
  plot(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr_UVM)[i]],main=as.vector(topGenes_NegCorr_UVM)[i],xlab = paste(expression(rho),'=',round(topNegCorr_UVM[i],digits = 4)),ylab = "log2(normalized_count+1)")
  lines(lowess(mycWideDF$MYC,gexpWideDF[,as.vector(topGenes_NegCorr_UVM)[i]]), col="red",lwd=3)
}

###GEXP##################
#correlation between MYC CNVR and MYC log2 GEXP
par(mfrow=c(1,1))
mycRho = cor(mycWideDF$MYC,gexpWideDF[,"MYC"],method = "spearman")
plot(mycWideDF$MYC,gexpWideDF[,"MYC"],main="MYC CNVR vs. GEXP",xlab = paste(expression(rho),'=',round(mycRho,digits = 4)),ylab = "log2(normalized_count+1)")
lines(lowess(mycWideDF$MYC,gexpWideDF[,"MYC"]), col="blue",lwd=3)
