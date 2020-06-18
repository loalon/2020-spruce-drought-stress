#' ---
#' title: "Spruce roots drought stress data analysis"
#' author: "Julia"
#' date: "20190611"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
library(DESeq2)

#' # Intro
#' this is the data preparation to generate seidr networks
#' 
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

projectDir <- "/mnt/picea/projects/spruce/vhurry/drought-stress-combined/"
dataDir <- file.path(projectDir, "data")
deDir <- file.path(projectDir, "de")

#' #Data load
#' Data was prepared by Julia before
#' The data has a shifted header, we use check.names=F
count_needles<-read.csv(file.path(dataDir,"counttable_needles.csv"), check.names=F)
count_roots<-read.csv(file.path(dataDir,"counttable_roots.csv"), check.names=F)



# count_combined <- cbind(count_needles,count_roots)
# samples_combined<-read.csv("droughtcombined.csv")
# samples_combined$Level <- factor(samples_combined$Level)
# samples_combined$Type<-factor(samples_combined$Type)
# samples_combined$Class<-factor(samples_combined$Class)
# group<-factor(paste(samples_combined$Type,samples_combined$Class))
# group2<-factor(paste(samples_combined$Level,samples_combined$Type))
# samples_combined<-cbind(samples_combined,group)
# samples_combined<-cbind(samples_combined,group2)
# 
# catalog.sel <- geneSelect(count_combined,samples_combined$group2,1,2)
# sum(catalog.sel)
# 
# nr_filtered <- count_combined[catalog.sel,]

#count_needles<-count_needles[,-15]
#count_combined <- cbind(count_needles,count_roots)
count_combined <- read.csv(file.path(dataDir,"count_table.csv"), check.names=F)

samples_combined<-read.csv(file.path(dataDir,"droughtcombined.csv"))
samples_combined<-samples_combined[-15,]
samples_combined$Level <- factor(samples_combined$Level)
samples_combined$Type<-factor(samples_combined$Type)
samples_combined$Class<-factor(samples_combined$Class)
samples_combined$group <- factor(paste(samples_combined$Level,samples_combined$Type))
samples_combined$Level<-factor(samples_combined$Level,levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
samples_combined$Type<-factor(samples_combined$Type,levels=c("Needle","Root"))

# group2<-factor(paste(samples_combined$Level,samples_combined$Type))
# samples_combined<-cbind(samples_combined,group)
# samples_combined<-cbind(samples_combined,group2)

catalog.sel <- geneSelect(count_combined,samples_combined$group,1,2)
sum(catalog.sel)

nr_filtered <- count_combined[catalog.sel,]

#save background for enrichment
write.table(rownames(nr_filtered), file=file.path(deDir,"background.txt"),row.names = F,col.names = F, quote = F )

dds <- DESeqDataSetFromMatrix(as.matrix(nr_filtered), samples_combined, ~group)
vsd.aware <- varianceStabilizingTransformation(dds, blind=FALSE)
#base level is not asjusted

combinedData <- t(assay(vsd.aware))
#' # Filter by zeroMAD
combinedData <- combinedData[,which(colMads(as.matrix(combinedData)) > 0)]

#write the files, 1 with names, other headless with data
write.table(combinedData, file =file.path(deDir,"vsdDroughtCombined.tsv"), sep='\t', row.names=FALSE, col.names=FALSE)
writeLines(colnames(combinedData), con=file.path(deDir,"vsdDroughtCombinedNames.txt"))
save(combinedData, file= file.path(deDir,"combinedData.RData"))
write.table(combinedData, file =file.path(deDir,"vsdDrought.tsv"), sep='\t', quote = F, col.names = NA)


