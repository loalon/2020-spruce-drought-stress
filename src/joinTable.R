#this file will take Julia's table and add the infomap cluster information
# 200120
# Alonso Serrano
# 
library(dplyr)
library(here)

projectDir <- "/mnt/picea/projects/spruce/vhurry/drought-stress-combined/"
dataDir <- file.path(projectDir, "data")
resultsDir <- file.path(projectDir, "results")
networkDir <- file.path(projectDir, "network")

table1 <- read.table(file.path(dataDir,"table1.tsv"), sep='\t', header = T,fill=T, quote="", stringsAsFactors = F)
infomapClusters <- read.table(file.path(networkDir,"cluster/InfomapClusters.tsv"), header=T,stringsAsFactors = F)
stats <- read.table(file.path(resultsDir,"seidrStats.txt"), header=T,stringsAsFactors = F)

x <- left_join(table1, infomapClusters, by = c("Congenie.Norway.spruce.ID" = "gene"))
head(x)
names(x)[15] <- "infomapCluster"
write.table(x, file=file.path(resultsDir,"newTable.tsv"), sep='\t', row.names = FALSE, quote = FALSE)


y <- left_join(x, stats, by = c("Congenie.Norway.spruce.ID" = "gene"))
head(y)
write.table(y, file=file.path(resultsDir, "table1_infomap_stats.tsv"), sep='\t', row.names = FALSE, quote = FALSE)
