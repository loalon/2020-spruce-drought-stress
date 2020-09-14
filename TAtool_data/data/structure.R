library(here)
library(RCurl)
library(dplyr)
library(stringr)

#structure
options(stringsAsFactors = FALSE)
#
expTitle <- ""

projectDir <- "/mnt/picea/projects/spruce/vhurry/drought-stress-combined"
appDir <- "/mnt/picea/shiny/sites/2020-spruce-drought-stress"

descriptionFile <- here("TAtool_data/data/intro.html")
description <-  readChar(descriptionFile, file.info(descriptionFile)$size)
#description <- str_replace_all(description, "(\\W)", "\\\\\\1")

#start fullData
fullData <- list()
speciesMap <- read.table(file.path(appDir, "config/speciesMap.tsv"),sep='\t', stringsAsFactors = F, header=T, quote = "", allowEscapes = T )

# main store general information about the experiment
fullData$main$title <- expTitle
fullData$main$description <- description
fullData$species <- "Picea abies"
speciesRow <- speciesMap[speciesMap$species == fullData$species,]
fullData$speciesNickname <- speciesRow$nickname
fullData$genePattern <- speciesRow$genePattern
fullData$enrOptions <- unlist(strsplit(speciesRow$enrichment,','))
# if next is "" no external source will be called
fullData$externalWeb <-  speciesRow$externalWeb
fullData$externalURL <-  speciesRow$externalURL

# conditions is a quick way to control how subdataset the experiment has, if only one, name it control
fullData$conditions <- c("combined")
fullData$profileColors[['combined']] <- c("#009444", "#009444", "#009444")

#prepare metadata
meta <- read.table(file.path(projectDir, "droughtcombined.csv"), sep=',', header=T, stringsAsFactors = F)

# delete sample 15 from meta, it was deleted in Julia's analysis
meta <- meta[-15,]

# keep only relevant meta
meta <- meta[,c("Level", "Type", "Class")]

# meta will store specific metadata for each dataset
fullData$meta[['combined']] <- meta

# expData store the expression data for each dataset
load(file.path(projectDir, "de/combinedData.RData"))
fullData$expData[['combined']] <- combinedData

# sample points, time, stage, % humidity, these will be represented in the x axis
fullData$samplePoints[['combined']] <- fullData$meta$combined$Level

# variables, tissues, mutants, control vs fertilised. Each one of this will have a line in the plots
fullData$variables[['combined']] <- fullData$meta$combined$Type

# Enrichment
load(file.path(projectDir, "network/cluster/clusterEnr.RData"))

# enr will store the enrichment results, better to have it precalculated for faster access
fullData$enr <- list() # initializing list prevents "more elements supplied than there are to replace" error
fullData$enr[['combined']] <- clusterEnr

# DE 
load(file.path(projectDir, "de/deResults.RData"))
fullData$de <- list()
fullData$de[['root']]$description <- "Roots only, control vs. each phase"
fullData$de[['root']]$results <- resultsRoots
fullData$de[['needle']]$description <- "Needles only, control vs. each phase"
fullData$de[['needle']]$results <- resultsNeedles
fullData$de[['combined']]$description<- "Roots vs Needles"
fullData$de[['combined']]$results <- resultsCombined
# 
# background
fullData$background[['combined']] <- read.table(file.path(projectDir, "de/background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]

# cluster
# 2 columns, gene and cluster
fullData$clusters[['combined']] <- read.table(file.path(projectDir, "network/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)

# #process images
netImage <- file.path(projectDir, "network/net_legend.PNG")
txt <- base64Encode(readBin(netImage, "raw", file.info(netImage)[1, "size"]), "txt")
fullData$networks[['combined']]$name <- "Combined network"
fullData$networks[['combined']]$description <- "Combined network generated with Seidr, visualization with Infomap"
fullData$networks[['combined']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', txt)

# Metrics 
temp <- read.table(file.path(projectDir, "network/results/backbone/stats.txt"), header=T, sep='\t')
temp2 <- temp %>% select(Target, PageRank_target, Betweenness_target, Strength_target, Eigenvector_target, Katz_target)
temp3 <- temp2 %>% distinct()
rownames(temp3) <- temp3$Target
temp3 <- temp3[ , !(names(temp3) %in% "Target")]
fullData$networks[['combined']]$stats <- temp3
fullData$networks[['combined']]$edgeList <- read.table(file.path(projectDir, "network/cluster", "edgeList.txt"), header=T, sep='\t')

# Save data
#save(fullData, file= file.path(appDir,"data/fullData.RData"))
save(fullData, file = here("TAtool_data/data/fullData.RData"))

# Additional for custom tab
metricsData <-  read.table("/mnt/picea/projects/spruce/vhurry/drought-stress-combined/results/table1_infomap_stats.tsv", 
                           sep = '\t',
                           header=T,
                           quote=""
)

#save(metricsData, file= file.path(appDir,"data/metricsData.RData"))
save(metricsData, file = here("TAtool_data/data/metricsData.RData"))

