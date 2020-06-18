# drought network seidr metric extraction
# 200120
# Alonso Serrano

library(tidyverse)
library(here)


projectDir <- "/mnt/picea/projects/spruce/vhurry/drought-stress-combined"
deDir <- file.path(projectDir, "de")
resultsDir <- file.path(projectDir, "results")
networkDir <- file.path(projectDir, "network")
statsFile <- file.path(networkDir, "results/backbone/stats.txt")

stats <- read.table(statsFile, header = T, stringsAsFactors = F,sep='\t')

statsDF1 <- data.frame(gene = stats$Source,
                      PageRank = stats$PageRank_Source,
                      Betweenness = stats$Betweenness_Source,
                      Strength = stats$Strength_Source,
                      Eigenvector = stats$Eigenvector_Source,
                      Katz = stats$Katz_Source)

statsDF2 <- data.frame(gene = stats$Target,
                       PageRank = stats$PageRank_target,
                       Betweenness = stats$Betweenness_target,
                       Strength = stats$Strength_target,
                       Eigenvector = stats$Eigenvector_target,
                       Katz = stats$Katz_target)

statsDF <- rbind(statsDF1, statsDF2)

stats.res <- statsDF[!duplicated(statsDF),]

write.table(stats.res, 
            file=file.path(resultsDir,"seidrStats.txt"),
            row.names = F, quote = F)
