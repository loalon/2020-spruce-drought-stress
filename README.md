# Drought transcriptome of Norway spruce
 
Julia C. Haas, Alexander Vergara, Alonso R. Serrano, Sanatkumar Mishra, Vaughan Hurry, Nathaniel R. Street.

## Data exploratory tool
http://terra.upsc.se/2020-spruce-drought-stress/

https://loalon.shinyapps.io/2020-spruce-drought-stress/ (Alternative)

## Repository content:

/doc - samples meta data.
* droughtcombined.csv: combined samples.
* samplesheet_roots.csv: only roots.
* samplesheet_needles.csv: only needles.

/src - scripts folder
* drought_pub_script1.R - RNA-Seq analysis of spruce-drought-stress.
* droughtNetworkAnalysis.Rmd - Clustering and enrichment of the drought network .
* droughtAUC.Rmd - Backbone analysis of the drought combined network.
* seidr.R - Expression data and gene name extraction for Seidr.
* joinTable.R - Create table with DE, Infomap clusters and Seidr stats.
* stats.R - Seidr stats extraction.
* ABAanalysis.Rmd - ABA and cluster analysis.

/TAtool_data -  Recipe to generate data for TAtool app.

