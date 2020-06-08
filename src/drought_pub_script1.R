#RNA-Seq analysis of spruce-drought-stress
##setwd("path_to_your_home_directory")
setwd("/mnt/picea/projects/spruce/vhurry/Drought_stress_CIS_element_analysis/")
# library(here)
# setwd(here("results"))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(plot3D))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(Vennerable))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(colorRamps))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

#' Register the default plot margin
mar <- par("mar")
pal <- brewer.pal(8,"Dark2")
##################################################################################
#needle data
samples_needles<- read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/samplesheet_needles.csv")
samples_needles$Class <- factor(samples_needles$Class)
samples_needles$Level <- factor(samples_needles$Level)
#count-table
count_needles<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/counttable_needles.csv", check.names=F)

#######################################################################
#root data
samples_roots<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/samplesheet_roots.csv")
samples_roots$Class <- factor(samples_roots$Class)
samples_roots$Level <- factor(samples_roots$Level)
count_roots<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/counttable_roots.csv", check.names=F)

#########################################################################
#' Select genes that are expressed above a treshold of counts 'exp' and in at least 'n' replicates
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}
##################################################################
#PCAs
#needles
catalog.sel <- geneSelect(count_needles,samples_needles$Level,1,2)
sum(catalog.sel)
n_filtered <- count_needles[catalog.sel,]
ddsFull <- DESeqDataSetFromMatrix(as.matrix(n_filtered), samples_needles,
                                  formula(~Level))
ddsFull
vsd <- varianceStabilizingTransformation(ddsFull, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_needles)
vst <- vst - min(vst)
vstn_filtered<-vst
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
samples_needles$Level<-factor(samples_needles$Level,levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pc3 <- cbind(pc$x[,1],pc$x[,2],pc$x[,3])
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_needles$Level)],
     pch=17,
     main="PCA of needles drought stress",sub="variance stabilized counts")
text2D(x,y,labels = row.names(pc$x), add= TRUE, colkey= FALSE, cex= 0.5)
legend("topright",inset=c(-0.2,0),pch=17,
       col=pal[1:8],
       legend=levels(samples_needles$Level),pt.cex=1,cex=0.4)
#identify sample 115 of 30% FC as outlier = remove
count_needles<-count_needles[,-15]
samples_needles<-samples_needles[-15,]
#redo PCA
catalog.sel <- geneSelect(count_needles,samples_needles$Level,1,2)
sum(catalog.sel)
n_filtered <- count_needles[catalog.sel,]
ddsFull <- DESeqDataSetFromMatrix(as.matrix(n_filtered), samples_needles,
                                  formula(~Level))
ddsFull
vsd <- varianceStabilizingTransformation(ddsFull, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_needles)
vst <- vst - min(vst)
vstn_filtered<-vst
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
samples_needles$Level<-factor(samples_needles$Level,levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pc3 <- cbind(pc$x[,1],pc$x[,2],pc$x[,3])
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_needles$Level)],
     pch=17,
     main="PCA of needles drought stress",sub="variance stabilized counts")
legend("topright",inset=c(-0.2,0),pch=17,
       col=pal[1:8],
       legend=levels(samples_needles$Level),pt.cex=1,cex=0.4)
############################################################################
#root PCA
catalog.sel <- geneSelect(count_roots,samples_roots$Level,1,2)
sum(catalog.sel)
r_filtered <- count_roots[catalog.sel,]
ddsFull <- DESeqDataSetFromMatrix(as.matrix(r_filtered), samples_roots,
                                  formula(~Level))
ddsFull
vsd <- varianceStabilizingTransformation(ddsFull, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_roots)
vst <- vst - min(vst)
vstr_filtered<-vst
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
samples_roots$Level<-factor(samples_roots$Level,levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pc3 <- cbind(pc$x[,1],pc$x[,2],pc$x[,3])
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_roots$Level)],
     pch=17,
     main="PCA of roots drought stress",sub="variance stabilized counts")
legend("topright",inset=c(-0.2,0),pch=17,
       col=pal[1:8],
       legend=levels(samples_roots$Level),pt.cex=1,cex=0.4)
###########################################################################
#combine data sets
count_combined <- cbind(count_needles,count_roots)
samples_combined<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/droughtcombined.csv")
samples_combined$Level <- factor(samples_combined$Level)
samples_combined$Type<-factor(samples_combined$Type)
samples_combined$Class<-factor(samples_combined$Class)
group<-factor(paste(samples_combined$Type,samples_combined$Class))
group2<-factor(paste(samples_combined$Level,samples_combined$Type))
samples_combined<-cbind(samples_combined,group)
samples_combined<-cbind(samples_combined,group2)

###########################################################################
#PCA of combined data
#genes need to be expressed at least once in 2 replicates of the same treatment, such as 60% FC, in the same sample type to be kept.
catalog.sel <- geneSelect(count_combined,samples_combined$group2,1,2)
sum(catalog.sel)
nr_filtered <- count_combined[catalog.sel,]#49682 genes of 70736 gene models remaining

#' Perform variance stabilisation before PCA analysis
ddsFull <- DESeqDataSetFromMatrix(as.matrix(nr_filtered), samples_combined,
                                  formula(~group))
ddsFull
vsd <- varianceStabilizingTransformation(ddsFull, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count_combined)
vstnr_filtered<-vst
###########################################################################
#calculate the average of the vst values per treatment for use in heatmaps
colnames(vst) <- factor(paste(samples_combined$Level,samples_combined$Type))
myMAmean_2<-sapply(unique(colnames(vst)), function(x) rowMeans( vst[ , grep(x, colnames(vst)), drop=FALSE]) )

#saving Julia vst object 
#write.table(myMAmean_2, file = "VST_means_data.csv", quote = F, sep = "\t", row.names = T)


###########################################################################
#2D plot as overview for both tissues
vst<-vstnr_filtered
vst <- vst - min(vst)
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
samples_combined$Level<-factor(samples_combined$Level,levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
samples_combined$Type<-factor(samples_combined$Type,levels=c("Needle","Root"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pc3 <- cbind(pc$x[,1],pc$x[,2],pc$x[,3])
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
shapes = c(17,16) 
shapes <- shapes[as.numeric(samples_combined$Type)]
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_combined$Level)],
     pch=shapes,
     main="PCA of drought stress",sub="variance stabilized counts")
legend("topright",inset=c(-0.2,0),pch=17,
       col=pal[1:8],
       legend=levels(samples_combined$Level),pt.cex=1,cex=0.4)
shapes = c(17, 16) 
legend("bottomright", inset=c(-0.2,0), pch=shapes,
       legend=levels(samples_combined$Type),cex=0.3,pt.cex=1)

#########################################################################
#DESEQ
#first with needle data only (excluding the outlier sample) and comparing the mean of the classes
#early response= (60,40 and 30%); late response (collapse, C2d), re-irrigation (R), control (80% FC); 
#30%7d not considered

ddsFull <- DESeqDataSetFromMatrix(as.matrix(n_filtered), samples_needles,formula(~Class))
ddsFull
class(ddsFull)

ddsC <- DESeq(ddsFull)
#########################################################################
#threshold for differentially expressed genes 
sig = 0.01
fold = 2 #down-regulated in comparison to control
foldd = -2 #up-regulated in comparison to control
#########################################################################
#first comparison control vs. mild drought (E)
resEn <- results(ddsC,contrast=c("Class","C","E"))
table(sign(resEn$log2FoldChange[resEn$padj <= .01 & !is.na(resEn$padj) & abs(resEn$log2FoldChange)>=2]))
resEn_sig = subset(resEn, padj <= sig  & !is.na(resEn$padj)&fold <= abs(log2FoldChange))

#second comparison control vs. severe drought (L)
resCL<- results(ddsC,contrast=c("Class","C","L"))
table(sign(resCL$log2FoldChange[resCL$padj <= .01 & !is.na(resCL$padj) & abs(resCL$log2FoldChange)>=2]))
resCL_sig = subset(resCL, padj <= sig  & !is.na(resCL$padj)&fold <= abs(log2FoldChange))

#third comparsion control vs. re-irrigation (R)
resCR<- results(ddsC,contrast=c("Class","C","R"))
table(sign(resCR$log2FoldChange[resCR$padj <= .01 & !is.na(resCR$padj) & abs(resCR$log2FoldChange)>=2]))
resCR_sig = subset(resCR, padj <= sig  & !is.na(resCR$padj)&fold <= abs(log2FoldChange))

rownames(resEn_sig)->CE
rownames(resCR_sig)->CR
rownames(resCL_sig)->CL
needleclass <- c(unique(CE),unique(CR),unique(CL))
needleclass <- needleclass[!duplicated(needleclass)]

#separate by up- and down regulated genes
resEn_down = subset(resEn, padj <= sig & !is.na(resEn$padj) & fold <= (log2FoldChange))
resEn_up = subset(resEn, padj <= sig & !is.na(resEn$padj) &foldd >= (log2FoldChange))
resCR_down = subset(resCR, padj <= sig & !is.na(resCR$padj) & fold <= (log2FoldChange))
resCR_up = subset(resCR, padj <= sig & !is.na(resCR$padj) &foldd >= (log2FoldChange))
resCL_down = subset(resCL, padj <= sig & !is.na(resCL$padj) & fold <= (log2FoldChange))
resCL_up = subset(resCL, padj <= sig & !is.na(resCL$padj) &foldd >= (log2FoldChange))
rownames(resEn_up)->Enup
rownames(resEn_down)->Endown
rownames(resCR_up)->CRup
rownames(resCR_down)->CRdown
rownames(resCL_up)->CLup
rownames(resCL_down)->CLdown
needlesup <- c(unique(Enup),unique(CRup),unique(CLup))
needlesup <- needlesup[!duplicated(needlesup)]
needlesupframe<-as.data.frame(needlesup)
needlesdown <- c(unique(Endown),unique(CRdown),unique(CLdown))
needlesdown <- needlesdown[!duplicated(needlesdown)]
needlesdownframe<-as.data.frame(needlesdown)
#checking sizes
dim(resEn_down) #7
dim(resEn_up) #91
dim(resCR_down) #9 
dim(resCR_up) #14
dim(resCL_down) #433
dim(resCL_up) #871

nchar(needlesup) #958
nchar(needlesdown) #446
#what is the overlap between treatment phases in numbers
length(intersect(CE,CL))#14
length(intersect(CE,CR))#3
length(intersect(CL,CR))#7
#Venn diagrams
plot.new()
grid.draw(venn.diagram(list(CE,CL,CR),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))

V <- Venn(SetNames = c("Mild", "Severe","Recovery"), Weight = c(`100`= 83,`010`= 1285,`001`= 15,`110`=12,`011`=5,`101`=1,`111`=2))
plot(V, doWeights = TRUE, type = "circles")

#what is the overlap in more detail, considering up- and down-regulation of genes at the different phases
length(intersect(Enup, CLup))
length(intersect(Enup,CRup))
length(intersect(CLup,CRup))
length(intersect(Endown, CLdown))
length(intersect(Endown, CRdown))
length(intersect(CRdown, CLdown))
length(intersect(Enup, CLdown))
#Venn diagrams 
plot.new()
grid.draw(venn.diagram(list(Enup,CLup,CRup),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))
plot.new()
grid.draw(venn.diagram(list(Endown,CLdown,CRdown),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))
V <- Venn(SetNames = c("Early", "Late","Recovery"), Weight = c(`100`= 79,`010`= 856,`001`= 7,`110`=9,`011`=4,`101`=1,`111`=2))
plot(V, doWeights = TRUE, type = "circles")
V <- Venn(SetNames = c("Early", "Late","Recovery"), Weight = c(`100`= 4,`010`= 430,`001`= 9,`110`=3,`011`=0,`101`=0,`111`=0))
plot(V, doWeights = TRUE, type = "circles")

#overlap not only in numbers but as data_frame with gene identifier
commonneedleup<-Reduce(intersect, list(Enup,CLup,CRup))
common1up<-intersect(Enup,CLup)
common2up<-intersect(Enup,CRup)
common3up<-intersect(CRup,CLup)
EnCLonly<-common1up[!common1up %in% commonneedleup]
EnCRonly<-common2up[!common2up %in% commonneedleup]
CLCRonly<-common3up[!common3up %in% commonneedleup]
commonup<-c(unique(common1up),unique(common2up),unique(common3up))
commonup <- commonup[!duplicated(commonup)]
Enupspecific<-resEn_up[!rownames(resEn_up) %in% commonup, ]
CLupspecific<-resCL_up[!rownames(resCL_up) %in% commonup, ]
CRupspecific<-resCR_up[!rownames(resCR_up) %in% commonup, ]

common1down<-intersect(Endown,CLdown)
common2down<-intersect(Endown,CRdown)
common3down<-intersect(CRdown,CLdown)
Endownspecific<-resEn_down[!rownames(resEn_down) %in% common1down, ]
CLdownspecific<-resCL_down[!rownames(resCL_down) %in% common1down, ]

###############################################################################
###############################################################################
#DESEQ roots
#remove the low expressed genes
ddsFull <- DESeqDataSetFromMatrix(as.matrix(r_filtered), samples_roots,formula(~Class))
ddsFull

ddsCR <- DESeq(ddsFull)

resrootE <- results(ddsCR,contrast=c("Class","C","E"))
table(sign(resrootE$log2FoldChange[resrootE$padj <= .01 & !is.na(resrootE$padj) & abs(resrootE$log2FoldChange)>=2]))
resrootE_sig = subset(resrootE, padj <= sig  & !is.na(resrootE$padj)&fold <= abs(log2FoldChange))

resrootCL<- results(ddsCR,contrast=c("Class","C","L"))
table(sign(resrootCL$log2FoldChange[resrootCL$padj <= .01 & !is.na(resrootCL$padj) & abs(resrootCL$log2FoldChange)>=2]))
resrootCL_sig = subset(resrootCL, padj <= sig  & !is.na(resrootCL$padj)&fold <= abs(log2FoldChange))

resrootCR<- results(ddsCR,contrast=c("Class","C","R"))
table(sign(resrootCR$log2FoldChange[resrootCR$padj <= .01 & !is.na(resrootCR$padj) & abs(resrootCR$log2FoldChange)>=2]))
resrootCR_sig = subset(resrootCR, padj <= sig  & !is.na(resrootCR$padj)&fold <= abs(log2FoldChange))

rownames(resrootE_sig)->CEroot
rownames(resrootCR_sig)->CRroot
rownames(resrootCL_sig)->CLroot
rootclass <- c(unique(CEroot),unique(CRroot),unique(CLroot))
rootclass <- rootclass[!duplicated(rootclass)]
#Venn diagrams
plot.new()
grid.draw(venn.diagram(list(CEroot,CLroot,CRroot),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))

V <- Venn(SetNames = c("Mild", "Severe","Recovery"), Weight = c(`100`= 7,`010`= 5277,`001`= 350,`110`=57,`011`=469,`101`=2,`111`=32))
plot(V, doWeights = TRUE, type = "circles")

#separate by up- and down regulated genes
resrootE_down = subset(resrootE, padj <= sig & !is.na(resrootE$padj) & fold <= (log2FoldChange))
resrootE_up = subset(resrootE, padj <= sig & !is.na(resrootE$padj) &foldd >= (log2FoldChange))
resrootCR_down = subset(resrootCR, padj <= sig & !is.na(resrootCR$padj) & fold <= (log2FoldChange))
resrootCR_up = subset(resrootCR, padj <= sig & !is.na(resrootCR$padj) &foldd >= (log2FoldChange))
resrootCL_down = subset(resrootCL, padj <= sig & !is.na(resrootCL$padj) & fold <= (log2FoldChange))
resrootCL_up = subset(resrootCL, padj <= sig & !is.na(resrootCL$padj) &foldd >= (log2FoldChange))
rownames(resrootE_up)->Erootup
rownames(resrootE_down)->Erootdown
rownames(resrootCR_up)->CRrootup
rownames(resrootCR_down)->CRrootdown
rownames(resrootCL_up)->CLrootup
rownames(resrootCL_down)->CLrootdown

#checking sizes
dim(resrootE_up) #90
dim(resrootE_down) #8
dim(resrootCR_up) #603
dim(resrootCR_down) #250
dim(resrootCL_up) #1926
dim(resrootCL_down) #3908


#Venn diagrams
plot.new()
grid.draw(venn.diagram(list(Erootup,CLrootup,CRrootup),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))
plot.new()
grid.draw(venn.diagram(list(Erootdown,CLrootdown,CRrootdown),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("yellow","red","blue"),category.names=c("Early","Late","Recovery")))
V <- Venn(SetNames = c("Early", "Late","Recovery"), Weight = c(`100`= 8,`010`= 1585,`001`= 309,`110`=50,`011`=262,`101`=3,`111`=29))
plot(V, doWeights = TRUE, type = "circles")
V <- Venn(SetNames = c("Early", "Late","Recovery"), Weight = c(`100`= 0,`010`= 3714,`001`= 61,`110`=6,`011`=187,`101`=0,`111`=2))
plot(V, doWeights = TRUE, type = "circles")

#gene names for DEGs overlapping between phases
commonrootup<-Reduce(intersect, list(Erootup,CLrootup,CRrootup))
common1up<-intersect(Erootup,CLrootup)
common2up<-intersect(Erootup,CRrootup)
common3up<-intersect(CRrootup,CLrootup)
EnCLrootonly<-common1up[!common1up %in% commonrootup]
EnCRrootonly<-common2up[!common2up %in% commonrootup]
CLCRrootonly<-common3up[!common3up %in% commonrootup]
commonup<-c(unique(common1up),unique(common2up),unique(common3up))
commonup <- commonup[!duplicated(commonup)]
Erootupspecific<-resrootE_up[!rownames(resrootE_up) %in% commonup, ]
CLrootupspecific<-resrootCL_up[!rownames(resrootCL_up) %in% commonup, ]
CRrootupspecific<-resrootCR_up[!rownames(resrootCR_up) %in% commonup, ]

commonrootdown<-Reduce(intersect, list(Erootdown,CLrootdown,CRrootdown))
common1down<-intersect(Erootdown,CLrootdown)
common2down<-intersect(Erootdown,CRrootdown)
common3down<-intersect(CRrootdown,CLrootdown)
EnCLrootdownonly<-common1down[!common1down %in% commonrootdown]
CLCRrootdownonly<-common3down[!common3down %in% commonrootdown]

commondown<-c(unique(common1down),unique(common2down),unique(common3down))
commondown <- commondown[!duplicated(commondown)]
CLrootdownspecific<-resrootCL_down[!rownames(resrootCL_down) %in% commondown, ]
CRrootdownspecific<-resrootCR_down[!rownames(resrootCR_down) %in% commondown, ]

#####################################
#overlap between tissues
rootclassup <- c(unique(Erootup),unique(CRrootup),unique(CLrootup))
rootclassup <- rootclassup[!duplicated(rootclassup)]
rootclassdown <- c(unique(Erootdown),unique(CRrootdown),unique(CLrootdown))
rootclassdown <- rootclassdown[!duplicated(rootclassdown)]
needleclassup <- c(unique(Enup),unique(CRup),unique(CLup))
needleclassup <- needleclassup[!duplicated(needleclassup)]
needleclassdown <- c(unique(Endown),unique(CRdown),unique(CLdown))
needleclassdown <- needleclassdown[!duplicated(needleclassdown)]
plot.new()
grid.draw(venn.diagram(list(needleclass,rootclass),scaled=TRUE,filename=NULL,lty="blank",
                       fill=c("green","brown"),category.names=c("needles","roots")))
plot.new()
grid.draw(venn.diagram(list(needleclassup,needleclassdown,rootclassup,rootclassdown), filename=NULL, category = c("needleup","needledown","rootup","rootdown"),  fill=c("green","darkgreen","#663300","#CC9933"),   ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed", margin = 0.2))

#####################################################################################################
#Saving DEGs list 
rootclassup
rootclassdown 
needleclassup 
needleclassdown

write.table(rootclassup, file = "DEGs_root_4FC_up.csv", quote = F, sep = "\t", row.names = F, col.names=F)
write.table(rootclassdown, file = "DEGs_root_4FC_down.csv", quote = F, sep = "\t", row.names = F, col.names=F)
write.table(needleclassup, file = "DEGs_needles_4FC_up.csv", quote = F, sep = "\t", row.names = F, col.names=F)
write.table(needleclassdown, file = "DEGs_needles_4FC_down.csv", quote = F, sep = "\t", row.names = F, col.names=F)

#####################################################################################################
#overlap between tissues and phases
length(intersect(CE,CEroot))#4 (98;98)
length(intersect(CL,CLroot))#595 (1304;5835)
length(intersect(CR,CRroot))#3 (23;853)

plot.new()
grid.draw(venn.diagram(list(Enup,Endown,Erootup,Erootdown), filename=NULL, category = c("needleup","needledown","rootup","rootdown"),  fill=c("green","darkgreen","#663300","#CC9933"),   ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed", margin = 0.2))

plot.new()
grid.draw(venn.diagram(list(CLup,CLdown,CLrootup,CLrootdown), filename=NULL, category = c("needleup","needledown","rootup","rootdown"),  fill=c("green","darkgreen","#663300","#CC9933"),   ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed", margin = 0.2))

plot.new()
grid.draw(venn.diagram(list(CRup,CRdown,CRrootup,CRrootdown), filename=NULL, category = c("needleup","needledown","rootup","rootdown"),  fill=c("green","darkgreen","#663300","#CC9933"),   ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed", margin = 0.2))


updown<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/venn_phasesupdown.csv",sep=";", head=T)
updown$treatment <- factor(updown$treatment, levels=c("mild","severe", "reirrigation"), ordered=T)
updown$origin <- factor(updown$origin,levels=c("needle","root","common"),ordered=T)
updown$direction <- factor(updown$direction,levels=c("up","down"),ordered=T)

p1<-ggplot(updown, aes(direction, count, fill = treatment, label=count)) + 
  geom_bar( stat="identity",position="stack") +geom_text( family="sans", size=3, position = position_stack(vjust = 0.5))+ facet_grid(.~origin)+scale_fill_manual(values = c("#669966","#FFCC66","#CCCCCC"))+
  theme(axis.text.y=element_text(size=8, angle=90, hjust=1)) +
  labs(x="Tissue",y="#DEG")+theme(strip.text.x = element_text(size=8, angle = 90))+
  theme(text=element_text(size=8,  family="Arial Unicode MS"))
p1<-p1+ theme_bw()+theme(strip.background =element_rect(fill="white"))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.ticks.x = element_blank())
p1
g <- ggplot_gtable(ggplot_build(p1))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#00FF00","#993300","#FFFFFF")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

phases<-read.csv("/mnt/picea/home/julia/Spruce/DATA/drought/needles/preprocessed/publication/venn_phases.csv",sep=";", head=T)
phases$treatment <- factor(phases$treatment, levels=c("mild","severe", "re-irrigation"), ordered=T)
phases$origin <- factor(phases$origin,levels=c("needle","root","common"),ordered=T)
p<-ggplot(phases, aes(origin, count, fill = treatment, label=count)) + 
  geom_bar( stat="identity",position="stack") +geom_text( family="sans", size=3, position = position_stack(vjust = 0.5))+ scale_fill_manual(values = c("#669966","#FFCC66","#CCCCCC"))+
  theme(axis.text.y=element_text(angle=90, hjust=1)) +
  labs(x="Tissue",y="#DEG")+theme(strip.text.x = element_text(angle = 90))+
  theme(text=element_text(size=8,  family="Arial Unicode MS"))
p+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.ticks.x = element_blank())

#######################################
#heatmaps
#up
exclude1<-c(unique(needleclassdown),unique(rootclassup),unique(rootclassdown))
exclude1 <- exclude1[!duplicated(exclude1)]
needleuponly<-needleclassup[!needleclassup %in% exclude1]
Eneedleuponly<-intersect(needleuponly,rownames(Enupspecific))
CLneedleuponly<-intersect(needleuponly,rownames(CLupspecific))
CRneedleuponly<-intersect(needleuponly,rownames(CRupspecific))
ECLneedleuponly<-intersect(needleuponly,EnCLonly)
ECRneedleuponly<-intersect(needleuponly,EnCRonly)
CLCRneedleuponly<-intersect(needleuponly,CLCRonly)
commonneedleuponly<-intersect(needleuponly,commonneedleup)

needlerootup<-intersect(needleclassup,rootclassup)

exclude2<-c(unique(needleclassdown),unique(needleclassup),unique(rootclassdown))
exclude2 <- exclude2[!duplicated(exclude2)]
rootuponly<-rootclassup[!rootclassup %in% exclude2]
Erootuponly<-intersect(rootuponly,rownames(Erootupspecific))
CLrootuponly<-intersect(rootuponly,rownames(CLrootupspecific))
CRrootuponly<-intersect(rootuponly,rownames(CRrootupspecific))
ECLrootuponly<-intersect(rootuponly,EnCLrootonly)
ECRrootuponly<-intersect(rootuponly,EnCRrootonly)
CLCRrootuponly<-intersect(rootuponly,CLCRrootonly)
commonrootuponly<-intersect(rootuponly,commonrootup)

Eneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% Eneedleuponly,]#56
CLneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLneedleuponly,]#400
CRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRneedleuponly,]#4
ECLneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLneedleuponly,]#3
ECRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECRneedleuponly,]#1
CLCRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLCRneedleuponly,]#1
commonneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% commonneedleuponly,]#2
Erootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% Erootuponly,]#7
CLrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLrootuponly,]#1236
CRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRrootuponly,]#278
ECLrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLrootuponly,]#29
ECRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECRrootuponly,]#2
CLCRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLCRrootuponly,]#192
commonrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% commonrootuponly,]#17
needlerootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% needlerootup,]#451

#######down
exclude3<-c(unique(needleclassup),unique(rootclassdown),unique(rootclassup))
exclude3 <- exclude3[!duplicated(exclude3)]
needledownonly<-needleclassdown[!needleclassdown %in% exclude3]
Eneedledownonly<-intersect(needledownonly,rownames(Endownspecific))
CLneedledownonly<-intersect(needledownonly,rownames(CLdownspecific))
CRneedledownonly<-intersect(needledownonly,rownames(resCR_down))
common1down<-intersect(Endown,CLdown)
ECLneedledownonly<-intersect(needledownonly,common1down)

needlerootdown<-intersect(needleclassdown,rootclassdown)
needlerootdown<-needlerootdown[!needlerootdown %in% rootclassup]

exclude4<-c(unique(needleclassdown),unique(needleclassup),unique(rootclassup))
exclude4 <- exclude4[!duplicated(exclude4)]
rootdownonly<-rootclassdown[!rootclassdown %in% exclude4]

CLrootdownonly<-intersect(rootdownonly,rownames(CLrootdownspecific))#3531
CRrootdownonly<-intersect(rootdownonly,rownames(CRrootdownspecific))#58
ECLrootdownonly<-intersect(rootdownonly,EnCLrootdownonly)#5
CLCRrootdownonly<-intersect(rootdownonly,CLCRrootdownonly)#183    
commonrootdownonly<-intersect(rootdownonly,commonrootdown)#2


Eneedledownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% Eneedledownonly,]#3
CLneedledownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLneedledownonly,]#290
CRneedledownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRneedledownonly,]#5
ECLneedledownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLneedledownonly,]#1

CLrootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLrootdownonly,]
CRrootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRrootdownonly,]
ECLrootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLrootdownonly,]
CLCRrootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLCRrootdownonly,]
commonrootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% commonrootdownonly,]
needlerootdownonlyvst<-myMAmean_2[rownames(myMAmean_2) %in% needlerootdown,]#130

#genes with variable expression in the same tissue or between tissues
opp<-intersect(needleclassup,rootclassdown)
opp_2<-intersect(needleclassdown,rootclassup)
opp_3<-c(unique(opp),unique(opp_2))
N_R_opp<-myMAmean_2[rownames(myMAmean_2) %in% opp_3,]#55

N_opp<-intersect(needleclassup,needleclassdown)
N_oppvst<-myMAmean_2[rownames(myMAmean_2) %in% N_opp,]  #1
dim(N_oppvst)

R_opp<-intersect(rootclassup,rootclassdown)
R_opp2<-R_opp[!R_opp %in% opp_3]
dim(R_opp2)

R_oppvst<-myMAmean_2[rownames(myMAmean_2) %in% R_opp2,]#18
dim(R_oppvst)
#heatmap up-regulated genes
vst<-rbind(Eneedleuponlyvst,ECLneedleuponlyvst,CLneedleuponlyvst,CLCRneedleuponlyvst,CRneedleuponlyvst,ECRneedleuponlyvst,commonneedleuponlyvst,needlerootuponlyvst,Erootuponlyvst,ECLrootuponlyvst,CLrootuponlyvst,CLCRrootuponlyvst,CRrootuponlyvst,ECRrootuponlyvst,commonrootuponlyvst)
annotation_row = data.frame( GeneClass = factor(rep(c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o"), c(56,3,400,1,4,1,2,451,7,29,1236,192,278,2,17))))
head(annotation_row)
head(vst)

dim(annotation_row)
dim(vst)

#cbind(annotation_row, vst)
#class(vst)
names <- colnames(vst) 

dim(annotation_row)
dim(vst)

#rownames(annotation_row) = rownames(vst)

#rowscaling separate for the tissues 
vst_1<-vst[,1:8]
vst_1<-vst_1-rowMeans(vst_1)
vst_2<-vst[,9:16]
vst_2<-vst_2-rowMeans(vst_2)
vst<-cbind(vst_1,vst_2)
summary(vst)


df = data.frame(type = factor(rep(c("control","M60","M40","M30","M307d","S","S48h","R"), 1,Time=1:11)),
                age = factor(rep(c(rep("Needle", 8)))))
ha_a = HeatmapAnnotation(df = df,
                         col = list(type=c("control" =  "#1B9E77", "M60" = "#D95F02","M40"="#7570B3","M30"="#E7298A","M307d"="#66A61E","S"="#FFCC00","S48h"="#996600","R"="#666666"),
                                    age =  c("Needle"="green")))

df = data.frame(type = factor(rep(c("control","M60","M40","M30","M307d","S","S48h","R"), 1,Time=1:11)),
                age = factor(rep(c(rep("Root", 8)))))
ha_b = HeatmapAnnotation(df = df,
                         col = list(type=c("control" =  "#1B9E77", "M60" = "#D95F02","M40"="#7570B3","M30"="#E7298A","M307d"="#66A61E","S"="#FFCC00","S48h"="#996600","R"="#666666"),
                                    age =  c("Root"="brown")))

#original code:
ht1 = Heatmap(vst_1,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,0,0,0,1,1,0,0,0,0,0,0),"mm"),top_annotation = ha_a,top_annotation_height = unit(3, "mm"))
ht2 = Heatmap(vst_2,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,0,0,0,1,1,0,0,0,0,0,0),"mm"),show_heatmap_legend = FALSE,top_annotation=ha_b,top_annotation_height = unit(5, "mm"))


#ht1 = Heatmap(vst_1,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,0,0,0,1,1,0,0,0,0,0,0),"mm"))
#ht2 = Heatmap(vst_2,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,0,0,0,1,1,0,0,0,0,0,0),"mm"))   #,show_heatmap_legend = FALSE,top_annotation=ha_b,top_annotation_height = unit(5, "mm"))


df = data.frame(type = factor(c(rep("a",56),rep("b",3),rep ("c",400),rep("d",1),rep("e",4),rep("f",1),rep("g",2),rep("h",451),rep("i",7),rep("j",29),rep("k",1236),rep("l",192),rep("m",278),rep("n",2),rep("o",17))))
colors<-c("#CCCC33","#FFFF00","#CCCC99" ,"#CCCCCC","#33FF00","#339900","#CCFF33","#FFFFCC","#FFCC00","#FF9900","#FFCC99","#996666","#CC6633","#CC3300","#CC9900")
ht_11<-Heatmap(df,col =colors  ,width = unit(5, "mm"))                             

ht_list = ht1 + ht2 
ht_list+ht_11
#down
vst_d<-rbind(Eneedledownonlyvst,ECLneedledownonlyvst,CLneedledownonlyvst,CRneedledownonlyvst,needlerootdownonlyvst,ECLrootdownonlyvst,CLrootdownonlyvst,CLCRrootdownonlyvst,CRrootdownonlyvst,commonrootdownonlyvst,N_R_opp,N_oppvst,R_oppvst)

annotation_row = data.frame( GeneClass = factor(rep(c("a", "b", "c","d","e","f","g","h","i","j","k","l","m"), c(3,1,290,5,130,5,3531,183,58,2,55,1,18))))
rownames(annotation_row) = rownames(vst_d)
#rowscaling if not log transformed, either separate for the tissues or together
vst_3<-vst_d[,1:8]
vst_3<-vst_3-rowMeans(vst_3)
vst_4<-vst_d[,9:16]
vst_4<-vst_4-rowMeans(vst_4)
vst_d<-cbind(vst_3,vst_4)

summary(vst_d)
ht3<-Heatmap(vst_3,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,1,1,0,0,0,0,1,0,0),"mm"),top_annotation = ha_a,top_annotation_height = unit(3, "mm"))
ht4<-Heatmap(vst_4,show_row_names=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,split = paste0("", annotation_row$GeneClass),col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),gap=unit(c(0,0,0,1,1,0,0,0,0,1,0,0),"mm"),show_heatmap_legend = FALSE,top_annotation=ha_b,top_annotation_height = unit(3, "mm"))

ht_list2 = ht3 + ht4

dfd = data.frame(type = factor(c(rep("a",3),rep("b",1),rep ("c",290),rep("d",5),rep("e",130),rep("f",5),rep("g",3531),rep("h",183),rep("i",58),rep("j",2),rep("k",55),rep("l",1),rep("m",18))))
colors<-c("#CCCC33","#FFFF00","#CCCC99" ,"#33FF00","#FFFFCC","#FF9900","#FFCC99","#996666","#CC6633","#CC9900","#660066","#00FF99","#663300")

ht_33<-Heatmap(dfd,col =colors  ,width = unit(5, "mm"))                             

ht_list2+ht_33

############### cluster profiles of gene expression
#Eneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% Eneedleuponly,]#56
dat_mean_needles <- melt(as.matrix(Eneedleuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(Eneedleuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean1<-dat_mean

p1<-ggplot(dat_mean1,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none")


#CLneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLneedleuponly,]#400
dat_mean_needles <- melt(as.matrix(CLneedleuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLneedleuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean2<-dat_mean
p2<-ggplot(dat_mean2,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun= mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#CRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRneedleuponly,]#4
dat_mean_needles <- melt(as.matrix(CRneedleuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CRneedleuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean3<-dat_mean
p3<-ggplot(dat_mean3,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#ECLneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLneedleuponly,]#3
dat_mean_needles <- melt(as.matrix(ECLneedleuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECLneedleuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean4<-dat_mean
p4<-ggplot(dat_mean4,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#ECRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECRneedleuponly,]#1
dat_mean_needles <- melt(as.matrix(ECRneedleuponlyvst[c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECRneedleuponlyvst[c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Gene <- sub("\\ .*","\\",dat_mean$Gene)
dat_mean$Gene <- factor(dat_mean$Gene, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean5<-dat_mean
p5<-ggplot(dat_mean5,aes(x=factor(Gene),y=VST, group= Tissue, color=Tissue)) + geom_line()+
  stat_summary(data= dat_mean, geom = "line",size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 


#CLCRneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLCRneedleuponly,]#1
dat_mean_needles <- melt(as.matrix(CLCRneedleuponlyvst[c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLCRneedleuponlyvst[c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Gene <- sub("\\ .*","\\",dat_mean$Gene)
dat_mean$Gene <- factor(dat_mean$Gene, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean6<-dat_mean
p6<-ggplot(dat_mean6,aes(x=factor(Gene),y=VST, group= Tissue, color=Tissue)) + geom_line()+
  stat_summary(data= dat_mean, geom = "line",size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 


#commonneedleuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% commonneedleuponly,]#2
dat_mean_needles <- melt(as.matrix(commonneedleuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(commonneedleuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean7<-dat_mean
p7<-ggplot(dat_mean7,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p7
#Erootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% Erootuponly,]#7
dat_mean_needles <- melt(as.matrix(Erootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(Erootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean8<-dat_mean
p8<-ggplot(dat_mean8,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p8
#CLrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLrootuponly,]#1236
dat_mean_needles <- melt(as.matrix(CLrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean9<-dat_mean
p9<-ggplot(dat_mean9,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#CRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CRrootuponly,]#278
dat_mean_needles <- melt(as.matrix(CRrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CRrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean10<-dat_mean
p10<-ggplot(dat_mean10,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun= mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#ECLrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECLrootuponly,]#29
dat_mean_needles <- melt(as.matrix(ECLrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECLrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean11<-dat_mean
p11<-ggplot(dat_mean11,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#ECRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% ECRrootuponly,]#2
dat_mean_needles <- melt(as.matrix(ECRrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECRrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean12<-dat_mean
p12<-ggplot(dat_mean12,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#CLCRrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% CLCRrootuponly,]#192
dat_mean_needles <- melt(as.matrix(CLCRrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLCRrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean13<-dat_mean
p13<-ggplot(dat_mean13,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun= mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#commonrootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% commonrootuponly,]#17
dat_mean_needles <- melt(as.matrix(commonrootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(commonrootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean14<-dat_mean
p14<-ggplot(dat_mean14,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

#needlerootuponlyvst<-myMAmean_2[rownames(myMAmean_2) %in% needlerootup,]#451
dat_mean_needles <- melt(as.matrix(needlerootuponlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(needlerootuponlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean15<-dat_mean
p15<-ggplot(dat_mean15,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

multiplot(p1, p2, p3, p4, p5,p6, p7,p8,p9,p10,p11,p12,p13,p14,p15, cols=3)

#profiles down
#vst_d<-rbind(Eneedledownonlyvst,ECLneedledownonlyvst,CLneedledownonlyvst,CRneedledownonlyvst,needlerootdownonlyvst,ECLrootdownonlyvst,CLrootdownonlyvst,CLCRrootdownonlyvst,CRrootdownonlyvst,commonrootdownonlyvst,N_R_oppositevst,N_updownvst,R_updownvst)

#Eneedledownonlyvst #3
dat_mean_needles <- melt(as.matrix(Eneedledownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(Eneedledownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean1d<-dat_mean
p1d<-ggplot(dat_mean1d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

p1d

#ECLneedledownonlyvst #1
dat_mean_needles <- melt(as.matrix(ECLneedledownonlyvst[c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECLneedledownonlyvst[c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Gene <- sub("\\ .*","\\",dat_mean$Gene)
dat_mean$Gene <- factor(dat_mean$Gene, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean2d<-dat_mean
p2d<-ggplot(dat_mean2d,aes(x=factor(Gene),y=VST, group= Tissue, color=Tissue)) + geom_line()+
  stat_summary(data= dat_mean, geom = "line",size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p2d
#CLneedledownonlyvst #290
dat_mean_needles <- melt(as.matrix(CLneedledownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLneedledownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean3d<-dat_mean
p3d<-ggplot(dat_mean3d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p3d
#CRneedledownonlyvst #5
dat_mean_needles <- melt(as.matrix(CRneedledownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CRneedledownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean4d<-dat_mean
p4d<-ggplot(dat_mean4d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

p4d
#needlerootdownonlyvst #130
dat_mean_needles <- melt(as.matrix(needlerootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(needlerootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean5d<-dat_mean
p5d<-ggplot(dat_mean5d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 

p5d
#ECLrootdownonlyvst #5
dat_mean_needles <- melt(as.matrix(ECLrootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(ECLrootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean6d<-dat_mean
p6d<-ggplot(dat_mean6d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p6d
#CLrootdownonlyvst #3531
dat_mean_needles <- melt(as.matrix(CLrootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLrootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean7d<-dat_mean
p7d<-ggplot(dat_mean7d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p7d
#CLCRrootdownonlyvst #183
dat_mean_needles <- melt(as.matrix(CLCRrootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CLCRrootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean8d<-dat_mean
p8d<-ggplot(dat_mean8d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p8d
#CRrootdownonlyvst#58
dat_mean_needles <- melt(as.matrix(CRrootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(CRrootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean9d<-dat_mean
p9d<-ggplot(dat_mean9d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p9d
#commonrootdownonlyvst#2
dat_mean_needles <- melt(as.matrix(commonrootdownonlyvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(commonrootdownonlyvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean10d<-dat_mean
p10d<-ggplot(dat_mean10d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p10d

#N_R_oppositevst #55
dat_mean_needles <- melt(as.matrix(N_R_opp[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(N_R_opp[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean11d<-dat_mean
p11d<-ggplot(dat_mean11d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p11d

#N_updownvst #1
dat_mean_needles <- melt(as.matrix(N_oppvst[c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(N_oppvst[c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Gene <- sub("\\ .*","\\",dat_mean$Gene)
dat_mean$Gene <- factor(dat_mean$Gene, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean12d<-dat_mean
p12d<-ggplot(dat_mean12d,aes(x=factor(Gene),y=VST, group= Tissue, color=Tissue)) + geom_line()+
  stat_summary(data= dat_mean, geom = "line",size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p12d
#R_updownvst #18
dat_mean_needles <- melt(as.matrix(R_oppvst[,c(1:8)]))
names(dat_mean_needles) <- c("Gene", "Treatment", "VST")
dat_mean_needles["Tissue"]<-"Needle"
dat_mean_needles2 <- melt(as.matrix(R_oppvst[,c(9:16)]))
names(dat_mean_needles2) <- c("Gene", "Treatment", "VST")
dat_mean_needles2["Tissue"]<-"Root"
dat_mean<-rbind(dat_mean_needles,dat_mean_needles2)
dat_mean$Treatment <- sub("\\ .*","\\",dat_mean$Treatment)
dat_mean$Treatment <- factor(dat_mean$Treatment, levels=c("80%","60%","40%","30%","30%7d","Collapse","C2d","Rehydrate"))
dat_mean$Tissue <- factor(dat_mean$Tissue, levels=c("Root","Needle"))
dat_mean13d<-dat_mean
p13d<-ggplot(dat_mean13d,aes(x=factor(Treatment),y=VST, group= Tissue, color=Tissue, geom="line")) + 
  stat_summary(fun = mean, geom = "line",size = 0.5)+
  stat_summary(fun.data="mean_cl_normal", geom='smooth', se=T)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) + theme(legend.position="none") 
p13d
multiplot(p1d, p2d, p3d, p4d, p5d,p6d, p7d,p8d,p9d,p10d,p11d,p12d,p13d, cols=3)

sessionInfo()
