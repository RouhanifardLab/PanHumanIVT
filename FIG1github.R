## ========================================== ##
##               Figure 1 Code                ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 03/28/2023
#@ version: 1.0.0

# README ------------------------------------------------------------------
# Software Requirement: R version >= 4.1 
# 
# Inputs were generated via command line, they are commented above the introduction of these inputs
# CL packages: samtools, NanoPlot

# R: Software Requirements ------------------------------------------------
list.of.packages <- c("ggplot2", "ggExtra","dplyr", "tidyr", "plyr", "readr", "stringr", "tibble", "stringr", "reshape2", "limma", "vtable", "caret", "venneuler")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) # from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them. Thanks!
library(ggplot2)
library(ggExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(reshape2)
library(limma)
library(caret)
library(venneuler)
library(vtable)

# DATA FRAME: Aligned Gene Counts  -----------------------------------------

##--COMMAND LINE--##
#columns: readID | gene name | readLENGTH | MAPQ
#samtools view -F 4 -F 256 -F 2048 <PATH TO SAM> | awk '{ split($3,a,"|"); print $1,a[6],length($10),$5}' > <PATH TO OUTPUT TXT FILE>

HeLa.rep1.MAP <- read.table(file = paste(dataPath,"HeLa_IVT_rep1.gencodev38.filtered.txt", sep = "/"), sep = '', header = FALSE)
HeLa.rep2.MAP <- read.table(file = paste(dataPath,"HeLa_IVT_rep2.gencodev38.filtered.txt", sep = "/"), sep = '', header = FALSE)
HeLa.rep3.MAP <- read.table(file = paste(dataPath,"HeLa_IVT_rep3.gencodev38.filtered.txt", sep = "/"), sep = '', header = FALSE)

colName <- c("readID", "annot", "pos", "length", "MAPQ")
colnames(HeLa.rep1.MAP) <- colName
colnames(HeLa.rep2.MAP) <- colName
colnames(HeLa.rep3.MAP) <- colName

HeLa.ALL <- rbind(HeLa.rep1.MAP, HeLa.rep2.MAP, HeLa.rep3.MAP)
HeLa.MAP <- data.frame(annot = unique(HeLa.ALL$annot))

#if (T) {
HeLa.MAP$HeLa.rep1 = -1
HeLa.MAP$HeLa.rep2 = -1
HeLa.MAP$HeLa.rep3 = -1

for (n in c(1:nrow(HeLa.MAP))) {
  if(n %% 1000 == 0) { print (n/nrow(HeLa.MAP)*100)}
  HeLa.MAP$HeLa.rep1[n] <- length(HeLa.rep1.MAP[which(HeLa.rep1.MAP$annot == HeLa.MAP$annot[n]),"annot"])
  HeLa.MAP$HeLa.rep2[n] <- length(HeLa.rep2.MAP[which(HeLa.rep2.MAP$annot == HeLa.MAP$annot[n]),"annot"])
  HeLa.MAP$HeLa.rep3[n] <- length(HeLa.rep3.MAP[which(HeLa.rep3.MAP$annot == HeLa.MAP$annot[n]),"annot"])
}

###---SUM REPLICATES---###
HeLa.MAP$HeLa.total = -1
HeLa.MAP$HeLa.total = HeLa.MAP$HeLa.rep1 + HeLa.MAP$HeLa.rep2 + HeLa.MAP$HeLa.rep3
HeLa.MAP$HeLa.total[HeLa.MAP$HeLa.total< 0] <- 0

HeLa.MAP <- HeLa.MAP[which(HeLa.MAP$HeLa.total != 0),]
#}  

# FIGURE 1A: %identity v length ------------------------------------------

##--COMMAND LINE--##
#NanoPlot -t 64 --raw -o <PATH TO OUTPUT DIRECTORY> --bam <PATH TO .bam FILE> 

HeLa.rep1.nanoP <- read_tsv(print(paste(dataPath, "HeLa_IVT_rep1_NanoPlot-data.tsv", sep = "/")))
HeLa.rep2.nanoP <- read_tsv(print(paste(dataPath, "HeLa_IVT_rep2_NanoPlot-data.tsv", sep = "/")))
HeLa.rep3.nanoP <- read_tsv(print(paste(dataPath, "HeLa_IVT_rep3_NanoPlot-data.tsv", sep = "/")))

##--Rep. 1--##
AL <- HeLa.rep1.nanoP$aligned_lengths
PI <- HeLa.rep1.nanoP$percentIdentity
smoothScatter(AL,PI,ylim=c(30,100), xlim=c(0,5000))


##--Rep. 2--##
AL <- HeLa.rep2.nanoP$aligned_lengths
PI <- HeLa.rep2.nanoP$percentIdentity
smoothScatter(AL,PI, ylim=c(30,100), xlim=c(0,5000))


##--Rep. 3--##
AL <- HeLa.rep3.nanoP$aligned_lengths
PI <- HeLa.rep3.nanoP$percentIdentity
smoothScatter(AL,PI, ylim=c(30,100), xlim=c(0,5000))


# FIGURE 1B: Heat Map Quartiles----------------------------------------------
###---DATA FRAME---###
minCount = 1
genes.quart <- HeLa.MAP[,c("annot", "HeLa.rep1", "HeLa.rep2", "HeLa.rep3")]
genes.quart$tot <- genes.quart$HeLa.rep1 + genes.quart$HeLa.rep2 + genes.quart$HeLa.rep3
genes.quart <- genes.quart[which(genes.quart$tot >= minCount),]

###---QUARTILES FROM TOTAL COUNT---###
tbl <- sumtable(genes.quart, out ='return')
as.numeric(tbl$Mean[which(tbl$Variable == "tot")])

quart1 <- genes.quart[which(genes.quart$tot <= 3), 1:4]
quart2 <- genes.quart[which(genes.quart$tot >= 4 & genes.quart$tot <= 24), 1:4]
quart3 <- genes.quart[which(genes.quart$tot >= 25 & genes.quart$tot <= 105), 1:4]
quart4 <- genes.quart[which(genes.quart$tot >= 106), 1:4]

###---RANDOM SAMPLING---###
set.seed(1)
quart1.rs <- sample_n(quart1, 25)
set.seed(2)
quart2.rs <- sample_n(quart2, 25)
set.seed(3)
quart3.rs <- sample_n(quart3, 25)
set.seed(4)
quart4.rs <- sample_n(quart4, 25)

###---BINARY---###
for (j in 1:nrow(quart1.rs)) {
  if (quart1.rs$HeLa.rep1[j] >= minCount) {quart1.rs$HeLa.rep1[j] = 1} else {quart1.rs$HeLa.rep1[j] = 0}
  if (quart1.rs$HeLa.rep2[j] >= minCount) {quart1.rs$HeLa.rep2[j] = 1} else {quart1.rs$HeLa.rep2[j] = 0}
  if (quart1.rs$HeLa.rep3[j] >= minCount) {quart1.rs$HeLa.rep3[j] = 1} else {quart1.rs$HeLa.rep3[j] = 0} }

for (j in 1:nrow(quart2.rs)) {
  if (quart2.rs$HeLa.rep1[j] >= minCount) {quart2.rs$HeLa.rep1[j] = 1} else {quart2.rs$HeLa.rep1[j] = 0}
  if (quart2.rs$HeLa.rep2[j] >= minCount) {quart2.rs$HeLa.rep2[j] = 1} else {quart2.rs$HeLa.rep2[j] = 0}
  if (quart2.rs$HeLa.rep3[j] >= minCount) {quart2.rs$HeLa.rep3[j] = 1} else {quart2.rs$HeLa.rep3[j] = 0} }

for (j in 1:nrow(quart3.rs)) {
  if (quart3.rs$HeLa.rep1[j] >= minCount) {quart3.rs$HeLa.rep1[j] = 1} else {quart3.rs$HeLa.rep1[j] = 0}
  if (quart3.rs$HeLa.rep2[j] >= minCount) {quart3.rs$HeLa.rep2[j] = 1} else {quart3.rs$HeLa.rep2[j] = 0}
  if (quart3.rs$HeLa.rep3[j] >= minCount) {quart3.rs$HeLa.rep3[j] = 1} else {quart3.rs$HeLa.rep3[j] = 0} }

for (j in 1:nrow(quart4.rs)) {
  if (quart4.rs$HeLa.rep1[j] >= minCount) {quart4.rs$HeLa.rep1[j] = 1} else {quart4.rs$HeLa.rep1[j] = 0}
  if (quart4.rs$HeLa.rep2[j] >= minCount) {quart4.rs$HeLa.rep2[j] = 1} else {quart4.rs$HeLa.rep2[j] = 0}
  if (quart4.rs$HeLa.rep3[j] >= minCount) {quart4.rs$HeLa.rep3[j] = 1} else {quart4.rs$HeLa.rep3[j] = 0} }

###---HEAT MAPS---###
colors <- c("#443983", "#35b779")
colors2 <- c("#35b779", "#443983")

##--QUARTILE 1--##
quart1.rs.melt <- melt(quart1.rs, id.vars = "annot")

q1 <- quart1.rs.melt %>% 
  ggplot(aes(x = variable , y = annot, fill = factor(value))) +
  geom_tile()+
  scale_fill_manual(values=colors)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        #       # panel.background=element_blank(),
        #       # panel.border=element_blank(),
        #       # panel.grid.major=element_blank(),
        #       # panel.grid.minor=element_blank(),
        plot.background=element_blank())
q1

##--QUARTILE 2--##
quart2.rs.melt <- melt(quart2.rs, id.vars = "annot")

q2 <- quart2.rs.melt %>% 
  ggplot(aes(x = variable , y = annot, fill = factor(value))) +
  geom_tile()+
  scale_fill_manual(values=colors)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        #       # panel.background=element_blank(),
        #       # panel.border=element_blank(),
        #       # panel.grid.major=element_blank(),
        #       # panel.grid.minor=element_blank(),
        plot.background=element_blank())
q2

##--QUARTILE 3--##

quart3.rs.melt <- melt(quart3.rs, id.vars = "annot")

q3 <- quart3.rs.melt %>% 
  ggplot(aes(x = variable , y = annot, fill = factor(value))) +
  geom_tile()+
  scale_fill_manual(values=colors2)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        #       # panel.background=element_blank(),
        #       # panel.border=element_blank(),
        #       # panel.grid.major=element_blank(),
        #       # panel.grid.minor=element_blank(),
        plot.background=element_blank())
q3

##--QUARTILE 4--##
quart4.rs.melt <- melt(quart4.rs, id.vars = "annot")

q4 <- quart4.rs.melt %>% 
  ggplot(aes(x = variable , y = annot, fill = factor(value))) +
  geom_tile()+
  scale_fill_manual(values=colors2)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        #       # panel.background=element_blank(),
        #       # panel.border=element_blank(),
        #       # panel.grid.major=element_blank(),
        #       # panel.grid.minor=element_blank(),
        plot.background=element_blank())
q4



# FIGURE 1C: Percent of Genes Covered by All ---------------------------------
if (T){
  numGenes <- data.frame("annot" = HeLa.MAP$annot,"HeLa.rep1" = HeLa.MAP$HeLa.rep1, "HeLa.rep2" = HeLa.MAP$HeLa.rep2, "HeLa.rep3" = HeLa.MAP$HeLa.rep3)
  numGenes$IVTtot = numGenes$HeLa.rep1 + numGenes$HeLa.rep2 + numGenes$HeLa.rep3
  
  numrows <- 10
  nGenes <- data.frame(matrix(ncol = 7, nrow = numrows))
  x <- c("minReads","HeLa.rep1", "HeLa.rep2", "HeLa.rep3", "Ngenes.all", "Ngenes.rep", "Ngenes.init")
  colnames(nGenes) <- x
  
  for (i in 1:10) {
    minReads = i
    nGenes$minReads[i] = i
    nGenes$Ngenes.init[i] = nrow(numGenes)
    
    numGenes.temp <- numGenes[which(numGenes$IVTtot >= minReads),]
    
    HeLa.rep1 <- numGenes.temp$annot[which(numGenes.temp$HeLa.rep1 > minReads)]
    HeLa.rep2 <- numGenes.temp$annot[which(numGenes.temp$HeLa.rep2 > minReads)]
    HeLa.rep3 <- numGenes.temp$annot[which(numGenes.temp$HeLa.rep3 > minReads)]
    
    # Assign Variables
    nGenes$HeLa.rep1[i] = length(HeLa.rep1)
    nGenes$HeLa.rep2[i] = length(HeLa.rep2)
    nGenes$HeLa.rep3[i] = length(HeLa.rep3)
    nGenes$Ngenes.rep[i] = nrow(numGenes.temp)
  }
  
  
  GC <- nGenes
  GC$HeLa.rep1perc <- (GC$HeLa.rep1/GC$Ngenes.init)*100
  GC$HeLa.rep2perc <- (GC$HeLa.rep2/GC$Ngenes.init)*100
  GC$HeLa.rep3perc <- (GC$HeLa.rep3/GC$Ngenes.init)*100
  GC$Ngenes.repperc <- (GC$Ngenes.rep/GC$Ngenes.init)*100
  
  GC.perc <- GC[,c(1,8:11)]
  
  coverageMelt <- melt(GC.perc, id.vars = "minReads")
  cols <- c("#440154", "#21918C", "#90d743", "#38588c")
}  

numGenesVcount <- coverageMelt %>% 
  ggplot(aes(x = minReads, y = value, color = variable)) +
  geom_line()+
  #geom_smooth(se = FALSE, method = "auto", linewidth = 0.5) + 
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(0, 10, by = 2)), labels = c(seq(0, 10, by = 2)))+
  ylim(30,100)+
  # scale_y_continuous(limits =c(40, 100),
  #                    breaks = seq(40, 100, by = 10),
  #                    labels = c(0, seq(50, 100, by = 10)),
  #                    expand = c(0,0,0.05,0)) +
  theme(axis.line=element_line(linewidth = 0.33),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        #       # panel.background=element_blank(),
        #       # panel.border=element_blank(),
        #       # panel.grid.major=element_blank(),
        #       # panel.grid.minor=element_blank(),
        plot.background=element_blank())

numGenesVcount

# FIGURE 1D: Replicate Contributions --------------------------------------
minReadCount <- data.frame("annot" = HeLa.MAP$annot,"HeLa.rep1" = HeLa.MAP$HeLa.rep1, "HeLa.rep2" = HeLa.MAP$HeLa.rep2, "HeLa.rep3" = HeLa.MAP$HeLa.rep3)
minReadCount$IVTtot = minReadCount$HeLa.rep1 + minReadCount$HeLa.rep2 + minReadCount$HeLa.rep3


numrows <- 20
GeneCountVD <- data.frame(matrix(ncol = 8, nrow = numrows))
x <- c("minReads","HeLa.rep1", "HeLa.rep2", "HeLa.rep3", "HeLa.rep12", "HeLa.rep13", "HeLa.rep23", "HeLa.total")
colnames(GeneCountVD) <- x


for (i in 1:20) {
  minReads = i
  GeneCountVD$minReads[i] = i
  minReadCount.temp <- minReadCount[which(minReadCount$IVTtot >= minReads),]
  
  data1 <- minReadCount.temp$annot[which(minReadCount.temp$HeLa.rep1 > 0)]
  data2 <- minReadCount.temp$annot[which(minReadCount.temp$HeLa.rep2 > 0)]
  data3 <- minReadCount.temp$annot[which(minReadCount.temp$HeLa.rep3 > 0)]
  
  # Find the overlap between all replicates
  overlap123 <- Reduce(intersect,(list(data1, data2, data3)))
  
  # Find the overlap between two replicates but not the third
  overlap12 <- intersect(data1, data2)
  overlap12 <- setdiff(overlap12, data3)
  
  overlap13 <- intersect(data1, data3)
  overlap13 <- setdiff(overlap13, data2)
  
  overlap23 <- intersect(data2, data3)
  overlap23 <- setdiff(overlap23, data1)
  
  # Find the genes in one replicate that are unique 
  unique1 <- setdiff(data1, union(data2, data3))
  unique2 <- setdiff(data2, union(data1, data3))
  unique3 <- setdiff(data3, union(data1, data2))
  
  # Assign Variables
  GeneCountVD$HeLa.rep1[i] = length(unique1)
  GeneCountVD$HeLa.rep2[i] = length(unique2)
  GeneCountVD$HeLa.rep3[i] = length(unique3)
  GeneCountVD$HeLa.rep12[i] = length(overlap12)
  GeneCountVD$HeLa.rep13[i] = length(overlap13)
  GeneCountVD$HeLa.rep23[i] = length(overlap23)
  GeneCountVD$HeLa.total[i] = length(overlap123)
  GeneCountVD$total[i] = nrow(minReadCount.temp)
}

